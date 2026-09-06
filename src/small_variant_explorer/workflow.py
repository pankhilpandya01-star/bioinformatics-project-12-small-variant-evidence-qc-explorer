"""Atomic duplicate marking and quality-filtered depth orchestration."""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import os
import platform
import re
import shutil
import subprocess
import tempfile
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any

from small_variant_explorer import __version__
from small_variant_explorer.depth import collect_depth, write_depth_outputs
from small_variant_explorer.validation import (
    ValidationError,
    parse_bam_header,
    read_fasta_lengths,
)


EXPECTED_SAMTOOLS_VERSION = "1.24"
PRIVATE_MARKERS = ("/home/", "/mnt/", "/users/", "c:" + "\\users\\")


class WorkflowError(RuntimeError):
    """Raised when BAM preparation or its verification fails."""


@dataclass(frozen=True, slots=True)
class WorkflowConfig:
    reference_fasta: Path
    input_bam: Path
    sample_id: str
    output_dir: Path
    ploidy: int = 1
    minimum_mapq: int = 20
    minimum_base_quality: int = 20
    minimum_depth: int = 10
    minimum_alt_depth: int = 5
    minimum_alt_fraction: float = 0.8
    minimum_variant_quality: float = 30.0
    minimum_alt_strand_depth: int = 1
    maximum_depth_factor: float = 3.0
    threads: int = 1


@dataclass(frozen=True, slots=True)
class ToolPaths:
    samtools: str


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_config(config: WorkflowConfig) -> WorkflowConfig:
    reference = config.reference_fasta.expanduser().resolve()
    input_bam = config.input_bam.expanduser().resolve()
    output = config.output_dir.expanduser().resolve()
    for label, path in (("reference FASTA", reference), ("input BAM", input_bam)):
        if not path.exists() or not path.is_file():
            raise WorkflowError(f"{label} is not a regular file: {path}")
    if output.exists():
        raise WorkflowError(f"output directory already exists: {output}")
    sample_id = config.sample_id.strip()
    if not sample_id:
        raise WorkflowError("sample ID must not be empty")
    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample_id) is None:
        raise WorkflowError(
            "sample ID must start with a letter or number and contain only "
            "letters, numbers, dots, underscores, or hyphens"
        )
    if config.ploidy not in (1, 2) or isinstance(config.ploidy, bool):
        raise WorkflowError("ploidy must be 1 or 2")
    for label, value, minimum, maximum in (
        ("minimum MAPQ", config.minimum_mapq, 0, 255),
        ("minimum base quality", config.minimum_base_quality, 0, 255),
    ):
        if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
            raise WorkflowError(f"{label} must be an integer from {minimum} through {maximum}")
    for label, value, minimum in (
        ("minimum depth", config.minimum_depth, 1),
        ("minimum alternate depth", config.minimum_alt_depth, 1),
        ("minimum alternate strand depth", config.minimum_alt_strand_depth, 0),
        ("threads", config.threads, 1),
    ):
        if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
            raise WorkflowError(f"{label} must be an integer of at least {minimum}")
    if (
        isinstance(config.minimum_alt_fraction, bool)
        or not isinstance(config.minimum_alt_fraction, (int, float))
        or not math.isfinite(config.minimum_alt_fraction)
        or not 0 < config.minimum_alt_fraction <= 1
    ):
        raise WorkflowError("minimum alternate fraction must be greater than 0 and at most 1")
    if (
        isinstance(config.minimum_variant_quality, bool)
        or not isinstance(config.minimum_variant_quality, (int, float))
        or not math.isfinite(config.minimum_variant_quality)
        or config.minimum_variant_quality < 0
    ):
        raise WorkflowError("minimum variant quality must be non-negative")
    if (
        isinstance(config.maximum_depth_factor, bool)
        or not isinstance(config.maximum_depth_factor, (int, float))
        or not math.isfinite(config.maximum_depth_factor)
        or config.maximum_depth_factor <= 1
    ):
        raise WorkflowError("maximum depth factor must be greater than 1")
    return replace(
        config,
        reference_fasta=reference,
        input_bam=input_bam,
        output_dir=output,
        sample_id=sample_id,
    )


def resolve_tools() -> ToolPaths:
    samtools = shutil.which("samtools")
    if samtools is None:
        raise WorkflowError("required executable not found on PATH: samtools")
    return ToolPaths(samtools)


def _run_capture(command: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    try:
        result = subprocess.run(
            command,
            cwd=cwd,
            text=True,
            encoding="utf-8",
            errors="replace",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
    except OSError as error:
        raise WorkflowError(f"could not launch {Path(command[0]).name}: {error}") from error
    if result.returncode != 0:
        detail = result.stderr.strip().splitlines()
        tail = detail[-1] if detail else "no diagnostic output"
        raise WorkflowError(
            f"{Path(command[0]).name} exited with code {result.returncode}: {tail}"
        )
    return result


def validate_tool_version(tools: ToolPaths) -> str:
    result = _run_capture([tools.samtools, "--version"], Path.cwd())
    match = re.search(r"^samtools\s+(\S+)", result.stdout, re.MULTILINE)
    version = match.group(1) if match else "unknown"
    if version != EXPECTED_SAMTOOLS_VERSION:
        raise WorkflowError(
            f"Samtools {EXPECTED_SAMTOOLS_VERSION} is required; found {version}"
        )
    return version


def _sanitize(text: str, replacements: dict[str, str]) -> str:
    sanitized = text
    for source, replacement in sorted(replacements.items(), key=lambda item: -len(item[0])):
        sanitized = sanitized.replace(source, replacement)
    return sanitized


def _portable(command: list[str], replacements: dict[str, str]) -> list[str]:
    portable = [_sanitize(value, replacements) for value in command]
    portable[0] = Path(command[0]).name
    return portable


def _run_logged(
    command: list[str],
    *,
    stage: str,
    cwd: Path,
    staging: Path,
    replacements: dict[str, str],
    commands: dict[str, list[str]],
) -> subprocess.CompletedProcess[str]:
    result = _run_capture(command, cwd)
    (staging / "logs" / f"{stage}.stdout.txt").write_text(
        _sanitize(result.stdout, replacements), encoding="utf-8", newline="\n"
    )
    (staging / "logs" / f"{stage}.stderr.txt").write_text(
        _sanitize(result.stderr, replacements), encoding="utf-8", newline="\n"
    )
    commands[stage] = _portable(command, replacements)
    return result


def _read_header(
    tools: ToolPaths,
    bam: Path,
    *,
    stage: str,
    staging: Path,
    replacements: dict[str, str],
    commands: dict[str, list[str]],
) -> str:
    """Read a BAM header while omitting provenance lines with native paths."""
    command = [tools.samtools, "view", "--no-PG", "-H", str(bam)]
    result = _run_capture(command, staging)
    header = "\n".join(
        line for line in result.stdout.splitlines() if not line.startswith("@PG\t")
    ) + "\n"
    (staging / "logs" / f"{stage}.stdout.txt").write_text(
        _sanitize(header, replacements), encoding="utf-8", newline="\n"
    )
    (staging / "logs" / f"{stage}.stderr.txt").write_text(
        _sanitize(result.stderr, replacements), encoding="utf-8", newline="\n"
    )
    commands[stage] = _portable(command, replacements)
    return header


def _count_records(
    tools: ToolPaths,
    bam: Path,
    *,
    flag: int | None,
    stage: str,
    staging: Path,
    replacements: dict[str, str],
    commands: dict[str, list[str]],
) -> int:
    command = [tools.samtools, "view", "-c"]
    if flag is not None:
        command.extend(("-f", str(flag)))
    command.append(str(bam))
    result = _run_logged(
        command,
        stage=stage,
        cwd=staging,
        staging=staging,
        replacements=replacements,
        commands=commands,
    )
    try:
        return int(result.stdout.strip())
    except ValueError as error:
        raise WorkflowError(f"{stage} did not return an integer record count") from error


def _write_sanitized_bam(
    tools: ToolPaths,
    source: Path,
    destination: Path,
    *,
    staging: Path,
    replacements: dict[str, str],
    commands: dict[str, list[str]],
) -> None:
    header = _read_header(
        tools,
        source,
        stage="marked-header",
        staging=staging,
        replacements=replacements,
        commands=commands,
    )
    header_path = staging / "intermediate" / "sanitized_header.sam"
    header_path.write_text(header, encoding="ascii", newline="\n")
    command = [tools.samtools, "reheader", "-P", str(header_path), str(source)]
    try:
        with destination.open("wb") as output:
            result = subprocess.run(
                command,
                cwd=staging,
                stdout=output,
                stderr=subprocess.PIPE,
                check=False,
            )
    except OSError as error:
        raise WorkflowError(f"could not launch samtools reheader: {error}") from error
    stderr = result.stderr.decode("utf-8", errors="replace")
    (staging / "logs" / "reheader.stderr.txt").write_text(
        _sanitize(stderr, replacements), encoding="utf-8", newline="\n"
    )
    commands["reheader"] = _portable(command, replacements)
    if result.returncode != 0:
        tail = stderr.strip().splitlines()
        raise WorkflowError(
            f"samtools reheader exited with code {result.returncode}: "
            f"{tail[-1] if tail else 'no diagnostic output'}"
        )


def _collect_depth(
    tools: ToolPaths,
    bam: Path,
    config: WorkflowConfig,
    contigs: tuple[tuple[str, int], ...],
    *,
    staging: Path,
    replacements: dict[str, str],
    commands: dict[str, list[str]],
):
    command = [
        tools.samtools,
        "depth",
        "-aa",
        "-q",
        str(config.minimum_base_quality),
        "-Q",
        str(config.minimum_mapq),
        str(bam),
    ]
    commands["depth"] = _portable(command, replacements)
    try:
        process = subprocess.Popen(
            command,
            cwd=staging,
            text=True,
            encoding="utf-8",
            errors="replace",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except OSError as error:
        raise WorkflowError(f"could not launch samtools depth: {error}") from error
    assert process.stdout is not None
    assert process.stderr is not None
    depth_path = staging / "depth" / "position_depth.tsv.gz"
    try:
        with gzip.open(depth_path, "wt", encoding="utf-8", newline="\n") as raw_output:
            dataset = collect_depth(
                process.stdout,
                contigs,
                minimum_depth=config.minimum_depth,
                maximum_depth_factor=config.maximum_depth_factor,
                raw_output=raw_output,
            )
        stderr = process.stderr.read()
        return_code = process.wait()
        process.stdout.close()
        process.stderr.close()
    except Exception:
        process.kill()
        process.wait()
        process.stdout.close()
        process.stderr.close()
        raise
    (staging / "logs" / "depth.stderr.txt").write_text(
        _sanitize(stderr, replacements), encoding="utf-8", newline="\n"
    )
    if return_code != 0:
        tail = stderr.strip().splitlines()
        raise WorkflowError(
            f"samtools depth exited with code {return_code}: "
            f"{tail[-1] if tail else 'no diagnostic output'}"
        )
    return dataset


def run_variant_ready_workflow(
    config: WorkflowConfig,
    *,
    tools: ToolPaths | None = None,
) -> Path:
    config = validate_config(config)
    tools = tools or resolve_tools()
    samtools_version = validate_tool_version(tools)
    try:
        reference_contigs = read_fasta_lengths(config.reference_fasta)
    except ValidationError as error:
        raise WorkflowError(str(error)) from error
    input_bam_hash = sha256_file(config.input_bam)
    reference_hash = sha256_file(config.reference_fasta)
    input_index = Path(str(config.input_bam) + ".bai")
    if not input_index.is_file():
        raise WorkflowError(f"input BAM index is missing: {input_index}")

    parent = config.output_dir.parent
    parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{config.output_dir.name}.staging-", dir=parent))
    for name in ("bam", "depth", "intermediate", "logs", "reports"):
        (staging / name).mkdir()
    replacements = {
        str(staging): ".",
        str(config.input_bam): "input.bam",
        str(config.reference_fasta): "reference.fasta",
        str(input_index): "input.bam.bai",
        tools.samtools: "samtools",
    }
    commands: dict[str, list[str]] = {}

    try:
        _run_logged(
            [tools.samtools, "quickcheck", "-v", str(config.input_bam)],
            stage="input-quickcheck",
            cwd=staging,
            staging=staging,
            replacements=replacements,
            commands=commands,
        )
        header_text = _read_header(
            tools,
            config.input_bam,
            stage="input-header",
            staging=staging,
            replacements=replacements,
            commands=commands,
        )
        header = parse_bam_header(header_text)
        if header.sort_order != "coordinate":
            raise WorkflowError("input BAM must declare coordinate sort order")
        if header.contigs != reference_contigs:
            raise WorkflowError("BAM @SQ records do not exactly match the reference FASTA")
        samples = {read_group["SM"] for read_group in header.read_groups}
        if samples != {config.sample_id}:
            raise WorkflowError(
                f"BAM read-group sample values {sorted(samples)!r} do not match "
                f"sample ID {config.sample_id!r}"
            )

        input_records = _count_records(
            tools,
            config.input_bam,
            flag=None,
            stage="input-record-count",
            staging=staging,
            replacements=replacements,
            commands=commands,
        )
        input_duplicate_records = _count_records(
            tools,
            config.input_bam,
            flag=1024,
            stage="input-duplicate-count",
            staging=staging,
            replacements=replacements,
            commands=commands,
        )

        name_collated = staging / "intermediate" / "name_collated.bam"
        fixmate = staging / "intermediate" / "fixmate.bam"
        position_sorted = staging / "intermediate" / "position_sorted.bam"
        marked = staging / "intermediate" / "marked.bam"
        variant_ready = staging / "bam" / "variant_ready.sorted.bam"
        # Markdup needs mate-score tags from a name-collated fixmate pass, then
        # coordinate order. It flags duplicate records but does not remove them.
        stages = (
            (
                "collate",
                [tools.samtools, "collate", "--no-PG", "-@", str(config.threads), "-o", str(name_collated), str(config.input_bam)],
            ),
            (
                "fixmate",
                [tools.samtools, "fixmate", "--no-PG", "-@", str(config.threads), "-m", str(name_collated), str(fixmate)],
            ),
            (
                "coordinate-sort",
                [tools.samtools, "sort", "--no-PG", "-@", str(config.threads), "-o", str(position_sorted), str(fixmate)],
            ),
            (
                "markdup",
                [
                    tools.samtools,
                    "markdup",
                    "--no-PG",
                    "-@",
                    str(config.threads),
                    "-c",
                    "-s",
                    "--json",
                    "-f",
                    str(staging / "reports" / "markdup.json"),
                    str(position_sorted),
                    str(marked),
                ],
            ),
        )
        for stage, command in stages:
            _run_logged(
                command,
                stage=stage,
                cwd=staging,
                staging=staging,
                replacements=replacements,
                commands=commands,
            )
        markdup_report = staging / "reports" / "markdup.json"
        markdup_report.write_text(
            _sanitize(markdup_report.read_text(encoding="utf-8"), replacements),
            encoding="utf-8",
            newline="\n",
        )
        # Preserve scientific header fields while removing native command paths.
        _write_sanitized_bam(
            tools,
            marked,
            variant_ready,
            staging=staging,
            replacements=replacements,
            commands=commands,
        )

        _run_logged(
            [tools.samtools, "index", "-@", str(config.threads), str(variant_ready)],
            stage="index",
            cwd=staging,
            staging=staging,
            replacements=replacements,
            commands=commands,
        )
        _run_logged(
            [tools.samtools, "quickcheck", "-v", str(variant_ready)],
            stage="output-quickcheck",
            cwd=staging,
            staging=staging,
            replacements=replacements,
            commands=commands,
        )
        output_records = _count_records(
            tools,
            variant_ready,
            flag=None,
            stage="output-record-count",
            staging=staging,
            replacements=replacements,
            commands=commands,
        )
        duplicate_records = _count_records(
            tools,
            variant_ready,
            flag=1024,
            stage="output-duplicate-count",
            staging=staging,
            replacements=replacements,
            commands=commands,
        )
        if output_records != input_records:
            raise WorkflowError("duplicate marking changed the BAM record count")

        for stage, command, report_name in (
            ("flagstat", [tools.samtools, "flagstat", "-@", str(config.threads), "-O", "json", str(variant_ready)], "flagstat.json"),
            ("stats", [tools.samtools, "stats", "-@", str(config.threads), str(variant_ready)], "stats.txt"),
            ("coverage", [tools.samtools, "coverage", str(variant_ready)], "coverage.tsv"),
            ("idxstats", [tools.samtools, "idxstats", str(variant_ready)], "idxstats.tsv"),
        ):
            result = _run_logged(
                command,
                stage=stage,
                cwd=staging,
                staging=staging,
                replacements=replacements,
                commands=commands,
            )
            (staging / "reports" / report_name).write_text(
                _sanitize(result.stdout, replacements), encoding="utf-8", newline="\n"
            )

        depth_dataset = _collect_depth(
            tools,
            variant_ready,
            config,
            reference_contigs,
            staging=staging,
            replacements=replacements,
            commands=commands,
        )
        write_depth_outputs(
            depth_dataset,
            staging / "depth",
            minimum_depth=config.minimum_depth,
        )
        # Every depth percentage uses the complete FASTA length as denominator.
        gate_passed = (
            depth_dataset.summary.mean_depth >= 15
            and depth_dataset.summary.percent_ge_minimum >= 90
        )

        if sha256_file(config.input_bam) != input_bam_hash:
            raise WorkflowError("input BAM checksum changed during processing")
        if sha256_file(config.reference_fasta) != reference_hash:
            raise WorkflowError("reference FASTA checksum changed during processing")

        manifest: dict[str, Any] = {
            "workflow_version": __version__,
            "milestone": "variant_ready_bam",
            "tools": {"python": platform.python_version(), "samtools": samtools_version},
            "parameters": {
                key: value
                for key, value in asdict(config).items()
                if key not in {"reference_fasta", "input_bam", "output_dir"}
            },
            "inputs": {
                "reference": {"filename": config.reference_fasta.name, "sha256": reference_hash},
                "bam": {"filename": config.input_bam.name, "sha256": input_bam_hash},
            },
            "counts": {
                "input_records": input_records,
                "input_duplicate_records": input_duplicate_records,
                "output_records": output_records,
                "duplicate_records": duplicate_records,
                "nonduplicate_records": output_records - duplicate_records,
            },
            "depth": asdict(depth_dataset.summary),
            "dataset_gate": {
                "minimum_mean_depth": 15.0,
                "minimum_percent_positions_at_minimum_depth": 90.0,
                "passed": gate_passed,
            },
            "outputs": {
                "variant_ready_bam": "bam/variant_ready.sorted.bam",
                "variant_ready_bai": "bam/variant_ready.sorted.bam.bai",
                "bam_sha256": sha256_file(variant_ready),
                "bai_sha256": sha256_file(Path(str(variant_ready) + ".bai")),
            },
            "commands": commands,
            "duplicate_policy": "Records are marked with the DUP flag and preserved; they are excluded by default from depth evidence.",
        }
        (staging / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8", newline="\n"
        )

        shutil.rmtree(staging / "intermediate")
        # A completed staging tree appears at the final path in one rename.
        os.replace(staging, config.output_dir)
        return config.output_dir
    except WorkflowError:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    except ValidationError as error:
        shutil.rmtree(staging, ignore_errors=True)
        raise WorkflowError(str(error)) from error
    except Exception as error:
        shutil.rmtree(staging, ignore_errors=True)
        raise WorkflowError(f"variant-ready workflow failed: {error}") from error
