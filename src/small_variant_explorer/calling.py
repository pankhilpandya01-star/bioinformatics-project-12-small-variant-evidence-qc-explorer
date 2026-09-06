"""Atomic haploid candidate calling, normalization, and evidence filtering."""

from __future__ import annotations

import json
import os
import platform
import re
import shutil
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from small_variant_explorer import __version__
from small_variant_explorer.dashboard import create_dashboard
from small_variant_explorer.depth import write_depth_outputs
from small_variant_explorer.validation import ValidationError, parse_bam_header, read_fasta_lengths
from small_variant_explorer.variants import EvidenceThresholds, evaluate_normalized_vcf
from small_variant_explorer.workflow import (
    PRIVATE_MARKERS,
    ToolPaths,
    WorkflowConfig,
    WorkflowError,
    _collect_depth,
    _read_header,
    _run_capture,
    _run_logged,
    _sanitize,
    sha256_file,
    validate_config,
    validate_tool_version,
)


EXPECTED_BCFTOOLS_VERSION = "1.24"


@dataclass(frozen=True, slots=True)
class VariantToolPaths:
    samtools: str
    bcftools: str


def resolve_variant_tools() -> VariantToolPaths:
    samtools = shutil.which("samtools")
    bcftools = shutil.which("bcftools")
    missing = [name for name, path in (("samtools", samtools), ("bcftools", bcftools)) if path is None]
    if missing:
        raise WorkflowError(f"required executable not found on PATH: {', '.join(missing)}")
    assert samtools is not None and bcftools is not None
    return VariantToolPaths(samtools, bcftools)


def validate_bcftools_version(tools: VariantToolPaths) -> str:
    result = _run_capture([tools.bcftools, "--version"], Path.cwd())
    match = re.search(r"^bcftools\s+(\S+)", result.stdout, re.MULTILINE)
    version = match.group(1) if match else "unknown"
    if version != EXPECTED_BCFTOOLS_VERSION:
        raise WorkflowError(
            f"BCFtools {EXPECTED_BCFTOOLS_VERSION} is required; found {version}"
        )
    return version


def run_variant_calling_workflow(
    config: WorkflowConfig,
    *,
    tools: VariantToolPaths | None = None,
) -> Path:
    """Call and evaluate normalized variants from a duplicate-marked BAM."""
    config = validate_config(config)
    tools = tools or resolve_variant_tools()
    samtools_version = validate_tool_version(ToolPaths(tools.samtools))
    bcftools_version = validate_bcftools_version(tools)
    try:
        reference_contigs = read_fasta_lengths(config.reference_fasta)
    except ValidationError as error:
        raise WorkflowError(str(error)) from error
    input_index = Path(str(config.input_bam) + ".bai")
    if not input_index.is_file():
        raise WorkflowError(f"input BAM index is missing: {input_index}")
    input_bam_hash = sha256_file(config.input_bam)
    reference_hash = sha256_file(config.reference_fasta)

    parent = config.output_dir.parent
    parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{config.output_dir.name}.staging-", dir=parent))
    for name in ("dashboard", "depth", "intermediate", "logs", "raw", "reports", "vcf"):
        (staging / name).mkdir()
    replacements = {
        str(staging): ".",
        str(config.input_bam): "input.bam",
        str(input_index): "input.bam.bai",
        str(config.reference_fasta): "reference.fasta",
        tools.samtools: "samtools",
        tools.bcftools: "bcftools",
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
            ToolPaths(tools.samtools),
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

        local_reference = staging / "intermediate" / "reference.fasta"
        shutil.copyfile(config.reference_fasta, local_reference)
        _run_logged(
            [tools.samtools, "faidx", "intermediate/reference.fasta"],
            stage="reference-index",
            cwd=staging,
            staging=staging,
            replacements=replacements,
            commands=commands,
        )

        depth_dataset = _collect_depth(
            ToolPaths(tools.samtools),
            config.input_bam,
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

        pileup = staging / "intermediate" / "pileup.bcf"
        called_bcf = staging / "intermediate" / "called.bcf"
        raw_bcf = staging / "raw" / "candidates.bcf"
        normalized = staging / "vcf" / "candidates.normalized.vcf.gz"
        _run_logged(
            [
                tools.bcftools,
                "mpileup",
                "--no-version",
                "--threads",
                str(config.threads),
                "--fasta-ref",
                "intermediate/reference.fasta",
                "--min-MQ",
                str(config.minimum_mapq),
                "--min-BQ",
                str(config.minimum_base_quality),
                "--max-depth",
                "1000000",
                "--annotate",
                "FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP",
                "--output-type",
                "u",
                "--output",
                str(pileup),
                str(config.input_bam),
            ],
            stage="mpileup",
            cwd=staging,
            staging=staging,
            replacements=replacements,
            commands=commands,
        )
        _run_logged(
            [
                tools.bcftools,
                "call",
                "--no-version",
                "--threads",
                str(config.threads),
                "--multiallelic-caller",
                "--variants-only",
                "--keep-alts",
                "--ploidy",
                str(config.ploidy),
                "--output-type",
                "b",
                "--output",
                str(called_bcf),
                str(pileup),
            ],
            stage="call",
            cwd=staging,
            staging=staging,
            replacements=replacements,
            commands=commands,
        )
        # BCFtools 1.24 writes INFO/MQ with a legacy integer declaration that
        # triggers its own VCF sanity warning. It is unused here, so omit it.
        _run_logged(
            [
                tools.bcftools,
                "annotate",
                "--no-version",
                "--remove",
                "INFO/MQ",
                "--output-type",
                "b",
                "--output",
                str(raw_bcf),
                str(called_bcf),
            ],
            stage="remove-unused-mq",
            cwd=staging,
            staging=staging,
            replacements=replacements,
            commands=commands,
        )
        _index_variant_file(
            tools,
            raw_bcf,
            "raw-index",
            staging,
            replacements,
            commands,
        )
        _run_logged(
            [
                tools.bcftools,
                "norm",
                "--no-version",
                "--threads",
                str(config.threads),
                "--fasta-ref",
                "intermediate/reference.fasta",
                "--multiallelics",
                "-any",
                "--rm-dup",
                "exact",
                "--output-type",
                "z",
                "--output",
                str(normalized),
                str(raw_bcf),
            ],
            stage="normalize",
            cwd=staging,
            staging=staging,
            replacements=replacements,
            commands=commands,
        )
        _index_variant_file(
            tools,
            normalized,
            "normalized-index",
            staging,
            replacements,
            commands,
        )

        # Haploid genotypes come from BCFtools; these filters grade read evidence
        # without changing the normalized allele representation.
        thresholds = EvidenceThresholds(
            minimum_depth=config.minimum_depth,
            minimum_alt_depth=config.minimum_alt_depth,
            minimum_alt_fraction=config.minimum_alt_fraction,
            minimum_variant_quality=config.minimum_variant_quality,
            minimum_alt_strand_depth=config.minimum_alt_strand_depth,
            high_depth_threshold=depth_dataset.summary.high_depth_threshold,
        )
        evaluation = evaluate_normalized_vcf(
            normalized,
            staging / "intermediate",
            staging / "reports",
            thresholds,
        )
        dashboard = staging / "dashboard" / "variant_evidence_dashboard.png"
        create_dashboard(depth_dataset, evaluation.evidence, thresholds, dashboard)
        for name in ("evaluated", "accepted", "rejected"):
            source = staging / "intermediate" / f"{name}.vcf"
            destination = staging / "vcf" / f"{name}.vcf.gz"
            _run_logged(
                [
                    tools.bcftools,
                    "view",
                    "--no-version",
                    "--output-type",
                    "z",
                    "--output",
                    str(destination),
                    str(source),
                ],
                stage=f"compress-{name}",
                cwd=staging,
                staging=staging,
                replacements=replacements,
                commands=commands,
            )
            _index_variant_file(
                tools,
                destination,
                f"index-{name}",
                staging,
                replacements,
                commands,
            )
            stats = _run_logged(
                [tools.bcftools, "stats", str(destination)],
                stage=f"stats-{name}",
                cwd=staging,
                staging=staging,
                replacements=replacements,
                commands=commands,
            )
            (staging / "reports" / f"{name}.bcftools.stats.txt").write_text(
                _sanitize(stats.stdout, replacements).replace("# \n", "#\n"),
                encoding="utf-8",
                newline="\n",
            )

        _audit_variant_headers(
            tools,
            (
                raw_bcf,
                normalized,
                staging / "vcf" / "evaluated.vcf.gz",
                staging / "vcf" / "accepted.vcf.gz",
                staging / "vcf" / "rejected.vcf.gz",
            ),
            staging,
        )
        if sha256_file(config.input_bam) != input_bam_hash:
            raise WorkflowError("input BAM checksum changed during variant calling")
        if sha256_file(config.reference_fasta) != reference_hash:
            raise WorkflowError("reference FASTA checksum changed during variant calling")

        output_files = (
            raw_bcf,
            Path(str(raw_bcf) + ".csi"),
            normalized,
            Path(str(normalized) + ".csi"),
            staging / "vcf" / "evaluated.vcf.gz",
            staging / "vcf" / "evaluated.vcf.gz.csi",
            staging / "vcf" / "accepted.vcf.gz",
            staging / "vcf" / "accepted.vcf.gz.csi",
            staging / "vcf" / "rejected.vcf.gz",
            staging / "vcf" / "rejected.vcf.gz.csi",
            dashboard,
        )
        manifest: dict[str, Any] = {
            "workflow_version": __version__,
            "milestone": "variant_calling_and_evidence_filtering",
            "tools": {
                "python": platform.python_version(),
                "samtools": samtools_version,
                "bcftools": bcftools_version,
            },
            "parameters": {
                key: value
                for key, value in asdict(config).items()
                if key not in {"reference_fasta", "input_bam", "output_dir"}
            },
            "inputs": {
                "reference": {"filename": config.reference_fasta.name, "sha256": reference_hash},
                "bam": {"filename": config.input_bam.name, "sha256": input_bam_hash},
            },
            "depth": asdict(depth_dataset.summary),
            "counts": {
                "normalized_candidates": evaluation.candidates,
                "accepted_candidates": evaluation.accepted,
                "rejected_candidates": evaluation.rejected,
            },
            "outputs": {
                str(path.relative_to(staging)).replace("\\", "/"): sha256_file(path)
                for path in output_files
            },
            "commands": commands,
            "partition_rule": "Every normalized candidate appears exactly once in accepted or rejected output.",
        }
        (staging / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8", newline="\n"
        )

        shutil.rmtree(staging / "intermediate")
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
        raise WorkflowError(f"variant-calling workflow failed: {error}") from error


def _index_variant_file(
    tools: VariantToolPaths,
    path: Path,
    stage: str,
    staging: Path,
    replacements: dict[str, str],
    commands: dict[str, list[str]],
) -> None:
    _run_logged(
        [tools.bcftools, "index", "--force", "--csi", str(path)],
        stage=stage,
        cwd=staging,
        staging=staging,
        replacements=replacements,
        commands=commands,
    )


def _audit_variant_headers(
    tools: VariantToolPaths,
    paths: tuple[Path, ...],
    staging: Path,
) -> None:
    for path in paths:
        result = _run_capture(
            [tools.bcftools, "view", "--no-version", "--header-only", str(path)],
            staging,
        )
        lowered = result.stdout.lower()
        for marker in PRIVATE_MARKERS:
            if marker in lowered:
                raise WorkflowError(f"native path marker found in {path.name} header")
