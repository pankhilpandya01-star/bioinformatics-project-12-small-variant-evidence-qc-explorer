"""Build the deterministic 2 kb controlled variant example."""

from __future__ import annotations

import csv
import hashlib
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class ControlledVariant:
    name: str
    position: int
    ref: str
    alt: str
    expected_status: str


CONTROLLED_VARIANTS = (
    ControlledVariant("snp_one", 300, "T", "G", "accepted"),
    ControlledVariant("snp_two", 700, "T", "A", "accepted"),
    ControlledVariant("insertion", 1000, "T", "TGG", "accepted"),
    ControlledVariant("deletion", 1300, "TAC", "T", "accepted"),
    ControlledVariant("weak_error", 1600, "T", "C", "rejected"),
)


def build_controlled_variant_fixture(output_dir: Path, samtools: str | None = None) -> Path:
    """Create a reference, truth table, SAM, sorted BAM, and BAI for tests."""
    if output_dir.exists():
        raise ValueError(f"controlled fixture directory already exists: {output_dir}")
    executable = samtools or shutil.which("samtools")
    if executable is None:
        raise RuntimeError("samtools is required to build the controlled fixture")
    output_dir.mkdir(parents=True)
    reference = _reference_sequence()
    reference_path = output_dir / "reference.fasta"
    reference_path.write_text(
        ">control\n" + "\n".join(
            reference[index : index + 80] for index in range(0, len(reference), 80)
        ) + "\n",
        encoding="ascii",
        newline="\n",
    )
    _write_truth(output_dir / "expected_variants.csv")
    sam_path = output_dir / "controlled_variants.sam"
    _write_sam(sam_path, reference)

    unsorted = output_dir / "controlled.unsorted.bam"
    sorted_bam = output_dir / "controlled.sorted.bam"
    for command in (
        [executable, "faidx", str(reference_path)],
        [executable, "view", "--no-PG", "-b", "-o", str(unsorted), str(sam_path)],
        [executable, "sort", "--no-PG", "-o", str(sorted_bam), str(unsorted)],
        [executable, "index", str(sorted_bam)],
    ):
        result = subprocess.run(
            command,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        if result.returncode != 0:
            raise RuntimeError(result.stderr.strip() or "controlled fixture command failed")
    unsorted.unlink()
    return sorted_bam


def _reference_sequence() -> str:
    alphabet = "ACGT"
    bases: list[str] = []
    block = 0
    while len(bases) < 2000:
        digest = hashlib.sha256(f"project12-control-{block}".encode("ascii")).digest()
        bases.extend(alphabet[value & 3] for value in digest)
        block += 1
    bases = bases[:2000]
    for position, base in (
        (300, "T"),
        (700, "T"),
        (1000, "T"),
        (1300, "T"),
        (1301, "A"),
        (1302, "C"),
        (1600, "T"),
    ):
        bases[position - 1] = base
    return "".join(bases)


def _write_truth(path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("name", "chrom", "position", "ref", "alt", "expected_status"))
        for variant in CONTROLLED_VARIANTS:
            writer.writerow(
                (
                    variant.name,
                    "control",
                    variant.position,
                    variant.ref,
                    variant.alt,
                    variant.expected_status,
                )
            )


def _write_sam(path: Path, reference: str) -> None:
    lines = (
        "@HD\tVN:1.6\tSO:unsorted",
        "@SQ\tSN:control\tLN:2000",
        "@RG\tID:controlled_rg\tSM:controlled\tLB:controlled_library\tPL:ILLUMINA",
    )
    records: list[str] = []
    for variant in CONTROLLED_VARIANTS:
        for read_index in range(10):
            start = variant.position - 49 - read_index
            use_alt = variant.name != "weak_error" or read_index < 6
            if variant.name == "weak_error":
                # The six weak alternate reads are deliberately forward-only.
                reverse = read_index >= 6
            else:
                reverse = read_index >= 5
            sequence, cigar = _read_sequence(reference, variant, start, use_alt)
            records.append(
                "\t".join(
                    (
                        f"{variant.name}_{read_index + 1}",
                        "16" if reverse else "0",
                        "control",
                        str(start),
                        "60",
                        cigar,
                        "*",
                        "0",
                        "0",
                        sequence,
                        "I" * len(sequence),
                        "RG:Z:controlled_rg",
                    )
                )
            )
    path.write_text(
        "\n".join((*lines, *records)) + "\n",
        encoding="ascii",
        newline="\n",
    )


def _read_sequence(
    reference: str,
    variant: ControlledVariant,
    start: int,
    use_alt: bool,
) -> tuple[str, str]:
    start_index = start - 1
    if not use_alt:
        return reference[start_index : start_index + 100], "100M"
    if len(variant.ref) == len(variant.alt) == 1:
        sequence = list(reference[start_index : start_index + 100])
        sequence[variant.position - start] = variant.alt
        return "".join(sequence), "100M"

    left_match = variant.position - start + 1
    right_match = 100 - left_match
    if len(variant.alt) > len(variant.ref):
        inserted = variant.alt[1:]
        sequence = (
            reference[start_index : variant.position]
            + inserted
            + reference[variant.position : variant.position + right_match]
        )
        return sequence, f"{left_match}M{len(inserted)}I{right_match}M"

    deleted_length = len(variant.ref) - len(variant.alt)
    sequence = (
        reference[start_index : variant.position]
        + reference[
            variant.position + deleted_length : variant.position + deleted_length + right_match
        ]
    )
    return sequence, f"{left_match}M{deleted_length}D{right_match}M"
