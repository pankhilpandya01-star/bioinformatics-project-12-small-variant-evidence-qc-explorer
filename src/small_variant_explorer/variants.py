"""Normalized VCF evidence parsing and deterministic filter assignment."""

from __future__ import annotations

import csv
import gzip
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from small_variant_explorer.validation import ValidationError


FILTER_DESCRIPTIONS = {
    "LowQual": "Variant quality is below the configured minimum",
    "LowDepth": "Reliable read depth is below the configured minimum",
    "HighDepth": "Read depth exceeds the genome-wide high-depth boundary",
    "LowAltDepth": "Alternate-allele read depth is below the configured minimum",
    "LowAltFraction": "Alternate-allele fraction is below the configured minimum",
    "OneStrandOnly": "Alternate support is absent or insufficient on one strand",
}


@dataclass(frozen=True, slots=True)
class EvidenceThresholds:
    minimum_depth: int
    minimum_alt_depth: int
    minimum_alt_fraction: float
    minimum_variant_quality: float
    minimum_alt_strand_depth: int
    high_depth_threshold: int


@dataclass(frozen=True, slots=True)
class CandidateEvidence:
    chrom: str
    position: int
    ref: str
    alt: str
    variant_type: str
    quality: float
    depth: int
    ref_depth: int
    alt_depth: int
    alt_fraction: float
    alt_forward_depth: int
    alt_reverse_depth: int
    genotype: str
    status: str
    filter_reasons: tuple[str, ...]

    @property
    def key(self) -> str:
        return f"{self.chrom}:{self.position}:{self.ref}:{self.alt}"


@dataclass(frozen=True, slots=True)
class EvaluationSummary:
    candidates: int
    accepted: int
    rejected: int
    evidence: tuple[CandidateEvidence, ...]


def evaluate_normalized_vcf(
    normalized_vcf: Path,
    intermediate_dir: Path,
    report_dir: Path,
    thresholds: EvidenceThresholds,
) -> EvaluationSummary:
    """Evaluate each biallelic normalized record and write plain VCF partitions."""
    header, records = _read_vcf(normalized_vcf)
    output_header = _evaluated_header(header)
    evaluated_lines: list[str] = []
    accepted_lines: list[str] = []
    rejected_lines: list[str] = []
    evidence: list[CandidateEvidence] = []
    seen_keys: set[str] = set()

    for fields in records:
        candidate = _evaluate_record(fields, thresholds)
        if candidate.key in seen_keys:
            raise ValidationError(f"duplicate normalized variant key: {candidate.key}")
        seen_keys.add(candidate.key)
        fields[6] = "PASS" if candidate.status == "accepted" else ";".join(
            candidate.filter_reasons
        )
        line = "\t".join(fields)
        evaluated_lines.append(line)
        if candidate.status == "accepted":
            accepted_lines.append(line)
        else:
            rejected_lines.append(line)
        evidence.append(candidate)

    for name, lines in (
        ("evaluated", evaluated_lines),
        ("accepted", accepted_lines),
        ("rejected", rejected_lines),
    ):
        (intermediate_dir / f"{name}.vcf").write_text(
            "\n".join((*output_header, *lines)) + "\n",
            encoding="utf-8",
            newline="\n",
        )

    _write_evidence_csv(report_dir / "variant_evidence.csv", evidence)
    _write_summaries(report_dir, evidence)
    accepted = sum(candidate.status == "accepted" for candidate in evidence)
    return EvaluationSummary(
        candidates=len(evidence),
        accepted=accepted,
        rejected=len(evidence) - accepted,
        evidence=tuple(evidence),
    )


def _read_vcf(path: Path) -> tuple[list[str], list[list[str]]]:
    opener = gzip.open if path.suffix == ".gz" else open
    header: list[str] = []
    records: list[list[str]] = []
    column_header: list[str] | None = None
    try:
        with opener(path, "rt", encoding="utf-8", newline=None) as handle:
            for line_number, line in enumerate(handle, start=1):
                value = line.rstrip("\r\n")
                if value.startswith("##"):
                    header.append(value)
                elif value.startswith("#CHROM\t"):
                    column_header = value.split("\t")
                    header.append(value)
                elif value:
                    fields = value.split("\t")
                    if len(fields) != 10:
                        raise ValidationError(
                            f"VCF line {line_number} must contain exactly one sample"
                        )
                    records.append(fields)
    except (OSError, UnicodeError) as error:
        raise ValidationError(f"could not read normalized VCF: {error}") from error
    if column_header is None:
        raise ValidationError("VCF has no #CHROM header")
    if len(column_header) != 10:
        raise ValidationError("VCF must contain exactly one sample")
    return header, records


def _evaluated_header(header: list[str]) -> list[str]:
    clean = [
        line
        for line in header
        if not (line.startswith("##bcftools_") and "Command=" in line)
    ]
    column_index = next(index for index, line in enumerate(clean) if line.startswith("#CHROM"))
    existing = {
        line.split("ID=", 1)[1].split(",", 1)[0]
        for line in clean
        if line.startswith("##FILTER=<ID=")
    }
    definitions = [
        f'##FILTER=<ID={name},Description="{description}">'
        for name, description in FILTER_DESCRIPTIONS.items()
        if name not in existing
    ]
    return clean[:column_index] + definitions + clean[column_index:]


def _evaluate_record(
    fields: list[str], thresholds: EvidenceThresholds
) -> CandidateEvidence:
    chrom, position_text, _, ref, alt, quality_text, _, _, format_text, sample_text = fields
    if not chrom or not ref or not alt or "," in alt or alt.startswith("<"):
        raise ValidationError("normalized VCF record is not a simple biallelic variant")
    try:
        position = int(position_text)
        quality = 0.0 if quality_text == "." else float(quality_text)
    except ValueError as error:
        raise ValidationError("VCF position or quality is not numeric") from error
    keys = format_text.split(":")
    values = sample_text.split(":")
    sample = dict(zip(keys, values, strict=False))
    for required in ("GT", "DP", "AD", "ADF", "ADR"):
        if required not in sample:
            raise ValidationError(f"VCF sample is missing FORMAT/{required}")
    depth = _parse_single_int(sample["DP"], "DP")
    allelic_depth = _parse_int_list(sample["AD"], "AD")
    forward_depth = _parse_int_list(sample["ADF"], "ADF")
    reverse_depth = _parse_int_list(sample["ADR"], "ADR")
    if min(len(allelic_depth), len(forward_depth), len(reverse_depth)) < 2:
        raise ValidationError("VCF allelic-depth fields must contain REF and ALT values")
    ref_depth = allelic_depth[0]
    alt_depth = allelic_depth[1]
    alt_forward = forward_depth[1]
    alt_reverse = reverse_depth[1]
    alt_fraction = alt_depth / depth if depth else 0.0

    # Reasons are independent: weak candidates retain every applicable warning.
    reasons: list[str] = []
    if quality < thresholds.minimum_variant_quality:
        reasons.append("LowQual")
    if depth < thresholds.minimum_depth:
        reasons.append("LowDepth")
    if depth > thresholds.high_depth_threshold:
        reasons.append("HighDepth")
    if alt_depth < thresholds.minimum_alt_depth:
        reasons.append("LowAltDepth")
    if alt_fraction < thresholds.minimum_alt_fraction:
        reasons.append("LowAltFraction")
    if (
        alt_forward < thresholds.minimum_alt_strand_depth
        or alt_reverse < thresholds.minimum_alt_strand_depth
    ):
        reasons.append("OneStrandOnly")

    return CandidateEvidence(
        chrom=chrom,
        position=position,
        ref=ref,
        alt=alt,
        variant_type=_variant_type(ref, alt),
        quality=quality,
        depth=depth,
        ref_depth=ref_depth,
        alt_depth=alt_depth,
        alt_fraction=alt_fraction,
        alt_forward_depth=alt_forward,
        alt_reverse_depth=alt_reverse,
        genotype=sample["GT"],
        status="rejected" if reasons else "accepted",
        filter_reasons=tuple(reasons),
    )


def _parse_single_int(value: str, name: str) -> int:
    values = _parse_int_list(value, name)
    if len(values) != 1:
        raise ValidationError(f"FORMAT/{name} must contain one integer")
    return values[0]


def _parse_int_list(value: str, name: str) -> list[int]:
    try:
        values = [0 if item == "." else int(item) for item in value.split(",")]
    except ValueError as error:
        raise ValidationError(f"FORMAT/{name} is not an integer list") from error
    if not values or any(item < 0 for item in values):
        raise ValidationError(f"FORMAT/{name} contains invalid depth")
    return values


def _variant_type(ref: str, alt: str) -> str:
    if len(ref) == len(alt) == 1:
        return "SNP"
    if len(alt) > len(ref):
        return "insertion"
    if len(ref) > len(alt):
        return "deletion"
    return "complex"


def _write_evidence_csv(path: Path, evidence: list[CandidateEvidence]) -> None:
    fields = (
        "chrom",
        "position",
        "ref",
        "alt",
        "variant_key",
        "variant_type",
        "quality",
        "depth",
        "ref_depth",
        "alt_depth",
        "alt_fraction",
        "alt_forward_depth",
        "alt_reverse_depth",
        "genotype",
        "status",
        "filter_reasons",
    )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for candidate in evidence:
            writer.writerow(
                {
                    "chrom": candidate.chrom,
                    "position": candidate.position,
                    "ref": candidate.ref,
                    "alt": candidate.alt,
                    "variant_key": candidate.key,
                    "variant_type": candidate.variant_type,
                    "quality": candidate.quality,
                    "depth": candidate.depth,
                    "ref_depth": candidate.ref_depth,
                    "alt_depth": candidate.alt_depth,
                    "alt_fraction": f"{candidate.alt_fraction:.6f}",
                    "alt_forward_depth": candidate.alt_forward_depth,
                    "alt_reverse_depth": candidate.alt_reverse_depth,
                    "genotype": candidate.genotype,
                    "status": candidate.status,
                    "filter_reasons": ";".join(candidate.filter_reasons) or "PASS",
                }
            )


def _write_summaries(report_dir: Path, evidence: list[CandidateEvidence]) -> None:
    types = Counter(candidate.variant_type for candidate in evidence)
    accepted_types = Counter(
        candidate.variant_type for candidate in evidence if candidate.status == "accepted"
    )
    with (report_dir / "variant_type_summary.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("variant_type", "total", "accepted", "rejected"))
        for variant_type in ("SNP", "insertion", "deletion", "complex"):
            writer.writerow(
                (
                    variant_type,
                    types[variant_type],
                    accepted_types[variant_type],
                    types[variant_type] - accepted_types[variant_type],
                )
            )

    reason_counts = Counter(
        reason for candidate in evidence for reason in candidate.filter_reasons
    )
    with (report_dir / "filter_reason_summary.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("filter_reason", "candidates"))
        writer.writerow(("PASS", sum(candidate.status == "accepted" for candidate in evidence)))
        for reason in FILTER_DESCRIPTIONS:
            writer.writerow((reason, reason_counts[reason]))

    substitutions = Counter(
        f"{candidate.ref}>{candidate.alt}"
        for candidate in evidence
        if candidate.variant_type == "SNP"
    )
    accepted_substitutions = Counter(
        f"{candidate.ref}>{candidate.alt}"
        for candidate in evidence
        if candidate.variant_type == "SNP" and candidate.status == "accepted"
    )
    with (report_dir / "substitution_summary.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("substitution", "accepted", "rejected", "total"))
        for substitution in sorted(substitutions):
            accepted = accepted_substitutions[substitution]
            writer.writerow(
                (
                    substitution,
                    accepted,
                    substitutions[substitution] - accepted,
                    substitutions[substitution],
                )
            )
