"""Reference and BAM-header validation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


FASTA_ALPHABET = frozenset("ACGTURYSWKMBDHVN")


class ValidationError(ValueError):
    """Raised when a workflow input is malformed or inconsistent."""


@dataclass(frozen=True, slots=True)
class BamHeader:
    sort_order: str
    contigs: tuple[tuple[str, int], ...]
    read_groups: tuple[dict[str, str], ...]


def read_fasta_lengths(path: Path) -> tuple[tuple[str, int], ...]:
    """Validate a FASTA and return ordered identifier/length pairs."""
    records: list[tuple[str, int]] = []
    identifiers: set[str] = set()
    current: str | None = None
    length = 0

    def finish() -> None:
        nonlocal current, length
        if current is None:
            return
        if length == 0:
            raise ValidationError(f"FASTA record {current!r} has no sequence")
        records.append((current, length))
        current = None
        length = 0

    try:
        with path.open("rt", encoding="ascii", newline=None) as handle:
            for line_number, line in enumerate(handle, start=1):
                value = line.rstrip("\r\n")
                if value.startswith(">"):
                    finish()
                    description = value[1:].strip()
                    identifier = description.split(maxsplit=1)[0] if description else ""
                    if not identifier:
                        raise ValidationError("FASTA record has an empty identifier")
                    if identifier in identifiers:
                        raise ValidationError(f"duplicate FASTA identifier: {identifier!r}")
                    identifiers.add(identifier)
                    current = identifier
                    continue
                if current is None:
                    raise ValidationError(
                        f"FASTA line {line_number} appears before the first header"
                    )
                if not value:
                    continue
                if any(character.isspace() for character in value):
                    raise ValidationError(
                        f"FASTA line {line_number} contains sequence whitespace"
                    )
                invalid = sorted(set(value.upper()) - FASTA_ALPHABET)
                if invalid:
                    raise ValidationError(
                        f"FASTA line {line_number} contains invalid symbols: {''.join(invalid)}"
                    )
                length += len(value)
    except (OSError, UnicodeError) as error:
        raise ValidationError(f"could not read reference FASTA: {error}") from error
    finish()
    if not records:
        raise ValidationError("reference FASTA contains no records")
    return tuple(records)


def parse_bam_header(text: str) -> BamHeader:
    """Parse only the SAM header fields required by the workflow contract."""
    sort_order = ""
    contigs: list[tuple[str, int]] = []
    read_groups: list[dict[str, str]] = []
    seen_contigs: set[str] = set()
    seen_read_groups: set[str] = set()

    for line in text.splitlines():
        fields = line.split("\t")
        if fields[0] == "@HD":
            values = _tag_values(fields[1:])
            sort_order = values.get("SO", "")
        elif fields[0] == "@SQ":
            values = _tag_values(fields[1:])
            name = values.get("SN", "")
            try:
                length = int(values.get("LN", ""))
            except ValueError as error:
                raise ValidationError("BAM @SQ header has an invalid LN value") from error
            if not name or length <= 0:
                raise ValidationError("BAM @SQ header requires positive LN and non-empty SN")
            if name in seen_contigs:
                raise ValidationError(f"duplicate BAM contig: {name!r}")
            seen_contigs.add(name)
            contigs.append((name, length))
        elif fields[0] == "@RG":
            values = _tag_values(fields[1:])
            read_group = values.get("ID", "")
            sample = values.get("SM", "")
            if not read_group or not sample:
                raise ValidationError("every BAM @RG line requires ID and SM values")
            if read_group in seen_read_groups:
                raise ValidationError(f"duplicate BAM read-group ID: {read_group!r}")
            seen_read_groups.add(read_group)
            read_groups.append(values)

    if not contigs:
        raise ValidationError("BAM header contains no @SQ records")
    if not read_groups:
        raise ValidationError("BAM header contains no @RG records")
    return BamHeader(sort_order, tuple(contigs), tuple(read_groups))


def _tag_values(fields: list[str]) -> dict[str, str]:
    values: dict[str, str] = {}
    for field in fields:
        key, separator, value = field.partition(":")
        if separator:
            values[key] = value
    return values
