"""Quality-filtered positional depth accounting."""

from __future__ import annotations

import csv
import math
from array import array
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, TextIO

from small_variant_explorer.validation import ValidationError


@dataclass(frozen=True, slots=True)
class DepthSummary:
    positions: int
    mean_depth: float
    median_depth: float
    covered_positions: int
    covered_percent: float
    positions_ge_minimum: int
    percent_ge_minimum: float
    high_depth_threshold: int
    callable_positions: int
    callable_percent: float
    maximum_depth: int


@dataclass(slots=True)
class DepthDataset:
    depths: dict[str, array]
    summary: DepthSummary


def collect_depth(
    lines: Iterable[str],
    expected_contigs: tuple[tuple[str, int], ...],
    *,
    minimum_depth: int,
    maximum_depth_factor: float,
    raw_output: TextIO | None = None,
) -> DepthDataset:
    """Validate a complete ``samtools depth -aa`` stream and summarize it."""
    depths = {name: array("I") for name, _ in expected_contigs}
    lengths = dict(expected_contigs)
    contig_order = {name: index for index, (name, _) in enumerate(expected_contigs)}
    current_contig_index = 0

    for line_number, line in enumerate(lines, start=1):
        if raw_output is not None:
            raw_output.write(line)
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) != 3:
            raise ValidationError(f"depth line {line_number} does not contain three fields")
        contig = fields[0]
        if contig not in depths:
            raise ValidationError(f"depth stream contains unknown contig {contig!r}")
        index = contig_order[contig]
        if index < current_contig_index:
            raise ValidationError("depth contigs are not in reference order")
        current_contig_index = index
        try:
            position = int(fields[1])
            depth = int(fields[2])
        except ValueError as error:
            raise ValidationError(f"depth line {line_number} is not numeric") from error
        expected_position = len(depths[contig]) + 1
        if position != expected_position:
            raise ValidationError(
                f"depth position for {contig!r} is {position}; expected {expected_position}"
            )
        if depth < 0 or position > lengths[contig]:
            raise ValidationError(f"depth line {line_number} is outside the reference")
        depths[contig].append(depth)

    for contig, expected_length in expected_contigs:
        if len(depths[contig]) != expected_length:
            raise ValidationError(
                f"depth stream has {len(depths[contig])} positions for {contig!r}; "
                f"expected {expected_length}"
            )

    histogram = Counter(value for values in depths.values() for value in values)
    positions = sum(histogram.values())
    depth_sum = sum(depth * count for depth, count in histogram.items())
    covered = positions - histogram.get(0, 0)
    positions_ge_minimum = sum(
        count for depth, count in histogram.items() if depth >= minimum_depth
    )
    median = _histogram_median(histogram, positions)
    high_depth_threshold = max(50, math.ceil(maximum_depth_factor * median))
    callable_positions = sum(
        count
        for depth, count in histogram.items()
        if minimum_depth <= depth <= high_depth_threshold
    )
    summary = DepthSummary(
        positions=positions,
        mean_depth=depth_sum / positions,
        median_depth=median,
        covered_positions=covered,
        covered_percent=100 * covered / positions,
        positions_ge_minimum=positions_ge_minimum,
        percent_ge_minimum=100 * positions_ge_minimum / positions,
        high_depth_threshold=high_depth_threshold,
        callable_positions=callable_positions,
        callable_percent=100 * callable_positions / positions,
        maximum_depth=max(histogram),
    )
    return DepthDataset(depths, summary)


def _histogram_median(histogram: Counter[int], count: int) -> float:
    lower_rank = (count - 1) // 2
    upper_rank = count // 2
    cumulative = 0
    lower: int | None = None
    for value in sorted(histogram):
        cumulative += histogram[value]
        if lower is None and cumulative > lower_rank:
            lower = value
        if cumulative > upper_rank:
            assert lower is not None
            return (lower + value) / 2
    raise AssertionError("non-empty depth histogram did not contain a median")


def write_depth_outputs(
    dataset: DepthDataset,
    output_dir: Path,
    *,
    minimum_depth: int,
    window_size: int = 10_000,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = dataset.summary
    with (output_dir / "coverage_summary.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("metric", "value", "unit"))
        units = {
            "positions": "positions",
            "mean_depth": "depth",
            "median_depth": "depth",
            "covered_positions": "positions",
            "covered_percent": "percent",
            "positions_ge_minimum": "positions",
            "percent_ge_minimum": "percent",
            "high_depth_threshold": "depth",
            "callable_positions": "positions",
            "callable_percent": "percent",
            "maximum_depth": "depth",
        }
        for metric, value in asdict(summary).items():
            writer.writerow((metric, value, units[metric]))

    with (output_dir / "genome_windows.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            ("contig", "start", "end", "positions", "mean_depth", "callable_percent")
        )
        for contig, values in dataset.depths.items():
            for start_index in range(0, len(values), window_size):
                window = values[start_index : start_index + window_size]
                callable_count = sum(
                    minimum_depth <= value <= summary.high_depth_threshold
                    for value in window
                )
                writer.writerow(
                    (
                        contig,
                        start_index + 1,
                        start_index + len(window),
                        len(window),
                        f"{sum(window) / len(window):.6f}",
                        f"{100 * callable_count / len(window):.6f}",
                    )
                )

    # BED uses zero-based, half-open intervals; adjacent callable bases are merged.
    with (output_dir / "callable_regions.bed").open(
        "w", encoding="utf-8", newline="\n"
    ) as handle:
        for contig, values in dataset.depths.items():
            start: int | None = None
            for index, value in enumerate(values):
                callable_base = minimum_depth <= value <= summary.high_depth_threshold
                if callable_base and start is None:
                    start = index
                elif not callable_base and start is not None:
                    handle.write(f"{contig}\t{start}\t{index}\n")
                    start = None
            if start is not None:
                handle.write(f"{contig}\t{start}\t{len(values)}\n")
