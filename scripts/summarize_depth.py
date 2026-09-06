"""Summarize a three-column ``samtools depth`` stream for the dataset gate."""

from __future__ import annotations

import sys


def main() -> int:
    positions = 0
    depth_sum = 0
    covered = 0
    depth_at_least_ten = 0
    maximum = 0

    for line_number, line in enumerate(sys.stdin, start=1):
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 3:
            raise ValueError(f"depth line {line_number} does not contain three fields")
        depth = int(fields[2])
        positions += 1
        depth_sum += depth
        covered += depth > 0
        depth_at_least_ten += depth >= 10
        maximum = max(maximum, depth)

    if positions == 0:
        raise ValueError("depth stream is empty")

    values = {
        "positions": positions,
        "mean_depth": depth_sum / positions,
        "covered_positions": covered,
        "covered_percent": 100 * covered / positions,
        "positions_ge_10": depth_at_least_ten,
        "percent_ge_10": 100 * depth_at_least_ten / positions,
        "max_depth": maximum,
    }
    for name, value in values.items():
        print(f"{name}={value:.6f}" if isinstance(value, float) else f"{name}={value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
