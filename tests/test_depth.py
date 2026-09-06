from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from small_variant_explorer.depth import collect_depth, write_depth_outputs
from small_variant_explorer.validation import ValidationError


class DepthTests(unittest.TestCase):
    def test_depth_summary_and_outputs_use_fixed_denominators(self) -> None:
        lines = ["chr1\t1\t0\n", "chr1\t2\t10\n", "chr1\t3\t10\n", "chr1\t4\t100\n"]
        dataset = collect_depth(
            lines,
            (("chr1", 4),),
            minimum_depth=10,
            maximum_depth_factor=3,
        )
        self.assertEqual(dataset.summary.positions, 4)
        self.assertEqual(dataset.summary.mean_depth, 30)
        self.assertEqual(dataset.summary.median_depth, 10)
        self.assertEqual(dataset.summary.high_depth_threshold, 50)
        self.assertEqual(dataset.summary.positions_ge_minimum, 3)
        self.assertEqual(dataset.summary.callable_positions, 2)
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            write_depth_outputs(dataset, output, minimum_depth=10, window_size=2)
            self.assertEqual(
                (output / "callable_regions.bed").read_text(encoding="utf-8"),
                "chr1\t1\t3\n",
            )
            self.assertEqual(len((output / "genome_windows.csv").read_text().splitlines()), 3)

    def test_missing_or_out_of_order_positions_are_rejected(self) -> None:
        with self.assertRaisesRegex(ValidationError, "expected 2"):
            collect_depth(
                ["chr1\t1\t1\n", "chr1\t3\t1\n"],
                (("chr1", 2),),
                minimum_depth=1,
                maximum_depth_factor=3,
            )


if __name__ == "__main__":
    unittest.main()
