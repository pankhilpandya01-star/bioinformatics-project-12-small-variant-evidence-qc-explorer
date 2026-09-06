from __future__ import annotations

import csv
import gzip
import tempfile
import unittest
from pathlib import Path

from small_variant_explorer.variants import EvidenceThresholds, evaluate_normalized_vcf
from small_variant_explorer.validation import ValidationError


class VariantEvidenceTests(unittest.TestCase):
    def test_independent_filter_reasons_and_partitioning(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "candidates.vcf.gz"
            with gzip.open(source, "wt", encoding="utf-8", newline="\n") as handle:
                handle.write(
                    "##fileformat=VCFv4.2\n"
                    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
                    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Depth\">\n"
                    "##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Forward\">\n"
                    "##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Reverse\">\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
                    "chr1\t10\t.\tA\tG\t60\t.\t.\tGT:DP:AD:ADF:ADR\t1:10:0,10:0,5:0,5\n"
                    "chr1\t20\t.\tA\tAT\t10\t.\t.\tGT:DP:AD:ADF:ADR\t1:8:4,4:2,4:2,0\n"
                )
            intermediate = root / "intermediate"
            reports = root / "reports"
            intermediate.mkdir()
            reports.mkdir()
            summary = evaluate_normalized_vcf(
                source,
                intermediate,
                reports,
                EvidenceThresholds(10, 5, 0.8, 30, 1, 50),
            )
            self.assertEqual((summary.candidates, summary.accepted, summary.rejected), (2, 1, 1))
            self.assertEqual(summary.evidence[0].variant_type, "SNP")
            self.assertEqual(summary.evidence[0].filter_reasons, ())
            self.assertEqual(
                summary.evidence[1].filter_reasons,
                ("LowQual", "LowDepth", "LowAltDepth", "LowAltFraction", "OneStrandOnly"),
            )
            with (reports / "variant_evidence.csv").open(
                "r", encoding="utf-8", newline=""
            ) as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual({row["status"] for row in rows}, {"accepted", "rejected"})

    def test_high_depth_is_reported_independently(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "candidates.vcf.gz"
            _write_vcf(
                source,
                (
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample",
                    "chr1\t30\t.\tC\tT\t60\t.\t.\tGT:DP:AD:ADF:ADR\t1:51:0,51:0,26:0,25",
                ),
            )
            intermediate = root / "intermediate"
            reports = root / "reports"
            intermediate.mkdir()
            reports.mkdir()
            summary = evaluate_normalized_vcf(
                source,
                intermediate,
                reports,
                EvidenceThresholds(10, 5, 0.8, 30, 1, 50),
            )
            self.assertEqual(summary.evidence[0].filter_reasons, ("HighDepth",))

    def test_malformed_multisample_missing_format_and_duplicate_vcfs_are_rejected(self) -> None:
        cases = (
            (("chr1\t1\t.\tA\tG\t60\t.\t.\tGT:DP:AD:ADF:ADR\t1:10:0,10:0,5:0,5",), "#CHROM"),
            (
                (
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb",
                    "chr1\t1\t.\tA\tG\t60\t.\t.\tGT:DP:AD:ADF:ADR\t1:10:0,10:0,5:0,5\t1:10:0,10:0,5:0,5",
                ),
                "exactly one sample",
            ),
            (
                (
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample",
                    "chr1\t1\t.\tA\tG\t60\t.\t.\tGT:DP:ADF:ADR\t1:10:0,5:0,5",
                ),
                "FORMAT/AD",
            ),
            (
                (
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample",
                    "chr1\t1\t.\tA\tG\t60\t.\t.\tGT:DP:AD:ADF:ADR\t1:10:0,10:0,5:0,5",
                    "chr1\t1\t.\tA\tG\t60\t.\t.\tGT:DP:AD:ADF:ADR\t1:10:0,10:0,5:0,5",
                ),
                "duplicate normalized variant key",
            ),
        )
        for lines, message in cases:
            with self.subTest(message=message), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                source = root / "candidates.vcf.gz"
                _write_vcf(source, lines)
                intermediate = root / "intermediate"
                reports = root / "reports"
                intermediate.mkdir()
                reports.mkdir()
                with self.assertRaisesRegex(ValidationError, message):
                    evaluate_normalized_vcf(
                        source,
                        intermediate,
                        reports,
                        EvidenceThresholds(10, 5, 0.8, 30, 1, 50),
                    )


def _write_vcf(path: Path, body_lines: tuple[str, ...]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="\n") as handle:
        handle.write(
            "##fileformat=VCFv4.2\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Depth\">\n"
            "##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Forward\">\n"
            "##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Reverse\">\n"
        )
        handle.write("\n".join(body_lines) + "\n")


if __name__ == "__main__":
    unittest.main()
