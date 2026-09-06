from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from small_variant_explorer.validation import ValidationError, parse_bam_header, read_fasta_lengths


class FastaValidationTests(unittest.TestCase):
    def test_valid_fasta_lengths_are_preserved_in_order(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reference.fasta"
            path.write_text(">a first\nACGT\n>b\nNN\n", encoding="ascii")
            self.assertEqual(read_fasta_lengths(path), (("a", 4), ("b", 2)))

    def test_malformed_fasta_is_rejected(self) -> None:
        cases = (
            ("ACGT\n", "before the first header"),
            (">x\n", "no sequence"),
            (">x\nAZ\n", "invalid symbols"),
            (">x\nA C\n", "whitespace"),
            (">x\nA\n>x\nC\n", "duplicate"),
        )
        for content, message in cases:
            with self.subTest(content=content), tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "reference.fasta"
                path.write_text(content, encoding="ascii")
                with self.assertRaisesRegex(ValidationError, message):
                    read_fasta_lengths(path)


class BamHeaderTests(unittest.TestCase):
    def test_required_header_fields_are_parsed(self) -> None:
        header = parse_bam_header(
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chr1\tLN:1000\n"
            "@RG\tID:rg1\tSM:sample\tLB:library\n"
        )
        self.assertEqual(header.sort_order, "coordinate")
        self.assertEqual(header.contigs, (("chr1", 1000),))
        self.assertEqual(header.read_groups[0]["SM"], "sample")

    def test_missing_contigs_and_read_groups_are_rejected(self) -> None:
        with self.assertRaisesRegex(ValidationError, "no @SQ"):
            parse_bam_header("@HD\tSO:coordinate\n@RG\tID:x\tSM:y\n")
        with self.assertRaisesRegex(ValidationError, "no @RG"):
            parse_bam_header("@HD\tSO:coordinate\n@SQ\tSN:x\tLN:1\n")


if __name__ == "__main__":
    unittest.main()
