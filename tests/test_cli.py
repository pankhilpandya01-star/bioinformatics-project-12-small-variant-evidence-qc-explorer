from __future__ import annotations

import unittest

from small_variant_explorer.cli import build_parser


class CliTests(unittest.TestCase):
    def test_public_options_are_exposed(self) -> None:
        parser = build_parser()
        options = {option for action in parser._actions for option in action.option_strings}
        expected = {
            "--reference-fasta",
            "--input-bam",
            "--sample-id",
            "--ploidy",
            "--minimum-mapq",
            "--minimum-base-quality",
            "--minimum-depth",
            "--minimum-alt-depth",
            "--minimum-alt-fraction",
            "--minimum-variant-quality",
            "--minimum-alt-strand-depth",
            "--maximum-depth-factor",
            "--threads",
            "--output-dir",
        }
        self.assertTrue(expected <= options)


if __name__ == "__main__":
    unittest.main()
