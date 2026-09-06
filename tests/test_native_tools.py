from __future__ import annotations

import subprocess
import tempfile
import unittest
from pathlib import Path

from small_variant_explorer.calling import resolve_variant_tools
from tests.helpers import create_bam_from_sam


class NativeToolIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.tools = resolve_variant_tools()

    def test_mapq_base_quality_and_duplicate_filters_control_depth(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bam = create_bam_from_sam(
                root,
                "@HD\tVN:1.6\tSO:coordinate\n"
                "@SQ\tSN:chr1\tLN:20\n"
                "@RG\tID:rg\tSM:sample\n"
                "good\t0\tchr1\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:rg\n"
                "low_mapq\t0\tchr1\t1\t19\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:rg\n"
                "low_baseq\t0\tchr1\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\t4444444444\tRG:Z:rg\n"
                "duplicate\t1024\tchr1\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:rg\n",
                "quality",
            )
            filtered = _depths(self.tools.samtools, bam, 20, 20)
            unfiltered = _depths(self.tools.samtools, bam, 0, 0)
            self.assertEqual(filtered, [1] * 10 + [0] * 10)
            self.assertEqual(unfiltered, [3] * 10 + [0] * 10)

    def test_bcftools_splits_multiallelics_left_aligns_and_deduplicates(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference = root / "reference.fasta"
            reference.write_text(">chr1\nAAAAACCCCCGGGGGTTTTT\n", encoding="ascii")
            subprocess.run(
                [self.tools.samtools, "faidx", str(reference)],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
            )
            source = root / "input.vcf"
            source.write_text(
                "##fileformat=VCFv4.2\n"
                "##contig=<ID=chr1,length=20>\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "chr1\t4\t.\tAA\tA\t60\t.\t.\n"
                "chr1\t4\t.\tAA\tA\t60\t.\t.\n"
                "chr1\t11\t.\tG\tA,T\t60\t.\t.\n",
                encoding="ascii",
                newline="\n",
            )
            output = root / "normalized.vcf.gz"
            result = subprocess.run(
                [
                    self.tools.bcftools,
                    "norm",
                    "--no-version",
                    "--fasta-ref",
                    str(reference),
                    "--multiallelics",
                    "-any",
                    "--rm-dup",
                    "exact",
                    "--output-type",
                    "z",
                    "--output",
                    str(output),
                    str(source),
                ],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            records = subprocess.run(
                [self.tools.bcftools, "query", "-f", "%POS\\t%REF\\t%ALT\\n", str(output)],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
            ).stdout.splitlines()
            self.assertEqual(records, ["1\tAA\tA", "11\tG\tA", "11\tG\tT"])


def _depths(samtools: str, bam: Path, base_quality: int, mapq: int) -> list[int]:
    result = subprocess.run(
        [
            samtools,
            "depth",
            "-aa",
            "-q",
            str(base_quality),
            "-Q",
            str(mapq),
            str(bam),
        ],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    return [int(line.split("\t")[2]) for line in result.stdout.splitlines()]


if __name__ == "__main__":
    unittest.main()
