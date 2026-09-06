from __future__ import annotations

import csv
import gzip
import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from small_variant_explorer.calling import (
    VariantToolPaths,
    resolve_variant_tools,
    run_variant_calling_workflow,
    validate_bcftools_version,
)
from small_variant_explorer.controlled import CONTROLLED_VARIANTS, build_controlled_variant_fixture
from small_variant_explorer.workflow import (
    PRIVATE_MARKERS,
    WorkflowConfig,
    WorkflowError,
    sha256_file,
)
from tests.helpers import CONTROLLED, create_bam_from_sam


class VariantCallingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.tools = resolve_variant_tools()

    def test_controlled_variants_are_called_normalized_and_partitioned(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture = root / "fixture"
            bam = build_controlled_variant_fixture(fixture, self.tools.samtools)
            reference = fixture / "reference.fasta"
            bam_hash = sha256_file(bam)
            reference_hash = sha256_file(reference)
            output = root / "result"
            run_variant_calling_workflow(
                WorkflowConfig(
                    reference_fasta=reference,
                    input_bam=bam,
                    sample_id="controlled",
                    output_dir=output,
                ),
                tools=self.tools,
            )

            manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["counts"]["normalized_candidates"], 5)
            self.assertEqual(manifest["counts"]["accepted_candidates"], 4)
            self.assertEqual(manifest["counts"]["rejected_candidates"], 1)
            self.assertEqual(sha256_file(bam), bam_hash)
            self.assertEqual(sha256_file(reference), reference_hash)

            with (output / "reports" / "variant_evidence.csv").open(
                "r", encoding="utf-8", newline=""
            ) as handle:
                evidence = list(csv.DictReader(handle))
            expected = {
                f"control:{variant.position}:{variant.ref}:{variant.alt}": variant.expected_status
                for variant in CONTROLLED_VARIANTS
            }
            observed = {row["variant_key"]: row["status"] for row in evidence}
            self.assertEqual(observed, expected)
            weak = next(row for row in evidence if row["position"] == "1600")
            self.assertIn("LowAltFraction", weak["filter_reasons"])
            self.assertIn("OneStrandOnly", weak["filter_reasons"])
            dashboard = output / "dashboard" / "variant_evidence_dashboard.png"
            self.assertTrue(dashboard.is_file())
            self.assertGreater(dashboard.stat().st_size, 20_000)
            self.assertEqual(dashboard.read_bytes()[:8], b"\x89PNG\r\n\x1a\n")

            normalized_keys = _vcf_keys(output / "vcf" / "candidates.normalized.vcf.gz")
            evaluated_keys = _vcf_keys(output / "vcf" / "evaluated.vcf.gz")
            accepted_keys = _vcf_keys(output / "vcf" / "accepted.vcf.gz")
            rejected_keys = _vcf_keys(output / "vcf" / "rejected.vcf.gz")
            self.assertEqual(normalized_keys, evaluated_keys)
            self.assertTrue(accepted_keys.isdisjoint(rejected_keys))
            self.assertEqual(normalized_keys, accepted_keys | rejected_keys)

            for name in (
                "candidates.normalized.vcf.gz",
                "evaluated.vcf.gz",
                "accepted.vcf.gz",
                "rejected.vcf.gz",
            ):
                path = output / "vcf" / name
                self.assertTrue(path.is_file())
                self.assertTrue(Path(str(path) + ".csi").is_file())
                result = subprocess.run(
                    [self.tools.bcftools, "view", "--no-version", "--header-only", str(path)],
                    text=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    check=False,
                )
                self.assertEqual(result.returncode, 0, result.stderr)

            for path in output.rglob("*"):
                if path.is_file() and path.suffix.lower() in {
                    ".txt",
                    ".json",
                    ".csv",
                    ".tsv",
                    ".bed",
                }:
                    text = path.read_text(encoding="utf-8", errors="replace").lower()
                    for marker in (
                        *PRIVATE_MARKERS,
                        str(root).lower(),
                        str(fixture).lower(),
                        str(output).lower(),
                    ):
                        self.assertNotIn(marker, text, f"{marker} in {path}")

    def test_empty_callset_is_valid_and_fully_published(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bam = create_bam_from_sam(
                root,
                "@HD\tVN:1.6\tSO:coordinate\n"
                "@SQ\tSN:chr1\tLN:1000\n"
                "@RG\tID:controlled_rg\tSM:controlled\tLB:library\tPL:ILLUMINA\n",
                "empty",
            )
            output = root / "result"
            run_variant_calling_workflow(
                WorkflowConfig(CONTROLLED / "reference.fasta", bam, "controlled", output),
                tools=self.tools,
            )
            manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(
                manifest["counts"],
                {
                    "normalized_candidates": 0,
                    "accepted_candidates": 0,
                    "rejected_candidates": 0,
                },
            )
            for name in ("candidates.normalized", "evaluated", "accepted", "rejected"):
                path = output / "vcf" / f"{name}.vcf.gz"
                self.assertEqual(_vcf_keys(path), set())
                self.assertTrue(Path(str(path) + ".csi").is_file())
            self.assertTrue((output / "dashboard" / "variant_evidence_dashboard.png").is_file())

    def test_bad_bam_and_multiple_samples_are_rejected_atomically(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            malformed = root / "malformed.bam"
            malformed.write_bytes(b"not a BAM")
            Path(str(malformed) + ".bai").write_bytes(b"not an index")
            multisample = create_bam_from_sam(
                root,
                "@HD\tVN:1.6\tSO:coordinate\n"
                "@SQ\tSN:chr1\tLN:1000\n"
                "@RG\tID:rg1\tSM:controlled\n"
                "@RG\tID:rg2\tSM:other\n",
                "multisample",
            )
            for bam, message in (
                (malformed, "not identified as sequence data"),
                (multisample, "do not match"),
            ):
                output = root / f"result-{bam.stem}"
                with self.subTest(bam=bam), self.assertRaisesRegex(WorkflowError, message):
                    run_variant_calling_workflow(
                        WorkflowConfig(
                            CONTROLLED / "reference.fasta", bam, "controlled", output
                        ),
                        tools=self.tools,
                    )
                self.assertFalse(output.exists())
                self.assertEqual(list(root.glob(f".{output.name}.staging-*")), [])

    def test_caller_subprocess_failure_cleans_staging(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture = root / "fixture"
            bam = build_controlled_variant_fixture(fixture, self.tools.samtools)
            output = root / "result"
            with patch(
                "small_variant_explorer.calling._run_logged",
                side_effect=WorkflowError("forced caller failure"),
            ), self.assertRaisesRegex(WorkflowError, "forced caller failure"):
                run_variant_calling_workflow(
                    WorkflowConfig(fixture / "reference.fasta", bam, "controlled", output),
                    tools=self.tools,
                )
            self.assertFalse(output.exists())
            self.assertEqual(list(root.glob(".result.staging-*")), [])

    def test_wrong_bcftools_version_is_rejected(self) -> None:
        wrong = subprocess.CompletedProcess(
            args=[], returncode=0, stdout="bcftools 1.23\n", stderr=""
        )
        with patch("small_variant_explorer.calling._run_capture", return_value=wrong):
            with self.assertRaisesRegex(WorkflowError, "BCFtools 1.24"):
                validate_bcftools_version(VariantToolPaths("samtools", "bcftools"))


def _vcf_keys(path: Path) -> set[str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return {
            f"{fields[0]}:{fields[1]}:{fields[3]}:{fields[4]}"
            for line in handle
            if not line.startswith("#")
            for fields in (line.rstrip("\n").split("\t"),)
        }


if __name__ == "__main__":
    unittest.main()
