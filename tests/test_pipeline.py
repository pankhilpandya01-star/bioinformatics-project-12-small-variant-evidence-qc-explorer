from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from small_variant_explorer.controlled import build_controlled_variant_fixture
from small_variant_explorer.pipeline import run_small_variant_workflow
from small_variant_explorer.workflow import WorkflowConfig, WorkflowError, sha256_file


class CompletePipelineTests(unittest.TestCase):
    def test_controlled_workflow_publishes_both_stages_atomically(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture = root / "fixture"
            bam = build_controlled_variant_fixture(fixture)
            reference = fixture / "reference.fasta"
            bam_hash = sha256_file(bam)
            reference_hash = sha256_file(reference)
            output = root / "complete"
            run_small_variant_workflow(
                WorkflowConfig(reference, bam, "controlled", output)
            )
            manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["counts"]["normalized_candidates"], 5)
            self.assertEqual(manifest["counts"]["accepted_candidates"], 4)
            self.assertEqual(manifest["counts"]["rejected_candidates"], 1)
            self.assertTrue((output / "preparation" / "bam" / "variant_ready.sorted.bam").is_file())
            self.assertTrue((output / "calling" / "vcf" / "accepted.vcf.gz").is_file())
            self.assertEqual(sha256_file(bam), bam_hash)
            self.assertEqual(sha256_file(reference), reference_hash)

    def test_calling_failure_removes_the_complete_workflow_staging_tree(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture = root / "fixture"
            bam = build_controlled_variant_fixture(fixture)
            reference = fixture / "reference.fasta"
            output = root / "complete"
            with patch(
                "small_variant_explorer.pipeline.run_variant_calling_workflow",
                side_effect=WorkflowError("forced complete-workflow failure"),
            ), self.assertRaisesRegex(WorkflowError, "forced complete-workflow failure"):
                run_small_variant_workflow(
                    WorkflowConfig(reference, bam, "controlled", output)
                )
            self.assertFalse(output.exists())
            self.assertEqual(list(root.glob(".complete.staging-*")), [])


if __name__ == "__main__":
    unittest.main()
