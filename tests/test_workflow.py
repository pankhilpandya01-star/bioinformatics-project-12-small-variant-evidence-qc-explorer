from __future__ import annotations

import dataclasses
import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from small_variant_explorer.workflow import (
    PRIVATE_MARKERS,
    ToolPaths,
    WorkflowConfig,
    WorkflowError,
    resolve_tools,
    run_variant_ready_workflow,
    sha256_file,
    validate_config,
    validate_tool_version,
)
from tests.helpers import CONTROLLED, create_controlled_bam


class WorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.tools = resolve_tools()

    def config(self, bam: Path, output: Path) -> WorkflowConfig:
        return WorkflowConfig(
            reference_fasta=CONTROLLED / "reference.fasta",
            input_bam=bam,
            sample_id="controlled",
            output_dir=output,
        )

    def test_controlled_duplicate_marking_preserves_records(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bam = create_controlled_bam(root)
            before = sha256_file(bam)
            output = root / "result"
            run_variant_ready_workflow(self.config(bam, output), tools=self.tools)

            manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["counts"]["input_records"], 6)
            self.assertEqual(manifest["counts"]["output_records"], 6)
            self.assertEqual(manifest["counts"]["input_duplicate_records"], 0)
            self.assertEqual(manifest["counts"]["duplicate_records"], 2)
            self.assertEqual(manifest["counts"]["nonduplicate_records"], 4)
            self.assertAlmostEqual(manifest["depth"]["mean_depth"], 0.2)
            self.assertFalse(manifest["dataset_gate"]["passed"])
            self.assertEqual(sha256_file(bam), before)

            prepared = output / "bam" / "variant_ready.sorted.bam"
            check = subprocess.run(
                [self.tools.samtools, "quickcheck", "-v", str(prepared)],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
            )
            self.assertEqual(check.returncode, 0, check.stderr)
            header = subprocess.run(
                [self.tools.samtools, "view", "--no-PG", "-H", str(prepared)],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
            ).stdout
            self.assertNotIn("@PG\t", header)
            self.assertTrue((output / "reports" / "markdup.json").is_file())
            self.assertTrue((output / "depth" / "position_depth.tsv.gz").is_file())
            self.assertFalse((output / "intermediate").exists())

            for path in output.rglob("*"):
                if path.is_file() and path.suffix.lower() in {".txt", ".json", ".csv", ".tsv", ".bed"}:
                    text = path.read_text(encoding="utf-8", errors="replace").lower()
                    for marker in PRIVATE_MARKERS:
                        self.assertNotIn(marker, text, f"{marker} in {path}")

    def test_sample_and_reference_disagreement_are_rejected_atomically(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bam = create_controlled_bam(root)
            for changes, message in (
                ({"sample_id": "wrong"}, "do not match"),
                ({"reference_fasta": root / "wrong.fasta"}, "do not exactly match"),
            ):
                if "reference_fasta" in changes:
                    changes["reference_fasta"].write_text(">other\nAAAA\n", encoding="ascii")
                output = root / f"result-{message.split()[0]}"
                with self.subTest(changes=changes), self.assertRaisesRegex(WorkflowError, message):
                    run_variant_ready_workflow(
                        dataclasses.replace(self.config(bam, output), **changes),
                        tools=self.tools,
                    )
                self.assertFalse(output.exists())
                self.assertEqual(list(root.glob(f".{output.name}.staging-*")), [])

    def test_forced_subprocess_failure_cleans_staging(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bam = create_controlled_bam(root)
            output = root / "result"
            with patch(
                "small_variant_explorer.workflow._run_logged",
                side_effect=WorkflowError("forced failure"),
            ), self.assertRaisesRegex(WorkflowError, "forced failure"):
                run_variant_ready_workflow(self.config(bam, output), tools=self.tools)
            self.assertFalse(output.exists())
            self.assertEqual(list(root.glob(".result.staging-*")), [])

    def test_invalid_values_and_existing_output_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bam = create_controlled_bam(root)
            output = root / "result"
            cases = (
                ({"ploidy": 3}, "ploidy"),
                ({"ploidy": True}, "ploidy"),
                ({"sample_id": "bad sample"}, "sample ID"),
                ({"sample_id": "-bad"}, "sample ID"),
                ({"minimum_mapq": 256}, "MAPQ"),
                ({"minimum_mapq": True}, "MAPQ"),
                ({"minimum_base_quality": -1}, "base quality"),
                ({"minimum_depth": 0}, "minimum depth"),
                ({"minimum_alt_fraction": 0}, "alternate fraction"),
                ({"minimum_alt_fraction": float("nan")}, "alternate fraction"),
                ({"minimum_alt_fraction": float("inf")}, "alternate fraction"),
                ({"minimum_alt_fraction": True}, "alternate fraction"),
                ({"minimum_variant_quality": -1}, "variant quality"),
                ({"minimum_variant_quality": float("nan")}, "variant quality"),
                ({"minimum_variant_quality": float("inf")}, "variant quality"),
                ({"maximum_depth_factor": 1}, "depth factor"),
                ({"maximum_depth_factor": float("nan")}, "depth factor"),
                ({"maximum_depth_factor": float("inf")}, "depth factor"),
                ({"threads": 0}, "threads"),
            )
            for changes, message in cases:
                with self.subTest(changes=changes), self.assertRaisesRegex(WorkflowError, message):
                    validate_config(dataclasses.replace(self.config(bam, output), **changes))
            output.mkdir()
            with self.assertRaisesRegex(WorkflowError, "already exists"):
                validate_config(self.config(bam, output))

    def test_missing_index_and_wrong_tool_version_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bam = create_controlled_bam(root)
            Path(str(bam) + ".bai").unlink()
            with self.assertRaisesRegex(WorkflowError, "index is missing"):
                run_variant_ready_workflow(self.config(bam, root / "result"), tools=self.tools)

        wrong = subprocess.CompletedProcess(
            args=[], returncode=0, stdout="samtools 1.23\n", stderr=""
        )
        with patch("small_variant_explorer.workflow._run_capture", return_value=wrong):
            with self.assertRaisesRegex(WorkflowError, "Samtools 1.24"):
                validate_tool_version(ToolPaths("samtools"))


if __name__ == "__main__":
    unittest.main()
