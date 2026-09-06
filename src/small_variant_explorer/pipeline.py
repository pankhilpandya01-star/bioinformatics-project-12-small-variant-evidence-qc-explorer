"""Atomic end-to-end preparation and small-variant evidence pipeline."""

from __future__ import annotations

import json
import os
import shutil
import tempfile
from dataclasses import replace
from pathlib import Path

from small_variant_explorer import __version__
from small_variant_explorer.calling import run_variant_calling_workflow
from small_variant_explorer.workflow import (
    WorkflowConfig,
    WorkflowError,
    run_variant_ready_workflow,
    sha256_file,
    validate_config,
)


def run_small_variant_workflow(config: WorkflowConfig) -> Path:
    """Prepare the BAM, call candidates, and publish both stages together."""
    config = validate_config(config)
    input_hash = sha256_file(config.input_bam)
    reference_hash = sha256_file(config.reference_fasta)
    parent = config.output_dir.parent
    parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{config.output_dir.name}.staging-", dir=parent))
    try:
        preparation = run_variant_ready_workflow(
            replace(config, output_dir=staging / "preparation")
        )
        calling = run_variant_calling_workflow(
            replace(
                config,
                input_bam=preparation / "bam" / "variant_ready.sorted.bam",
                output_dir=staging / "calling",
            )
        )
        if sha256_file(config.input_bam) != input_hash:
            raise WorkflowError("source BAM checksum changed during the complete workflow")
        if sha256_file(config.reference_fasta) != reference_hash:
            raise WorkflowError("reference checksum changed during the complete workflow")
        calling_manifest = json.loads(
            (calling / "run_manifest.json").read_text(encoding="utf-8")
        )
        manifest = {
            "workflow_version": __version__,
            "milestone": "small_variant_workflow",
            "inputs": {
                "reference_sha256": reference_hash,
                "bam_sha256": input_hash,
            },
            "stages": {
                "preparation": "preparation/run_manifest.json",
                "calling": "calling/run_manifest.json",
            },
            "counts": calling_manifest["counts"],
            "publication": "Both stages were published together after integrity checks passed.",
        }
        (staging / "run_manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8", newline="\n"
        )
        # The public output appears only after preparation and calling both pass.
        os.replace(staging, config.output_dir)
        return config.output_dir
    except Exception as error:
        shutil.rmtree(staging, ignore_errors=True)
        if isinstance(error, WorkflowError):
            raise
        raise WorkflowError(f"complete small-variant workflow failed: {error}") from error
