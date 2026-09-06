from __future__ import annotations

import shutil
import subprocess
import unittest
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
CONTROLLED = PROJECT / "examples" / "controlled"


def samtools_path() -> str:
    path = shutil.which("samtools")
    if path is None:
        raise unittest.SkipTest("samtools is not available")
    return path


def create_controlled_bam(directory: Path) -> Path:
    samtools = shutil.which("samtools")
    if samtools is None:
        raise RuntimeError("samtools is required for controlled BAM tests")
    bam = directory / "controlled.sorted.bam"
    subprocess.run(
        [samtools, "view", "--no-PG", "-b", "-o", str(bam), str(CONTROLLED / "duplicate_example.sam")],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    subprocess.run(
        [samtools, "index", "-@", "1", str(bam)],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    return bam


def create_bam_from_sam(directory: Path, sam_text: str, name: str = "input") -> Path:
    """Create and index a small coordinate-sorted BAM from controlled SAM text."""
    samtools = shutil.which("samtools")
    if samtools is None:
        raise RuntimeError("samtools is required for controlled BAM tests")
    sam = directory / f"{name}.sam"
    bam = directory / f"{name}.sorted.bam"
    sam.write_text(sam_text, encoding="ascii", newline="\n")
    subprocess.run(
        [samtools, "view", "--no-PG", "-b", "-o", str(bam), str(sam)],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    subprocess.run(
        [samtools, "index", "-@", "1", str(bam)],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    return bam
