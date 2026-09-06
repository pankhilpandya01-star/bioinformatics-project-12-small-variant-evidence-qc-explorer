"""Command-line interface for small-variant calling and evidence QC."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from small_variant_explorer import __version__
from small_variant_explorer.pipeline import run_small_variant_workflow
from small_variant_explorer.validation import ValidationError
from small_variant_explorer.workflow import WorkflowConfig, WorkflowError


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="small-variant-explorer",
        description="Prepare a BAM, call normalized small variants, and evaluate their evidence.",
    )
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--input-bam", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--ploidy", type=int, choices=(1, 2), default=1)
    parser.add_argument("--minimum-mapq", type=int, default=20)
    parser.add_argument("--minimum-base-quality", type=int, default=20)
    parser.add_argument("--minimum-depth", type=int, default=10)
    parser.add_argument("--minimum-alt-depth", type=int, default=5)
    parser.add_argument("--minimum-alt-fraction", type=float, default=0.8)
    parser.add_argument("--minimum-variant-quality", type=float, default=30.0)
    parser.add_argument("--minimum-alt-strand-depth", type=int, default=1)
    parser.add_argument("--maximum-depth-factor", type=float, default=3.0)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        output = run_small_variant_workflow(
            WorkflowConfig(
                reference_fasta=args.reference_fasta,
                input_bam=args.input_bam,
                sample_id=args.sample_id,
                output_dir=args.output_dir,
                ploidy=args.ploidy,
                minimum_mapq=args.minimum_mapq,
                minimum_base_quality=args.minimum_base_quality,
                minimum_depth=args.minimum_depth,
                minimum_alt_depth=args.minimum_alt_depth,
                minimum_alt_fraction=args.minimum_alt_fraction,
                minimum_variant_quality=args.minimum_variant_quality,
                minimum_alt_strand_depth=args.minimum_alt_strand_depth,
                maximum_depth_factor=args.maximum_depth_factor,
                threads=args.threads,
            )
        )
    except (WorkflowError, ValidationError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    print(f"Completed small-variant evidence workflow: {output}")
    return 0
