from __future__ import annotations

import csv
import gzip
import json
import unittest
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT / "results" / "variant_ready"
CONTROLLED_RESULTS = PROJECT / "results" / "controlled_variant_calling"
CONTROLLED_TRUTH = PROJECT / "examples" / "controlled" / "variant_fixture" / "expected_variants.csv"
AUTHENTIC_RESULTS = PROJECT / "results" / "authentic_analysis"


class PublishedResultTests(unittest.TestCase):
    def test_duplicate_and_depth_totals_reconcile(self) -> None:
        manifest = json.loads((RESULTS / "run_manifest.json").read_text(encoding="utf-8"))
        markdup = json.loads((RESULTS / "markdup.json").read_text(encoding="utf-8"))
        flagstat = json.loads((RESULTS / "flagstat.json").read_text(encoding="utf-8"))
        with (RESULTS / "coverage_summary.csv").open(
            "r", encoding="utf-8", newline=""
        ) as handle:
            coverage = {row["metric"]: row["value"] for row in csv.DictReader(handle)}

        counts = manifest["counts"]
        self.assertEqual(counts["input_records"], counts["output_records"])
        self.assertEqual(
            counts["duplicate_records"] + counts["nonduplicate_records"],
            counts["output_records"],
        )
        self.assertEqual(markdup["READ"], counts["input_records"])
        self.assertEqual(markdup["WRITTEN"], counts["output_records"])
        self.assertEqual(markdup["DUPLICATE TOTAL"], counts["duplicate_records"])
        self.assertEqual(flagstat["QC-passed reads"]["total"], counts["output_records"])
        self.assertEqual(
            flagstat["QC-passed reads"]["duplicates"], counts["duplicate_records"]
        )
        self.assertEqual(int(coverage["positions"]), manifest["depth"]["positions"])
        self.assertAlmostEqual(
            float(coverage["mean_depth"]), manifest["depth"]["mean_depth"]
        )
        self.assertAlmostEqual(
            float(coverage["callable_percent"]), manifest["depth"]["callable_percent"]
        )
        self.assertTrue(manifest["dataset_gate"]["passed"])

    def test_controlled_truth_and_vcf_partitions_reconcile(self) -> None:
        manifest = json.loads(
            (CONTROLLED_RESULTS / "run_manifest.json").read_text(encoding="utf-8")
        )
        with CONTROLLED_TRUTH.open("r", encoding="utf-8", newline="") as handle:
            truth = {
                f"{row['chrom']}:{row['position']}:{row['ref']}:{row['alt']}": row[
                    "expected_status"
                ]
                for row in csv.DictReader(handle)
            }
        with (CONTROLLED_RESULTS / "reports" / "variant_evidence.csv").open(
            "r", encoding="utf-8", newline=""
        ) as handle:
            evidence = {row["variant_key"]: row["status"] for row in csv.DictReader(handle)}
        normalized = _vcf_keys(
            CONTROLLED_RESULTS / "vcf" / "candidates.normalized.vcf.gz"
        )
        accepted = _vcf_keys(CONTROLLED_RESULTS / "vcf" / "accepted.vcf.gz")
        rejected = _vcf_keys(CONTROLLED_RESULTS / "vcf" / "rejected.vcf.gz")
        self.assertEqual(evidence, truth)
        self.assertEqual(manifest["counts"]["normalized_candidates"], len(normalized))
        self.assertEqual(manifest["counts"]["accepted_candidates"], len(accepted))
        self.assertEqual(manifest["counts"]["rejected_candidates"], len(rejected))
        self.assertTrue(accepted.isdisjoint(rejected))
        self.assertEqual(normalized, accepted | rejected)

    def test_authentic_vcf_csv_dashboard_and_source_totals_reconcile(self) -> None:
        manifest = json.loads(
            (AUTHENTIC_RESULTS / "run_manifest.json").read_text(encoding="utf-8")
        )
        preparation_manifest = json.loads(
            (RESULTS / "run_manifest.json").read_text(encoding="utf-8")
        )
        dataset_manifest = json.loads(
            (PROJECT / "results" / "dataset_review" / "run_manifest.json").read_text(
                encoding="utf-8"
            )
        )
        with (AUTHENTIC_RESULTS / "reports" / "variant_evidence.csv").open(
            "r", encoding="utf-8", newline=""
        ) as handle:
            evidence_rows = list(csv.DictReader(handle))
        evidence_keys = [row["variant_key"] for row in evidence_rows]
        normalized_list = _vcf_key_list(
            AUTHENTIC_RESULTS / "vcf" / "candidates.normalized.vcf.gz"
        )
        evaluated = _vcf_keys(AUTHENTIC_RESULTS / "vcf" / "evaluated.vcf.gz")
        accepted = _vcf_keys(AUTHENTIC_RESULTS / "vcf" / "accepted.vcf.gz")
        rejected = _vcf_keys(AUTHENTIC_RESULTS / "vcf" / "rejected.vcf.gz")
        counts = manifest["counts"]

        self.assertEqual(len(normalized_list), len(set(normalized_list)))
        self.assertEqual(len(evidence_keys), len(set(evidence_keys)))
        self.assertEqual(set(normalized_list), set(evidence_keys))
        self.assertEqual(set(normalized_list), evaluated)
        self.assertTrue(accepted.isdisjoint(rejected))
        self.assertEqual(set(normalized_list), accepted | rejected)
        self.assertEqual(counts["normalized_candidates"], len(normalized_list))
        self.assertEqual(counts["accepted_candidates"], len(accepted))
        self.assertEqual(counts["rejected_candidates"], len(rejected))
        self.assertEqual(
            manifest["inputs"]["bam"]["sha256"],
            preparation_manifest["outputs"]["bam_sha256"],
        )
        self.assertEqual(
            manifest["inputs"]["reference"]["sha256"],
            dataset_manifest["checksums"]["reference_sha256"],
        )
        dashboard = AUTHENTIC_RESULTS / "dashboard" / "variant_evidence_dashboard.png"
        self.assertGreater(dashboard.stat().st_size, 100_000)
        self.assertEqual(dashboard.read_bytes()[:8], b"\x89PNG\r\n\x1a\n")


def _vcf_keys(path: Path) -> set[str]:
    return set(_vcf_key_list(path))


def _vcf_key_list(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [
            f"{fields[0]}:{fields[1]}:{fields[3]}:{fields[4]}"
            for line in handle
            if not line.startswith("#")
            for fields in (line.rstrip("\n").split("\t"),)
        ]


if __name__ == "__main__":
    unittest.main()
