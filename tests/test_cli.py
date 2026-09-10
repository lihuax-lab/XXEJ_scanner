from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from XXEJ_scanner.cli import build_parser, run_scan
from XXEJ_scanner.models import CandidateRegion, ScannerConfig


class SkipChrMTest(unittest.TestCase):
    def test_skip_chrm_flag(self) -> None:
        args = build_parser().parse_args(
            [
                "scan",
                "--treated-bam",
                "treated.bam",
                "--reference-fasta",
                "reference.fa",
                "--output-dir",
                "out",
                "--skip-chrm",
            ]
        )

        self.assertTrue(args.skip_chrm)

    def test_skip_chrm_filters_regions_before_evidence_extraction(self) -> None:
        regions = [
            CandidateRegion("chr1", 10, 20, "region_1"),
            CandidateRegion("chrM", 10, 20, "region_2"),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            config = ScannerConfig(
                "treated.bam",
                "reference.fa",
                tmpdir,
                skip_chrm=True,
                evidence_backend="python",
            )
            with (
                patch("XXEJ_scanner.cli.validate_inputs"),
                patch(
                    "XXEJ_scanner.cli.call_candidate_regions", return_value=regions
                ),
                patch(
                    "XXEJ_scanner.cli.call_structural_evidence_regions",
                    return_value=[],
                ),
                patch(
                    "XXEJ_scanner.cli.annotate_region_coverage",
                    side_effect=lambda candidate_regions, *_: candidate_regions,
                ),
                patch("XXEJ_scanner.cli.iter_region_evidence", return_value=iter(())),
                patch("XXEJ_scanner.cli.ReferenceGenome"),
            ):
                summary = run_scan(config)

            self.assertEqual(summary["candidate_regions"], 1)
            self.assertEqual(
                Path(tmpdir, "candidate_regions.bed").read_text().split("\t", 1)[0],
                "chr1",
            )


if __name__ == "__main__":
    unittest.main()
