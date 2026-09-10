from __future__ import annotations

import sys
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from XXEJ_scanner.coverage import call_candidate_regions
from XXEJ_scanner.models import ScannerConfig


class CoverageThresholdTest(unittest.TestCase):
    def test_explicit_minimum_replaces_percentile(self) -> None:
        bins = {
            ("chr1", 0, 10): 1.0,
            ("chr1", 20, 30): 6.0,
            ("chr1", 40, 50): 10.0,
        }
        config = ScannerConfig(
            "control-only.bam",
            "reference.fa",
            "out",
            top_percentile=100,
            min_treated_coverage=5,
            merge_distance=0,
        )

        with patch(
            "XXEJ_scanner.coverage.compute_binned_coverage", return_value=bins
        ):
            regions = call_candidate_regions("control-only.bam", config)

        self.assertEqual([region.treated_coverage for region in regions], [6.0, 10.0])


if __name__ == "__main__":
    unittest.main()
