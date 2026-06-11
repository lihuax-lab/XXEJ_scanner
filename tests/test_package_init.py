from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import XXEJ_scanner


class PackageInitTest(unittest.TestCase):
    def test_exports_version_string(self) -> None:
        self.assertIsInstance(XXEJ_scanner.__version__, str)
        self.assertTrue(XXEJ_scanner.__version__)

    def test_public_api_exports_stable_entrypoints(self) -> None:
        expected_names = {
            "ScannerConfig",
            "CandidateRegion",
            "RegionEvidence",
            "BreakpointCluster",
            "RepairEvent",
            "EventEvidence",
            "ReferenceGenome",
            "run_scan",
            "parse_bed_regions",
            "call_candidate_regions",
            "collect_region_evidence",
            "cluster_clip_sites",
            "cluster_evidence_graph",
            "classify_local_events",
            "classify_bnd_events",
            "find_microhomology",
            "write_events_tsv",
        }

        self.assertTrue(expected_names.issubset(set(XXEJ_scanner.__all__)))
        for name in expected_names:
            self.assertTrue(hasattr(XXEJ_scanner, name), name)

    def test_public_api_does_not_export_private_names(self) -> None:
        private_names = [
            name
            for name in XXEJ_scanner.__all__
            if name.startswith("_") and name != "__version__"
        ]

        self.assertEqual(private_names, [])


if __name__ == "__main__":
    unittest.main()
