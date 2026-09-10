from __future__ import annotations

import sys
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from XXEJ_scanner.genotype import count_ref_like_breakends, update_event_fraction
from XXEJ_scanner.models import RepairEvent, ScannerConfig


def config() -> ScannerConfig:
    return ScannerConfig("treated.bam", "reference.fa", "out")


class BreakendReferenceSupportTest(unittest.TestCase):
    @patch("XXEJ_scanner.genotype.count_spanning_reads", side_effect=[7, 3])
    def test_deletion_counts_each_breakend_independently(self, count) -> None:
        event = RepairEvent(
            "event1",
            "LOCAL_DEL",
            "chr1",
            100,
            10_000,
            "chr1",
            100,
            "junction",
            bkp_B_chrom="chr1",
            bkp_B_pos=10_000,
            bkp_B_side="junction",
        )

        support = count_ref_like_breakends("treated.bam", event, config())

        self.assertEqual(support, (7, 3))
        self.assertEqual(count.call_count, 2)

    def test_unresolved_event_fraction_remains_unavailable(self) -> None:
        event = RepairEvent(
            "event1", "BND_INTER", "chr1", 99, 101, "chr1", 100, "pair_only"
        )
        event.ref_spanning_support = 10
        event.support_read_names = {"r1"}

        update_event_fraction(event)

        self.assertEqual(event.repair_evidence_fraction, "NA")


if __name__ == "__main__":
    unittest.main()
