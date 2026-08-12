from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from XXEJ_scanner.models import RepairEvent, ScannerConfig
from XXEJ_scanner.validation import assign_event_filter, second_pass_validate_event


def scanner_config(**overrides: object) -> ScannerConfig:
    values = {
        "treated_bam": "treated.bam",
        "reference_fasta": "reference.fa",
        "output_dir": "out",
        "min_alt_support": 3,
        "min_bnd_support": 3,
    }
    values.update(overrides)
    return ScannerConfig(**values)


class PairOnlyBndFilterTest(unittest.TestCase):
    def test_pair_only_bnd_filter_survives_second_pass(self) -> None:
        event = RepairEvent(
            event_id="TMP_BND_PAIRONLY_region1_chr1_100_chr2_500",
            event_type="BND_INTER",
            chrom="chr1",
            start=99,
            end=101,
            bkp_A_chrom="chr1",
            bkp_A_pos=100,
            bkp_A_side="pair_only",
            bkp_B_chrom="chr2",
            bkp_B_pos=500,
            bkp_B_side="remote",
            remote_chrom="chr2",
            remote_pos=500,
            alt_discordant_pair_support=1,
            support_read_names={"pair1"},
        )
        config = scanner_config()

        self.assertEqual(assign_event_filter(event, config), "PairOnlyBnd")

        validated = second_pass_validate_event(event, "missing.bam", config)

        self.assertEqual(validated.filter, "PairOnlyBnd")


if __name__ == "__main__":
    unittest.main()
