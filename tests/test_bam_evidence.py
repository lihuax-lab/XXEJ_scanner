"""Focused tests for BAM evidence classification."""

import sys
from types import SimpleNamespace
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from XXEJ_scanner.bam_evidence import (
    extract_discordant_pair_from_read,
    is_discordant_pair,
)
from XXEJ_scanner.models import CandidateRegion, ScannerConfig


def make_config(**overrides):
    values = {
        "treated_bam": "sample.bam",
        "reference_fasta": "genome.fa",
        "output_dir": "out",
    }
    values.update(overrides)
    return ScannerConfig(**values)


def make_read(**overrides):
    values = {
        "query_name": "read-1",
        "is_unmapped": False,
        "is_secondary": False,
        "is_supplementary": False,
        "is_duplicate": False,
        "is_paired": True,
        "mate_is_unmapped": False,
        "reference_name": "chr1",
        "next_reference_name": "chr1",
        "reference_start": 100,
        "next_reference_start": 150,
        "template_length": 50,
        "is_reverse": False,
        "mate_is_reverse": True,
        "mapping_quality": 60,
        "cigartuples": [(0, 50)],
    }
    values.update(overrides)
    return SimpleNamespace(**values)


class DiscordantPairTests(unittest.TestCase):
    def test_different_chromosome_pair_is_discordant(self):
        config = make_config()
        read = make_read(next_reference_name="chr2")

        is_discordant, reason = is_discordant_pair(read, config)

        self.assertTrue(is_discordant)
        self.assertEqual(reason, "different_chrom")

    def test_proper_pair_outside_region_is_not_discordant(self):
        config = make_config()
        region = CandidateRegion("chr11", 68034323, 68034478, "peak-1")
        read = make_read(
            query_name="E250150437L1C004R0291164482",
            reference_name="chr11",
            next_reference_name="chr11",
            reference_start=68034010,
            next_reference_start=68034010,
            template_length=124,
            cigartuples=[(0, 124)],
        )

        pair = extract_discordant_pair_from_read(read, config, region)

        self.assertIsNone(pair)

    def test_distant_mate_keeps_outside_region_reason(self):
        config = make_config(discordant_min_distance=100, max_insert_size=10_000)
        region = CandidateRegion("chr1", 50, 200, "peak-1")
        read = make_read(next_reference_start=700, template_length=650)

        pair = extract_discordant_pair_from_read(read, config, region)

        self.assertIsNotNone(pair)
        assert pair is not None
        self.assertEqual(pair.reason, "distant_mate,mate_outside_candidate_region")


if __name__ == "__main__":
    unittest.main()
