"""Focused tests for BAM evidence classification."""

import sys
from types import SimpleNamespace
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from XXEJ_scanner.bam_evidence import (
    extract_cigar_indels_from_read,
    extract_clip_sites_from_read,
    extract_discordant_pair_from_read,
    extract_split_reads_from_sa_tag,
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
    query_length = sum(
        length for op, length in values["cigartuples"] if op in {0, 1, 4, 7, 8}
    )
    if "query_sequence" not in overrides:
        values["query_sequence"] = "A" * query_length
    if "query_qualities" not in overrides:
        values["query_qualities"] = [40] * query_length
    return SimpleNamespace(**values)


class BreakpointQualityTests(unittest.TestCase):
    def test_soft_clip_requires_eighty_percent_q20(self):
        qualities = [40] * 50
        qualities[5:7] = [10, 10]
        passing = make_read(cigartuples=[(4, 10), (0, 40)], query_qualities=qualities)

        self.assertEqual(len(extract_clip_sites_from_read(passing, make_config())), 1)

        qualities[7] = 10
        failing = make_read(cigartuples=[(4, 10), (0, 40)], query_qualities=qualities)
        self.assertEqual(extract_clip_sites_from_read(failing, make_config()), [])

    def test_missing_qualities_and_hard_clips_fail_when_filtering_is_enabled(self):
        missing = make_read(cigartuples=[(4, 10), (0, 40)], query_qualities=None)
        hard_clipped = make_read(cigartuples=[(5, 10), (0, 40)])

        self.assertEqual(extract_clip_sites_from_read(missing, make_config()), [])
        self.assertEqual(extract_clip_sites_from_read(hard_clipped, make_config()), [])
        self.assertEqual(
            len(
                extract_clip_sites_from_read(
                    hard_clipped, make_config(min_breakpoint_baseq=0)
                )
            ),
            1,
        )

    def test_indels_require_quality_at_each_resolved_boundary(self):
        deletion = make_read(cigartuples=[(0, 20), (2, 5), (0, 20)])
        insertion = make_read(cigartuples=[(0, 20), (1, 10), (0, 20)])

        self.assertEqual(len(extract_cigar_indels_from_read(deletion, make_config())), 1)
        self.assertEqual(len(extract_cigar_indels_from_read(insertion, make_config())), 1)

        insertion.query_qualities[30:33] = [10, 10, 10]
        self.assertEqual(extract_cigar_indels_from_read(insertion, make_config()), [])

    def test_strict_baseq_override_is_used(self):
        read = make_read(
            cigartuples=[(4, 10), (0, 40)], query_qualities=[22] * 50
        )

        self.assertEqual(len(extract_clip_sites_from_read(read, make_config())), 1)
        self.assertEqual(
            extract_clip_sites_from_read(read, make_config(), min_baseq=25), []
        )


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


class SplitReadTests(unittest.TestCase):
    def test_sa_breakpoints_use_query_order_and_cigar_ends(self):
        read = make_read(
            is_paired=False,
            reference_start=100,
            cigartuples=[(0, 50), (4, 50)],
            has_tag=lambda tag: tag == "SA",
            get_tag=lambda tag: "chr2,501,+,50S50M,60,2;",
        )

        splits = extract_split_reads_from_sa_tag(read, make_config())

        self.assertEqual(len(splits), 1)
        self.assertEqual((splits[0].pos, splits[0].side), (150, "right_clip"))
        self.assertEqual(splits[0].remote_pos, 500)
        self.assertEqual(splits[0].remote_nm, 2)

    def test_reverse_sa_uses_reference_end_for_query_start(self):
        read = make_read(
            is_paired=False,
            reference_start=100,
            cigartuples=[(0, 50), (4, 50)],
            has_tag=lambda tag: tag == "SA",
            get_tag=lambda tag: "chr2,501,-,50M50S,60,1;",
        )

        split = extract_split_reads_from_sa_tag(read, make_config())[0]

        self.assertEqual(split.pos, 150)
        self.assertEqual(split.remote_pos, 550)

    def test_sa_remote_quality_filters_are_enforced(self):
        read = make_read(
            is_paired=False,
            cigartuples=[(0, 50), (4, 50)],
            has_tag=lambda tag: tag == "SA",
            get_tag=lambda tag: "chr2,501,+,50S50M,40,11;",
        )

        splits = extract_split_reads_from_sa_tag(
            read, make_config(max_sa_nm=10), min_mapq=50
        )

        self.assertEqual(splits, [])

    def test_sa_junction_requires_local_base_quality(self):
        read = make_read(
            is_paired=False,
            reference_start=100,
            cigartuples=[(0, 50), (4, 50)],
            query_qualities=[40] * 45 + [10] * 10 + [40] * 45,
            has_tag=lambda tag: tag == "SA",
            get_tag=lambda tag: "chr2,501,+,50S50M,60,2;",
        )

        self.assertEqual(extract_split_reads_from_sa_tag(read, make_config()), [])


if __name__ == "__main__":
    unittest.main()
