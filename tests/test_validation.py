from __future__ import annotations

import sys
import unittest
from unittest.mock import patch
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from XXEJ_scanner.models import (
    CandidateRegion,
    CigarIndel,
    DiscordantPair,
    RegionEvidence,
    RepairEvent,
    ScannerConfig,
    SplitReadEvidence,
)
from XXEJ_scanner.validation import (
    assign_event_filter,
    matching_event_read_names,
    second_pass_validate_event,
)


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


def local_event(event_type: str = "LOCAL_INS") -> RepairEvent:
    return RepairEvent(
        event_id="event1",
        event_type=event_type,
        chrom="chr1",
        start=99,
        end=101 if event_type == "LOCAL_INS" else 120,
        bkp_A_chrom="chr1",
        bkp_A_pos=100,
        bkp_A_side="junction",
        bkp_B_chrom="chr1" if event_type == "LOCAL_DEL" else "NA",
        bkp_B_pos=120 if event_type == "LOCAL_DEL" else "NA",
        bkp_B_side="junction" if event_type == "LOCAL_DEL" else "NA",
        inserted_sequence="AC" if event_type == "LOCAL_INS" else "NA",
        inserted_length=2 if event_type == "LOCAL_INS" else "NA",
        deleted_length=20 if event_type == "LOCAL_DEL" else "NA",
        score=1,
        evidence_level="RESOLVED",
        junction_resolved=True,
        support_read_names={"r1", "r2", "r3"},
    )


def indel(start: int, end: int, sequence: str, read_name: str) -> CigarIndel:
    operation = "INS" if start == end else "DEL"
    return CigarIndel(
        chrom="chr1",
        start=start,
        end=end,
        operation=operation,
        length=len(sequence) if operation == "INS" else end - start,
        sequence=sequence,
        read_name=read_name,
        mapq=60,
        cigar="test",
    )


class EventSpecificEvidenceTest(unittest.TestCase):
    def test_insertion_requires_exact_position_and_sequence(self) -> None:
        event = local_event()
        evidence = RegionEvidence(
            CandidateRegion("chr1", 0, 200, "r"),
            indels=[
                indel(100, 100, "AC", "match"),
                indel(100, 100, "GT", "wrong_sequence"),
                indel(101, 101, "AC", "wrong_position"),
            ],
        )

        self.assertEqual(
            matching_event_read_names(event, evidence, scanner_config()), {"match"}
        )

    def test_deletion_requires_both_exact_endpoints(self) -> None:
        event = local_event("LOCAL_DEL")
        evidence = RegionEvidence(
            CandidateRegion("chr1", 0, 200, "r"),
            indels=[
                indel(100, 120, "NA", "match"),
                indel(100, 121, "NA", "wrong_end"),
            ],
        )

        self.assertEqual(
            matching_event_read_names(event, evidence, scanner_config()), {"match"}
        )

    def test_bnd_matches_both_directions_and_split_orientation(self) -> None:
        event = RepairEvent(
            event_id="bnd1",
            event_type="BND_INTER",
            chrom="chr1",
            start=99,
            end=101,
            bkp_A_chrom="chr1",
            bkp_A_pos=100,
            bkp_A_side="right_clip",
            bkp_B_chrom="chr2",
            bkp_B_pos=500,
            bkp_B_side="remote",
            remote_chrom="chr2",
            remote_pos=500,
            orientation="+-",
        )
        evidence = RegionEvidence(
            CandidateRegion("chr1", 0, 200, "r"),
            discordant_pairs=[
                DiscordantPair(
                    "pair", "chr2", 500, "chr1", 100, "-+", 60, True, False, "50M", "different_chrom"
                )
            ],
            split_reads=[
                SplitReadEvidence(
                    "split", "chr1", 100, "right_clip", "chr2", 500, "-", "50M", 60, 0, "+-", 60, "50M", "SA"
                ),
                SplitReadEvidence(
                    "wrong_orientation", "chr1", 100, "right_clip", "chr2", 500, "+", "50M", 60, 0, "++", 60, "50M", "SA"
                ),
            ],
        )

        self.assertEqual(
            matching_event_read_names(event, evidence, scanner_config()),
            {"pair", "split"},
        )

    def test_control_event_and_missing_coverage_are_filtered(self) -> None:
        event = local_event()
        event.control_assessed = True
        event.control_alt_support = 2
        self.assertEqual(assign_event_filter(event, scanner_config()), "PresentInControl")

        event.control_alt_support = 0
        event.control_ref_support = 0
        event.control_depth = 0
        self.assertEqual(assign_event_filter(event, scanner_config()), "NoControlCoverage")

    @patch("XXEJ_scanner.validation.collect_region_evidence")
    def test_unrelated_strict_indels_cannot_rescue_event(self, collect) -> None:
        event = local_event()
        collect.return_value = RegionEvidence(
            CandidateRegion("chr1", 0, 200, "strict"),
            indels=[indel(100, 100, "GT", "wrong1"), indel(101, 101, "AC", "wrong2")],
        )

        validated = second_pass_validate_event(event, "treated.bam", scanner_config())

        self.assertEqual(validated.filter, "LowSupport")


if __name__ == "__main__":
    unittest.main()
