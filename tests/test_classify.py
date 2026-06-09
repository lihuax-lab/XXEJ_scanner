from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from XXEJ_scanner.classify import _classify_mmej_del, assign_final_event_ids
from XXEJ_scanner.models import (
    BreakpointCluster,
    CandidateRegion,
    CigarIndel,
    ClipSite,
    EventEvidence,
    RegionEvidence,
    RepairEvent,
    ScannerConfig,
    SplitReadEvidence,
)


class FakeReference:
    def __init__(self, sequence: str) -> None:
        self.sequence = sequence

    def fetch(self, _chrom: str, start: int, end: int) -> str:
        start = max(0, min(start, len(self.sequence)))
        end = max(0, min(end, len(self.sequence)))
        if end <= start:
            return ""
        return self.sequence[start:end].upper()


def sequence_with_matches(matches: dict[int, str], length: int = 180) -> str:
    filler = "0123456789BDEFHIJKLMNOPQRSTUVWXYZ!#$%&()*+,-./:;<=>?@[]^_{|}~"
    sequence = list((filler * ((length // len(filler)) + 1))[:length])
    for start, bases in matches.items():
        sequence[start : start + len(bases)] = bases
    return "".join(sequence)


def scanner_config(**overrides: object) -> ScannerConfig:
    values = {
        "treated_bam": "treated.bam",
        "reference_fasta": "reference.fa",
        "output_dir": "out",
        "min_alt_support": 1,
        "clip_cluster_window": 5,
        "max_local_event_distance": 1000,
    }
    values.update(overrides)
    return ScannerConfig(**values)


def region() -> CandidateRegion:
    return CandidateRegion("chr1", 0, 200, "region1")


def cluster(pos: int, side: str, count: int = 1) -> BreakpointCluster:
    left_count = count if side in {"left_clip", "both"} else 0
    right_count = count if side in {"right_clip", "both"} else 0
    if side == "both":
        left_count = max(1, count // 2)
        right_count = count - left_count
    return BreakpointCluster(
        region_id="region1",
        chrom="chr1",
        cluster_start=pos,
        cluster_end=pos + 1,
        peak_pos=pos,
        clip_side=side,
        clip_count=count,
        left_clip_count=left_count,
        right_clip_count=right_count,
        treated_depth=20,
    )


def clip_site(pos: int, side: str, read_name: str, sequence: str = "NNNN") -> ClipSite:
    return ClipSite(
        chrom="chr1",
        pos=pos,
        side=side,
        clip_length=len(sequence),
        clip_sequence=sequence,
        read_name=read_name,
        strand="+",
        mapq=60,
        cigar="10M4S" if side == "right_clip" else "4S10M",
        is_reverse=False,
        reference_start=max(0, pos - 10),
        reference_end=pos,
    )


def deletion_indel(start: int, end: int, read_name: str) -> CigarIndel:
    return CigarIndel(
        chrom="chr1",
        start=start,
        end=end,
        operation="DEL",
        length=end - start,
        sequence="NA",
        read_name=read_name,
        mapq=60,
        cigar=f"10M{end - start}D10M",
    )


def split_read(pos: int, remote_pos: int, read_name: str) -> SplitReadEvidence:
    return SplitReadEvidence(
        read_name=read_name,
        chrom="chr1",
        pos=pos,
        side="SA",
        remote_chrom="chr1",
        remote_pos=remote_pos,
        remote_strand="+",
        remote_cigar="10M",
        remote_mapq=60,
        orientation="++",
        mapq=60,
        cigar="10M",
        sa_tag="chr1,1,+,10M,60,0;",
    )


def repair_event(event_id: str, pos: int = 10) -> RepairEvent:
    return RepairEvent(
        event_id=event_id,
        event_type="NHEJ_INS",
        chrom="chr1",
        start=pos,
        end=pos + 1,
        bkp_A_chrom="chr1",
        bkp_A_pos=pos,
        bkp_A_side="left_clip",
    )


def event_evidence(event_id: str, read_name: str, pos: int = 10) -> EventEvidence:
    return EventEvidence(
        event_id=event_id,
        read_name=read_name,
        evidence_type="soft_clip",
        chrom="chr1",
        pos=pos,
    )


class FinalEventIdAssignmentTest(unittest.TestCase):
    def test_assigns_global_ids_and_updates_evidence_links(self) -> None:
        events = [
            repair_event("TMP_INS_region1_chr1_10_left_clip", 10),
            repair_event("TMP_INS_region2_chr1_20_left_clip", 20),
        ]
        evidence = [
            event_evidence("TMP_INS_region1_chr1_10_left_clip", "read1", 10),
            event_evidence("TMP_INS_region2_chr1_20_left_clip", "read2", 20),
        ]

        assign_final_event_ids(events, evidence)

        self.assertEqual([event.event_id for event in events], ["XEJ_000001", "XEJ_000002"])
        self.assertEqual(
            [row.event_id for row in evidence],
            ["XEJ_000001", "XEJ_000002"],
        )

    def test_rejects_duplicate_temporary_ids_before_relinking_evidence(self) -> None:
        events = [
            repair_event("TMP_DUPLICATE", 10),
            repair_event("TMP_DUPLICATE", 20),
        ]
        evidence = [event_evidence("TMP_DUPLICATE", "read1", 10)]

        with self.assertRaisesRegex(ValueError, "Temporary event IDs are not unique"):
            assign_final_event_ids(events, evidence)


class ClassifyMmejDeletionTest(unittest.TestCase):
    def test_mmej_pair_is_sorted_internally(self) -> None:
        left = cluster(10, "right_clip")
        right = cluster(20, "left_clip")
        evidence = RegionEvidence(
            region=region(),
            clip_sites=[
                clip_site(10, "right_clip", "left_clip_read"),
                clip_site(20, "left_clip", "right_clip_read"),
            ],
        )

        events, _event_evidence, used = _classify_mmej_del(
            region(),
            [right, left],
            evidence,
            FakeReference(sequence_with_matches({})),
            scanner_config(),
        )

        self.assertEqual(len(events), 1)
        self.assertEqual(events[0].bkp_A_pos, 10)
        self.assertEqual(events[0].bkp_A_side, "right_clip")
        self.assertEqual(events[0].bkp_B_pos, 20)
        self.assertEqual(events[0].bkp_B_side, "left_clip")
        self.assertIn(("chr1", 10, "right_clip"), used)
        self.assertIn(("chr1", 20, "left_clip"), used)

    def test_incompatible_clip_orientation_is_not_mmej(self) -> None:
        left = cluster(10, "right_clip")
        right = cluster(20, "right_clip")
        evidence = RegionEvidence(
            region=region(),
            clip_sites=[
                clip_site(10, "right_clip", "left_read"),
                clip_site(20, "right_clip", "right_read"),
            ],
        )

        events, _event_evidence, used = _classify_mmej_del(
            region(),
            [left, right],
            evidence,
            FakeReference(sequence_with_matches({})),
            scanner_config(),
        )

        self.assertEqual(events, [])
        self.assertEqual(used, set())

    def test_pair_scoring_can_prefer_lower_clip_support_with_mh_and_junction(self) -> None:
        high_clip_left = cluster(10, "right_clip", count=5)
        high_clip_right = cluster(30, "left_clip", count=5)
        mh_left = cluster(100, "right_clip", count=3)
        mh_right = cluster(120, "left_clip", count=3)
        evidence = RegionEvidence(
            region=region(),
            clip_sites=[
                clip_site(10, "right_clip", "a1"),
                clip_site(30, "left_clip", "a2"),
                clip_site(100, "right_clip", "b1"),
                clip_site(120, "left_clip", "b2"),
            ],
            indels=[
                deletion_indel(100, 126, "indel1"),
                deletion_indel(100, 126, "indel2"),
            ],
        )

        events, _event_evidence, _used = _classify_mmej_del(
            region(),
            [high_clip_left, high_clip_right, mh_left, mh_right],
            evidence,
            FakeReference(sequence_with_matches({94: "GATTAC", 120: "GATTAC"})),
            scanner_config(max_local_event_distance=50, max_microhomology_length=6),
        )

        self.assertEqual(len(events), 1)
        self.assertEqual(events[0].bkp_A_pos, 100)
        self.assertEqual(events[0].bkp_B_pos, 120)
        self.assertEqual(events[0].microhomology, "GATTAC")
        self.assertEqual(events[0].alt_indel_support, 2)
        self.assertIn("cigar_del", events[0].junction_evidence_types)

    def test_adjusted_microhomology_span_matches_cigar_deletion(self) -> None:
        left = cluster(100, "right_clip")
        right = cluster(120, "left_clip")
        evidence = RegionEvidence(
            region=region(),
            clip_sites=[
                clip_site(100, "right_clip", "left_read"),
                clip_site(120, "left_clip", "right_read"),
            ],
            indels=[deletion_indel(102, 123, "indel1")],
        )

        events, _event_evidence, _used = _classify_mmej_del(
            region(),
            [left, right],
            evidence,
            FakeReference(sequence_with_matches({98: "TTGA", 119: "TTGA"})),
            scanner_config(
                clip_cluster_window=1,
                min_microhomology_length=4,
                max_microhomology_length=4,
                microhomology_search_window=3,
            ),
        )

        self.assertEqual(len(events), 1)
        event = events[0]
        self.assertEqual(event.microhomology_left_end, 102)
        self.assertEqual(event.microhomology_right_start, 119)
        self.assertEqual(event.microhomology_deletion_start, 102)
        self.assertEqual(event.microhomology_deletion_end, 123)
        self.assertEqual(event.microhomology_deletion_length, 21)
        self.assertEqual(event.alt_indel_support, 1)
        self.assertEqual(event.junction_evidence_support, 1)
        self.assertEqual(event.junction_evidence_types, {"cigar_del"})

    def test_split_read_and_soft_clip_remap_count_as_junction_evidence(self) -> None:
        left = cluster(50, "right_clip")
        right = cluster(80, "left_clip")
        evidence = RegionEvidence(
            region=region(),
            clip_sites=[
                clip_site(50, "right_clip", "left_remap", "RGHT"),
                clip_site(80, "left_clip", "right_remap", "LEFT"),
            ],
            split_reads=[split_read(50, 80, "split1")],
        )

        events, event_evidence, _used = _classify_mmej_del(
            region(),
            [left, right],
            evidence,
            FakeReference(sequence_with_matches({46: "LEFT", 80: "RGHT"})),
            scanner_config(
                clip_cluster_window=4,
                min_clip_length=4,
                microhomology_search_window=0,
            ),
        )

        self.assertEqual(len(events), 1)
        event = events[0]
        self.assertEqual(event.junction_evidence_support, 3)
        self.assertEqual(
            event.junction_evidence_types,
            {"soft_clip_remap", "split_read_sa"},
        )
        self.assertIn("split_read_sa", {row.evidence_type for row in event_evidence})


if __name__ == "__main__":
    unittest.main()
