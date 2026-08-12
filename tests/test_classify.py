from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from XXEJ_scanner.classify import (
    _classify_local_del,
    assign_final_event_ids,
    classify_bnd_events,
    classify_local_events,
)
from XXEJ_scanner.models import (
    BreakpointCluster,
    CandidateRegion,
    CigarIndel,
    ClipSite,
    DiscordantPair,
    EventEvidence,
    RegionEvidence,
    RepairEvent,
    ScannerConfig,
    SplitReadEvidence,
)
from XXEJ_scanner.validation import assign_event_filter


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
        remote_nm=0,
        orientation="++",
        mapq=60,
        cigar="10M",
        sa_tag="chr1,1,+,10M,60,0;",
    )


def discordant_pair(
    pos: int, remote_chrom: str, remote_pos: int, read_name: str
) -> DiscordantPair:
    return DiscordantPair(
        read_name=read_name,
        chrom="chr1",
        pos=pos,
        mate_chrom=remote_chrom,
        mate_pos=remote_pos,
        orientation="+-",
        mapq=60,
        is_reverse=False,
        mate_is_reverse=True,
        cigar="100M",
        reason="different_chrom",
    )


def repair_event(event_id: str, pos: int = 10) -> RepairEvent:
    return RepairEvent(
        event_id=event_id,
        event_type="LOCAL_INS",
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


class ClassifyLocalDeletionTest(unittest.TestCase):
    def test_clip_pair_without_a_junction_is_not_a_deletion(self) -> None:
        left = cluster(10, "right_clip")
        right = cluster(20, "left_clip")
        evidence = RegionEvidence(
            region=region(),
            clip_sites=[
                clip_site(10, "right_clip", "left_clip_read"),
                clip_site(20, "left_clip", "right_clip_read"),
            ],
        )

        events, event_evidence = _classify_local_del(
            region(),
            [right, left],
            evidence,
            FakeReference(sequence_with_matches({})),
            scanner_config(),
        )

        self.assertEqual(events, [])
        self.assertEqual(event_evidence, [])

    def test_reports_each_read_resolved_deletion_allele(self) -> None:
        evidence = RegionEvidence(
            region=region(),
            indels=[
                deletion_indel(10, 20, "a1"),
                deletion_indel(10, 20, "a2"),
                deletion_indel(30, 40, "b1"),
                deletion_indel(30, 40, "b2"),
            ],
        )

        events, _event_evidence = _classify_local_del(
            region(),
            [],
            evidence,
            FakeReference(sequence_with_matches({})),
            scanner_config(min_alt_support=2),
        )

        self.assertEqual([(event.start, event.end) for event in events], [(10, 20), (30, 40)])
        self.assertTrue(all(event.event_type == "LOCAL_DEL" for event in events))
        self.assertTrue(all(event.junction_resolved for event in events))

    def test_microhomology_is_annotated_only_at_resolved_coordinates(self) -> None:
        left = cluster(100, "right_clip")
        right = cluster(123, "left_clip")
        evidence = RegionEvidence(
            region=region(),
            indels=[deletion_indel(102, 123, "indel1")],
        )

        events, _event_evidence = _classify_local_del(
            region(),
            [left, right],
            evidence,
            FakeReference(sequence_with_matches({98: "TTGA", 123: "TTGA"})),
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
        self.assertEqual(event.microhomology_right_start, 123)
        self.assertEqual(event.microhomology_offset_a, 0)
        self.assertEqual(event.microhomology_offset_b, 0)
        self.assertEqual(event.microhomology_deletion_start, 102)
        self.assertEqual(event.microhomology_deletion_end, 127)
        self.assertEqual(event.microhomology_deletion_length, 25)
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

        events, event_evidence = _classify_local_del(
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

    def test_local_deletion_and_remote_bnd_are_both_reported(self) -> None:
        left = cluster(50, "right_clip")
        right = cluster(80, "left_clip")
        evidence = RegionEvidence(
            region=region(),
            clip_sites=[
                clip_site(50, "right_clip", "left_read"),
                clip_site(80, "left_clip", "right_read"),
            ],
            indels=[deletion_indel(50, 80, "del1")],
            discordant_pairs=[discordant_pair(50, "chr2", 500, "pair1")],
        )

        events, event_evidence = classify_local_events(
            region(),
            [left, right],
            evidence,
            FakeReference(sequence_with_matches({})),
            scanner_config(min_bnd_support=1),
        )

        event_types = [event.event_type for event in events]
        self.assertEqual(event_types, ["LOCAL_DEL", "BND_INTER"])
        self.assertEqual(events[1].bkp_A_pos, 50)
        self.assertEqual(events[1].remote_chrom, "chr2")
        self.assertIn(
            "discordant_pair", {row.evidence_type for row in event_evidence}
        )


class ClassifyBndEventsTest(unittest.TestCase):
    def test_pair_only_discordant_pair_emits_bnd_without_cluster(self) -> None:
        config = scanner_config(min_alt_support=3, min_bnd_support=3)
        evidence = RegionEvidence(
            region=region(),
            discordant_pairs=[discordant_pair(50, "chr2", 500, "pair1")],
        )

        events, event_evidence = classify_local_events(
            region(),
            [],
            evidence,
            FakeReference(sequence_with_matches({})),
            config,
        )

        self.assertEqual(len(events), 1)
        event = events[0]
        self.assertEqual(event.event_type, "BND_INTER")
        self.assertEqual(event.bkp_A_pos, 50)
        self.assertEqual(event.bkp_A_side, "pair_only")
        self.assertEqual(event.remote_chrom, "chr2")
        self.assertEqual(event.alt_clip_support, 0)
        self.assertEqual(event.alt_discordant_pair_support, 1)
        self.assertEqual(assign_event_filter(event, config), "PairOnlyBnd")
        self.assertIn(
            "discordant_pair", {row.evidence_type for row in event_evidence}
        )

    def test_discordant_pair_does_not_need_to_be_near_peak(self) -> None:
        local_cluster = cluster(100, "left_clip")
        evidence = RegionEvidence(
            region=region(),
            clip_sites=[clip_site(100, "left_clip", "clip1")],
            discordant_pairs=[discordant_pair(20, "chr2", 500, "pair1")],
        )

        events, event_evidence = classify_bnd_events(
            region(),
            [local_cluster],
            evidence,
            scanner_config(min_bnd_support=1),
        )

        self.assertEqual(len(events), 1)
        self.assertEqual(events[0].event_type, "BND_INTER")
        self.assertEqual(events[0].bkp_A_pos, 100)
        self.assertIn(
            "discordant_pair", {row.evidence_type for row in event_evidence}
        )

    def test_discordant_pair_is_assigned_to_nearest_cluster_once(self) -> None:
        left_cluster = cluster(100, "left_clip")
        right_cluster = cluster(160, "right_clip")
        evidence = RegionEvidence(
            region=region(),
            clip_sites=[
                clip_site(100, "left_clip", "clip1"),
                clip_site(160, "right_clip", "clip2"),
            ],
            discordant_pairs=[discordant_pair(150, "chr2", 500, "pair1")],
        )

        events, _event_evidence = classify_bnd_events(
            region(),
            [left_cluster, right_cluster],
            evidence,
            scanner_config(min_bnd_support=1),
        )

        self.assertEqual(len(events), 1)
        self.assertEqual(events[0].bkp_A_pos, 160)

    def test_nearby_same_chromosome_split_is_not_bnd(self) -> None:
        local_cluster = cluster(100, "left_clip")
        evidence = RegionEvidence(
            region=region(),
            clip_sites=[clip_site(100, "left_clip", "clip1")],
            split_reads=[split_read(100, 120, "split1")],
        )

        events, event_evidence = classify_bnd_events(
            region(),
            [local_cluster],
            evidence,
            scanner_config(min_bnd_support=1, max_local_event_distance=1000),
        )

        self.assertEqual(events, [])
        self.assertEqual(event_evidence, [])

    def test_distant_same_chromosome_split_is_bnd_intra(self) -> None:
        local_cluster = cluster(100, "left_clip")
        evidence = RegionEvidence(
            region=region(),
            clip_sites=[clip_site(100, "left_clip", "clip1")],
            split_reads=[split_read(100, 2000, "split1")],
        )

        events, event_evidence = classify_bnd_events(
            region(),
            [local_cluster],
            evidence,
            scanner_config(min_bnd_support=1, max_local_event_distance=1000),
        )

        self.assertEqual(len(events), 1)
        self.assertEqual(events[0].event_type, "BND_INTRA")
        self.assertEqual(events[0].bkp_B_pos, 2000)
        self.assertIn("split_read_sa", {row.evidence_type for row in event_evidence})

    def test_split_anchor_absorbs_nearby_pair_across_bin_boundary(self) -> None:
        local_cluster = cluster(100, "left_clip")
        split = split_read(100, 1999, "split1")
        pair = discordant_pair(100, "chr1", 2001, "pair1")
        evidence = RegionEvidence(
            region=region(),
            split_reads=[split],
            discordant_pairs=[pair],
        )

        events, _rows = classify_bnd_events(
            region(),
            [local_cluster],
            evidence,
            scanner_config(
                min_bnd_support=1,
                max_local_event_distance=1000,
                coverage_bin_size=10,
            ),
        )

        self.assertEqual(len(events), 1)
        self.assertEqual(events[0].alt_split_support, 1)
        self.assertEqual(events[0].alt_discordant_pair_support, 1)
        self.assertEqual(events[0].remote_pos, 2000)


if __name__ == "__main__":
    unittest.main()
