from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from XXEJ_scanner.breakpoints import cluster_evidence_graph
from XXEJ_scanner.models import (
    CandidateRegion,
    CigarIndel,
    ClipSite,
    RegionEvidence,
    ScannerConfig,
)


def scanner_config(**overrides: object) -> ScannerConfig:
    values = {
        "treated_bam": "treated.bam",
        "reference_fasta": "reference.fa",
        "output_dir": "out",
        "clip_cluster_window": 5,
    }
    values.update(overrides)
    return ScannerConfig(**values)


def region() -> CandidateRegion:
    return CandidateRegion("chr1", 0, 100, "region1")


def clip_site(pos: int, side: str, read_name: str) -> ClipSite:
    return ClipSite(
        chrom="chr1",
        pos=pos,
        side=side,
        clip_length=10,
        clip_sequence="AAAAAAAAAA",
        read_name=read_name,
        strand="+",
        mapq=60,
        cigar="10M10S" if side == "right_clip" else "10S10M",
        is_reverse=False,
        reference_start=max(0, pos - 10),
        reference_end=pos,
    )


def insertion(pos: int, read_name: str) -> CigarIndel:
    return CigarIndel(
        chrom="chr1",
        start=pos,
        end=pos,
        operation="INS",
        length=2,
        sequence="AA",
        read_name=read_name,
        mapq=60,
        cigar="10M2I10M",
    )


def deletion(start: int, end: int, read_name: str) -> CigarIndel:
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


class EvidenceGraphClusteringTest(unittest.TestCase):
    def test_shared_insertion_anchor_links_clips_beyond_window(self) -> None:
        evidence = RegionEvidence(
            region=region(),
            clip_sites=[
                clip_site(10, "right_clip", "clip1"),
                clip_site(20, "right_clip", "clip2"),
            ],
            indels=[insertion(15, "ins1")],
        )

        clusters = cluster_evidence_graph(evidence, scanner_config(), region=region())

        self.assertEqual(len(clusters), 1)
        self.assertEqual(clusters[0].cluster_start, 10)
        self.assertEqual(clusters[0].cluster_end, 21)
        self.assertEqual(clusters[0].clip_count, 2)

    def test_deletion_endpoints_remain_separate_anchors(self) -> None:
        evidence = RegionEvidence(
            region=region(),
            clip_sites=[
                clip_site(10, "right_clip", "left_clip"),
                clip_site(30, "left_clip", "right_clip"),
            ],
            indels=[deletion(10, 30, "del1")],
        )

        clusters = cluster_evidence_graph(evidence, scanner_config(), region=region())

        self.assertEqual([cluster.peak_pos for cluster in clusters], [10, 30])

    def test_region_filtering_matches_window_cluster_scope(self) -> None:
        evidence = RegionEvidence(
            region=region(),
            clip_sites=[
                clip_site(10, "right_clip", "inside"),
                ClipSite(
                    chrom="chr2",
                    pos=10,
                    side="right_clip",
                    clip_length=10,
                    clip_sequence="AAAAAAAAAA",
                    read_name="outside_chrom",
                    strand="+",
                    mapq=60,
                    cigar="10M10S",
                    is_reverse=False,
                    reference_start=0,
                    reference_end=10,
                ),
            ],
        )

        clusters = cluster_evidence_graph(evidence, scanner_config(), region=region())

        self.assertEqual(len(clusters), 1)
        self.assertEqual(clusters[0].chrom, "chr1")


if __name__ == "__main__":
    unittest.main()
