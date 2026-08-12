from __future__ import annotations

import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from XXEJ_scanner.cli import _deduplicate_bnd_events
from XXEJ_scanner.coverage import call_structural_evidence_regions
from XXEJ_scanner.models import CigarIndel, DiscordantPair, RepairEvent, ScannerConfig


class FakeBam:
    def __init__(self, reads: list[object]):
        self.reads = reads

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return None

    def fetch(self):
        return iter(self.reads)


class StructuralDiscoveryTest(unittest.TestCase):
    def test_discovers_exact_indels_and_bounded_remote_clusters(self) -> None:
        reads = [SimpleNamespace(kind="ins", query_name=f"i{i}") for i in range(3)]
        reads += [SimpleNamespace(kind="pair", query_name=f"p{i}") for i in range(3)]
        config = ScannerConfig(
            treated_bam="treated.bam",
            reference_fasta="ref.fa",
            output_dir="out",
            min_alt_support=3,
            min_bnd_support=3,
            coverage_bin_size=10,
            scan_padding=10,
        )

        def indels(read, _config):
            if read.kind != "ins":
                return []
            return [
                CigarIndel(
                    "chr1", 100, 100, "INS", 2, "AA", read.query_name, 60, "20M2I20M"
                )
            ]

        def pairs(read, _config):
            if read.kind != "pair":
                return None
            offset = int(read.query_name[-1])
            return DiscordantPair(
                read.query_name,
                "chr1",
                995 + offset * 5,
                "chr2",
                5000 + offset,
                "+-",
                60,
                False,
                True,
                "50M",
                "different_chrom",
            )

        with (
            patch(
                "XXEJ_scanner.coverage.pysam.AlignmentFile",
                return_value=FakeBam(reads),
            ),
            patch("XXEJ_scanner.coverage.passes_read_filters", return_value=True),
            patch(
                "XXEJ_scanner.coverage.extract_cigar_indels_from_read",
                side_effect=indels,
            ),
            patch(
                "XXEJ_scanner.coverage.extract_discordant_pair_from_read",
                side_effect=pairs,
            ),
            patch(
                "XXEJ_scanner.coverage.extract_split_reads_from_sa_tag",
                return_value=[],
            ),
        ):
            regions = call_structural_evidence_regions("treated.bam", config)

        self.assertEqual(
            [
                (region.chrom, region.start, region.end, region.score)
                for region in regions
            ],
            [("chr1", 90, 111, 3.0), ("chr1", 985, 1016, 3.0)],
        )

    def test_collapses_mirrored_bnd_calls(self) -> None:
        left = RepairEvent(
            "left",
            "BND_INTER",
            "chr1",
            99,
            101,
            "chr1",
            100,
            "pair_only",
            bkp_B_chrom="chr2",
            bkp_B_pos=5000,
            orientation="+-",
            alt_discordant_pair_support=3,
            support_read_names={"r1", "r2", "r3"},
        )
        right = RepairEvent(
            "right",
            "BND_INTER",
            "chr2",
            4999,
            5001,
            "chr2",
            5000,
            "left_clip",
            bkp_B_chrom="chr1",
            bkp_B_pos=100,
            orientation="-+",
            alt_split_support=3,
            junction_resolved=True,
            support_read_names={"r1", "r2", "r3"},
        )

        events, evidence = _deduplicate_bnd_events([left, right], [], 100)

        self.assertEqual([event.event_id for event in events], ["right"])
        self.assertEqual(evidence, [])


if __name__ == "__main__":
    unittest.main()
