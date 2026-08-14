from __future__ import annotations

from copy import deepcopy
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from XXEJ_scanner.cli import _deduplicate_bnd_events
from XXEJ_scanner.coverage import call_structural_evidence_regions
from XXEJ_scanner.models import (
    CigarIndel,
    DiscordantPair,
    EventEvidence,
    RepairEvent,
    ScannerConfig,
)


def _reference_deduplicate_bnd_events(events, evidence, window):
    """Previous exhaustive implementation used as an equivalence oracle."""
    def canonical_breakends(event):
        endpoints = [
            (event.bkp_A_chrom, int(event.bkp_A_pos)),
            (event.bkp_B_chrom, int(event.bkp_B_pos)),
        ]
        orientation = event.orientation
        if endpoints[1] < endpoints[0]:
            endpoints.reverse()
            if orientation != "NA" and len(orientation) == 2:
                orientation = orientation[::-1]
        return endpoints, orientation

    kept = []
    aliases = {}
    for event in events:
        if event.event_type not in {"BND_INTRA", "BND_INTER"}:
            kept.append(event)
            continue
        endpoints, orientation = canonical_breakends(event)
        matching = next(
            (
                candidate
                for candidate in kept
                if candidate.event_type == event.event_type
                and (
                    orientation == "NA"
                    or canonical_breakends(candidate)[1] == "NA"
                    or canonical_breakends(candidate)[1] == orientation
                )
                and all(
                    left[0] == right[0]
                    and abs(left[1] - right[1]) <= max(1, window)
                    for left, right in zip(
                        canonical_breakends(candidate)[0], endpoints
                    )
                )
            ),
            None,
        )
        if matching is None:
            kept.append(event)
            continue
        current_rank = (
            matching.filter == "PASS",
            matching.junction_resolved,
            matching.bkp_A_side != "pair_only",
            matching.alt_support,
        )
        new_rank = (
            event.filter == "PASS",
            event.junction_resolved,
            event.bkp_A_side != "pair_only",
            event.alt_support,
        )
        if new_rank > current_rank:
            kept[kept.index(matching)] = event
            aliases[matching.event_id] = event.event_id
            event.support_read_names.update(matching.support_read_names)
        else:
            aliases[event.event_id] = matching.event_id
            matching.support_read_names.update(event.support_read_names)

    for row in evidence:
        while row.event_id in aliases:
            row.event_id = aliases[row.event_id]
    kept_ids = {event.event_id for event in kept}
    return kept, [row for row in evidence if row.event_id in kept_ids]


def _bnd_event(index: int, *, unique: bool = False) -> RepairEvent:
    group = index if unique else (index * 17) % 40
    jitter_a = 0 if unique else (index * 29) % 241 - 120
    jitter_b = 0 if unique else (index * 43) % 241 - 120
    chrom_a, chrom_b = (("chr1", "chr1") if index % 3 == 0 else ("chr1", "chr2"))
    pos_a = group * 1000 + 1000 + jitter_a
    pos_b = group * 1500 + 1_000_000 + jitter_b
    orientation = ("+-", "-+", "NA")[index % 3]
    if index % 4 == 0:
        chrom_a, chrom_b = chrom_b, chrom_a
        pos_a, pos_b = pos_b, pos_a
        if orientation != "NA":
            orientation = orientation[::-1]
    return RepairEvent(
        f"event_{index}",
        "BND_INTRA" if chrom_a == chrom_b else "BND_INTER",
        chrom_a,
        pos_a,
        pos_a + 1,
        chrom_a,
        pos_a,
        "pair_only" if index % 2 else "left_clip",
        bkp_B_chrom=chrom_b,
        bkp_B_pos=pos_b,
        orientation=orientation,
        junction_resolved=index % 5 == 0,
        filter="PASS" if index % 7 == 0 else "LowSupport",
        support_read_names={f"read_{index}_{n}" for n in range(index % 5 + 1)},
    )


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

    def test_deduplicates_bnds_with_the_same_temporary_id(self) -> None:
        event = _bnd_event(0)
        duplicate = deepcopy(event)
        evidence = [EventEvidence(event.event_id, "read", "split_read", "chr1", 1)]

        events, evidence = _deduplicate_bnd_events(
            [event, duplicate], evidence, 100
        )

        self.assertEqual([item.event_id for item in events], [event.event_id])
        self.assertEqual([row.event_id for row in evidence], [event.event_id])

    def test_rejects_cyclic_event_aliases(self) -> None:
        low = _bnd_event(0)
        low.event_id = "A"
        low.support_read_names = {"r1"}
        middle = deepcopy(low)
        middle.event_id = "B"
        middle.support_read_names.add("r2")
        high = deepcopy(middle)
        high.event_id = "A"
        high.support_read_names.add("r3")

        with self.assertRaisesRegex(ValueError, "Cyclic event ID alias"):
            _deduplicate_bnd_events(
                [low, middle, high],
                [EventEvidence("A", "read", "split_read", "chr1", 1)],
                100,
            )

    def test_indexed_bnd_deduplication_matches_exhaustive_results(self) -> None:
        events = []
        for index in range(600):
            events.append(_bnd_event(index))
            if index % 53 == 0:
                events.append(
                    RepairEvent(
                        f"local_{index}",
                        "LOCAL_INS",
                        "chr1",
                        index,
                        index + 1,
                        "chr1",
                        index,
                        "junction",
                    )
                )
        evidence = [
            EventEvidence(
                event.event_id, f"evidence_{index}", "split_read", "chr1", index
            )
            for index, event in enumerate(events)
            if event.event_type.startswith("BND_")
        ]

        expected = _reference_deduplicate_bnd_events(
            deepcopy(events), deepcopy(evidence), 100
        )
        actual = _deduplicate_bnd_events(deepcopy(events), deepcopy(evidence), 100)

        self.assertEqual(actual, expected)

    def test_deduplicates_many_sparse_bnds_without_quadratic_scan(self) -> None:
        events = [_bnd_event(index, unique=True) for index in range(10_000)]

        kept, evidence = _deduplicate_bnd_events(events, [], 100)

        self.assertEqual(len(kept), len(events))
        self.assertEqual(evidence, [])


if __name__ == "__main__":
    unittest.main()
