from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from XXEJ_scanner.io import (
    RawEvidenceWriter,
    write_events_tsv,
    write_raw_clip_sites_tsv,
    write_raw_discordant_pairs_tsv,
    write_raw_split_reads_tsv,
)
from XXEJ_scanner.models import (
    CandidateRegion,
    ClipSite,
    DiscordantPair,
    RegionEvidence,
    RepairEvent,
    SplitReadEvidence,
)


class EventsTsvTest(unittest.TestCase):
    def test_microhomology_fields_are_appended_after_existing_event_fields(self) -> None:
        event = RepairEvent(
            event_id="XEJ_000001",
            event_type="LOCAL_DEL",
            chrom="chr1",
            start=10,
            end=20,
            bkp_A_chrom="chr1",
            bkp_A_pos=10,
            bkp_A_side="right_clip",
            bkp_B_chrom="chr1",
            bkp_B_pos=20,
            bkp_B_side="left_clip",
            notes="test",
            microhomology_left_end=12,
            microhomology_right_start=19,
            junction_evidence_support=1,
            junction_evidence_types={"cigar_del"},
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "events.tsv"
            write_events_tsv(str(path), [event])
            header = path.read_text().splitlines()[0].split("\t")

        notes_index = header.index("notes")
        self.assertEqual(header[notes_index + 1 :], [
            "microhomology_left_end",
            "microhomology_right_start",
            "microhomology_offset_a",
            "microhomology_offset_b",
            "microhomology_deletion_start",
            "microhomology_deletion_end",
            "microhomology_deletion_length",
            "microhomology_ambiguity_bases",
            "microhomology_equivalent_hits",
            "microhomology_low_complexity",
            "junction_evidence_support",
            "junction_evidence_types",
            "evidence_level",
            "junction_resolved",
        ])

    def test_streamed_raw_evidence_matches_batch_writers(self) -> None:
        evidence = RegionEvidence(
            region=CandidateRegion("chr1", 10, 20, "region-1"),
            clip_sites=[
                ClipSite(
                    "chr1", 10, "left_clip", 5, "AAAAA", "read-1", "+", 60,
                    "5S45M", False, 10, 55,
                )
            ],
            discordant_pairs=[
                DiscordantPair(
                    "read-2", "chr1", 12, "chr2", 30, "+-", 50, False,
                    True, "50M", "different_chrom",
                )
            ],
            split_reads=[
                SplitReadEvidence(
                    "read-3", "chr1", 15, "right_clip", "chr2", 40, "+",
                    "10S40M", 55, 1, "++", 60, "40M10S",
                    "chr2,41,+,10S40M,55,1;",
                )
            ],
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            expected = [root / f"expected-{index}.tsv" for index in range(3)]
            actual = [root / f"actual-{index}.tsv" for index in range(3)]
            write_raw_clip_sites_tsv(str(expected[0]), evidence.clip_sites)
            write_raw_discordant_pairs_tsv(
                str(expected[1]), evidence.discordant_pairs
            )
            write_raw_split_reads_tsv(str(expected[2]), evidence.split_reads)
            with RawEvidenceWriter(*(str(path) for path in actual)) as writer:
                writer.write(evidence)

            self.assertEqual(
                [path.read_bytes() for path in actual],
                [path.read_bytes() for path in expected],
            )


if __name__ == "__main__":
    unittest.main()
