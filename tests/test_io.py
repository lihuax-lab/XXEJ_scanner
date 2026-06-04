from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from XXEJ_scanner.io import write_events_tsv
from XXEJ_scanner.models import RepairEvent


class EventsTsvTest(unittest.TestCase):
    def test_microhomology_fields_are_appended_after_existing_event_fields(self) -> None:
        event = RepairEvent(
            event_id="XEJ_000001",
            event_type="MMEJ_DEL",
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
        ])


if __name__ == "__main__":
    unittest.main()
