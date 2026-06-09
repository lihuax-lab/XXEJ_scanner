from __future__ import annotations

import unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from XXEJ_scanner.reference import detect_microhomology, find_microhomology


class FakeReference:
    def __init__(self, sequence: str) -> None:
        self.sequence = sequence

    def fetch(self, _chrom: str, start: int, end: int) -> str:
        start = max(0, min(start, len(self.sequence)))
        end = max(0, min(end, len(self.sequence)))
        if end <= start:
            return ""
        return self.sequence[start:end].upper()


def sequence_with_matches(matches: dict[int, str], length: int = 64) -> str:
    filler = "0123456789BDEFHIJKLMNOPQRSTUVWXYZ!#$%&()*+,-./:;<=>?@[]^_{|}~"
    sequence = list(filler[:length])
    for start, bases in matches.items():
        sequence[start : start + len(bases)] = bases
    return "".join(sequence)


class DetectMicrohomologyTest(unittest.TestCase):
    def test_detects_exact_coordinate_microhomology(self) -> None:
        reference = FakeReference(
            sequence_with_matches(
                {
                    7: "ACG",
                    20: "ACG",
                }
            )
        )

        sequence, length = detect_microhomology(
            reference,
            "chr1",
            bkp_a=10,
            bkp_b=20,
            min_len=1,
            max_len=5,
        )

        self.assertEqual(sequence, "ACG")
        self.assertEqual(length, 3)

    def test_detects_microhomology_near_shifted_breakpoints(self) -> None:
        reference = FakeReference(
            sequence_with_matches(
                {
                    8: "TTGA",
                    24: "TTGA",
                }
            )
        )

        exact_sequence, exact_length = detect_microhomology(
            reference,
            "chr1",
            bkp_a=10,
            bkp_b=25,
            min_len=4,
            max_len=4,
            search_window=0,
        )
        shifted_sequence, shifted_length = detect_microhomology(
            reference,
            "chr1",
            bkp_a=10,
            bkp_b=25,
            min_len=4,
            max_len=4,
            search_window=2,
        )

        self.assertEqual((exact_sequence, exact_length), ("NA", 0))
        self.assertEqual((shifted_sequence, shifted_length), ("TTGA", 4))

    def test_find_microhomology_returns_adjusted_coordinates(self) -> None:
        reference = FakeReference(
            sequence_with_matches(
                {
                    8: "TTGA",
                    24: "TTGA",
                }
            )
        )

        hit = find_microhomology(
            reference,
            "chr1",
            bkp_a=10,
            bkp_b=25,
            min_len=4,
            max_len=4,
            search_window=2,
        )

        self.assertEqual(hit.sequence, "TTGA")
        self.assertEqual(hit.length, 4)
        self.assertEqual(hit.left_end, 12)
        self.assertEqual(hit.right_start, 24)
        self.assertEqual(hit.offset_a, 2)
        self.assertEqual(hit.offset_b, -1)
        self.assertEqual(hit.deletion_start, 12)
        self.assertEqual(hit.deletion_end, 28)
        self.assertEqual(hit.deletion_length, 16)
        self.assertEqual(hit.ambiguity_bases, 4)

    def test_prefers_longest_nearby_microhomology(self) -> None:
        reference = FakeReference(
            sequence_with_matches(
                {
                    7: "ACG",
                    20: "ACG",
                    10: "GATTA",
                    25: "GATTA",
                }
            )
        )

        sequence, length = detect_microhomology(
            reference,
            "chr1",
            bkp_a=10,
            bkp_b=20,
            min_len=1,
            max_len=5,
            search_window=5,
        )

        self.assertEqual(sequence, "GATTA")
        self.assertEqual(length, 5)

    def test_counts_equivalent_hits_and_marks_poly_base_low_complexity(self) -> None:
        reference = FakeReference("A" * 40)

        hit = find_microhomology(
            reference,
            "chr1",
            bkp_a=10,
            bkp_b=20,
            min_len=2,
            max_len=2,
            search_window=1,
        )

        self.assertEqual(hit.sequence, "AA")
        self.assertEqual(hit.length, 2)
        self.assertGreater(hit.equivalent_hit_count, 1)
        self.assertTrue(hit.low_complexity)

    def test_marks_simple_repeat_low_complexity(self) -> None:
        reference = FakeReference(
            sequence_with_matches(
                {
                    6: "ATAT",
                    20: "ATAT",
                }
            )
        )

        hit = find_microhomology(
            reference,
            "chr1",
            bkp_a=10,
            bkp_b=20,
            min_len=4,
            max_len=4,
            search_window=0,
        )

        self.assertEqual(hit.sequence, "ATAT")
        self.assertTrue(hit.low_complexity)


if __name__ == "__main__":
    unittest.main()
