from __future__ import annotations

import runpy
import unittest
from pathlib import Path


EVALUATOR = runpy.run_path(
    str(Path(__file__).resolve().parents[1] / "scripts/evaluate_xxej_benchmark.py")
)


class BenchmarkEvaluatorTest(unittest.TestCase):
    def test_matches_bnd_with_reversed_reported_endpoints(self) -> None:
        truth = {
            "event_type": "BND_INTER",
            "chrom": "chr1",
            "bkp_A_pos": "100",
            "remote_chrom": "chr2",
            "remote_pos": "5000",
        }
        event = {
            "event_type": "BND_INTER",
            "chrom": "chr2",
            "bkp_A_chrom": "chr2",
            "bkp_A_pos": "5000",
            "bkp_B_chrom": "chr1",
            "bkp_B_pos": "100",
        }

        self.assertTrue(EVALUATOR["_is_compatible"](truth, event, 25))
        self.assertEqual(EVALUATOR["_match_distance"](truth, event, 25), 0)

    def test_insertion_match_requires_the_same_sequence(self) -> None:
        truth = {
            "event_type": "LOCAL_INS",
            "chrom": "chr1",
            "bkp_A_pos": "100",
            "inserted_sequence": "AA",
        }
        event = {
            "event_type": "LOCAL_INS",
            "chrom": "chr1",
            "bkp_A_pos": "100",
            "inserted_sequence": "TT",
        }

        self.assertFalse(EVALUATOR["_is_compatible"](truth, event, 25))


if __name__ == "__main__":
    unittest.main()
