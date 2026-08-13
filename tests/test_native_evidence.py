"""Native evidence extraction must remain identical to the Python reference."""

from dataclasses import asdict
from pathlib import Path
import sys
import tempfile
import unittest

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from XXEJ_scanner.bam_evidence import collect_region_evidence
from XXEJ_scanner.models import CandidateRegion, ScannerConfig
from XXEJ_scanner.native_evidence import collect_regions_native, native_available


@unittest.skipUnless(native_available(), "HTSlib extension was not built")
class NativeEvidenceTests(unittest.TestCase):
    def test_native_matches_python_for_all_evidence_types(self):
        with tempfile.TemporaryDirectory() as tmp:
            bam_path = str(Path(tmp) / "evidence.bam")
            header = pysam.AlignmentHeader.from_dict(
                {
                    "HD": {"VN": "1.6", "SO": "coordinate"},
                    "SQ": [
                        {"SN": "chr1", "LN": 10_000},
                        {"SN": "chr2", "LN": 10_000},
                    ],
                }
            )
            with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
                bam.write(
                    self._read(
                        header,
                        "clip-sa-pair",
                        100,
                        "40M10S",
                        flag=1 | 32,
                        mate_tid=1,
                        mate_pos=500,
                        sa="chr2,501,+,10S40M,60,2;",
                    )
                )
                bam.write(self._read(header, "insertion", 200, "20M5I25M"))
                bam.write(self._read(header, "deletion", 300, "20M5D25M"))
                bam.write(
                    self._read(header, "duplicate", 350, "10S40M", flag=1024)
                )
            pysam.index(bam_path)

            region = CandidateRegion("chr1", 50, 400, "region-1")
            config = ScannerConfig(
                treated_bam=bam_path,
                reference_fasta="unused.fa",
                output_dir=tmp,
                scan_padding=0,
            )

            expected = collect_region_evidence(bam_path, region, config)
            actual = collect_regions_native(bam_path, [region], config, padding=0)[0]

            self.assertEqual(asdict(actual), asdict(expected))

    @staticmethod
    def _read(
        header: pysam.AlignmentHeader,
        name: str,
        start: int,
        cigar: str,
        *,
        flag: int = 0,
        mate_tid: int = -1,
        mate_pos: int = -1,
        sa: str | None = None,
    ) -> pysam.AlignedSegment:
        read = pysam.AlignedSegment(header)
        read.query_name = name
        read.flag = flag
        read.reference_id = 0
        read.reference_start = start
        read.mapping_quality = 60
        read.cigarstring = cigar
        query_length = read.infer_query_length(always=False)
        assert query_length is not None
        read.query_sequence = "A" * query_length
        read.query_qualities = pysam.qualitystring_to_array("I" * query_length)
        read.next_reference_id = mate_tid
        read.next_reference_start = mate_pos
        if sa:
            read.set_tag("SA", sa)
        return read


if __name__ == "__main__":
    unittest.main()
