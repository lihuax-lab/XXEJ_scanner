"""Reference FASTA helpers."""

from __future__ import annotations

from dataclasses import dataclass

import pysam


@dataclass
class ReferenceGenome:
    path: str

    def __post_init__(self) -> None:
        # Keep a single FASTA handle open across the scan; repeated open/close
        # calls are expensive when many candidate regions are processed.
        self._fasta = pysam.FastaFile(self.path)

    def close(self) -> None:
        self._fasta.close()

    def fetch(self, chrom: str, start: int, end: int) -> str:
        # Clamp intervals to contig bounds so callers can request small flanks
        # around breakpoints at chromosome edges without special handling.
        chrom_len = self._fasta.get_reference_length(chrom)
        start = max(0, min(start, chrom_len))
        end = max(0, min(end, chrom_len))
        if end <= start:
            return ""
        return self._fasta.fetch(chrom, start, end).upper()

    def __enter__(self) -> "ReferenceGenome":
        return self

    def __exit__(self, *_exc: object) -> None:
        self.close()


def fetch_reference_sequence(reference: ReferenceGenome, chrom: str, start: int, end: int) -> str:
    return reference.fetch(chrom, start, end)


def detect_microhomology(
    reference: ReferenceGenome,
    chrom: str,
    bkp_a: int,
    bkp_b: int,
    min_len: int = 1,
    max_len: int = 20,
) -> tuple[str, int]:
    """Return the longest short sequence shared immediately before both breakpoints."""
    if bkp_a > bkp_b:
        bkp_a, bkp_b = bkp_b, bkp_a
    for length in range(max_len, min_len - 1, -1):
        # Detect microhomology for a canonical L-m-X-m-R -> L-m-R deletion model.
        left = reference.fetch(chrom, bkp_a - length, bkp_a)
        right = reference.fetch(chrom, bkp_b, bkp_b + length)
        if len(left) == length and left == right:
            return left, length
    return "NA", 0


def check_clipped_sequence_against_reference(
    reference: ReferenceGenome,
    chrom: str,
    pos: int,
    clipped_sequence: str,
    window: int = 30,
) -> bool:
    if not clipped_sequence or clipped_sequence == "NA":
        return False
    # This helper is a coarse artifact/context check: an exact local match means
    # the clip may simply be nearby reference sequence rather than filler DNA.
    context = reference.fetch(chrom, pos - window, pos + window)
    return clipped_sequence.upper() in context
