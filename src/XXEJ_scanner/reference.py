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


@dataclass(frozen=True, slots=True)
class MicrohomologyHit:
    sequence: str
    length: int
    left_end: int
    right_start: int
    offset_a: int
    offset_b: int
    equivalent_hit_count: int = 0
    low_complexity: bool = False

    @property
    def found(self) -> bool:
        return self.length > 0 and self.sequence != "NA"

    @property
    def deletion_start(self) -> int:
        return self.left_end

    @property
    def deletion_end(self) -> int:
        return self.right_start + self.length

    @property
    def deletion_length(self) -> int:
        return self.deletion_end - self.deletion_start if self.found else 0

    @property
    def ambiguity_bases(self) -> int:
        return self.length if self.found else 0


def _is_low_complexity(sequence: str) -> bool:
    if not sequence or sequence == "NA":
        return False
    if len(set(sequence)) == 1:
        return True
    for unit_length in (1, 2):
        unit = sequence[:unit_length]
        repeated = (unit * ((len(sequence) // unit_length) + 1))[: len(sequence)]
        if repeated == sequence:
            return True
    return False


def find_microhomology(
    reference: ReferenceGenome,
    chrom: str,
    bkp_a: int,
    bkp_b: int,
    min_len: int = 1,
    max_len: int = 20,
    search_window: int = 5,
) -> MicrohomologyHit:
    """Return the best microhomology hit near both breakpoints."""
    if bkp_a > bkp_b:
        bkp_a, bkp_b = bkp_b, bkp_a

    search_window = max(0, search_window)
    offsets = sorted(
        range(-search_window, search_window + 1),
        key=lambda offset: (abs(offset), offset),
    )

    for length in range(max_len, min_len - 1, -1):
        hits: list[MicrohomologyHit] = []
        for offset_a in offsets:
            adjusted_a = bkp_a + offset_a
            left = reference.fetch(chrom, adjusted_a - length, adjusted_a)
            if len(left) != length:
                continue
            for offset_b in offsets:
                adjusted_b = bkp_b + offset_b
                if adjusted_a > adjusted_b:
                    continue
                # Detect microhomology for a canonical L-m-X-m-R -> L-m-R
                # deletion model, allowing small coordinate uncertainty.
                right = reference.fetch(chrom, adjusted_b, adjusted_b + length)
                if len(right) == length and left == right:
                    hits.append(
                        MicrohomologyHit(
                            sequence=left,
                            length=length,
                            left_end=adjusted_a,
                            right_start=adjusted_b,
                            offset_a=offset_a,
                            offset_b=offset_b,
                            low_complexity=_is_low_complexity(left),
                        )
                    )
        if hits:
            best = min(
                hits,
                key=lambda hit: (
                    abs(hit.offset_a) + abs(hit.offset_b),
                    max(abs(hit.offset_a), abs(hit.offset_b)),
                    abs(hit.offset_a),
                    abs(hit.offset_b),
                    hit.offset_a,
                    hit.offset_b,
                ),
            )
            return MicrohomologyHit(
                sequence=best.sequence,
                length=best.length,
                left_end=best.left_end,
                right_start=best.right_start,
                offset_a=best.offset_a,
                offset_b=best.offset_b,
                equivalent_hit_count=len(hits),
                low_complexity=best.low_complexity,
            )
    return MicrohomologyHit(
        sequence="NA",
        length=0,
        left_end=bkp_a,
        right_start=bkp_b,
        offset_a=0,
        offset_b=0,
        equivalent_hit_count=0,
        low_complexity=False,
    )


def detect_microhomology(
    reference: ReferenceGenome,
    chrom: str,
    bkp_a: int,
    bkp_b: int,
    min_len: int = 1,
    max_len: int = 20,
    search_window: int = 5,
) -> tuple[str, int]:
    """Return the longest short sequence shared near both breakpoints."""
    hit = find_microhomology(
        reference,
        chrom,
        bkp_a,
        bkp_b,
        min_len,
        max_len,
        search_window,
    )
    return hit.sequence, hit.length


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
