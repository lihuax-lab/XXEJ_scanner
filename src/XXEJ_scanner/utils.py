"""Small utility helpers used across the scanner."""

from __future__ import annotations

import math
import os
import sys
from pathlib import Path
from typing import Iterable, Sequence, TypeVar

import pysam

from .models import ScannerConfig

T = TypeVar("T")

CIGAR_MATCH = 0
CIGAR_INS = 1
CIGAR_DEL = 2
CIGAR_SKIP = 3
CIGAR_SOFT_CLIP = 4
CIGAR_HARD_CLIP = 5
CIGAR_PAD = 6
CIGAR_EQUAL = 7
CIGAR_DIFF = 8

REF_CONSUMING_OPS = {CIGAR_MATCH, CIGAR_DEL, CIGAR_SKIP, CIGAR_EQUAL, CIGAR_DIFF}
QUERY_CONSUMING_OPS = {CIGAR_MATCH, CIGAR_INS, CIGAR_SOFT_CLIP, CIGAR_EQUAL, CIGAR_DIFF}
ALIGNED_QUERY_OPS = {CIGAR_MATCH, CIGAR_EQUAL, CIGAR_DIFF}
# pysam uses the SAM numeric CIGAR op codes above. These op sets define how
# coordinate cursors move during manual CIGAR parsing.


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def safe_mkdir(path: str | os.PathLike[str]) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)


def clamp_start(start: int) -> int:
    return max(0, int(start))


def unique_preserve_order(values: Iterable[T]) -> list[T]:
    return list(dict.fromkeys(values))


def format_float(value: float | int | None, digits: int = 4) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return "NA"
        return f"{value:.{digits}f}"
    return str(value)


def percentile(values: Sequence[float], pct: float) -> float:
    # Linear interpolation keeps the threshold stable for small numbers of bins.
    if not values:
        return 0.0
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    rank = (pct / 100.0) * (len(ordered) - 1)
    low = int(math.floor(rank))
    high = int(math.ceil(rank))
    if low == high:
        return ordered[low]
    frac = rank - low
    return ordered[low] * (1.0 - frac) + ordered[high] * frac


def cigar_to_string(
    cigartuples: list[tuple[int, int]] | tuple[tuple[int, int], ...] | None,
) -> str:
    if not cigartuples:
        return "*"
    op_map = {
        CIGAR_MATCH: "M",
        CIGAR_INS: "I",
        CIGAR_DEL: "D",
        CIGAR_SKIP: "N",
        CIGAR_SOFT_CLIP: "S",
        CIGAR_HARD_CLIP: "H",
        CIGAR_PAD: "P",
        CIGAR_EQUAL: "=",
        CIGAR_DIFF: "X",
    }
    return "".join(f"{length}{op_map.get(op, '?')}" for op, length in cigartuples)


def parse_cigar_string(cigar: str) -> list[tuple[int, int]]:
    # Tiny parser used by tests and future callers that need pysam-style tuples
    # without constructing an AlignedSegment.
    op_map = {"M": 0, "I": 1, "D": 2, "N": 3, "S": 4, "H": 5, "P": 6, "=": 7, "X": 8}
    number = ""
    tuples: list[tuple[int, int]] = []
    for char in cigar:
        if char.isdigit():
            number += char
            continue
        if not number or char not in op_map:
            raise ValueError(f"Invalid CIGAR string: {cigar}")
        tuples.append((op_map[char], int(number)))
        number = ""
    if number:
        raise ValueError(f"Invalid CIGAR string: {cigar}")
    return tuples


def reference_consumed_length(cigartuples: Sequence[tuple[int, int]] | None) -> int:
    if not cigartuples:
        return 0
    return sum(length for op, length in cigartuples if op in REF_CONSUMING_OPS)


def aligned_query_length(cigartuples: Sequence[tuple[int, int]] | None) -> int:
    if not cigartuples:
        return 0
    return sum(length for op, length in cigartuples if op in ALIGNED_QUERY_OPS)


def pair_orientation(is_reverse: bool, mate_is_reverse: bool) -> str:
    return ("-" if is_reverse else "+") + ("-" if mate_is_reverse else "+")


def ensure_bam_index(path: str) -> None:
    # Fail early and clearly. pysam.fetch requires an index, and a missing index
    # otherwise appears later as a less helpful region-iteration error.
    candidate_indexes = [
        f"{path}.bai",
        f"{path}.csi",
        str(Path(path).with_suffix(".bai")),
        str(Path(path).with_suffix(".csi")),
    ]
    if not any(Path(index).exists() for index in candidate_indexes):
        raise FileNotFoundError(
            f"BAM index not found for {path}. Expected .bai or .csi next to the BAM."
        )
    with pysam.AlignmentFile(path, "rb") as bam:
        if not bam.has_index():
            raise FileNotFoundError(f"pysam could not load a BAM index for {path}.")


def ensure_fasta_index(path: str) -> None:
    if not Path(f"{path}.fai").exists():
        raise FileNotFoundError(
            f"FASTA index not found for {path}. Run samtools faidx first."
        )
    with pysam.FastaFile(path):
        return


def validate_inputs(config: ScannerConfig) -> None:
    # Centralized validation keeps CLI startup errors deterministic before any
    # output files are written.
    ensure_bam_index(config.treated_bam)
    if config.control_bam:
        ensure_bam_index(config.control_bam)
    ensure_fasta_index(config.reference_fasta)
    for bed in (config.candidate_bed, config.peak_bed):
        if bed and not Path(bed).exists():
            raise FileNotFoundError(f"BED file not found: {bed}")
