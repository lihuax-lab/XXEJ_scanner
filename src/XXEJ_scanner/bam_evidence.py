"""BAM parsing and evidence extraction for XXEJ scanning."""

from __future__ import annotations

from collections import Counter
from math import ceil
from typing import Iterable, Iterator

import pysam

from .models import (
    CandidateRegion,
    CigarIndel,
    ClipSite,
    DiscordantPair,
    RegionEvidence,
    ScannerConfig,
    SplitReadEvidence,
)
from .utils import (
    CIGAR_DEL,
    CIGAR_HARD_CLIP,
    CIGAR_INS,
    CIGAR_SOFT_CLIP,
    QUERY_CONSUMING_OPS,
    REF_CONSUMING_OPS,
    aligned_query_length,
    cigar_to_string,
    clamp_start,
    pair_orientation,
    parse_cigar_string,
    reference_consumed_length,
)


def get_reference_end(read: object) -> int:
    # pysam usually provides reference_end, but tests and light-weight mocks may
    # only expose reference_start + CIGAR. Falling back here keeps coordinate
    # logic centralized.
    reference_end = getattr(read, "reference_end", None)
    if reference_end is not None:
        return int(reference_end)
    return int(getattr(read, "reference_start")) + reference_consumed_length(
        getattr(read, "cigartuples", None)
    )


def has_tag(read: object, tag: str) -> bool:
    if hasattr(read, "has_tag"):
        return bool(read.has_tag(tag))
    try:
        read.get_tag(tag)
        return True
    except Exception:
        return False


def get_tag(read: object, tag: str) -> str:
    return str(read.get_tag(tag))


def read_reference_name(read: object) -> str:
    name = getattr(read, "reference_name", None)
    if name is None:
        raise ValueError("Read object does not expose reference_name")
    return str(name)


def read_mate_reference_name(read: object) -> str:
    # Some SAM representations store "=" for "same as reference_name"; normalize
    # it so discordant-pair grouping never has to special-case that sentinel.
    name = getattr(read, "next_reference_name", None)
    if name is None or name == "=":
        name = getattr(read, "reference_name", None)
    return str(name)


def passes_read_filters(
    read: object,
    config: ScannerConfig,
    *,
    min_mapq: int | None = None,
) -> bool:
    # These are discovery filters. Supplementary reads are skipped here because
    # SA-tag parsing captures split evidence from the primary alignment without
    # double-counting the same molecule.
    if getattr(read, "is_unmapped", False):
        return False
    if getattr(read, "is_secondary", False):
        return False
    if getattr(read, "is_supplementary", False) and not config.include_supplementary:
        return False
    if getattr(read, "is_duplicate", False) and not config.allow_duplicates:
        return False
    if getattr(read, "cigartuples", None) is None:
        return False
    if int(getattr(read, "mapping_quality", 0)) < (
        config.min_mapq if min_mapq is None else min_mapq
    ):
        return False
    if (
        aligned_query_length(getattr(read, "cigartuples", None))
        < config.min_aligned_length
    ):
        return False
    return True


def passes_breakpoint_quality(
    read: object,
    query_pos: int | None,
    config: ScannerConfig,
    *,
    min_baseq: int | None = None,
) -> bool:
    """Check a query window centered on one resolved breakpoint boundary."""
    threshold = config.min_breakpoint_baseq if min_baseq is None else min_baseq
    if threshold <= 0:
        return True
    qualities = getattr(read, "query_qualities", None)
    radius = config.breakpoint_quality_window
    if (
        query_pos is None
        or qualities is None
        or query_pos < radius
        or query_pos + radius > len(qualities)
    ):
        return False
    window = qualities[query_pos - radius : query_pos + radius]
    required = ceil(len(window) * config.min_breakpoint_quality_fraction)
    return sum(quality >= threshold for quality in window) >= required


def clip_lengths(
    cigartuples: list[tuple[int, int]] | tuple[tuple[int, int], ...] | None,
) -> tuple[int, int, str, str]:
    # Only end clipping is breakpoint-informative for this first version. Internal
    # clipping is rare in standard CIGARs and would need separate interpretation.
    if not cigartuples:
        return 0, 0, "", ""
    first_op, first_len = cigartuples[0]
    last_op, last_len = cigartuples[-1]
    left = first_len if first_op in {CIGAR_SOFT_CLIP, CIGAR_HARD_CLIP} else 0
    right = last_len if last_op in {CIGAR_SOFT_CLIP, CIGAR_HARD_CLIP} else 0
    left_type = (
        "S"
        if first_op == CIGAR_SOFT_CLIP
        else ("H" if first_op == CIGAR_HARD_CLIP else "")
    )
    right_type = (
        "S"
        if last_op == CIGAR_SOFT_CLIP
        else ("H" if last_op == CIGAR_HARD_CLIP else "")
    )
    return left, right, left_type, right_type


def is_strongly_clipped(read: object, min_clip_length: int) -> bool:
    left, right, _, _ = clip_lengths(getattr(read, "cigartuples", None))
    return left >= min_clip_length or right >= min_clip_length


def extract_clip_sites_from_read(
    read: object,
    config: ScannerConfig,
    *,
    min_baseq: int | None = None,
) -> list[ClipSite]:
    if not passes_read_filters(read, config):
        return []

    cigartuples = getattr(read, "cigartuples", None)
    left_len, right_len, left_type, right_type = clip_lengths(cigartuples)
    query_sequence = getattr(read, "query_sequence", None) or ""
    chrom = read_reference_name(read)
    cigar = cigar_to_string(cigartuples)
    reference_start = int(getattr(read, "reference_start"))
    reference_end = get_reference_end(read)
    strand = "-" if getattr(read, "is_reverse", False) else "+"
    sites: list[ClipSite] = []

    if left_len >= config.min_clip_length and passes_breakpoint_quality(
        read,
        left_len if left_type == "S" else None,
        config,
        min_baseq=min_baseq,
    ):
        # Left soft clipping means the clipped bases occur before the aligned
        # portion of the query, so the candidate break lies at reference_start.
        sequence = (
            query_sequence[:left_len] if left_type == "S" and query_sequence else "NA"
        )
        sites.append(
            ClipSite(
                chrom=chrom,
                pos=reference_start,
                side="left_clip",
                clip_length=left_len,
                clip_sequence=sequence,
                read_name=str(getattr(read, "query_name", "")),
                strand=strand,
                mapq=int(getattr(read, "mapping_quality", 0)),
                cigar=cigar,
                is_reverse=bool(getattr(read, "is_reverse", False)),
                reference_start=reference_start,
                reference_end=reference_end,
                clip_type=left_type,
            )
        )

    if right_len >= config.min_clip_length and passes_breakpoint_quality(
        read,
        len(query_sequence) - right_len if right_type == "S" else None,
        config,
        min_baseq=min_baseq,
    ):
        # Right soft clipping points to the reference end of the alignment. This
        # is a 0-based half-open end coordinate from pysam.
        sequence = (
            query_sequence[-right_len:]
            if right_type == "S" and query_sequence
            else "NA"
        )
        sites.append(
            ClipSite(
                chrom=chrom,
                pos=reference_end,
                side="right_clip",
                clip_length=right_len,
                clip_sequence=sequence,
                read_name=str(getattr(read, "query_name", "")),
                strand=strand,
                mapq=int(getattr(read, "mapping_quality", 0)),
                cigar=cigar,
                is_reverse=bool(getattr(read, "is_reverse", False)),
                reference_start=reference_start,
                reference_end=reference_end,
                clip_type=right_type,
            )
        )

    return sites


def extract_cigar_indels_from_read(
    read: object,
    config: ScannerConfig,
    *,
    min_baseq: int | None = None,
) -> list[CigarIndel]:
    if not passes_read_filters(read, config):
        return []
    chrom = read_reference_name(read)
    ref_pos = int(getattr(read, "reference_start"))
    query_pos = 0
    query_sequence = getattr(read, "query_sequence", None) or ""
    cigar = cigar_to_string(getattr(read, "cigartuples", None))
    indels: list[CigarIndel] = []

    for op, length in getattr(read, "cigartuples", []) or []:
        # ref_pos and query_pos are advanced independently because insertions
        # consume query only, while deletions consume reference only.
        if op == CIGAR_INS:
            if length >= config.min_indel_length and all(
                passes_breakpoint_quality(
                    read, boundary, config, min_baseq=min_baseq
                )
                for boundary in (query_pos, query_pos + length)
            ):
                sequence = (
                    query_sequence[query_pos : query_pos + length]
                    if query_sequence
                    else "NA"
                )
                indels.append(
                    CigarIndel(
                        chrom=chrom,
                        start=ref_pos,
                        end=ref_pos,
                        operation="INS",
                        length=length,
                        sequence=sequence or "NA",
                        read_name=str(getattr(read, "query_name", "")),
                        mapq=int(getattr(read, "mapping_quality", 0)),
                        cigar=cigar,
                    )
                )
        elif op == CIGAR_DEL:
            if length >= config.min_indel_length and passes_breakpoint_quality(
                read, query_pos, config, min_baseq=min_baseq
            ):
                indels.append(
                    CigarIndel(
                        chrom=chrom,
                        start=ref_pos,
                        end=ref_pos + length,
                        operation="DEL",
                        length=length,
                        sequence="NA",
                        read_name=str(getattr(read, "query_name", "")),
                        mapq=int(getattr(read, "mapping_quality", 0)),
                        cigar=cigar,
                    )
                )
        if op in REF_CONSUMING_OPS:
            ref_pos += length
        if op in QUERY_CONSUMING_OPS:
            query_pos += length

    return indels


def is_discordant_pair(
    read: object, config: ScannerConfig, region: CandidateRegion | None = None
) -> tuple[bool, str]:
    if not passes_read_filters(read, config):
        return False, ""
    if not getattr(read, "is_paired", False):
        return False, ""
    if getattr(read, "mate_is_unmapped", False):
        return False, ""

    chrom = read_reference_name(read)
    mate_chrom = read_mate_reference_name(read)
    read_pos = int(getattr(read, "reference_start"))
    mate_pos = int(getattr(read, "next_reference_start", -1))
    is_diff_chrom = chrom != mate_chrom
    reasons: list[str] = []

    # Discordance is intentionally permissive and reason-coded. CUT&Tag repair
    # evidence may not follow a single library-size/orientation model perfectly.
    if is_diff_chrom:
        reasons.append("different_chrom")
    else:
        distance = abs(mate_pos - read_pos)
        if distance > config.discordant_min_distance:
            reasons.append("distant_mate")
        if abs(int(getattr(read, "template_length", 0))) > config.max_insert_size:
            reasons.append("large_insert")

    orientation = config.library_orientation.lower()
    if orientation != "any":
        same_strand = bool(getattr(read, "is_reverse", False)) == bool(
            getattr(read, "mate_is_reverse", False)
        )
        if orientation in {"fr", "rf"} and same_strand:
            reasons.append("same_strand_pair")
        elif orientation in {"ff", "rr"} and not same_strand:
            reasons.append("opposite_strand_pair")

    mate_outside_candidate_region = False
    if region is not None and chrom == mate_chrom:
        # A same-chromosome mate can still support wrong-end joining if it lands
        # outside the enriched local search interval plus a small merge buffer,
        # but that alone should not turn an otherwise proper pair into
        # discordant evidence.
        padded_start = region.start - config.merge_distance
        padded_end = region.end + config.merge_distance
        mate_outside_candidate_region = mate_pos < padded_start or mate_pos > padded_end

    if mate_outside_candidate_region and reasons:
        reasons.append("mate_outside_candidate_region")

    return bool(reasons), ",".join(reasons)


def extract_discordant_pair_from_read(
    read: object, config: ScannerConfig, region: CandidateRegion | None = None
) -> DiscordantPair | None:
    is_discordant, reason = is_discordant_pair(read, config, region)
    if not is_discordant:
        return None
    return DiscordantPair(
        read_name=str(getattr(read, "query_name", "")),
        chrom=read_reference_name(read),
        pos=int(getattr(read, "reference_start")),
        mate_chrom=read_mate_reference_name(read),
        mate_pos=int(getattr(read, "next_reference_start", -1)),
        orientation=pair_orientation(
            bool(getattr(read, "is_reverse", False)),
            bool(getattr(read, "mate_is_reverse", False)),
        ),
        mapq=int(getattr(read, "mapping_quality", 0)),
        is_reverse=bool(getattr(read, "is_reverse", False)),
        mate_is_reverse=bool(getattr(read, "mate_is_reverse", False)),
        cigar=cigar_to_string(getattr(read, "cigartuples", None)),
        reason=reason,
    )


def _alignment_geometry(
    reference_start: int,
    cigartuples: list[tuple[int, int]],
    strand: str,
) -> tuple[int, int, int, int]:
    """Return query start/end and reference start/end for one alignment segment."""
    left_clip, right_clip, _, _ = clip_lengths(cigartuples)
    query_span = sum(
        length
        for op, length in cigartuples
        if op in QUERY_CONSUMING_OPS or op == CIGAR_HARD_CLIP
    )
    if strand == "+":
        query_start, query_end = left_clip, query_span - right_clip
    else:
        query_start, query_end = right_clip, query_span - left_clip
    return (
        query_start,
        query_end,
        reference_start,
        reference_start + reference_consumed_length(cigartuples),
    )


def _breakpoint_for_query_edge(
    geometry: tuple[int, int, int, int], strand: str, edge: str
) -> tuple[int, str]:
    _query_start, _query_end, reference_start, reference_end = geometry
    uses_reference_start = (edge == "start") == (strand == "+")
    if uses_reference_start:
        return reference_start, "left_clip"
    return reference_end, "right_clip"


def extract_split_reads_from_sa_tag(
    read: object,
    config: ScannerConfig,
    *,
    min_mapq: int | None = None,
    min_baseq: int | None = None,
) -> list[SplitReadEvidence]:
    if not passes_read_filters(read, config):
        return []
    if not has_tag(read, "SA"):
        return []
    # SA tags encode supplementary alignments as
    # rname,pos,strand,CIGAR,mapQ,NM;... with 1-based positions.
    sa_tag = get_tag(read, "SA")
    chrom = read_reference_name(read)
    orientation_local = "-" if getattr(read, "is_reverse", False) else "+"
    cigar = cigar_to_string(getattr(read, "cigartuples", None))
    local_geometry = _alignment_geometry(
        int(getattr(read, "reference_start")),
        list(getattr(read, "cigartuples", None) or []),
        orientation_local,
    )
    left_clip, right_clip, left_type, right_type = clip_lengths(
        getattr(read, "cigartuples", None)
    )
    splits: list[SplitReadEvidence] = []

    for item in sa_tag.rstrip(";").split(";"):
        if not item:
            continue
        fields = item.split(",")
        if len(fields) < 6:
            continue
        remote_chrom, remote_pos_1, remote_strand, remote_cigar, remote_mapq, nm = (
            fields[:6]
        )
        try:
            remote_start = int(remote_pos_1) - 1
            remote_mapq_int = int(remote_mapq)
            remote_nm_int = int(nm)
            remote_cigartuples = parse_cigar_string(remote_cigar)
        except ValueError:
            continue
        threshold = config.min_mapq if min_mapq is None else min_mapq
        if remote_mapq_int < threshold or remote_nm_int > config.max_sa_nm:
            continue
        remote_geometry = _alignment_geometry(
            remote_start, remote_cigartuples, remote_strand
        )
        local_first = (
            local_geometry[0] + local_geometry[1]
            <= remote_geometry[0] + remote_geometry[1]
        )
        local_edge, remote_edge = (
            ("end", "start") if local_first else ("start", "end")
        )
        local_pos, side = _breakpoint_for_query_edge(
            local_geometry, orientation_local, local_edge
        )
        query_boundary = (
            left_clip
            if side == "left_clip" and left_type == "S"
            else (
                len(getattr(read, "query_sequence", None) or "") - right_clip
                if side == "right_clip" and right_type == "S"
                else None
            )
        )
        if not passes_breakpoint_quality(
            read, query_boundary, config, min_baseq=min_baseq
        ):
            continue
        remote_pos, _remote_side = _breakpoint_for_query_edge(
            remote_geometry, remote_strand, remote_edge
        )
        splits.append(
            SplitReadEvidence(
                read_name=str(getattr(read, "query_name", "")),
                chrom=chrom,
                pos=local_pos,
                side=side,
                remote_chrom=remote_chrom,
                remote_pos=remote_pos,
                remote_strand=remote_strand,
                remote_cigar=remote_cigar,
                remote_mapq=remote_mapq_int,
                remote_nm=remote_nm_int,
                orientation=orientation_local + remote_strand,
                mapq=int(getattr(read, "mapping_quality", 0)),
                cigar=cigar,
                sa_tag=sa_tag,
            )
        )
    return splits


def iter_bam_records(
    bam: pysam.AlignmentFile,
    chrom: str,
    start: int,
    end: int,
    config: ScannerConfig,
    *,
    include_supplementary: bool = False,
    min_mapq: int | None = None,
) -> Iterator[pysam.AlignedSegment]:
    for read in bam.fetch(chrom, clamp_start(start), end):
        if passes_read_filters(
            read,
            config,
            min_mapq=min_mapq,
        ):
            yield read


def collect_region_evidence(
    bam_path: str,
    region: CandidateRegion,
    config: ScannerConfig,
    *,
    padding: int | None = None,
    min_mapq: int | None = None,
    min_baseq: int | None = None,
) -> RegionEvidence:
    pad = config.scan_padding if padding is None else padding
    evidence = RegionEvidence(region=region)
    # Avoid counting the same read name multiple times when both mates or
    # multiple SA records point to the same remote locus.
    seen_discordant: set[str] = set()
    seen_split: set[tuple[object, ...]] = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in iter_bam_records(
            bam,
            region.chrom,
            clamp_start(region.start - pad),
            region.end + pad,
            config,
            min_mapq=min_mapq,
        ):
            evidence.clip_sites.extend(
                extract_clip_sites_from_read(read, config, min_baseq=min_baseq)
            )
            evidence.indels.extend(
                extract_cigar_indels_from_read(read, config, min_baseq=min_baseq)
            )
            discordant = extract_discordant_pair_from_read(read, config, region)
            if discordant and discordant.read_name not in seen_discordant:
                evidence.discordant_pairs.append(discordant)
                seen_discordant.add(discordant.read_name)
            for split in extract_split_reads_from_sa_tag(
                read, config, min_mapq=min_mapq, min_baseq=min_baseq
            ):
                key = (
                    split.read_name,
                    split.chrom,
                    split.pos,
                    split.remote_chrom,
                    split.remote_pos,
                    split.orientation,
                )
                if key not in seen_split:
                    evidence.split_reads.append(split)
                    seen_split.add(key)
    return evidence


def _count_read_name_depth(
    bam_path: str,
    chrom: str,
    start: int,
    end: int,
    config: ScannerConfig,
    *,
    min_mapq: int | None = None,
) -> int:
    """Count unique read names overlapping a genomic region."""
    read_names: set[str] = set()

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in iter_bam_records(bam, chrom, start, end, config, min_mapq=min_mapq):
            if get_reference_end(read) <= start or int(read.reference_start) >= end:
                continue

            read_names.add(str(read.query_name))

    return len(read_names)


def _count_mean_pileup_depth(
    bam_path: str,
    chrom: str,
    start: int,
    end: int,
    config: ScannerConfig,
    *,
    min_mapq: int | None = None,
) -> float:
    """Calculate mean pileup depth across a genomic region."""
    total_depth = 0
    region_length = end - start
    mapq_threshold = config.min_mapq if min_mapq is None else min_mapq

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for column in bam.pileup(
            chrom,
            start,
            end,
            truncate=True,
            stepper="all",
            min_mapping_quality=mapq_threshold,
        ):
            column_depth = 0

            for pileup_read in column.pileups:
                read = pileup_read.alignment

                if read.is_unmapped:
                    continue

                if read.is_secondary or read.is_supplementary:
                    continue

                if read.is_duplicate and not config.allow_duplicates:
                    continue

                if read.mapping_quality < mapq_threshold:
                    continue

                if pileup_read.is_del:
                    continue

                if pileup_read.is_refskip:
                    continue

                column_depth += 1

            total_depth += column_depth

    return total_depth / region_length


def count_depth(
    bam_path: str,
    chrom: str,
    start: int,
    end: int,
    config: ScannerConfig,
    *,
    min_mapq: int | None = None,
) -> float:
    """Count local depth in a genomic region.

    Parameters
    ----------
    bam_path
        Input BAM path.
    chrom
        Reference contig name.
    start
        Region start coordinate.
    end
        Region end coordinate.
    config
        Scanner configuration.
        - config.depth_count_method = "region": count unique read names overlapping the region.
        - config.depth_count_method = "pileup": calculate mean pileup depth across the region.
    min_mapq
        Optional mapping quality threshold. If None, config.min_mapq is used.

    Returns
    -------
    float
        Local depth value. The "region" method returns an integer-like
        float, while the "pileup" method returns mean base-level depth.
    """
    start = clamp_start(start)

    if end <= start:
        return 0.0

    method = config.depth_count_method
    if method == "region":
        return float(
            _count_read_name_depth(
                bam_path,
                chrom,
                start,
                end,
                config,
                min_mapq=min_mapq,
            )
        )
    if method == "pileup":
        return _count_mean_pileup_depth(
            bam_path,
            chrom,
            start,
            end,
            config,
            min_mapq=min_mapq,
        )

    raise ValueError(f"Unsupported depth counting method: {method!r}")


def count_spanning_reads(
    bam_path: str,
    chrom: str,
    start: int,
    end: int,
    config: ScannerConfig,
    *,
    min_mapq: int | None = None,
) -> int:
    # REF-like reads must bridge the queried interval without strong clipping and
    # without local pair discordance. This is an evidence fraction denominator,
    # not a diploid genotype model.
    read_names: set[str] = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in iter_bam_records(
            bam, chrom, clamp_start(start - 1), end + 1, config, min_mapq=min_mapq
        ):
            if is_strongly_clipped(read, config.min_clip_length):
                continue
            if getattr(read, "is_paired", False) and not getattr(
                read, "is_proper_pair", False
            ):
                continue
            if int(read.reference_start) < start and get_reference_end(read) > end:
                read_names.add(str(read.query_name))
    return len(read_names)


def clip_site_count_near(
    sites: Iterable[ClipSite], chrom: str, pos: int, window: int
) -> int:
    return sum(
        1 for site in sites if site.chrom == chrom and abs(site.pos - pos) <= window
    )


def dominant_sequence(sequences: Iterable[str]) -> str:
    usable = [seq for seq in sequences if seq and seq != "NA"]
    if not usable:
        return "NA"
    return Counter(usable).most_common(1)[0][0]
