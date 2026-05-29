"""BAM parsing and evidence extraction for XXEJ scanning."""

from __future__ import annotations

from collections import Counter
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
    include_supplementary: bool = False,
    min_mapq: int | None = None,
) -> bool:
    # These are discovery filters. Supplementary reads are skipped here because
    # SA-tag parsing captures split evidence from the primary alignment without
    # double-counting the same molecule.
    if getattr(read, "is_unmapped", False):
        return False
    if getattr(read, "is_secondary", False):
        return False
    if getattr(read, "is_supplementary", False) and not include_supplementary:
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


def extract_clip_sites_from_read(read: object, config: ScannerConfig) -> list[ClipSite]:
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

    if left_len >= config.min_clip_length:
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

    if right_len >= config.min_clip_length:
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
    read: object, config: ScannerConfig
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
            if length >= config.min_indel_length:
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
            if length >= config.min_indel_length:
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

    if region is not None and chrom == mate_chrom:
        # A same-chromosome mate can still support wrong-end joining if it lands
        # outside the enriched local search interval plus a small merge buffer.
        padded_start = region.start - config.merge_distance
        padded_end = region.end + config.merge_distance
        if mate_pos < padded_start or mate_pos > padded_end:
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


def extract_split_reads_from_sa_tag(
    read: object, config: ScannerConfig
) -> list[SplitReadEvidence]:
    if not passes_read_filters(read, config):
        return []
    if not has_tag(read, "SA"):
        return []
    # SA tags encode supplementary alignments as
    # rname,pos,strand,CIGAR,mapQ,NM;... with 1-based positions.
    sa_tag = get_tag(read, "SA")
    chrom = read_reference_name(read)
    local_pos = int(getattr(read, "reference_start"))
    orientation_local = "-" if getattr(read, "is_reverse", False) else "+"
    cigar = cigar_to_string(getattr(read, "cigartuples", None))
    splits: list[SplitReadEvidence] = []

    left_len, right_len, _, _ = clip_lengths(getattr(read, "cigartuples", None))
    if left_len >= config.min_clip_length and right_len >= config.min_clip_length:
        side = "both"
    elif left_len >= config.min_clip_length:
        side = "left_clip"
    elif right_len >= config.min_clip_length:
        side = "right_clip"
    else:
        side = "SA"

    for item in sa_tag.rstrip(";").split(";"):
        if not item:
            continue
        fields = item.split(",")
        if len(fields) < 6:
            continue
        remote_chrom, remote_pos_1, remote_strand, remote_cigar, remote_mapq, _nm = (
            fields[:6]
        )
        try:
            # Convert SAM's 1-based SA position to the 0-based convention used
            # everywhere else in this package.
            remote_pos = int(remote_pos_1) - 1
            remote_mapq_int = int(remote_mapq)
        except ValueError:
            continue
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
            include_supplementary=include_supplementary,
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
) -> RegionEvidence:
    pad = config.scan_padding if padding is None else padding
    evidence = RegionEvidence(region=region)
    # Avoid counting the same read name multiple times when both mates or
    # multiple SA records point to the same remote locus.
    seen_discordant: set[str] = set()
    seen_split: set[tuple[str, str, int]] = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in iter_bam_records(
            bam,
            region.chrom,
            clamp_start(region.start - pad),
            region.end + pad,
            config,
            min_mapq=min_mapq,
        ):
            evidence.clip_sites.extend(extract_clip_sites_from_read(read, config))
            evidence.indels.extend(extract_cigar_indels_from_read(read, config))
            discordant = extract_discordant_pair_from_read(read, config, region)
            if discordant and discordant.read_name not in seen_discordant:
                evidence.discordant_pairs.append(discordant)
                seen_discordant.add(discordant.read_name)
            for split in extract_split_reads_from_sa_tag(read, config):
                key = (split.read_name, split.remote_chrom, split.remote_pos)
                if key not in seen_split:
                    evidence.split_reads.append(split)
                    seen_split.add(key)
    return evidence


def count_depth(
    bam_path: str,
    chrom: str,
    start: int,
    end: int,
    config: ScannerConfig,
    *,
    min_mapq: int | None = None,
) -> int:
    # Count unique read names rather than pileup bases. This is a simple local
    # molecule/read support measure and is less inflated by read length.
    read_names: set[str] = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in iter_bam_records(
            bam, chrom, clamp_start(start), end, config, min_mapq=min_mapq
        ):
            if get_reference_end(read) <= start or int(read.reference_start) >= end:
                continue
            read_names.add(str(read.query_name))
    return len(read_names)


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
