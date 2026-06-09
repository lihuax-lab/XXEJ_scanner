"""Rule-based repair event classification."""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from dataclasses import dataclass

from .bam_evidence import dominant_sequence
from .models import (
    BreakpointCluster,
    CandidateRegion,
    ClipSite,
    CigarIndel,
    DiscordantPair,
    EventEvidence,
    RegionEvidence,
    RepairEvent,
    ScannerConfig,
    SplitReadEvidence,
)
from .reference import MicrohomologyHit, ReferenceGenome, find_microhomology


def _near(pos_a: int, pos_b: int, window: int) -> bool:
    return abs(pos_a - pos_b) <= window


def _cluster_side(cluster: BreakpointCluster) -> str:
    return cluster.clip_side


@dataclass(slots=True)
class _MmejPairCandidate:
    left: BreakpointCluster
    right: BreakpointCluster
    microhomology: MicrohomologyHit
    deletion_indels: list[CigarIndel]
    split_reads: list[SplitReadEvidence]
    remap_clips: list[ClipSite]
    score: float
    junction_evidence_types: set[str]

    @property
    def indel_support(self) -> int:
        return len({indel.read_name for indel in self.deletion_indels})

    @property
    def junction_read_names(self) -> set[str]:
        reads = {indel.read_name for indel in self.deletion_indels}
        reads.update(split.read_name for split in self.split_reads)
        reads.update(site.read_name for site in self.remap_clips)
        return reads


def _microhomology_score_bonus(length: int, low_complexity: bool = False) -> float:
    if length <= 1:
        return 0.0
    if length <= 3:
        bonus = 0.5
    elif length <= 5:
        bonus = 1.0
    else:
        bonus = 2.0
    return min(bonus, 0.5) if low_complexity else bonus


def _event_score(event: RepairEvent) -> float:
    # First-version scoring is deliberately transparent. It is a ranking aid for
    # review/IGV triage, not a calibrated probability or genotype likelihood.
    microhomology_bonus = _microhomology_score_bonus(
        event.microhomology_length,
        event.microhomology_low_complexity,
    )
    junction_bonus = math.log2(1 + event.junction_evidence_support)
    control_noise_penalty = event.normal_noise
    mapping_quality_penalty = 0.0
    return (
        1.0 * math.log2(1 + max(event.treated_depth, 0))
        + 2.0 * math.log2(1 + event.alt_clip_support)
        + 2.5 * math.log2(1 + event.alt_split_support)
        + 2.5 * math.log2(1 + event.alt_discordant_pair_support)
        + 1.5 * math.log2(1 + event.alt_indel_support)
        + 1.0 * microhomology_bonus
        + 1.0 * junction_bonus
        - 2.0 * control_noise_penalty
        - mapping_quality_penalty
    )


def classify_local_events(
    region: CandidateRegion,
    clusters: list[BreakpointCluster],
    evidence: RegionEvidence,
    reference: ReferenceGenome,
    config: ScannerConfig,
) -> tuple[list[RepairEvent], list[EventEvidence]]:
    # Priority order matters. A compatible two-breakpoint MMEJ deletion consumes
    # its clusters first; remaining clusters are tested for remote BND evidence,
    # then for local insertion/filler evidence.
    events: list[RepairEvent] = []
    event_evidence: list[EventEvidence] = []
    sorted_clusters = sorted(clusters, key=lambda cluster: cluster.peak_pos)

    mmej_events, mmej_evidence, used_cluster_keys = _classify_mmej_del(
        region, sorted_clusters, evidence, reference, config
    )
    events.extend(mmej_events)
    event_evidence.extend(mmej_evidence)

    for cluster in sorted_clusters:
        key = (cluster.chrom, cluster.peak_pos, cluster.clip_side)
        if key in used_cluster_keys:
            continue
        bnd_events, bnd_evidence = classify_bnd_events(
            region, [cluster], evidence, config
        )
        if bnd_events:
            events.extend(bnd_events)
            event_evidence.extend(bnd_evidence)
            continue

        insertion_event, insertion_evidence = _classify_nhej_ins(
            region, cluster, evidence, config
        )
        if insertion_event:
            events.append(insertion_event)
            event_evidence.extend(insertion_evidence)

    return events, event_evidence


def assign_final_event_ids(
    events: list[RepairEvent],
    evidence: list[EventEvidence],
) -> tuple[list[RepairEvent], list[EventEvidence]]:
    # Classifiers emit temporary IDs so their read-level evidence can be joined
    # locally. Final IDs must be assigned once globally after all regions have
    # been scanned, otherwise each region starts again at XEJ_000001.
    temporary_counts = Counter(event.event_id for event in events)
    duplicate_temporary_ids = sorted(
        event_id for event_id, count in temporary_counts.items() if count > 1
    )
    if duplicate_temporary_ids:
        duplicates = ", ".join(duplicate_temporary_ids[:5])
        raise ValueError(
            "Temporary event IDs are not unique; cannot safely assign final IDs. "
            f"Duplicate IDs include: {duplicates}"
        )

    id_map: dict[str, str] = {}
    for idx, event in enumerate(events, 1):
        old_event_id = event.event_id
        new_event_id = f"XEJ_{idx:06d}"
        id_map[old_event_id] = new_event_id
        event.event_id = new_event_id

    for ev in evidence:
        if ev.event_id in id_map:
            ev.event_id = id_map[ev.event_id]

    return events, evidence


def classify_bnd_events(
    region: CandidateRegion,
    clusters: list[BreakpointCluster],
    evidence: RegionEvidence,
    config: ScannerConfig,
) -> tuple[list[RepairEvent], list[EventEvidence]]:
    events: list[RepairEvent] = []
    event_evidence: list[EventEvidence] = []
    for cluster in clusters:
        linked_discordants = [
            pair
            for pair in evidence.discordant_pairs
            if pair.chrom == cluster.chrom
            and _near(pair.pos, cluster.peak_pos, config.clip_cluster_window * 2)
        ]
        linked_splits = [
            split
            for split in evidence.split_reads
            if split.chrom == cluster.chrom
            and _near(split.pos, cluster.peak_pos, config.clip_cluster_window * 3)
        ]
        # Remote support is grouped into coarse bins. This avoids splitting a
        # real remote locus into many tiny calls due to mate-position jitter.
        grouped: dict[tuple[str, int], dict[str, list[object]]] = defaultdict(
            lambda: {"pairs": [], "splits": []}
        )
        for pair in linked_discordants:
            remote_bin = pair.mate_pos - (
                pair.mate_pos % max(1, config.coverage_bin_size)
            )
            grouped[(pair.mate_chrom, remote_bin)]["pairs"].append(pair)
        for split in linked_splits:
            remote_bin = split.remote_pos - (
                split.remote_pos % max(1, config.coverage_bin_size)
            )
            grouped[(split.remote_chrom, remote_bin)]["splits"].append(split)

        for (remote_chrom, _remote_bin), group in grouped.items():
            pairs = group["pairs"]
            splits = group["splits"]
            support = len({item.read_name for item in pairs + splits})
            if support < config.min_bnd_support:
                continue
            # Report the average remote coordinate as an imprecise BND anchor.
            # This remains a candidate junction until validated downstream.
            remote_positions = [item.mate_pos for item in pairs] + [
                item.remote_pos for item in splits
            ]
            remote_pos = (
                int(round(sum(remote_positions) / len(remote_positions)))
                if remote_positions
                else "NA"
            )
            orientation_values = [
                item.orientation
                for item in pairs + splits
                if getattr(item, "orientation", "NA") != "NA"
            ]
            orientation = (
                Counter(orientation_values).most_common(1)[0][0]
                if orientation_values
                else "NA"
            )
            event_type = (
                "NHEJ_BND_INS_INTER"
                if remote_chrom != cluster.chrom
                else "NHEJ_BND_INS_INTRA"
            )
            local_clip_reads = {
                site.read_name
                for site in evidence.clip_sites
                if site.chrom == cluster.chrom
                and _near(site.pos, cluster.peak_pos, config.clip_cluster_window)
            }
            temp_event_id = (
                f"TMP_BND_{region.region_id}_{cluster.chrom}_{cluster.peak_pos}_"
                f"{remote_chrom}_{remote_pos}"
            )
            event = RepairEvent(
                event_id=temp_event_id,
                event_type=event_type,
                chrom=cluster.chrom,
                start=max(region.start, cluster.peak_pos - 1),
                end=min(region.end, cluster.peak_pos + 1),
                bkp_A_chrom=cluster.chrom,
                bkp_A_pos=cluster.peak_pos,
                bkp_A_side=_cluster_side(cluster),
                bkp_B_chrom=remote_chrom,
                bkp_B_pos=remote_pos,
                bkp_B_side="remote",
                remote_chrom=remote_chrom,
                remote_pos=remote_pos,
                orientation=orientation,
                alt_clip_support=cluster.clip_count,
                alt_split_support=len({split.read_name for split in splits}),
                alt_discordant_pair_support=len({pair.read_name for pair in pairs}),
                treated_depth=cluster.treated_depth,
                control_depth=cluster.control_depth,
                normal_noise=cluster.normal_noise,
                notes="Candidate wrong-end joining supported by clipped and remote evidence.",
                support_read_names=local_clip_reads
                | {item.read_name for item in pairs + splits},
            )
            event.score = _event_score(event)
            events.append(event)
            for site in evidence.clip_sites:
                if site.chrom == cluster.chrom and _near(
                    site.pos, cluster.peak_pos, config.clip_cluster_window
                ):
                    event_evidence.append(_clip_evidence(temp_event_id, site))
            for pair in pairs:
                event_evidence.append(_pair_evidence(temp_event_id, pair))
            for split in splits:
                event_evidence.append(_split_evidence(temp_event_id, split))
    return events, event_evidence


def _classify_mmej_del(
    region: CandidateRegion,
    clusters: list[BreakpointCluster],
    evidence: RegionEvidence,
    reference: ReferenceGenome,
    config: ScannerConfig,
) -> tuple[list[RepairEvent], list[EventEvidence], set[tuple[str, int, str]]]:
    events: list[RepairEvent] = []
    event_evidence: list[EventEvidence] = []
    used: set[tuple[str, int, str]] = set()
    if len(clusters) < 2:
        return events, event_evidence, used

    # Score deletion-compatible pairs after normalizing genomic left/right order.
    # Microhomology remains contextual evidence, but it now helps select the most
    # MMEJ-like pair when several local breakpoint hypotheses overlap.
    best_candidate: _MmejPairCandidate | None = None
    ordered_clusters = sorted(
        clusters,
        key=lambda cluster: (cluster.chrom, cluster.peak_pos),
    )
    for idx, cluster_a in enumerate(ordered_clusters):
        for cluster_b in ordered_clusters[idx + 1 :]:
            left, right = sorted(
                (cluster_a, cluster_b),
                key=lambda cluster: cluster.peak_pos,
            )
            if left.chrom != right.chrom:
                continue
            if not _is_deletion_compatible_pair(left, right):
                continue
            distance = abs(right.peak_pos - left.peak_pos)
            if distance <= 0 or distance > config.max_local_event_distance:
                continue

            mh_hit = find_microhomology(
                reference,
                left.chrom,
                left.peak_pos,
                right.peak_pos,
                config.min_microhomology_length,
                config.max_microhomology_length,
                config.microhomology_search_window,
            )
            deletion_start, deletion_end = _expected_deletion_span(
                left,
                right,
                mh_hit,
            )
            deletion_indels = _indels_matching_deletion_span(
                evidence.indels,
                left.chrom,
                deletion_start,
                deletion_end,
                config,
            )
            split_reads = _split_reads_matching_deletion_span(
                evidence.split_reads,
                left.chrom,
                deletion_start,
                deletion_end,
                config,
            )
            remap_clips = _soft_clips_matching_opposite_flanks(
                evidence.clip_sites,
                reference,
                left,
                right,
                deletion_start,
                deletion_end,
                config,
            )
            junction_types = _junction_evidence_types(
                deletion_indels,
                split_reads,
                remap_clips,
            )
            candidate = _MmejPairCandidate(
                left=left,
                right=right,
                microhomology=mh_hit,
                deletion_indels=deletion_indels,
                split_reads=split_reads,
                remap_clips=remap_clips,
                score=_mmej_pair_score(
                    left,
                    right,
                    mh_hit,
                    deletion_indels,
                    split_reads,
                    remap_clips,
                ),
                junction_evidence_types=junction_types,
            )
            if best_candidate is None or _mmej_candidate_key(
                candidate
            ) > _mmej_candidate_key(best_candidate):
                best_candidate = candidate

    if best_candidate is None:
        return events, event_evidence, used

    left = best_candidate.left
    right = best_candidate.right
    mh_hit = best_candidate.microhomology
    local_clips = [
        site
        for site in evidence.clip_sites
        if site.chrom == left.chrom
        and (
            _near(site.pos, left.peak_pos, config.clip_cluster_window)
            or _near(site.pos, right.peak_pos, config.clip_cluster_window)
        )
    ]
    support_reads = {
        site.read_name
        for site in local_clips
    }
    support_reads.update(indel.read_name for indel in best_candidate.deletion_indels)
    support_reads.update(split.read_name for split in best_candidate.split_reads)
    support_reads.update(site.read_name for site in best_candidate.remap_clips)
    if len(support_reads) < config.min_alt_support:
        return events, event_evidence, used

    event = RepairEvent(
        event_id=f"TMP_MMEJ_{region.region_id}_{left.chrom}_{left.peak_pos}_{right.peak_pos}",
        event_type="MMEJ_DEL",
        chrom=left.chrom,
        start=min(left.peak_pos, right.peak_pos),
        end=max(left.peak_pos, right.peak_pos),
        bkp_A_chrom=left.chrom,
        bkp_A_pos=left.peak_pos,
        bkp_A_side=_cluster_side(left),
        bkp_B_chrom=right.chrom,
        bkp_B_pos=right.peak_pos,
        bkp_B_side=_cluster_side(right),
        deleted_length=abs(right.peak_pos - left.peak_pos),
        microhomology=mh_hit.sequence,
        microhomology_length=mh_hit.length,
        alt_clip_support=left.clip_count + right.clip_count,
        alt_split_support=len({split.read_name for split in best_candidate.split_reads}),
        alt_indel_support=best_candidate.indel_support,
        treated_depth=max(left.treated_depth, right.treated_depth),
        control_depth=max(left.control_depth, right.control_depth),
        normal_noise=max(left.normal_noise, right.normal_noise),
        notes=_mmej_notes(mh_hit, best_candidate.junction_evidence_types),
        microhomology_left_end=mh_hit.left_end if mh_hit.found else "NA",
        microhomology_right_start=mh_hit.right_start if mh_hit.found else "NA",
        microhomology_offset_a=mh_hit.offset_a if mh_hit.found else "NA",
        microhomology_offset_b=mh_hit.offset_b if mh_hit.found else "NA",
        microhomology_deletion_start=mh_hit.deletion_start if mh_hit.found else "NA",
        microhomology_deletion_end=mh_hit.deletion_end if mh_hit.found else "NA",
        microhomology_deletion_length=mh_hit.deletion_length if mh_hit.found else "NA",
        microhomology_ambiguity_bases=mh_hit.ambiguity_bases,
        microhomology_equivalent_hits=mh_hit.equivalent_hit_count,
        microhomology_low_complexity=mh_hit.low_complexity,
        junction_evidence_support=len(best_candidate.junction_read_names),
        junction_evidence_types=best_candidate.junction_evidence_types,
        support_read_names=support_reads,
    )
    event.score = _event_score(event)
    events.append(event)
    for site in local_clips:
        event_evidence.append(_clip_evidence(event.event_id, site))
    for indel in best_candidate.deletion_indels:
        event_evidence.append(_indel_evidence(event.event_id, indel))
    for split in best_candidate.split_reads:
        event_evidence.append(_split_evidence(event.event_id, split))
    used.add((left.chrom, left.peak_pos, left.clip_side))
    used.add((right.chrom, right.peak_pos, right.clip_side))
    return events, event_evidence, used


def _is_deletion_compatible_pair(
    left: BreakpointCluster,
    right: BreakpointCluster,
) -> bool:
    return left.clip_side in {"right_clip", "both"} and right.clip_side in {
        "left_clip",
        "both",
    }


def _orientation_bonus(left: BreakpointCluster, right: BreakpointCluster) -> float:
    if left.clip_side == "right_clip" and right.clip_side == "left_clip":
        return 2.0
    return 1.0


def _expected_deletion_span(
    left: BreakpointCluster,
    right: BreakpointCluster,
    hit: MicrohomologyHit,
) -> tuple[int, int]:
    if hit.found:
        return hit.deletion_start, hit.deletion_end
    return left.peak_pos, right.peak_pos


def _mmej_pair_score(
    left: BreakpointCluster,
    right: BreakpointCluster,
    hit: MicrohomologyHit,
    deletion_indels: list[CigarIndel],
    split_reads: list[SplitReadEvidence],
    remap_clips: list[ClipSite],
) -> float:
    clip_support = left.clip_count + right.clip_count
    indel_support = len({indel.read_name for indel in deletion_indels})
    junction_reads = {indel.read_name for indel in deletion_indels}
    junction_reads.update(split.read_name for split in split_reads)
    junction_reads.update(site.read_name for site in remap_clips)
    return (
        clip_support
        + 1.5 * indel_support
        + 2.0 * len(junction_reads)
        + _microhomology_score_bonus(hit.length, hit.low_complexity)
        + _orientation_bonus(left, right)
        - 2.0 * max(left.normal_noise, right.normal_noise)
    )


def _mmej_candidate_key(
    candidate: _MmejPairCandidate,
) -> tuple[float, int, int, int, int]:
    return (
        candidate.score,
        len(candidate.junction_read_names),
        candidate.microhomology.length,
        candidate.left.clip_count + candidate.right.clip_count,
        -abs(candidate.right.peak_pos - candidate.left.peak_pos),
    )


def _junction_evidence_types(
    deletion_indels: list[CigarIndel],
    split_reads: list[SplitReadEvidence],
    remap_clips: list[ClipSite],
) -> set[str]:
    evidence_types: set[str] = set()
    if deletion_indels:
        evidence_types.add("cigar_del")
    if split_reads:
        evidence_types.add("split_read_sa")
    if remap_clips:
        evidence_types.add("soft_clip_remap")
    return evidence_types


def _mmej_notes(hit: MicrohomologyHit, junction_evidence_types: set[str]) -> str:
    notes = ["Candidate local deletion-like event from paired breakpoint clusters."]
    if hit.found:
        notes.append("Reference microhomology context detected.")
    else:
        notes.append("No reference microhomology detected near the nominal breakpoints.")
    if junction_evidence_types:
        evidence = ",".join(sorted(junction_evidence_types))
        notes.append(f"Junction-level evidence types: {evidence}.")
    else:
        notes.append("No junction-level read evidence detected.")
    if hit.equivalent_hit_count > 1:
        notes.append(f"{hit.equivalent_hit_count} equivalent microhomology placements.")
    if hit.low_complexity:
        notes.append("Microhomology sequence is low-complexity.")
    return " ".join(notes)


def _classify_nhej_ins(
    region: CandidateRegion,
    cluster: BreakpointCluster,
    evidence: RegionEvidence,
    config: ScannerConfig,
) -> tuple[RepairEvent | None, list[EventEvidence]]:
    # NHEJ_INS is intentionally local: it uses nearby short CIGAR insertions and
    # optional clipped filler sequence, but avoids turning every unsupported clip
    # cluster into a local insertion.
    local_insertions = [
        indel
        for indel in evidence.indels
        if indel.operation == "INS"
        and indel.length <= config.max_insertion_length
        and indel.chrom == cluster.chrom
        and _near(indel.start, cluster.peak_pos, config.clip_cluster_window * 2)
    ]
    local_clips = [
        site
        for site in evidence.clip_sites
        if site.chrom == cluster.chrom
        and _near(site.pos, cluster.peak_pos, config.clip_cluster_window)
    ]
    short_clip_proxies = [
        site
        for site in local_clips
        if site.clip_length <= config.max_insertion_length
        and site.clip_sequence
        and site.clip_sequence != "NA"
    ]
    insertion_reads = {indel.read_name for indel in local_insertions}

    if len(insertion_reads) < config.min_nhej_ins_indel_support:
        if not config.allow_clip_only_nhej_ins:
            return None, []
        support_reads = {site.read_name for site in short_clip_proxies}
    else:
        support_reads = {site.read_name for site in local_clips} | insertion_reads

    if len(support_reads) < config.min_alt_support:
        return None, []

    inserted_sequence = dominant_sequence(indel.sequence for indel in local_insertions)
    if inserted_sequence == "NA":
        # If explicit insertion sequence is unavailable, a short terminal soft
        # clip can still act as a filler-sequence proxy. Long clips are not used
        # for NHEJ_INS length inference because they often mark unresolved split
        # or mapping evidence.
        inserted_sequence = dominant_sequence(
            site.clip_sequence for site in short_clip_proxies
        )
    if inserted_sequence == "NA":
        return None, []

    inserted_length: int | str = len(inserted_sequence)
    notes = (
        "Candidate local NHEJ-like insertion with local CIGAR insertion evidence."
        if insertion_reads
        else "Candidate local NHEJ-like insertion from clipped filler-sequence evidence only."
    )
    event = RepairEvent(
        event_id=f"TMP_INS_{region.region_id}_{cluster.chrom}_{cluster.peak_pos}_{cluster.clip_side}",
        event_type="NHEJ_INS",
        chrom=cluster.chrom,
        start=max(region.start, cluster.peak_pos - 1),
        end=min(region.end, cluster.peak_pos + 1),
        bkp_A_chrom=cluster.chrom,
        bkp_A_pos=cluster.peak_pos,
        bkp_A_side=_cluster_side(cluster),
        inserted_sequence=inserted_sequence,
        inserted_length=inserted_length,
        alt_clip_support=cluster.clip_count,
        alt_indel_support=len({indel.read_name for indel in local_insertions}),
        treated_depth=cluster.treated_depth,
        control_depth=cluster.control_depth,
        normal_noise=cluster.normal_noise,
        notes=notes,
        support_read_names=support_reads,
    )
    event.score = _event_score(event)
    ev_rows = [_clip_evidence(event.event_id, site) for site in local_clips]
    ev_rows.extend(_indel_evidence(event.event_id, indel) for indel in local_insertions)
    return event, ev_rows


def _indels_matching_deletion_span(
    indels: list[CigarIndel],
    chrom: str,
    deletion_start: int,
    deletion_end: int,
    config: ScannerConfig,
) -> list[CigarIndel]:
    return [
        indel
        for indel in indels
        if indel.chrom == chrom
        and indel.operation == "DEL"
        and _near(indel.start, deletion_start, config.clip_cluster_window)
        and _near(indel.end, deletion_end, config.clip_cluster_window)
    ]


def _split_reads_matching_deletion_span(
    splits: list[SplitReadEvidence],
    chrom: str,
    deletion_start: int,
    deletion_end: int,
    config: ScannerConfig,
) -> list[SplitReadEvidence]:
    window = config.clip_cluster_window * 2
    return [
        split
        for split in splits
        if split.chrom == chrom
        and split.remote_chrom == chrom
        and (
            (
                _near(split.pos, deletion_start, window)
                and _near(split.remote_pos, deletion_end, window)
            )
            or (
                _near(split.pos, deletion_end, window)
                and _near(split.remote_pos, deletion_start, window)
            )
        )
    ]


def _soft_clips_matching_opposite_flanks(
    sites: list[ClipSite],
    reference: ReferenceGenome,
    left: BreakpointCluster,
    right: BreakpointCluster,
    deletion_start: int,
    deletion_end: int,
    config: ScannerConfig,
) -> list[ClipSite]:
    matches: list[ClipSite] = []
    for site in sites:
        if site.chrom != left.chrom:
            continue
        if (
            site.side == "right_clip"
            and _near(site.pos, left.peak_pos, config.clip_cluster_window)
            and _right_clip_matches_flank(
                reference,
                site.chrom,
                deletion_end,
                site.clip_sequence,
                config,
            )
        ):
            matches.append(site)
            continue
        if (
            site.side == "left_clip"
            and _near(site.pos, right.peak_pos, config.clip_cluster_window)
            and _left_clip_matches_flank(
                reference,
                site.chrom,
                deletion_start,
                site.clip_sequence,
                config,
            )
        ):
            matches.append(site)
    return matches


def _right_clip_matches_flank(
    reference: ReferenceGenome,
    chrom: str,
    flank_start: int,
    clipped_sequence: str,
    config: ScannerConfig,
) -> bool:
    sequence = clipped_sequence.upper()
    if not sequence or sequence == "NA":
        return False
    match_len = min(len(sequence), config.clip_cluster_window)
    if match_len < min(config.min_clip_length, len(sequence)):
        return False
    flank = reference.fetch(chrom, flank_start, flank_start + match_len)
    return len(flank) == match_len and sequence[:match_len] == flank


def _left_clip_matches_flank(
    reference: ReferenceGenome,
    chrom: str,
    flank_end: int,
    clipped_sequence: str,
    config: ScannerConfig,
) -> bool:
    sequence = clipped_sequence.upper()
    if not sequence or sequence == "NA":
        return False
    match_len = min(len(sequence), config.clip_cluster_window)
    if match_len < min(config.min_clip_length, len(sequence)):
        return False
    flank = reference.fetch(chrom, flank_end - match_len, flank_end)
    return len(flank) == match_len and sequence[-match_len:] == flank


def _clip_evidence(event_id: str, site: object) -> EventEvidence:
    return EventEvidence(
        event_id=event_id,
        read_name=site.read_name,
        evidence_type="soft_clip",
        chrom=site.chrom,
        pos=site.pos,
        cigar=site.cigar,
        mapq=site.mapq,
        is_reverse=site.is_reverse,
        clip_side=site.side,
        clip_length=site.clip_length,
        clip_sequence=site.clip_sequence,
    )


def _indel_evidence(event_id: str, indel: CigarIndel) -> EventEvidence:
    return EventEvidence(
        event_id=event_id,
        read_name=indel.read_name,
        evidence_type=f"cigar_{indel.operation.lower()}",
        chrom=indel.chrom,
        pos=indel.start,
        cigar=indel.cigar,
        mapq=indel.mapq,
        clip_sequence=indel.sequence,
    )


def _pair_evidence(event_id: str, pair: DiscordantPair) -> EventEvidence:
    return EventEvidence(
        event_id=event_id,
        read_name=pair.read_name,
        evidence_type="discordant_pair",
        chrom=pair.chrom,
        pos=pair.pos,
        mate_chrom=pair.mate_chrom,
        mate_pos=pair.mate_pos,
        cigar=pair.cigar,
        mapq=pair.mapq,
        is_reverse=pair.is_reverse,
        mate_is_reverse=pair.mate_is_reverse,
    )


def _split_evidence(event_id: str, split: SplitReadEvidence) -> EventEvidence:
    return EventEvidence(
        event_id=event_id,
        read_name=split.read_name,
        evidence_type="split_read_sa",
        chrom=split.chrom,
        pos=split.pos,
        mate_chrom=split.remote_chrom,
        mate_pos=split.remote_pos,
        cigar=split.cigar,
        mapq=split.mapq,
        clip_side=split.side,
        sa_tag=split.sa_tag,
    )


def classify_event_type(event: RepairEvent) -> str:
    return event.event_type
