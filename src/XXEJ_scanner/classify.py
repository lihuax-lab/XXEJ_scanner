"""Rule-based repair event classification."""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from statistics import median

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
from .reference import ReferenceGenome, find_microhomology


def _near(pos_a: int, pos_b: int, window: int) -> bool:
    return abs(pos_a - pos_b) <= window


def _cluster_side(cluster: BreakpointCluster) -> str:
    return cluster.clip_side


def _cluster_key(cluster: BreakpointCluster) -> tuple[str, int, str]:
    return (cluster.chrom, cluster.peak_pos, cluster.clip_side)


def _microhomology_score_bonus(length: int, low_complexity: bool = False) -> float:
    if length <= 3:
        return 0.0
    if length <= 5:
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
    sorted_clusters = sorted(clusters, key=lambda cluster: cluster.peak_pos)
    deletion_events, deletion_evidence = _classify_local_del(
        region, sorted_clusters, evidence, reference, config
    )
    bnd_events, bnd_evidence = classify_bnd_events(
        region, sorted_clusters, evidence, config
    )
    insertion_events, insertion_evidence = _classify_local_ins(
        region, sorted_clusters, evidence, config
    )
    return (
        deletion_events + bnd_events + insertion_events,
        deletion_evidence + bnd_evidence + insertion_evidence,
    )


def assign_final_event_ids(
    events: list[RepairEvent],
    evidence: list[EventEvidence],
) -> tuple[list[RepairEvent], list[EventEvidence]]:
    id_aliases: dict[str, str] = {}
    unique: dict[tuple[object, ...], RepairEvent] = {}
    for event in events:
        canonical = unique.get(event.allele_key)
        if canonical is None:
            unique[event.allele_key] = event
            continue
        id_aliases[event.event_id] = canonical.event_id
        canonical.support_read_names.update(event.support_read_names)
        canonical.junction_evidence_types.update(event.junction_evidence_types)
        canonical.junction_evidence_support = len(canonical.support_read_names)
        for field in (
            "alt_clip_support",
            "alt_split_support",
            "alt_discordant_pair_support",
            "alt_indel_support",
        ):
            setattr(canonical, field, max(getattr(canonical, field), getattr(event, field)))
    events[:] = list(unique.values())
    for row in evidence:
        row.event_id = id_aliases.get(row.event_id, row.event_id)

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


def _nearest_remote_cluster_for_pair(
    pair: DiscordantPair,
    clusters: list[BreakpointCluster],
    config: ScannerConfig,
) -> BreakpointCluster | None:
    candidates = [
        cluster
        for cluster in clusters
        if pair.chrom == cluster.chrom
        and _is_remote_bnd_anchor(
            cluster.chrom,
            cluster.peak_pos,
            pair.mate_chrom,
            pair.mate_pos,
            config,
        )
    ]
    if not candidates:
        return None
    return min(
        candidates,
        key=lambda cluster: (
            abs(pair.pos - cluster.peak_pos),
            -cluster.clip_count,
            cluster.peak_pos,
            cluster.clip_side,
        ),
    )


def classify_bnd_events(
    region: CandidateRegion,
    clusters: list[BreakpointCluster],
    evidence: RegionEvidence,
    config: ScannerConfig,
) -> tuple[list[RepairEvent], list[EventEvidence]]:
    events: list[RepairEvent] = []
    event_evidence: list[EventEvidence] = []
    discordants_by_cluster: dict[tuple[str, int, str], list[DiscordantPair]] = (
        defaultdict(list)
    )
    pair_only_discordants: list[DiscordantPair] = []
    for pair in evidence.discordant_pairs:
        cluster = _nearest_remote_cluster_for_pair(pair, clusters, config)
        if cluster is None:
            pair_only_discordants.append(pair)
            continue
        discordants_by_cluster[_cluster_key(cluster)].append(pair)

    for cluster in clusters:
        linked_discordants = discordants_by_cluster[_cluster_key(cluster)]
        linked_splits = [
            split
            for split in evidence.split_reads
            if split.chrom == cluster.chrom
            and _near(split.pos, cluster.peak_pos, config.clip_cluster_window * 3)
            and _is_remote_bnd_anchor(
                cluster.chrom,
                cluster.peak_pos,
                split.remote_chrom,
                split.remote_pos,
                config,
            )
        ]
        for group in _group_remote_evidence(
            linked_discordants, linked_splits, config.coverage_bin_size
        ):
            pairs = [item for item in group if isinstance(item, DiscordantPair)]
            splits = [item for item in group if isinstance(item, SplitReadEvidence)]
            first = group[0]
            remote_chrom = (
                first.mate_chrom
                if isinstance(first, DiscordantPair)
                else first.remote_chrom
            )
            support = len({item.read_name for item in pairs + splits})
            if support < config.min_bnd_support:
                continue
            remote_positions = [item.mate_pos for item in pairs] + [
                item.remote_pos for item in splits
            ]
            remote_pos = int(median(remote_positions))
            orientation_values = [
                item.orientation
                for item in (splits or pairs)
                if getattr(item, "orientation", "NA") != "NA"
            ]
            orientation = (
                Counter(orientation_values).most_common(1)[0][0]
                if orientation_values
                else "NA"
            )
            event_type = (
                "BND_INTER"
                if remote_chrom != cluster.chrom
                else "BND_INTRA"
            )
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
                notes="Candidate breakend supported by remote alignment evidence.",
                evidence_level="RESOLVED" if splits else "MULTI_SIGNAL",
                junction_resolved=bool(splits),
                support_read_names={item.read_name for item in pairs + splits},
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

    pair_only_events, pair_only_evidence = _classify_pair_only_bnd_events(
        region, pair_only_discordants, config
    )
    events.extend(pair_only_events)
    event_evidence.extend(pair_only_evidence)
    return events, event_evidence


def _classify_pair_only_bnd_events(
    region: CandidateRegion,
    pairs: list[DiscordantPair],
    config: ScannerConfig,
) -> tuple[list[RepairEvent], list[EventEvidence]]:
    events: list[RepairEvent] = []
    event_evidence: list[EventEvidence] = []
    eligible = [
        pair
        for pair in pairs
        if _is_remote_bnd_anchor(
            pair.chrom, pair.pos, pair.mate_chrom, pair.mate_pos, config
        )
    ]
    for group in _group_pair_evidence(eligible, config.coverage_bin_size):
        local_chrom = group[0].chrom
        remote_chrom = group[0].mate_chrom
        local_positions = [pair.pos for pair in group]
        remote_positions = [pair.mate_pos for pair in group]
        local_pos = int(median(local_positions))
        remote_pos = int(median(remote_positions))
        orientation_values = [
            pair.orientation for pair in group if pair.orientation != "NA"
        ]
        orientation = (
            Counter(orientation_values).most_common(1)[0][0]
            if orientation_values
            else "NA"
        )
        event_type = (
            "BND_INTER"
            if remote_chrom != local_chrom
            else "BND_INTRA"
        )
        temp_event_id = (
            f"TMP_BND_PAIRONLY_{region.region_id}_{local_chrom}_{local_pos}_"
            f"{remote_chrom}_{remote_pos}"
        )
        event = RepairEvent(
            event_id=temp_event_id,
            event_type=event_type,
            chrom=local_chrom,
            start=max(0, local_pos - 1),
            end=local_pos + 1,
            bkp_A_chrom=local_chrom,
            bkp_A_pos=local_pos,
            bkp_A_side="pair_only",
            bkp_B_chrom=remote_chrom,
            bkp_B_pos=remote_pos,
            bkp_B_side="remote",
            remote_chrom=remote_chrom,
            remote_pos=remote_pos,
            orientation=orientation,
            alt_discordant_pair_support=len({pair.read_name for pair in group}),
            treated_depth=region.treated_coverage,
            control_depth=region.control_coverage,
            notes=(
                "Pair-only candidate breakend from discordant mate evidence."
            ),
            evidence_level="PAIR_ONLY",
            support_read_names={pair.read_name for pair in group},
        )
        event.score = _event_score(event)
        events.append(event)
        for pair in group:
            event_evidence.append(_pair_evidence(temp_event_id, pair))

    return events, event_evidence


def _group_remote_evidence(
    pairs: list[DiscordantPair],
    splits: list[SplitReadEvidence],
    window: int,
) -> list[list[object]]:
    def values(item: object) -> tuple[str, str, int]:
        if isinstance(item, DiscordantPair):
            return item.mate_chrom, item.orientation, item.mate_pos
        assert isinstance(item, SplitReadEvidence)
        return item.remote_chrom, item.orientation, item.remote_pos

    groups: list[list[object]] = []
    for item in sorted(splits, key=values):
        key = values(item)
        if groups:
            first_key = values(groups[-1][0])
            if key[:2] == first_key[:2] and key[2] - first_key[2] <= max(1, window):
                groups[-1].append(item)
                continue
        groups.append([item])

    for pair in sorted(pairs, key=values):
        pair_chrom, pair_orientation, pair_pos = values(pair)
        compatible = [
            group
            for group in groups
            if values(group[0])[0] == pair_chrom
            and all(abs(values(item)[2] - pair_pos) <= max(1, window) for item in group)
        ]
        if compatible:
            min(
                compatible,
                key=lambda group: abs(
                    median(values(item)[2] for item in group) - pair_pos
                ),
            ).append(pair)
            continue
        pair_group = next(
            (
                group
                for group in groups
                if isinstance(group[0], DiscordantPair)
                and values(group[0])[:2] == (pair_chrom, pair_orientation)
                and all(
                    abs(values(item)[2] - pair_pos) <= max(1, window)
                    for item in group
                )
            ),
            None,
        )
        if pair_group is None:
            groups.append([pair])
        else:
            pair_group.append(pair)
    return groups


def _group_pair_evidence(
    pairs: list[DiscordantPair], window: int
) -> list[list[DiscordantPair]]:
    groups: list[list[DiscordantPair]] = []
    for pair in sorted(
        pairs,
        key=lambda item: (
            item.chrom,
            item.mate_chrom,
            item.orientation,
            item.pos,
            item.mate_pos,
        ),
    ):
        if groups:
            first = groups[-1][0]
            if (
                (pair.chrom, pair.mate_chrom, pair.orientation)
                == (first.chrom, first.mate_chrom, first.orientation)
                and all(
                    abs(pair.pos - member.pos) <= max(1, window)
                    and abs(pair.mate_pos - member.mate_pos) <= max(1, window)
                    for member in groups[-1]
                )
            ):
                groups[-1].append(pair)
                continue
        groups.append([pair])
    return groups


def _is_remote_bnd_anchor(
    local_chrom: str,
    local_pos: int,
    remote_chrom: str,
    remote_pos: int,
    config: ScannerConfig,
) -> bool:
    if remote_chrom != local_chrom:
        return True
    return abs(remote_pos - local_pos) > config.max_local_event_distance


def _nearest_cluster(
    clusters: list[BreakpointCluster],
    chrom: str,
    pos: int,
    sides: set[str],
    window: int,
) -> BreakpointCluster | None:
    candidates = [
        cluster
        for cluster in clusters
        if cluster.chrom == chrom
        and cluster.clip_side in sides
        and _near(cluster.peak_pos, pos, window)
    ]
    return min(candidates, key=lambda cluster: abs(cluster.peak_pos - pos), default=None)


def _classify_local_del(
    region: CandidateRegion,
    clusters: list[BreakpointCluster],
    evidence: RegionEvidence,
    reference: ReferenceGenome,
    config: ScannerConfig,
) -> tuple[list[RepairEvent], list[EventEvidence]]:
    grouped: dict[tuple[str, int, int], dict[str, list[object]]] = defaultdict(
        lambda: {"indels": [], "splits": []}
    )
    for indel in evidence.indels:
        if indel.operation == "DEL" and 0 < indel.length <= config.max_local_event_distance:
            grouped[(indel.chrom, indel.start, indel.end)]["indels"].append(indel)
    for split in evidence.split_reads:
        if split.chrom != split.remote_chrom:
            continue
        start, end = sorted((split.pos, split.remote_pos))
        if 0 < end - start <= config.max_local_event_distance:
            grouped[(split.chrom, start, end)]["splits"].append(split)

    events: list[RepairEvent] = []
    rows: list[EventEvidence] = []
    for (chrom, start, end), group in sorted(grouped.items()):
        indels = group["indels"]
        splits = group["splits"]
        support_reads = {item.read_name for item in indels + splits}
        if len(support_reads) < config.min_alt_support:
            continue

        left = _nearest_cluster(
            clusters, chrom, start, {"right_clip", "both"}, config.clip_cluster_window
        )
        right = _nearest_cluster(
            clusters, chrom, end, {"left_clip", "both"}, config.clip_cluster_window
        )
        local_clips = [
            site
            for site in evidence.clip_sites
            if site.chrom == chrom
            and (
                _near(site.pos, start, config.clip_cluster_window)
                or _near(site.pos, end, config.clip_cluster_window)
            )
        ]
        remap_clips = (
            _soft_clips_matching_opposite_flanks(
                evidence.clip_sites, reference, left, right, start, end, config
            )
            if left and right
            else []
        )
        support_reads.update(site.read_name for site in remap_clips)
        mh_hit = find_microhomology(
            reference,
            chrom,
            start,
            end,
            config.min_microhomology_length,
            config.max_microhomology_length,
            0,
        )
        junction_types = _junction_evidence_types(indels, splits, remap_clips)
        event_id = f"TMP_DEL_{region.region_id}_{chrom}_{start}_{end}"
        event = RepairEvent(
            event_id=event_id,
            event_type="LOCAL_DEL",
            chrom=chrom,
            start=start,
            end=end,
            bkp_A_chrom=chrom,
            bkp_A_pos=start,
            bkp_A_side=_cluster_side(left) if left else "junction",
            bkp_B_chrom=chrom,
            bkp_B_pos=end,
            bkp_B_side=_cluster_side(right) if right else "junction",
            deleted_length=end - start,
            microhomology=mh_hit.sequence,
            microhomology_length=mh_hit.length,
            alt_clip_support=len({site.read_name for site in local_clips}),
            alt_split_support=len({split.read_name for split in splits}),
            alt_indel_support=len({indel.read_name for indel in indels}),
            treated_depth=max(
                left.treated_depth if left else 0,
                right.treated_depth if right else 0,
            ),
            control_depth=max(
                left.control_depth if left else 0,
                right.control_depth if right else 0,
            ),
            normal_noise=max(
                left.normal_noise if left else 0,
                right.normal_noise if right else 0,
            ),
            notes="Resolved local deletion; microhomology is sequence context only.",
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
            junction_evidence_support=len(support_reads),
            junction_evidence_types=junction_types,
            evidence_level="RESOLVED",
            junction_resolved=True,
            support_read_names=support_reads,
        )
        event.score = _event_score(event)
        events.append(event)
        rows.extend(_clip_evidence(event_id, site) for site in local_clips)
        rows.extend(_indel_evidence(event_id, indel) for indel in indels)
        rows.extend(_split_evidence(event_id, split) for split in splits)
    return events, rows


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


def _classify_local_ins(
    region: CandidateRegion,
    clusters: list[BreakpointCluster],
    evidence: RegionEvidence,
    config: ScannerConfig,
) -> tuple[list[RepairEvent], list[EventEvidence]]:
    grouped: dict[tuple[str, int, str], list[CigarIndel]] = defaultdict(list)
    for indel in evidence.indels:
        if (
            indel.operation == "INS"
            and indel.length <= config.max_insertion_length
            and indel.sequence != "NA"
        ):
            grouped[(indel.chrom, indel.start, indel.sequence)].append(indel)

    events: list[RepairEvent] = []
    rows: list[EventEvidence] = []
    for (chrom, pos, sequence), insertions in sorted(grouped.items()):
        support_reads = {indel.read_name for indel in insertions}
        if len(support_reads) < max(
            config.min_alt_support, config.min_nhej_ins_indel_support
        ):
            continue
        cluster = _nearest_cluster(
            clusters,
            chrom,
            pos,
            {"left_clip", "right_clip", "both"},
            config.clip_cluster_window * 2,
        )
        local_clips = [
            site
            for site in evidence.clip_sites
            if site.chrom == chrom and _near(site.pos, pos, config.clip_cluster_window)
        ]
        event_id = f"TMP_INS_{region.region_id}_{chrom}_{pos}_{sequence}"
        event = RepairEvent(
            event_id=event_id,
            event_type="LOCAL_INS",
            chrom=chrom,
            start=max(0, pos - 1),
            end=pos + 1,
            bkp_A_chrom=chrom,
            bkp_A_pos=pos,
            bkp_A_side=_cluster_side(cluster) if cluster else "junction",
            inserted_sequence=sequence,
            inserted_length=len(sequence),
            alt_clip_support=len({site.read_name for site in local_clips}),
            alt_indel_support=len(support_reads),
            treated_depth=cluster.treated_depth if cluster else region.treated_coverage,
            control_depth=cluster.control_depth if cluster else region.control_coverage,
            normal_noise=cluster.normal_noise if cluster else 0,
            notes="Resolved local insertion grouped by position and inserted sequence.",
            junction_evidence_support=len(support_reads),
            junction_evidence_types={"cigar_ins"},
            evidence_level="RESOLVED",
            junction_resolved=True,
            support_read_names=support_reads,
        )
        event.score = _event_score(event)
        events.append(event)
        rows.extend(_clip_evidence(event_id, site) for site in local_clips)
        rows.extend(_indel_evidence(event_id, indel) for indel in insertions)
    return events, rows


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
