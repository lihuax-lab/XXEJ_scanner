"""Rule-based repair event classification."""

from __future__ import annotations

import math
from collections import Counter, defaultdict

from .bam_evidence import dominant_sequence
from .models import (
    BreakpointCluster,
    CandidateRegion,
    CigarIndel,
    DiscordantPair,
    EventEvidence,
    RegionEvidence,
    RepairEvent,
    ScannerConfig,
    SplitReadEvidence,
)
from .reference import ReferenceGenome, detect_microhomology


def _near(pos_a: int, pos_b: int, window: int) -> bool:
    return abs(pos_a - pos_b) <= window


def _cluster_side(cluster: BreakpointCluster) -> str:
    return cluster.clip_side


def _event_score(event: RepairEvent) -> float:
    # First-version scoring is deliberately transparent. It is a ranking aid for
    # review/IGV triage, not a calibrated probability or genotype likelihood.
    microhomology_bonus = 1.0 if event.microhomology_length else 0.0
    control_noise_penalty = event.normal_noise
    mapping_quality_penalty = 0.0
    return (
        1.0 * math.log2(1 + max(event.treated_depth, 0))
        + 2.0 * math.log2(1 + event.alt_clip_support)
        + 2.5 * math.log2(1 + event.alt_split_support)
        + 2.5 * math.log2(1 + event.alt_discordant_pair_support)
        + 1.5 * math.log2(1 + event.alt_indel_support)
        + 1.0 * microhomology_bonus
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

    for idx, event in enumerate(events, 1):
        # Classifiers use temporary IDs so their evidence rows can be joined
        # before final deterministic event IDs are assigned.
        old_event_id = event.event_id
        event.event_id = f"XEJ_{idx:06d}"
        for ev in event_evidence:
            if ev.event_id == old_event_id:
                ev.event_id = event.event_id
    return events, event_evidence


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
            temp_event_id = f"TMP_BND_{cluster.chrom}_{cluster.peak_pos}_{remote_chrom}_{remote_pos}"
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

    # Select the strongest same-chromosome cluster pair within the configured
    # local-event distance. This keeps the first version simple and avoids
    # enumerating many overlapping deletion hypotheses from one noisy region.
    best_pair: tuple[BreakpointCluster, BreakpointCluster, int] | None = None
    best_support = -1
    for idx, left in enumerate(clusters):
        for right in clusters[idx + 1 :]:
            if left.chrom != right.chrom:
                continue
            distance = abs(right.peak_pos - left.peak_pos)
            if distance <= 0 or distance > config.max_local_event_distance:
                continue
            deletion_support = _indels_near_pair(
                evidence.indels, left.peak_pos, right.peak_pos, config
            )
            support = (
                left.clip_count
                + right.clip_count
                + len({indel.read_name for indel in deletion_support})
            )
            if support > best_support:
                best_pair = (
                    left,
                    right,
                    len({indel.read_name for indel in deletion_support}),
                )
                best_support = support

    if best_pair is None:
        return events, event_evidence, used

    left, right, indel_support = best_pair
    # Microhomology is reference-context evidence. It boosts interpretation but
    # is not required on its own to create a candidate event.
    mh_seq, mh_len = detect_microhomology(
        reference,
        left.chrom,
        left.peak_pos,
        right.peak_pos,
        config.min_microhomology_length,
        config.max_microhomology_length,
    )
    support_reads = {
        site.read_name
        for site in evidence.clip_sites
        if site.chrom == left.chrom
        and (
            _near(site.pos, left.peak_pos, config.clip_cluster_window)
            or _near(site.pos, right.peak_pos, config.clip_cluster_window)
        )
    }
    deletion_indels = _indels_near_pair(
        evidence.indels, left.peak_pos, right.peak_pos, config
    )
    support_reads.update(indel.read_name for indel in deletion_indels)
    if len(support_reads) < config.min_alt_support:
        return events, event_evidence, used

    event = RepairEvent(
        event_id=f"TMP_MMEJ_{left.chrom}_{left.peak_pos}_{right.peak_pos}",
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
        microhomology=mh_seq,
        microhomology_length=mh_len,
        alt_clip_support=left.clip_count + right.clip_count,
        alt_indel_support=indel_support,
        treated_depth=max(left.treated_depth, right.treated_depth),
        control_depth=max(left.control_depth, right.control_depth),
        normal_noise=max(left.normal_noise, right.normal_noise),
        notes="Candidate local deletion with breakpoint clusters and microhomology context.",
        support_read_names=support_reads,
    )
    event.score = _event_score(event)
    events.append(event)
    for cluster_pos in (left.peak_pos, right.peak_pos):
        for site in evidence.clip_sites:
            if site.chrom == left.chrom and _near(
                site.pos, cluster_pos, config.clip_cluster_window
            ):
                event_evidence.append(_clip_evidence(event.event_id, site))
    for indel in deletion_indels:
        event_evidence.append(_indel_evidence(event.event_id, indel))
    used.add((left.chrom, left.peak_pos, left.clip_side))
    used.add((right.chrom, right.peak_pos, right.clip_side))
    return events, event_evidence, used


def _classify_nhej_ins(
    region: CandidateRegion,
    cluster: BreakpointCluster,
    evidence: RegionEvidence,
    config: ScannerConfig,
) -> tuple[RepairEvent | None, list[EventEvidence]]:
    # NHEJ_INS is intentionally local: it uses nearby short CIGAR insertions and
    # clipped sequence but avoids creating a local insertion if remote support is
    # strong enough to classify the cluster as a BND first.
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
    support_reads = {site.read_name for site in local_clips}
    support_reads.update(indel.read_name for indel in local_insertions)
    if len(support_reads) < config.min_alt_support:
        return None, []

    inserted_sequence = dominant_sequence(indel.sequence for indel in local_insertions)
    if inserted_sequence == "NA":
        # If no explicit CIGAR insertion exists, use short clipped sequence as a
        # conservative filler-sequence proxy.
        inserted_sequence = dominant_sequence(
            site.clip_sequence
            for site in local_clips
            if site.clip_length <= config.max_insertion_length
        )
    inserted_length: int | str = (
        len(inserted_sequence)
        if inserted_sequence != "NA"
        else (
            round(sum(site.clip_length for site in local_clips) / len(local_clips))
            if local_clips
            else "NA"
        )
    )
    event = RepairEvent(
        event_id=f"TMP_INS_{cluster.chrom}_{cluster.peak_pos}_{cluster.clip_side}",
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
        notes="Candidate local NHEJ-like insertion/filler sequence near a clipped breakpoint.",
        support_read_names=support_reads,
    )
    event.score = _event_score(event)
    ev_rows = [_clip_evidence(event.event_id, site) for site in local_clips]
    ev_rows.extend(_indel_evidence(event.event_id, indel) for indel in local_insertions)
    return event, ev_rows


def _indels_near_pair(
    indels: list[CigarIndel],
    left_pos: int,
    right_pos: int,
    config: ScannerConfig,
) -> list[CigarIndel]:
    # A deletion CIGAR may not exactly match clipped peak positions, so use the
    # cluster window as tolerance on both sides.
    lo, hi = sorted((left_pos, right_pos))
    return [
        indel
        for indel in indels
        if indel.operation == "DEL"
        and indel.start <= hi + config.clip_cluster_window
        and indel.end >= lo - config.clip_cluster_window
    ]


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
