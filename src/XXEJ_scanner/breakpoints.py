"""Soft-clipping based breakpoint clustering."""

from __future__ import annotations

from collections import Counter, defaultdict
from typing import Iterable

from .models import (
    BreakpointCluster,
    CandidateRegion,
    ClipSite,
    RegionEvidence,
    ScannerConfig,
)


def cluster_clip_sites(
    sites: Iterable[ClipSite],
    config: ScannerConfig,
    *,
    region: CandidateRegion | None = None,
) -> list[BreakpointCluster]:
    # Cluster by genomic proximity only. Strand and side are summarized after
    # clustering because both sides can mark the same local DSB-prone position.
    site_list = sorted(
        [site for site in sites if region is None or (site.chrom == region.chrom and region.start <= site.pos <= region.end)],
        key=lambda site: (site.chrom, site.pos),
    )
    if not site_list:
        return []

    clusters: list[list[ClipSite]] = []
    current: list[ClipSite] = [site_list[0]]
    for site in site_list[1:]:
        last = current[-1]
        if site.chrom == last.chrom and site.pos <= last.pos + config.clip_cluster_window:
            current.append(site)
            continue
        clusters.append(current)
        current = [site]
    clusters.append(current)

    return [_build_cluster(cluster, region.region_id if region else "NA") for cluster in clusters]


def cluster_evidence_graph(
    evidence: RegionEvidence,
    config: ScannerConfig,
    *,
    region: CandidateRegion | None = None,
) -> list[BreakpointCluster]:
    """Cluster clipped breakpoints with lightweight evidence-graph communities.

    The output remains a soft-clip-derived BreakpointCluster. Non-clip evidence
    only strengthens links between nearby clipped observations; deletion
    endpoints are represented as separate anchors so a CIGAR deletion does not
    collapse both breakpoints into one cluster.
    """
    site_list = sorted(
        [
            site
            for site in evidence.clip_sites
            if region is None
            or (
                site.chrom == region.chrom
                and region.start <= site.pos <= region.end
            )
        ],
        key=lambda site: (site.chrom, site.pos, site.read_name, site.side),
    )
    if not site_list:
        return []

    graph = _build_clip_evidence_graph(site_list, evidence, config, region)
    communities = _label_propagation_communities(graph)
    region_id = region.region_id if region else evidence.region.region_id
    clusters = [
        _build_cluster([site_list[idx] for idx in community], region_id)
        for community in communities
    ]
    return sorted(
        clusters,
        key=lambda cluster: (cluster.chrom, cluster.peak_pos, cluster.cluster_start),
    )


def _build_cluster(sites: list[ClipSite], region_id: str) -> BreakpointCluster:
    # Use the modal clipped coordinate as the peak. This is more stable than the
    # midpoint when a true breakpoint has a small tail of nearby noisy clips.
    positions = [site.pos for site in sites]
    position_counts = Counter(positions)
    peak_pos = position_counts.most_common(1)[0][0]
    left_count = sum(1 for site in sites if site.side == "left_clip")
    right_count = sum(1 for site in sites if site.side == "right_clip")
    plus_count = sum(1 for site in sites if site.strand == "+")
    minus_count = sum(1 for site in sites if site.strand == "-")
    if left_count and right_count:
        side = "both"
    elif left_count:
        side = "left_clip"
    else:
        side = "right_clip"
    return BreakpointCluster(
        region_id=region_id,
        chrom=sites[0].chrom,
        cluster_start=min(positions),
        cluster_end=max(positions) + 1,
        peak_pos=peak_pos,
        clip_side=side,
        clip_count=len(sites),
        left_clip_count=left_count,
        right_clip_count=right_count,
        strand_plus_count=plus_count,
        strand_minus_count=minus_count,
    )


def _build_clip_evidence_graph(
    sites: list[ClipSite],
    evidence: RegionEvidence,
    config: ScannerConfig,
    region: CandidateRegion | None,
) -> dict[int, dict[int, float]]:
    graph: dict[int, dict[int, float]] = {idx: {} for idx in range(len(sites))}
    indices_by_chrom: dict[str, list[int]] = defaultdict(list)
    for idx, site in enumerate(sites):
        indices_by_chrom[site.chrom].append(idx)

    window = max(1, config.clip_cluster_window)
    graph_link_window = window * 2

    for indices in indices_by_chrom.values():
        for offset, left_idx in enumerate(indices):
            left = sites[left_idx]
            for right_idx in indices[offset + 1 :]:
                right = sites[right_idx]
                distance = right.pos - left.pos
                if distance > graph_link_window:
                    break
                if distance <= window:
                    _add_edge(
                        graph,
                        left_idx,
                        right_idx,
                        _proximity_weight(distance, window),
                    )
                if left.read_name and left.read_name == right.read_name:
                    _add_edge(graph, left_idx, right_idx, 1.5)

    evidence_window = graph_link_window
    for chrom, anchors in _evidence_anchors(evidence, region).items():
        indices = indices_by_chrom.get(chrom, [])
        for anchor_pos, anchor_weight in anchors:
            nearby = [
                idx
                for idx in indices
                if abs(sites[idx].pos - anchor_pos) <= evidence_window
            ]
            for offset, left_idx in enumerate(nearby):
                left = sites[left_idx]
                for right_idx in nearby[offset + 1 :]:
                    right = sites[right_idx]
                    if abs(right.pos - left.pos) > graph_link_window:
                        continue
                    anchor_distance = max(
                        abs(left.pos - anchor_pos),
                        abs(right.pos - anchor_pos),
                    )
                    decay = max(
                        0.25,
                        1.0 - anchor_distance / max(1, evidence_window + 1),
                    )
                    _add_edge(
                        graph,
                        left_idx,
                        right_idx,
                        anchor_weight * decay,
                    )
    return graph


def _evidence_anchors(
    evidence: RegionEvidence,
    region: CandidateRegion | None,
) -> dict[str, list[tuple[int, float]]]:
    anchors: dict[str, list[tuple[int, float]]] = defaultdict(list)

    def add(chrom: str, pos: int, weight: float) -> None:
        if region is None or (
            chrom == region.chrom and region.start <= pos <= region.end
        ):
            anchors[chrom].append((pos, weight))

    for indel in evidence.indels:
        if indel.operation == "INS":
            add(indel.chrom, indel.start, 2.5)
        elif indel.operation == "DEL":
            add(indel.chrom, indel.start, 2.0)
            add(indel.chrom, indel.end, 2.0)

    for split in evidence.split_reads:
        add(split.chrom, split.pos, 2.0)
        if split.remote_chrom == split.chrom:
            add(split.remote_chrom, split.remote_pos, 1.5)

    for pair in evidence.discordant_pairs:
        add(pair.chrom, pair.pos, 1.5)
        if pair.mate_chrom == pair.chrom:
            add(pair.mate_chrom, pair.mate_pos, 1.0)

    return anchors


def _proximity_weight(distance: int, window: int) -> float:
    return 1.0 + 2.0 * (1.0 - min(distance, window) / window)


def _add_edge(
    graph: dict[int, dict[int, float]],
    left: int,
    right: int,
    weight: float,
) -> None:
    if left == right or weight <= 0:
        return
    graph[left][right] = graph[left].get(right, 0.0) + weight
    graph[right][left] = graph[right].get(left, 0.0) + weight


def _label_propagation_communities(
    graph: dict[int, dict[int, float]],
    *,
    max_iterations: int = 20,
) -> list[list[int]]:
    labels = {idx: idx for idx in graph}
    for _iteration in range(max_iterations):
        changed = False
        for idx in sorted(graph):
            if not graph[idx]:
                continue
            weights_by_label: dict[int, float] = defaultdict(float)
            for neighbor, weight in graph[idx].items():
                weights_by_label[labels[neighbor]] += weight
            current_label = labels[idx]
            current_weight = weights_by_label.get(current_label, 0.0)
            best_weight = max(weights_by_label.values())
            best_labels = [
                label
                for label, weight in weights_by_label.items()
                if weight == best_weight
            ]
            best_label = min(best_labels)
            if best_label != current_label and best_weight >= current_weight:
                labels[idx] = best_label
                changed = True
        if not changed:
            break

    communities_by_label: dict[int, list[int]] = defaultdict(list)
    for idx, label in labels.items():
        communities_by_label[label].append(idx)
    communities = [sorted(indices) for indices in communities_by_label.values()]
    communities.sort(key=lambda indices: indices[0])
    return communities


def score_breakpoint_cluster(cluster: BreakpointCluster, config: ScannerConfig) -> BreakpointCluster:
    # Transparent heuristic score: more support and higher local clip rate help,
    # while a high control clip rate penalizes the cluster.
    cluster.clip_rate = cluster.clip_count / cluster.treated_depth if cluster.treated_depth else 0.0
    control_rate = cluster.normal_noise
    support_component = min(1.0, cluster.clip_count / max(1, config.min_alt_support * 2))
    rate_component = min(1.0, cluster.clip_rate * 4.0)
    noise_penalty = min(1.0, control_rate / max(config.max_normal_clip_rate, 1e-6))
    cluster.confidence = max(0.0, 0.55 * support_component + 0.45 * rate_component - 0.40 * noise_penalty)
    return cluster


def filter_breakpoint_clusters(
    clusters: Iterable[BreakpointCluster],
    config: ScannerConfig,
) -> list[BreakpointCluster]:
    # Keep filtering conservative. A cluster can be biologically interesting only
    # after it has enough treated support and is not dominated by control noise.
    kept: list[BreakpointCluster] = []
    for cluster in clusters:
        if cluster.clip_count < config.min_alt_support:
            continue
        if cluster.normal_noise > config.max_normal_clip_rate:
            continue
        kept.append(cluster)
    return kept
