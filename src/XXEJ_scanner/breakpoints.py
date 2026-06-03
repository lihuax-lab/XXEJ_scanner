"""Soft-clipping based breakpoint clustering."""

from __future__ import annotations

from collections import Counter
from typing import Iterable

from .models import BreakpointCluster, CandidateRegion, ClipSite, ScannerConfig


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
