"""Simple second-pass validation and event filtering."""

from __future__ import annotations

from .bam_evidence import collect_region_evidence
from .models import CandidateRegion, RepairEvent, ScannerConfig


def assign_event_filter(event: RepairEvent, config: ScannerConfig) -> str:
    # Filters are ordered from most direct failure mode to most general. The
    # first matching label is reported to keep events.tsv easy to scan.
    if event.alt_support < config.min_alt_support:
        return "LowSupport"
    if event.normal_noise > config.max_normal_clip_rate:
        return "HighControlNoise"
    if event.alt_clip_support < config.min_alt_support and event.event_type in {"NHEJ_INS", "MMEJ_DEL"}:
        return "WeakClipCluster"
    if event.event_type.startswith("NHEJ_BND") and (
        event.alt_split_support + event.alt_discordant_pair_support < config.min_bnd_support
    ):
        return "NoRemoteSupport"
    if event.score <= 0:
        return "AmbiguousEventType"
    return "PASS"


def second_pass_validate_event(
    event: RepairEvent,
    treated_bam: str,
    config: ScannerConfig,
) -> RepairEvent:
    # Re-scan a tight local window with stricter MAPQ. This catches candidates
    # whose broad-region support disappears under stricter evidence criteria.
    region = CandidateRegion(
        chrom=event.chrom,
        start=max(0, int(event.bkp_A_pos) - config.second_pass_window),
        end=int(event.bkp_A_pos) + config.second_pass_window,
        region_id=f"{event.event_id}_second_pass",
    )
    strict_evidence = collect_region_evidence(
        treated_bam,
        region,
        config,
        padding=0,
        min_mapq=config.strict_min_mapq,
    )
    strict_clip_support = {
        site.read_name
        for site in strict_evidence.clip_sites
        if abs(site.pos - int(event.bkp_A_pos)) <= config.clip_cluster_window
    }
    strict_indel_support = {
        indel.read_name
        for indel in strict_evidence.indels
        if event.start - config.clip_cluster_window <= indel.start <= event.end + config.clip_cluster_window
    }
    strict_remote_support = {pair.read_name for pair in strict_evidence.discordant_pairs} | {
        split.read_name for split in strict_evidence.split_reads
    }

    # Any strict evidence class can rescue the event, but it must reach the same
    # minimum support threshold used during discovery.
    strict_total = len(strict_clip_support | strict_indel_support | strict_remote_support)
    if strict_total < config.min_alt_support:
        event.filter = "LowSupport"
        if event.notes:
            event.notes += " "
        event.notes += "Second-pass strict support was below threshold."
        return event

    event.filter = assign_event_filter(event, config)
    return event
