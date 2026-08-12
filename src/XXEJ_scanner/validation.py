"""Simple second-pass validation and event filtering."""

from __future__ import annotations

from .bam_evidence import collect_region_evidence
from .models import CandidateRegion, RepairEvent, ScannerConfig


def is_pair_only_bnd(event: RepairEvent) -> bool:
    return (
        event.event_type.startswith("BND_")
        and event.bkp_A_side == "pair_only"
        and event.alt_clip_support == 0
        and event.alt_discordant_pair_support > 0
    )


def assign_event_filter(event: RepairEvent, config: ScannerConfig) -> str:
    # Filters are ordered from most direct failure mode to most general. The
    # first matching label is reported to keep events.tsv easy to scan.
    if is_pair_only_bnd(event):
        return "PairOnlyBnd"
    if event.evidence_level == "CLIP_ONLY":
        return "ClipOnly"
    if event.alt_support < config.min_alt_support:
        return "LowSupport"
    if event.normal_noise > config.max_normal_clip_rate:
        return "HighControlNoise"
    if (
        event.event_type == "LOCAL_INS"
        and not config.allow_clip_only_nhej_ins
        and event.alt_indel_support < config.min_nhej_ins_indel_support
    ):
        return "NoInsertionEvidence"
    if event.event_type == "LOCAL_DEL" and not event.junction_resolved:
        return "NoJunctionEvidence"
    if event.event_type.startswith("BND_") and (
        event.alt_split_support + event.alt_discordant_pair_support
        < config.min_bnd_support
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
    if is_pair_only_bnd(event):
        event.filter = "PairOnlyBnd"
        return event

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
        if event.start - config.clip_cluster_window
        <= indel.start
        <= event.end + config.clip_cluster_window
    }
    strict_remote_support = {
        pair.read_name for pair in strict_evidence.discordant_pairs
    } | {split.read_name for split in strict_evidence.split_reads}

    if (
        event.event_type == "LOCAL_INS"
        and not config.allow_clip_only_nhej_ins
        and len(strict_indel_support) < config.min_nhej_ins_indel_support
    ):
        event.filter = "NoInsertionEvidence"
        if event.notes:
            event.notes += " "
        event.notes += "Second-pass strict insertion support was below threshold."
        return event

    # Any strict evidence class can rescue the event, but it must reach the same
    # minimum support threshold used during discovery.
    strict_total = len(
        strict_clip_support | strict_indel_support | strict_remote_support
    )
    if strict_total < config.min_alt_support:
        event.filter = "LowSupport"
        if event.notes:
            event.notes += " "
        event.notes += "Second-pass strict support was below threshold."
        return event

    event.filter = assign_event_filter(event, config)
    return event
