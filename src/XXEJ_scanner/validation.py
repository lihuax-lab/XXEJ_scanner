"""Simple second-pass validation and event filtering."""

from __future__ import annotations

from .bam_evidence import collect_region_evidence
from .models import CandidateRegion, RegionEvidence, RepairEvent, ScannerConfig


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
    if event.control_assessed:
        if event.control_alt_support > config.max_control_alt_support:
            return "PresentInControl"
        if (
            event.control_alt_support == 0
            and event.control_ref_support == 0
            and event.control_depth <= 0
        ):
            return "NoControlCoverage"
    if event.normal_noise > config.max_normal_clip_rate:
        return "HighControlNoise"
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

    loci = [(event.bkp_A_chrom, int(event.bkp_A_pos))]
    if event.bkp_B_chrom != "NA" and event.bkp_B_pos != "NA":
        loci.append((event.bkp_B_chrom, int(event.bkp_B_pos)))
    strict_support: set[str] = set()
    for idx, (chrom, pos) in enumerate(dict.fromkeys(loci), 1):
        region = CandidateRegion(
            chrom=chrom,
            start=max(0, pos - config.second_pass_window),
            end=pos + config.second_pass_window,
            region_id=f"{event.event_id}_second_pass_{idx}",
        )
        strict_evidence = collect_region_evidence(
            treated_bam,
            region,
            config,
            padding=0,
            min_mapq=config.strict_min_mapq,
            min_baseq=(
                config.strict_min_breakpoint_baseq
                if config.min_breakpoint_baseq > 0
                else 0
            ),
        )
        strict_support.update(
            matching_event_read_names(event, strict_evidence, config)
        )

    threshold = (
        config.min_bnd_support
        if event.event_type.startswith("BND_")
        else config.min_alt_support
    )
    if len(strict_support) < threshold:
        event.filter = "LowSupport"
        if event.notes:
            event.notes += " "
        event.notes += "Second-pass strict support was below threshold."
        return event

    event.filter = assign_event_filter(event, config)
    return event


def matching_event_read_names(
    event: RepairEvent,
    evidence: RegionEvidence,
    config: ScannerConfig,
) -> set[str]:
    if event.event_type == "LOCAL_INS":
        return {
            indel.read_name
            for indel in evidence.indels
            if indel.operation == "INS"
            and indel.chrom == event.bkp_A_chrom
            and indel.start == int(event.bkp_A_pos)
            and indel.sequence == event.inserted_sequence
        }
    if event.event_type == "LOCAL_DEL":
        support = {
            indel.read_name
            for indel in evidence.indels
            if indel.operation == "DEL"
            and indel.chrom == event.bkp_A_chrom
            and indel.start == int(event.bkp_A_pos)
            and indel.end == int(event.bkp_B_pos)
        }
        for split in evidence.split_reads:
            if (
                split.chrom != event.bkp_A_chrom
                or split.remote_chrom != event.bkp_B_chrom
            ):
                continue
            if sorted((split.pos, split.remote_pos)) == [
                int(event.bkp_A_pos),
                int(event.bkp_B_pos),
            ]:
                support.add(split.read_name)
        return support

    support: set[str] = set()
    for pair in evidence.discordant_pairs:
        if _relation_matches(
            event,
            pair.chrom,
            pair.pos,
            pair.mate_chrom,
            pair.mate_pos,
            config,
        ):
            support.add(pair.read_name)
    for split in evidence.split_reads:
        relation = _relation_matches(
            event,
            split.chrom,
            split.pos,
            split.remote_chrom,
            split.remote_pos,
            config,
        )
        if relation == "direct" and split.orientation == event.orientation:
            support.add(split.read_name)
        elif relation == "swapped" and split.orientation[::-1] == event.orientation:
            support.add(split.read_name)
    return support


def _relation_matches(
    event: RepairEvent,
    local_chrom: str,
    local_pos: int,
    remote_chrom: str,
    remote_pos: int,
    config: ScannerConfig,
) -> str:
    local_window = config.clip_cluster_window * 3
    remote_window = config.coverage_bin_size
    if (
        local_chrom == event.bkp_A_chrom
        and remote_chrom == event.bkp_B_chrom
        and abs(local_pos - int(event.bkp_A_pos)) <= local_window
        and abs(remote_pos - int(event.bkp_B_pos)) <= remote_window
    ):
        return "direct"
    if (
        local_chrom == event.bkp_B_chrom
        and remote_chrom == event.bkp_A_chrom
        and abs(local_pos - int(event.bkp_B_pos)) <= remote_window
        and abs(remote_pos - int(event.bkp_A_pos)) <= local_window
    ):
        return "swapped"
    return ""
