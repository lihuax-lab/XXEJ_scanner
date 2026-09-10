"""REF-like and ALT-like evidence counting."""

from __future__ import annotations

from .bam_evidence import count_spanning_reads
from .models import RepairEvent, ScannerConfig


def compute_repair_evidence_fraction(alt_support: int, ref_support: int) -> float:
    # This fraction is an enrichment-data evidence ratio. It deliberately avoids
    # the term VAF because CUT&Tag coverage is not genome-uniform.
    denominator = alt_support + ref_support
    return alt_support / denominator if denominator else 0.0


def count_ref_like_reads(
    bam_path: str,
    event: RepairEvent,
    config: ScannerConfig,
    *,
    min_mapq: int | None = None,
) -> int:
    support_a, support_b = count_ref_like_breakends(
        bam_path, event, config, min_mapq=min_mapq
    )
    return min(support_a, int(support_b)) if support_b != "NA" else support_a


def count_ref_like_breakends(
    bam_path: str,
    event: RepairEvent,
    config: ScannerConfig,
    *,
    min_mapq: int | None = None,
) -> tuple[int, int | str]:
    def count(chrom: str, pos: int) -> int:
        return count_spanning_reads(
            bam_path,
            chrom,
            max(0, pos - 1),
            pos + 1,
            config,
            min_mapq=min_mapq,
        )

    support_a = count(event.bkp_A_chrom, int(event.bkp_A_pos))
    support_b = (
        count(event.bkp_B_chrom, int(event.bkp_B_pos))
        if event.bkp_B_chrom != "NA" and event.bkp_B_pos != "NA"
        else "NA"
    )
    return support_a, support_b


def update_event_fraction(event: RepairEvent) -> RepairEvent:
    # Mutate the event in place so the final writer can serialize all support
    # fields without recalculating them.
    if event.junction_resolved:
        event.repair_evidence_fraction = compute_repair_evidence_fraction(
            event.alt_support,
            event.ref_spanning_support,
        )
        if event.control_assessed:
            event.control_repair_evidence_fraction = compute_repair_evidence_fraction(
                event.control_alt_support,
                event.control_ref_support,
            )
    return event
