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
    # For BNDs, assess REF-like reads only at the local breakend.
    # For MMEJ_DEL, the two breakpoints can be kilobases apart, so short paired
    # reads can never span the full interval — counting spanning reads across the
    # whole deletion yields a denominator of zero and a spurious fraction of 1.0.
    # Instead, count independently at each breakpoint and take the minimum, which
    # is the most conservative estimate of undeleted molecules in this window.
    if event.event_type.startswith("NHEJ_BND"):
        start = max(0, int(event.bkp_A_pos) - 1)
        end = int(event.bkp_A_pos) + 1
        return count_spanning_reads(bam_path, event.chrom, start, end, config, min_mapq=min_mapq)
    if event.event_type == "MMEJ_DEL":
        bkp_a = int(event.bkp_A_pos)
        bkp_b = int(event.bkp_B_pos)
        count_a = count_spanning_reads(
            bam_path, event.chrom,
            max(0, bkp_a - 1), bkp_a + 1,
            config, min_mapq=min_mapq,
        )
        count_b = count_spanning_reads(
            bam_path, event.chrom,
            max(0, bkp_b - 1), bkp_b + 1,
            config, min_mapq=min_mapq,
        )
        return min(count_a, count_b)
    start = max(0, int(event.start))
    end = max(start + 1, int(event.end))
    return count_spanning_reads(bam_path, event.chrom, start, end, config, min_mapq=min_mapq)


def update_event_fraction(event: RepairEvent) -> RepairEvent:
    # Mutate the event in place so the final writer can serialize all support
    # fields without recalculating them.
    event.repair_evidence_fraction = compute_repair_evidence_fraction(
        event.alt_support,
        event.ref_spanning_support,
    )
    event.control_repair_evidence_fraction = compute_repair_evidence_fraction(
        event.control_alt_support,
        event.control_ref_support,
    )
    return event
