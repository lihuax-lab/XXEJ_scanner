"""Output writers for scanner result tables."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path
from typing import Iterable

from .models import (
    BreakpointCluster,
    CandidateRegion,
    ClipSite,
    DiscordantPair,
    EventEvidence,
    RepairEvent,
    ScannerConfig,
    SplitReadEvidence,
)
from .utils import format_float, safe_mkdir


def _value(value: object) -> str:
    # Keep all result tables string-safe and spreadsheet-friendly. Missing values
    # are written as NA rather than empty fields.
    if value is None:
        return "NA"
    if isinstance(value, float):
        return format_float(value)
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, set):
        return ",".join(sorted(str(item) for item in value)) if value else "NA"
    return str(value)


def write_candidate_regions_bed(path: str, regions: Iterable[CandidateRegion]) -> None:
    # BED output intentionally has no header so it can be loaded directly in IGV
    # or intersected with bedtools.
    with Path(path).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for region in regions:
            writer.writerow(
                [
                    region.chrom,
                    region.start,
                    region.end,
                    region.region_id,
                    format_float(region.score),
                    format_float(region.treated_coverage),
                    format_float(region.control_coverage),
                    format_float(region.log2fc),
                    format_float(region.clip_rate),
                ]
            )


def write_breakpoint_clusters_tsv(
    path: str, clusters: Iterable[BreakpointCluster]
) -> None:
    fields = [
        "region_id",
        "chrom",
        "cluster_start",
        "cluster_end",
        "peak_pos",
        "clip_side",
        "clip_count",
        "left_clip_count",
        "right_clip_count",
        "treated_depth",
        "control_depth",
        "clip_rate",
        "strand_plus_count",
        "strand_minus_count",
        "normal_noise",
        "confidence",
    ]
    _write_dataclass_tsv(path, fields, clusters)


def write_events_tsv(path: str, events: Iterable[RepairEvent]) -> None:
    fields = [
        "event_id",
        "event_type",
        "chrom",
        "start",
        "end",
        "bkp_A_chrom",
        "bkp_A_pos",
        "bkp_A_side",
        "bkp_B_chrom",
        "bkp_B_pos",
        "bkp_B_side",
        "remote_chrom",
        "remote_pos",
        "orientation",
        "inserted_sequence",
        "inserted_length",
        "deleted_length",
        "microhomology",
        "microhomology_length",
        "alt_clip_support",
        "alt_split_support",
        "alt_discordant_pair_support",
        "alt_indel_support",
        "ref_spanning_support",
        "treated_depth",
        "control_depth",
        "repair_evidence_fraction",
        "control_alt_support",
        "control_ref_support",
        "control_repair_evidence_fraction",
        "normal_noise",
        "score",
        "filter",
        "notes",
        "microhomology_left_end",
        "microhomology_right_start",
        "microhomology_offset_a",
        "microhomology_offset_b",
        "microhomology_deletion_start",
        "microhomology_deletion_end",
        "microhomology_deletion_length",
        "microhomology_ambiguity_bases",
        "microhomology_equivalent_hits",
        "microhomology_low_complexity",
        "junction_evidence_support",
        "junction_evidence_types",
    ]
    _write_dataclass_tsv(path, fields, events)


def write_event_evidence_tsv(path: str, evidence: Iterable[EventEvidence]) -> None:
    fields = [
        "event_id",
        "read_name",
        "evidence_type",
        "chrom",
        "pos",
        "mate_chrom",
        "mate_pos",
        "cigar",
        "mapq",
        "is_reverse",
        "mate_is_reverse",
        "clip_side",
        "clip_length",
        "clip_sequence",
        "sa_tag",
        "classification",
    ]
    _write_dataclass_tsv(path, fields, evidence)


def write_raw_clip_sites_tsv(path: str, sites: Iterable[ClipSite]) -> None:
    fields = [
        "chrom",
        "pos",
        "side",
        "clip_length",
        "clip_sequence",
        "read_name",
        "strand",
        "mapq",
        "cigar",
        "is_reverse",
        "reference_start",
        "reference_end",
        "clip_type",
    ]
    _write_dataclass_tsv(path, fields, sites)


def write_raw_discordant_pairs_tsv(path: str, pairs: Iterable[DiscordantPair]) -> None:
    fields = [
        "read_name",
        "chrom",
        "pos",
        "mate_chrom",
        "mate_pos",
        "orientation",
        "mapq",
        "is_reverse",
        "mate_is_reverse",
        "cigar",
        "reason",
    ]
    _write_dataclass_tsv(path, fields, pairs)


def write_raw_split_reads_tsv(path: str, splits: Iterable[SplitReadEvidence]) -> None:
    fields = [
        "read_name",
        "chrom",
        "pos",
        "side",
        "remote_chrom",
        "remote_pos",
        "remote_strand",
        "remote_cigar",
        "remote_mapq",
        "orientation",
        "mapq",
        "cigar",
        "sa_tag",
    ]
    _write_dataclass_tsv(path, fields, splits)


def write_igv_loci_bed(
    path: str, events: Iterable[RepairEvent], flank: int = 100
) -> None:
    with Path(path).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for event in events:
            writer.writerow(
                [
                    event.chrom,
                    max(0, int(event.start) - flank),
                    int(event.end) + flank,
                    event.event_id,
                ]
            )
            if event.remote_chrom != "NA" and event.remote_pos != "NA":
                remote = int(event.remote_pos)
                writer.writerow(
                    [
                        event.remote_chrom,
                        max(0, remote - flank),
                        remote + flank,
                        f"{event.event_id}_remote",
                    ]
                )


def write_run_summary_json(path: str, summary: dict[str, object]) -> None:
    with Path(path).open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")


def _write_dataclass_tsv(path: str, fields: list[str], rows: Iterable[object]) -> None:
    # Shared TSV writer keeps column order explicit at each call site while still
    # using dataclass serialization for the row values.
    with Path(path).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        for row in rows:
            values = asdict(row)
            writer.writerow([_value(values.get(field)) for field in fields])


def prepare_output_dir(output_dir: str) -> dict[str, str]:
    # Return named paths instead of constructing filenames throughout cli.py.
    # This keeps the output contract in one place.
    safe_mkdir(output_dir)
    files = {
        "candidate_regions": "candidate_regions.bed",
        "breakpoint_clusters": "breakpoint_clusters.tsv",
        "events": "events.tsv",
        "event_evidence": "event_evidence.tsv",
        "run_summary": "run_summary.json",
        "raw_clip_sites": "raw_clip_sites.tsv",
        "raw_discordant_pairs": "raw_discordant_pairs.tsv",
        "raw_split_reads": "raw_split_reads.tsv",
        "igv_loci": "igv_loci.bed",
    }
    return {key: str(Path(output_dir) / name) for key, name in files.items()}
