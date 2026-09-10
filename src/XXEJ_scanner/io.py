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
    RegionEvidence,
    RepairEvent,
    ScannerConfig,
    SplitReadEvidence,
)
from .utils import format_float, safe_mkdir

RAW_CLIP_FIELDS = [
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
RAW_DISCORDANT_PAIR_FIELDS = [
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
RAW_SPLIT_READ_FIELDS = [
    "read_name",
    "chrom",
    "pos",
    "side",
    "remote_chrom",
    "remote_pos",
    "remote_strand",
    "remote_cigar",
    "remote_mapq",
    "remote_nm",
    "orientation",
    "mapq",
    "cigar",
    "sa_tag",
]


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
        "ref_support_A",
        "ref_support_B",
        "treated_depth",
        "control_depth",
        "repair_evidence_fraction",
        "control_alt_support",
        "control_ref_support",
        "control_ref_support_A",
        "control_ref_support_B",
        "control_repair_evidence_fraction",
        "control_assessed",
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
        "evidence_level",
        "junction_resolved",
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
    _write_dataclass_tsv(path, RAW_CLIP_FIELDS, sites)


def write_raw_discordant_pairs_tsv(path: str, pairs: Iterable[DiscordantPair]) -> None:
    _write_dataclass_tsv(path, RAW_DISCORDANT_PAIR_FIELDS, pairs)


def write_raw_split_reads_tsv(path: str, splits: Iterable[SplitReadEvidence]) -> None:
    _write_dataclass_tsv(path, RAW_SPLIT_READ_FIELDS, splits)


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


class _DataclassTsvWriter:
    def __init__(self, path: str, fields: list[str]) -> None:
        self.path = path
        self.fields = fields
        self.handle = None
        self.writer = None

    def open(self) -> None:
        self.handle = Path(self.path).open("w", newline="")
        self.writer = csv.writer(self.handle, delimiter="\t", lineterminator="\n")
        self.writer.writerow(self.fields)

    def write(self, rows: Iterable[object]) -> None:
        assert self.writer is not None
        for row in rows:
            values = asdict(row)
            self.writer.writerow([_value(values.get(field)) for field in self.fields])

    def close(self) -> None:
        if self.handle is not None:
            self.handle.close()

    def __enter__(self) -> "_DataclassTsvWriter":
        self.open()
        return self

    def __exit__(self, *_exc: object) -> None:
        self.close()


def _write_dataclass_tsv(path: str, fields: list[str], rows: Iterable[object]) -> None:
    # Shared TSV writer keeps column order explicit at each call site while still
    # using dataclass serialization for the row values.
    with _DataclassTsvWriter(path, fields) as writer:
        writer.write(rows)


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


class RawEvidenceWriter:
    """Write raw evidence incrementally so large scans do not retain every row."""

    def __init__(self, clip_path: str, pair_path: str, split_path: str) -> None:
        self._writers = (
            _DataclassTsvWriter(clip_path, RAW_CLIP_FIELDS),
            _DataclassTsvWriter(pair_path, RAW_DISCORDANT_PAIR_FIELDS),
            _DataclassTsvWriter(split_path, RAW_SPLIT_READ_FIELDS),
        )

    def __enter__(self) -> "RawEvidenceWriter":
        for writer in self._writers:
            writer.open()
        return self

    def write(self, evidence: RegionEvidence) -> None:
        self._writers[0].write(evidence.clip_sites)
        self._writers[1].write(evidence.discordant_pairs)
        self._writers[2].write(evidence.split_reads)

    def __exit__(self, *_exc: object) -> None:
        for writer in self._writers:
            writer.close()
