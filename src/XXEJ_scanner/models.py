"""Shared dataclasses for the XXEJ scanner."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

ClipSide = Literal["left_clip", "right_clip"]
EventType = Literal[
    "NHEJ_INS",
    "MMEJ_DEL",
    "NHEJ_BND_INS_INTRA",
    "NHEJ_BND_INS_INTER",
]


@dataclass(slots=True)
class ScannerConfig:
    # Runtime parameters collected from the CLI. Keeping thresholds here makes it
    # straightforward to pass the same configuration across all pipeline stages.
    treated_bam: str
    reference_fasta: str
    output_dir: str
    control_bam: str | None = None
    candidate_bed: str | None = None
    peak_bed: str | None = None
    sample_name: str = "treated"
    control_name: str = "control"
    min_mapq: int = 20
    strict_min_mapq: int = 30
    min_clip_length: int = 10
    clip_cluster_window: int = 20
    coverage_bin_size: int = 100
    merge_distance: int = 300
    max_normal_clip_rate: float = 0.05
    min_alt_support: int = 3
    min_bnd_support: int = 3
    min_treated_coverage: float = 5.0
    min_log2fc: float = 1.0
    top_percentile: float = 95.0
    pseudo_count: float = 1.0
    max_insert_size: int = 1000
    discordant_min_distance: int = 1000
    library_orientation: str = "fr"
    allow_duplicates: bool = False
    include_supplementary: bool = False
    min_aligned_length: int = 20
    scan_padding: int = 200
    max_local_event_distance: int = 10000
    max_insertion_length: int = 50
    min_indel_length: int = 1
    min_microhomology_length: int = 1
    max_microhomology_length: int = 20
    microhomology_search_window: int = 5
    second_pass_window: int = 150
    depth_count_method: str = "pileup"


@dataclass(slots=True)
class CandidateRegion:
    # A coverage-enriched or user-provided search interval. This is not itself a
    # repair call; it only scopes evidence extraction.
    chrom: str
    start: int
    end: int
    region_id: str
    score: float = 0.0
    treated_coverage: float = 0.0
    control_coverage: float = 0.0
    log2fc: float = 0.0
    clip_rate: float = 0.0


@dataclass(slots=True)
class ClipSite:
    # One clipped read end converted into a breakpoint-side observation.
    # Coordinates are 0-based; right clips use the half-open alignment end.
    chrom: str
    pos: int
    side: ClipSide
    clip_length: int
    clip_sequence: str
    read_name: str
    strand: str
    mapq: int
    cigar: str
    is_reverse: bool
    reference_start: int
    reference_end: int
    clip_type: str = "S"


@dataclass(slots=True)
class CigarIndel:
    chrom: str
    start: int
    end: int
    operation: Literal["INS", "DEL"]
    length: int
    sequence: str
    read_name: str
    mapq: int
    cigar: str


@dataclass(slots=True)
class DiscordantPair:
    read_name: str
    chrom: str
    pos: int
    mate_chrom: str
    mate_pos: int
    orientation: str
    mapq: int
    is_reverse: bool
    mate_is_reverse: bool
    cigar: str
    reason: str


@dataclass(slots=True)
class SplitReadEvidence:
    read_name: str
    chrom: str
    pos: int
    side: str
    remote_chrom: str
    remote_pos: int
    remote_strand: str
    remote_cigar: str
    remote_mapq: int
    orientation: str
    mapq: int
    cigar: str
    sa_tag: str


@dataclass(slots=True)
class RegionEvidence:
    # All discovery evidence extracted from one candidate region and one BAM.
    region: CandidateRegion
    clip_sites: list[ClipSite] = field(default_factory=list)
    indels: list[CigarIndel] = field(default_factory=list)
    discordant_pairs: list[DiscordantPair] = field(default_factory=list)
    split_reads: list[SplitReadEvidence] = field(default_factory=list)


@dataclass(slots=True)
class BreakpointCluster:
    # A local pileup of clipped read ends, summarized before event classification.
    region_id: str
    chrom: str
    cluster_start: int
    cluster_end: int
    peak_pos: int
    clip_side: str
    clip_count: int
    left_clip_count: int
    right_clip_count: int
    treated_depth: float = 0.0
    control_depth: float = 0.0
    clip_rate: float = 0.0
    strand_plus_count: int = 0
    strand_minus_count: int = 0
    normal_noise: float = 0.0
    confidence: float = 0.0


@dataclass(slots=True)
class RepairEvent:
    # Final candidate event record. It preserves separate evidence counts so
    # downstream users can filter by their preferred evidence combination.
    event_id: str
    event_type: EventType
    chrom: str
    start: int
    end: int
    bkp_A_chrom: str
    bkp_A_pos: int
    bkp_A_side: str
    bkp_B_chrom: str = "NA"
    bkp_B_pos: int | str = "NA"
    bkp_B_side: str = "NA"
    remote_chrom: str = "NA"
    remote_pos: int | str = "NA"
    orientation: str = "NA"
    inserted_sequence: str = "NA"
    inserted_length: int | str = "NA"
    deleted_length: int | str = "NA"
    microhomology: str = "NA"
    microhomology_length: int = 0
    alt_clip_support: int = 0
    alt_split_support: int = 0
    alt_discordant_pair_support: int = 0
    alt_indel_support: int = 0
    ref_spanning_support: int = 0
    treated_depth: int = 0
    control_depth: int = 0
    repair_evidence_fraction: float = 0.0
    control_alt_support: int = 0
    control_ref_support: int = 0
    control_repair_evidence_fraction: float = 0.0
    score: float = 0.0
    filter: str = "NA"
    notes: str = ""
    normal_noise: float = 0.0
    microhomology_left_end: int | str = "NA"
    microhomology_right_start: int | str = "NA"
    microhomology_offset_a: int | str = "NA"
    microhomology_offset_b: int | str = "NA"
    microhomology_deletion_start: int | str = "NA"
    microhomology_deletion_end: int | str = "NA"
    microhomology_deletion_length: int | str = "NA"
    microhomology_ambiguity_bases: int = 0
    microhomology_equivalent_hits: int = 0
    microhomology_low_complexity: bool = False
    junction_evidence_support: int = 0
    junction_evidence_types: set[str] = field(default_factory=set)
    support_read_names: set[str] = field(default_factory=set, repr=False)

    @property
    def alt_support(self) -> int:
        # Prefer unique read-name support when available so reads contributing
        # multiple evidence classes do not inflate the evidence fraction.
        if self.support_read_names:
            return len(self.support_read_names)
        return (
            self.alt_clip_support
            + self.alt_split_support
            + self.alt_discordant_pair_support
            + self.alt_indel_support
        )


@dataclass(slots=True)
class EventEvidence:
    # Debuggable read-level support row for event_evidence.tsv.
    event_id: str
    read_name: str
    evidence_type: str
    chrom: str
    pos: int
    mate_chrom: str = "NA"
    mate_pos: int | str = "NA"
    cigar: str = "NA"
    mapq: int = 0
    is_reverse: bool = False
    mate_is_reverse: bool = False
    clip_side: str = "NA"
    clip_length: int | str = "NA"
    clip_sequence: str = "NA"
    sa_tag: str = "NA"
    classification: str = "ALT-like"
