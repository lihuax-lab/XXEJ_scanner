"""Public Python API for the CUT&Tag-enriched XXEJ repair scanner."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

from .bam_evidence import (
    collect_region_evidence,
    count_depth,
    count_spanning_reads,
    is_discordant_pair,
)
from .breakpoints import (
    cluster_clip_sites,
    cluster_evidence_graph,
    filter_breakpoint_clusters,
    score_breakpoint_cluster,
)
from .classify import (
    assign_final_event_ids,
    classify_bnd_events,
    classify_event_type,
    classify_local_events,
)
from .cli import run_scan
from .coverage import (
    annotate_region_coverage,
    call_candidate_regions,
    parse_bed_regions,
)
from .genotype import (
    compute_repair_evidence_fraction,
    count_ref_like_reads,
    update_event_fraction,
)
from .io import (
    prepare_output_dir,
    write_breakpoint_clusters_tsv,
    write_candidate_regions_bed,
    write_event_evidence_tsv,
    write_events_tsv,
    write_igv_loci_bed,
    write_raw_clip_sites_tsv,
    write_raw_discordant_pairs_tsv,
    write_raw_split_reads_tsv,
    write_run_summary_json,
)
from .models import (
    BreakpointCluster,
    CandidateRegion,
    CigarIndel,
    ClipSide,
    ClipSite,
    ClusterMethod,
    DiscordantPair,
    EventEvidence,
    EventType,
    RegionEvidence,
    RepairEvent,
    ScannerConfig,
    SplitReadEvidence,
)
from .reference import (
    MicrohomologyHit,
    ReferenceGenome,
    check_clipped_sequence_against_reference,
    detect_microhomology,
    fetch_reference_sequence,
    find_microhomology,
)
from .validation import assign_event_filter, second_pass_validate_event

__all__ = [
    "__version__",
    "annotate_region_coverage",
    "assign_event_filter",
    "assign_final_event_ids",
    "BreakpointCluster",
    "call_candidate_regions",
    "CandidateRegion",
    "check_clipped_sequence_against_reference",
    "CigarIndel",
    "classify_bnd_events",
    "classify_event_type",
    "classify_local_events",
    "ClipSide",
    "ClipSite",
    "cluster_clip_sites",
    "cluster_evidence_graph",
    "ClusterMethod",
    "collect_region_evidence",
    "compute_repair_evidence_fraction",
    "count_depth",
    "count_ref_like_reads",
    "count_spanning_reads",
    "detect_microhomology",
    "DiscordantPair",
    "EventEvidence",
    "EventType",
    "fetch_reference_sequence",
    "filter_breakpoint_clusters",
    "find_microhomology",
    "is_discordant_pair",
    "MicrohomologyHit",
    "parse_bed_regions",
    "prepare_output_dir",
    "ReferenceGenome",
    "RegionEvidence",
    "RepairEvent",
    "run_scan",
    "ScannerConfig",
    "score_breakpoint_cluster",
    "second_pass_validate_event",
    "SplitReadEvidence",
    "update_event_fraction",
    "write_breakpoint_clusters_tsv",
    "write_candidate_regions_bed",
    "write_event_evidence_tsv",
    "write_events_tsv",
    "write_igv_loci_bed",
    "write_raw_clip_sites_tsv",
    "write_raw_discordant_pairs_tsv",
    "write_raw_split_reads_tsv",
    "write_run_summary_json",
]

try:
    __version__ = version("xxej-scanner")
except PackageNotFoundError:
    __version__ = "0.1.0"
