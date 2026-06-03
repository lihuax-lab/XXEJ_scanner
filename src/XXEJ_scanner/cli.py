"""Command-line interface for XXEJ_scanner."""

from __future__ import annotations

import argparse
from pathlib import Path

from .bam_evidence import clip_site_count_near, collect_region_evidence, count_depth
from .breakpoints import (
    cluster_clip_sites,
    filter_breakpoint_clusters,
    score_breakpoint_cluster,
)
from .classify import classify_local_events
from .coverage import (
    annotate_region_coverage,
    call_candidate_regions,
    parse_bed_regions,
)
from .genotype import count_ref_like_reads, update_event_fraction
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
    EventEvidence,
    RepairEvent,
    ScannerConfig,
)
from .reference import ReferenceGenome
from .utils import log, validate_inputs
from .validation import assign_event_filter, second_pass_validate_event


def build_parser() -> argparse.ArgumentParser:
    # Keep all first-version thresholds exposed here so runs are reproducible
    # from the command line and downstream notebooks/scripts do not need hidden
    # defaults.
    parser = argparse.ArgumentParser(prog="XXEJ_scanner")
    subparsers = parser.add_subparsers(dest="command", required=True)
    scan = subparsers.add_parser(
        "scan", help="Scan a paired-end CUT&Tag BAM for candidate XXEJ repair outcomes"
    )
    scan.add_argument("--treated-bam", "--tumor-bam", dest="treated_bam", required=True)
    scan.add_argument("--control-bam", default=None)
    scan.add_argument("--reference-fasta", required=True)
    scan.add_argument("--output-dir", required=True)
    scan.add_argument("--candidate-bed", default=None)
    scan.add_argument("--peak-bed", default=None)
    scan.add_argument("--sample-name", default="treated")
    scan.add_argument("--control-name", default="control")
    scan.add_argument("--min-mapq", type=int, default=20)
    scan.add_argument("--strict-min-mapq", type=int, default=30)
    scan.add_argument("--min-clip-length", type=int, default=10)
    scan.add_argument("--clip-cluster-window", type=int, default=20)
    scan.add_argument("--coverage-bin-size", type=int, default=100)
    scan.add_argument("--merge-distance", type=int, default=300)
    scan.add_argument("--max-normal-clip-rate", type=float, default=0.05)
    scan.add_argument("--min-alt-support", type=int, default=3)
    scan.add_argument("--min-bnd-support", type=int, default=3)
    scan.add_argument("--min-treated-coverage", type=float, default=5.0)
    scan.add_argument("--min-log2fc", type=float, default=1.0)
    scan.add_argument("--top-percentile", type=float, default=95.0)
    scan.add_argument("--pseudo-count", type=float, default=1.0)
    scan.add_argument("--max-insert-size", type=int, default=1000)
    scan.add_argument("--discordant-min-distance", type=int, default=1000)
    scan.add_argument(
        "--library-orientation", choices=["fr", "rf", "ff", "rr", "any"], default="fr"
    )
    scan.add_argument("--allow-duplicates", action="store_true")
    scan.add_argument("--min-aligned-length", type=int, default=20)
    scan.add_argument("--scan-padding", type=int, default=200)
    scan.add_argument("--max-local-event-distance", type=int, default=10000)
    scan.add_argument("--max-insertion-length", type=int, default=50)
    scan.add_argument("--min-indel-length", type=int, default=1)
    scan.add_argument("--min-microhomology-length", type=int, default=1)
    scan.add_argument("--max-microhomology-length", type=int, default=20)
    scan.add_argument("--second-pass-window", type=int, default=150)
    return parser


def _config_from_args(args: argparse.Namespace) -> ScannerConfig:
    # Convert argparse's loose namespace into a typed config object used by the
    # rest of the pipeline. This keeps non-CLI code independent of argparse.
    return ScannerConfig(
        treated_bam=args.treated_bam,
        control_bam=args.control_bam,
        reference_fasta=args.reference_fasta,
        output_dir=args.output_dir,
        candidate_bed=args.candidate_bed,
        peak_bed=args.peak_bed,
        sample_name=args.sample_name,
        control_name=args.control_name,
        min_mapq=args.min_mapq,
        strict_min_mapq=args.strict_min_mapq,
        min_clip_length=args.min_clip_length,
        clip_cluster_window=args.clip_cluster_window,
        coverage_bin_size=args.coverage_bin_size,
        merge_distance=args.merge_distance,
        max_normal_clip_rate=args.max_normal_clip_rate,
        min_alt_support=args.min_alt_support,
        min_bnd_support=args.min_bnd_support,
        min_treated_coverage=args.min_treated_coverage,
        min_log2fc=args.min_log2fc,
        top_percentile=args.top_percentile,
        pseudo_count=args.pseudo_count,
        max_insert_size=args.max_insert_size,
        discordant_min_distance=args.discordant_min_distance,
        library_orientation=args.library_orientation,
        allow_duplicates=args.allow_duplicates,
        min_aligned_length=args.min_aligned_length,
        scan_padding=args.scan_padding,
        max_local_event_distance=args.max_local_event_distance,
        max_insertion_length=args.max_insertion_length,
        min_indel_length=args.min_indel_length,
        min_microhomology_length=args.min_microhomology_length,
        max_microhomology_length=args.max_microhomology_length,
        second_pass_window=args.second_pass_window,
    )


def run_scan(config: ScannerConfig) -> dict[str, object]:
    log("[1/6] Loading configuration")
    validate_inputs(config)
    output_paths = prepare_output_dir(config.output_dir)

    log("[2/6] Calling candidate regions")
    # BED input is treated as a search-space hint, not as a called repair event.
    # Coverage statistics are still annotated so candidate_regions.bed remains
    # comparable between BED-driven and de novo scans.
    if config.candidate_bed or config.peak_bed:
        input_bed = config.candidate_bed or config.peak_bed
        assert input_bed is not None
        regions = parse_bed_regions(input_bed)
        regions = annotate_region_coverage(
            regions, config.treated_bam, config, config.control_bam
        )
    else:
        regions = call_candidate_regions(config.treated_bam, config, config.control_bam)

    all_clusters: list[BreakpointCluster] = []
    all_events: list[RepairEvent] = []
    all_event_evidence: list[EventEvidence] = []
    raw_clip_sites = []
    raw_discordant_pairs = []
    raw_split_reads = []

    log("[3/6] Extracting clipped reads")
    with ReferenceGenome(config.reference_fasta) as reference:
        for region in regions:
            # Evidence is collected per candidate region to avoid assuming WGS-
            # like uniform coverage. Each region can have its own local depth,
            # control noise, and breakpoint structure.
            treated_evidence = collect_region_evidence(
                config.treated_bam, region, config
            )
            control_evidence = (
                collect_region_evidence(config.control_bam, region, config)
                if config.control_bam
                else None
            )
            raw_clip_sites.extend(treated_evidence.clip_sites)
            raw_discordant_pairs.extend(treated_evidence.discordant_pairs)
            raw_split_reads.extend(treated_evidence.split_reads)

            log("[4/6] Clustering breakpoints")
            clusters = cluster_clip_sites(
                treated_evidence.clip_sites, config, region=region
            )
            for cluster in clusters:
                # Depth is counted around the cluster, not across the whole
                # enriched region, because CUT&Tag peaks can be highly uneven.
                cluster.treated_depth = count_depth(
                    config.treated_bam,
                    cluster.chrom,
                    cluster.cluster_start - config.clip_cluster_window,
                    cluster.cluster_end + config.clip_cluster_window,
                    config,
                )
                if config.control_bam and control_evidence:
                    cluster.control_depth = count_depth(
                        config.control_bam,
                        cluster.chrom,
                        cluster.cluster_start - config.clip_cluster_window,
                        cluster.cluster_end + config.clip_cluster_window,
                        config,
                    )
                    control_clips = clip_site_count_near(
                        control_evidence.clip_sites,
                        cluster.chrom,
                        cluster.peak_pos,
                        config.clip_cluster_window,
                    )
                    # normal_noise is a local clip rate in the control sample.
                    # It is used as an artifact penalty, not as a hard proof
                    # that the treated signal is false.
                    cluster.normal_noise = (
                        control_clips / cluster.control_depth
                        if cluster.control_depth
                        else 0.0
                    )
                score_breakpoint_cluster(cluster, config)
            kept_clusters = filter_breakpoint_clusters(clusters, config)
            # Region clip_rate is a debugging summary for the BED output; final
            # event calls still require breakpoint and repair-product evidence.
            region.clip_rate = (
                sum(cluster.clip_count for cluster in kept_clusters)
                / max(1, sum(cluster.treated_depth for cluster in kept_clusters))
                if kept_clusters
                else 0.0
            )
            all_clusters.extend(kept_clusters)

            log("[5/6] Classifying repair events")
            events, event_evidence = classify_local_events(
                region, kept_clusters, treated_evidence, reference, config
            )
            for event in events:
                # REF-like counting is deliberately separate from event
                # classification: event evidence asks "is there abnormal
                # structure?", while ref_spanning_support asks "how much intact
                # local sequence is still visible around this breakpoint?".
                event.ref_spanning_support = count_ref_like_reads(
                    config.treated_bam, event, config
                )
                if config.control_bam:
                    event.control_ref_support = count_ref_like_reads(
                        config.control_bam, event, config
                    )
                    event.control_alt_support = (
                        _control_alt_support(event, control_evidence, config)
                        if control_evidence
                        else 0
                    )
                update_event_fraction(event)
                event.filter = assign_event_filter(event, config)
                # The second pass re-queries a narrower interval with stricter
                # MAPQ. It is a lightweight validation step, not local assembly.
                event = second_pass_validate_event(event, config.treated_bam, config)
                all_events.append(event)
            all_event_evidence.extend(event_evidence)

    log("[6/6] Writing outputs")
    write_candidate_regions_bed(output_paths["candidate_regions"], regions)
    write_breakpoint_clusters_tsv(output_paths["breakpoint_clusters"], all_clusters)
    write_events_tsv(output_paths["events"], all_events)
    write_event_evidence_tsv(output_paths["event_evidence"], all_event_evidence)
    write_raw_clip_sites_tsv(output_paths["raw_clip_sites"], raw_clip_sites)
    write_raw_discordant_pairs_tsv(
        output_paths["raw_discordant_pairs"], raw_discordant_pairs
    )
    write_raw_split_reads_tsv(output_paths["raw_split_reads"], raw_split_reads)
    write_igv_loci_bed(output_paths["igv_loci"], all_events)
    summary = {
        "sample_name": config.sample_name,
        "control_name": config.control_name if config.control_bam else None,
        "treated_bam": str(Path(config.treated_bam)),
        "control_bam": str(Path(config.control_bam)) if config.control_bam else None,
        "reference_fasta": str(Path(config.reference_fasta)),
        "candidate_regions": len(regions),
        "breakpoint_clusters": len(all_clusters),
        "events": len(all_events),
        "pass_events": sum(1 for event in all_events if event.filter == "PASS"),
        "repair_evidence_fraction_note": (
            "repair_evidence_fraction is ALT-like support divided by ALT-like plus REF-like support "
            "in CUT&Tag-enriched data; it is not a true WGS allele fraction."
        ),
    }
    write_run_summary_json(output_paths["run_summary"], summary)
    return summary


def _control_alt_support(
    event: RepairEvent, control_evidence: object, config: ScannerConfig
) -> int:
    # Mirror the treated ALT-like evidence definition in the control sample so
    # recurrent mapping artifacts can be downweighted in filters and summaries.
    support = set()
    for site in control_evidence.clip_sites:
        if (
            site.chrom == event.chrom
            and abs(site.pos - int(event.bkp_A_pos)) <= config.clip_cluster_window
        ):
            support.add(site.read_name)
    for indel in control_evidence.indels:
        if (
            indel.chrom == event.chrom
            and int(event.start) - config.clip_cluster_window
            <= indel.start
            <= int(event.end) + config.clip_cluster_window
        ):
            support.add(indel.read_name)
    if event.event_type.startswith("NHEJ_BND"):
        for pair in control_evidence.discordant_pairs:
            if (
                pair.chrom == event.chrom
                and abs(pair.pos - int(event.bkp_A_pos)) <= config.clip_cluster_window
            ):
                support.add(pair.read_name)
        for split in control_evidence.split_reads:
            if (
                split.chrom == event.chrom
                and abs(split.pos - int(event.bkp_A_pos)) <= config.clip_cluster_window
            ):
                support.add(split.read_name)
    return len(support)


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "scan":
        run_scan(_config_from_args(args))
        return 0
    parser.error(f"Unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
