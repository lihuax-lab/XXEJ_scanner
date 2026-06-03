"""Coverage enrichment based candidate region detection."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import pysam

from .bam_evidence import passes_read_filters
from .models import CandidateRegion, ScannerConfig
from .utils import clamp_start, percentile


@dataclass(slots=True)
class BinnedCoverage:
    # Retained as a named shape for future extensions; current code returns a
    # dictionary for faster overlap lookups.
    chrom: str
    start: int
    end: int
    coverage: float


def parse_bed_regions(path: str, *, id_prefix: str = "region") -> list[CandidateRegion]:
    # BED files provide candidate search intervals. A score column is preserved
    # if present, but event calling still requires breakpoint evidence later.
    regions: list[CandidateRegion] = []
    with Path(path).open() as handle:
        for idx, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(
                    f"BED line {idx} has fewer than 3 columns: {line.rstrip()}"
                )
            chrom, start, end = fields[:3]
            name = (
                fields[3]
                if len(fields) >= 4 and fields[3]
                else f"{id_prefix}_{len(regions) + 1}"
            )
            score = (
                float(fields[4])
                if len(fields) >= 5 and fields[4] not in {"", "."}
                else 0.0
            )
            regions.append(
                CandidateRegion(
                    chrom=chrom,
                    start=int(start),
                    end=int(end),
                    region_id=name,
                    score=score,
                )
            )
    return regions


def _iter_scan_windows(
    bam: pysam.AlignmentFile, regions: list[CandidateRegion] | None
) -> list[tuple[str, int, int]]:
    # Without a BED, scan every contig advertised by the BAM header. This avoids
    # hard-coded chromosome naming assumptions such as chr1/chr2.
    if regions:
        return [(region.chrom, region.start, region.end) for region in regions]
    return [(chrom, 0, length) for chrom, length in zip(bam.references, bam.lengths)]


# TODO: optimize algorithm
def compute_binned_coverage(
    bam_path: str,
    config: ScannerConfig,
    regions: list[CandidateRegion] | None = None,
) -> dict[tuple[str, int, int], float]:
    # Coverage is stored as average aligned bases per bin. For CUT&Tag this is
    # only used to define candidate search space, not to infer copy number.
    bin_size = config.coverage_bin_size
    coverage_bases: dict[tuple[str, int], int] = {}
    bin_ends: dict[tuple[str, int], int] = {}

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for chrom, scan_start, scan_end in _iter_scan_windows(bam, regions):
            for bin_start in range(
                scan_start - (scan_start % bin_size), scan_end, bin_size
            ):
                bin_start = clamp_start(bin_start)
                bin_end = min(bin_start + bin_size, scan_end)
                coverage_bases.setdefault((chrom, bin_start), 0)
                bin_ends[(chrom, bin_start)] = bin_end
            for read in bam.fetch(chrom, clamp_start(scan_start), scan_end):
                if not passes_read_filters(read, config):
                    continue
                read_start = max(int(read.reference_start), scan_start)
                read_end = min(int(read.reference_end), scan_end)
                if read_end <= read_start:
                    continue
                first_bin = read_start - (read_start % bin_size)
                for bin_start in range(first_bin, read_end, bin_size):
                    # Split each read's contribution across all overlapped bins
                    # so long reads or edge-crossing reads do not get assigned
                    # wholly to their start bin.
                    bin_start = clamp_start(bin_start)
                    bin_end = min(bin_start + bin_size, scan_end)
                    overlap = max(
                        0, min(read_end, bin_end) - max(read_start, bin_start)
                    )
                    if overlap:
                        coverage_bases[(chrom, bin_start)] = (
                            coverage_bases.get((chrom, bin_start), 0) + overlap
                        )
                        bin_ends[(chrom, bin_start)] = bin_end

    binned: dict[tuple[str, int, int], float] = {}
    for (chrom, start), bases in coverage_bases.items():
        end = bin_ends[(chrom, start)]
        width = max(1, end - start)
        binned[(chrom, start, end)] = bases / width
    return binned


def _mapped_read_count(bam_path: str) -> int:
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        try:
            return int(bam.mapped)
        except Exception:
            return 0


def _log2fc(
    treated: float, control: float, pseudo: float, control_scale: float = 1.0
) -> float:
    return math.log2((treated + pseudo) / (control * control_scale + pseudo))


def call_candidate_regions(
    treated_bam: str,
    config: ScannerConfig,
    control_bam: str | None = None,
) -> list[CandidateRegion]:
    treated_bins = compute_binned_coverage(treated_bam, config)
    control_bins = compute_binned_coverage(control_bam, config) if control_bam else {}
    if control_bam:
        coverage_threshold = config.min_treated_coverage
    # If no control is supplied, use a top-percentile floor to avoid reporting
    # every low-coverage bin. With a control, rely on min coverage + log2FC.
    else:
        dynamic_threshold = percentile(treated_bins.values(), config.top_percentile)
        coverage_threshold = max(config.min_treated_coverage, dynamic_threshold)

    control_scale = 1.0
    if control_bam:
        # A small library-size normalization keeps log2FC interpretable when the
        # treated and control BAMs have noticeably different mapped read counts.
        treated_mapped = _mapped_read_count(treated_bam)
        control_mapped = _mapped_read_count(control_bam)
        if treated_mapped > 0 and control_mapped > 0:
            control_scale = treated_mapped / control_mapped

    selected: list[CandidateRegion] = []
    for idx, ((chrom, start, end), treated_cov) in enumerate(
        sorted(treated_bins.items()), 1
    ):
        control_cov = control_bins.get((chrom, start, end), 0.0)
        log2fc = (
            _log2fc(treated_cov, control_cov, config.pseudo_count, control_scale)
            if control_bam
            else 0.0
        )
        if treated_cov < coverage_threshold:
            continue
        if control_bam and log2fc < config.min_log2fc:
            continue
        selected.append(
            CandidateRegion(
                chrom=chrom,
                start=start,
                end=end,
                region_id=f"candidate_bin_{idx}",
                score=treated_cov + max(0.0, log2fc),
                treated_coverage=treated_cov,
                control_coverage=control_cov,
                log2fc=log2fc,
            )
        )

    return merge_candidate_bins(selected, config.merge_distance)


def merge_candidate_bins(
    regions: list[CandidateRegion], merge_distance: int
) -> list[CandidateRegion]:
    if not regions:
        return []
    merged: list[CandidateRegion] = []
    sorted_regions = sorted(
        regions, key=lambda region: (region.chrom, region.start, region.end)
    )
    current = sorted_regions[0]
    members = [current]

    def flush(items: list[CandidateRegion]) -> CandidateRegion:
        # Merge summary values as width-weighted means. The region score remains
        # the best member score so a sharp peak is not diluted by nearby bins.
        chrom = items[0].chrom
        start = min(item.start for item in items)
        end = max(item.end for item in items)
        score = max(item.score for item in items)
        width = sum(max(1, item.end - item.start) for item in items)
        treated = (
            sum(item.treated_coverage * max(1, item.end - item.start) for item in items)
            / width
        )
        control = (
            sum(item.control_coverage * max(1, item.end - item.start) for item in items)
            / width
        )
        log2fc = (
            sum(item.log2fc * max(1, item.end - item.start) for item in items) / width
        )
        return CandidateRegion(
            chrom=chrom,
            start=start,
            end=end,
            region_id=f"region_{len(merged) + 1}",
            score=score,
            treated_coverage=treated,
            control_coverage=control,
            log2fc=log2fc,
        )

    for region in sorted_regions[1:]:
        if (
            region.chrom == current.chrom
            and region.start <= current.end + merge_distance
        ):
            members.append(region)
            if region.end > current.end:
                current = region
            continue
        merged.append(flush(members))
        current = region
        members = [region]
    merged.append(flush(members))
    return merged


def annotate_region_coverage(
    regions: list[CandidateRegion],
    treated_bam: str,
    config: ScannerConfig,
    control_bam: str | None = None,
) -> list[CandidateRegion]:
    # BED-driven scans skip de novo peak calling, but still annotate the same
    # coverage fields expected by candidate_regions.bed.
    treated_bins = compute_binned_coverage(treated_bam, config, regions)
    control_bins = (
        compute_binned_coverage(control_bam, config, regions) if control_bam else {}
    )

    annotated: list[CandidateRegion] = []
    for idx, region in enumerate(regions, 1):
        treated = _mean_overlap_coverage(region, treated_bins)
        control = _mean_overlap_coverage(region, control_bins) if control_bam else 0.0
        log2fc = _log2fc(treated, control, config.pseudo_count) if control_bam else 0.0
        annotated.append(
            CandidateRegion(
                chrom=region.chrom,
                start=region.start,
                end=region.end,
                region_id=region.region_id or f"region_{idx}",
                score=max(region.score, treated + max(0.0, log2fc)),
                treated_coverage=treated,
                control_coverage=control,
                log2fc=log2fc,
                clip_rate=region.clip_rate,
            )
        )
    return annotated


def _mean_overlap_coverage(
    region: CandidateRegion, bins: dict[tuple[str, int, int], float]
) -> float:
    weighted = 0.0
    width = 0
    for (chrom, start, end), cov in bins.items():
        if chrom != region.chrom:
            continue
        overlap = max(0, min(region.end, end) - max(region.start, start))
        if overlap:
            weighted += cov * overlap
            width += overlap
    return weighted / width if width else 0.0
