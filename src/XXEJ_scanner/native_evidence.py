"""Optional HTSlib-backed evidence extraction."""

from __future__ import annotations

from collections.abc import Iterator, Sequence

from .bam_evidence import collect_region_evidence
from .models import (
    CandidateRegion,
    CigarIndel,
    ClipSite,
    DiscordantPair,
    RegionEvidence,
    ScannerConfig,
    SplitReadEvidence,
)
from .utils import log

try:
    from . import _evidence_native
except ImportError:
    _evidence_native = None


def native_available() -> bool:
    return _evidence_native is not None


def selected_backend(config: ScannerConfig) -> str:
    return (
        "native"
        if config.evidence_backend != "python" and native_available()
        else "python"
    )


def _native_config(config: ScannerConfig) -> dict[str, object]:
    return {
        "breakpoint_quality_window": config.breakpoint_quality_window,
        "min_breakpoint_quality_fraction": config.min_breakpoint_quality_fraction,
        "min_clip_length": config.min_clip_length,
        "min_indel_length": config.min_indel_length,
        "min_aligned_length": config.min_aligned_length,
        "max_sa_nm": config.max_sa_nm,
        "discordant_min_distance": config.discordant_min_distance,
        "max_insert_size": config.max_insert_size,
        "merge_distance": config.merge_distance,
        "allow_duplicates": config.allow_duplicates,
        "include_supplementary": config.include_supplementary,
        "library_orientation": config.library_orientation,
    }


def _from_native(
    region: CandidateRegion, raw: dict[str, list[tuple[object, ...]]]
) -> RegionEvidence:
    return RegionEvidence(
        region=region,
        clip_sites=[ClipSite(*row) for row in raw["clip_sites"]],
        indels=[CigarIndel(*row) for row in raw["indels"]],
        discordant_pairs=[
            DiscordantPair(*row) for row in raw["discordant_pairs"]
        ],
        split_reads=[SplitReadEvidence(*row) for row in raw["split_reads"]],
    )


def collect_regions_native(
    bam_path: str,
    regions: Sequence[CandidateRegion],
    config: ScannerConfig,
    *,
    padding: int | None = None,
    min_mapq: int | None = None,
    min_baseq: int | None = None,
) -> list[RegionEvidence]:
    if _evidence_native is None:
        raise RuntimeError(
            "Native evidence backend is unavailable; install HTSlib and rebuild the package."
        )
    raw_regions = _evidence_native.collect_regions(
        bam_path=bam_path,
        regions=[(region.chrom, region.start, region.end) for region in regions],
        config=_native_config(config),
        padding=config.scan_padding if padding is None else padding,
        min_mapq=config.min_mapq if min_mapq is None else min_mapq,
        min_baseq=config.min_breakpoint_baseq if min_baseq is None else min_baseq,
    )
    return [
        _from_native(region, raw)
        for region, raw in zip(regions, raw_regions, strict=True)
    ]


def iter_region_evidence(
    regions: Sequence[CandidateRegion], config: ScannerConfig
) -> Iterator[tuple[CandidateRegion, RegionEvidence, RegionEvidence | None]]:
    use_native = selected_backend(config) == "native"
    if config.evidence_backend == "native" and not use_native:
        raise RuntimeError(
            "--evidence-backend native requested, but the HTSlib extension was not built."
        )
    if not use_native:
        if config.evidence_backend == "auto":
            log("[3/6] Native evidence backend unavailable; using pysam")
        for region in regions:
            treated = collect_region_evidence(config.treated_bam, region, config)
            control = (
                collect_region_evidence(config.control_bam, region, config)
                if config.control_bam
                else None
            )
            yield region, treated, control
        return

    log("[3/6] Using HTSlib C++ evidence backend")
    batch_size = config.evidence_batch_size
    for start in range(0, len(regions), batch_size):
        batch = regions[start : start + batch_size]
        treated_batch = collect_regions_native(config.treated_bam, batch, config)
        control_batch = (
            collect_regions_native(config.control_bam, batch, config)
            if config.control_bam
            else [None] * len(batch)
        )
        yield from zip(batch, treated_batch, control_batch, strict=True)
