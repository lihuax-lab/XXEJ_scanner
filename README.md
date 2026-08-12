# XXEJ_scanner

`XXEJ_scanner` is a first-version, rule-based scanner for LIG3/LIG4 CUT&Tag-enriched paired-end BAM files. It reports candidate DSB repair-associated events rather than definitive structural variants.

The pipeline keeps three concepts separate:

1. Candidate regions from CUT&Tag coverage, user BED input, or strong structural evidence.
2. Breakpoint clusters from local soft-clipped read ends.
3. Candidate repair events supported by clipped reads, CIGAR indels, SA-tag split reads, and discordant pairs.

## Install

This repository is a `uv` project. From the project root:

```bash
uv sync --no-editable
```

Run the scanner through the installed console script:

```bash
uv run --no-editable XXEJ_scanner --help
uv run --no-editable XXEJ_scanner scan --help
```

The lowercase alias is also available:

```bash
uv run --no-editable xxej-scanner scan --help
```

## Python API

The package exposes a small public API from the package root:

```python
from XXEJ_scanner import ScannerConfig, run_scan

config = ScannerConfig(
    treated_bam="sample.sorted.bam",
    reference_fasta="genome.fa",
    output_dir="results/XXEJ_scanner",
    peak_bed="peaks/union.narrow.bed",
)
summary = run_scan(config)
```

Common programmatic entry points include `parse_bed_regions`,
`call_candidate_regions`, `call_structural_evidence_regions`,
`collect_region_evidence`, `cluster_clip_sites`,
`cluster_evidence_graph`, `classify_local_events`, `classify_bnd_events`, and
the TSV writer helpers such as `write_events_tsv`.

## Inputs

Required:

```bash
--treated-bam sample.sorted.bam
--reference-fasta genome.fa
--output-dir results/XXEJ_scanner
```

Optional:

```bash
--control-bam control.sorted.bam
--candidate-bed candidates.bed
--peak-bed peaks.bed
--sample-name treated
--control-name control
--depth-count-method pileup
--cluster-method window
--max-sa-nm 10
--max-control-alt-support 1
--microhomology-search-window 5
--include-supplementary
```

BAM and FASTA indexes are required. Missing `.bai`/`.csi` or `.fai` files are reported as errors.

## Example

```bash
uv run --no-editable XXEJ_scanner scan \
  --treated-bam bam/090-ETO-L189.sorted.bam \
  --control-bam bam/090-control.sorted.bam \
  --reference-fasta ref/genome.fa \
  --output-dir results/XXEJ_scanner \
  --min-mapq 20 \
  --min-clip-length 10 \
  --clip-cluster-window 20 \
  --cluster-method window \
  --coverage-bin-size 100 \
  --merge-distance 300 \
  --depth-count-method pileup \
  --min-alt-support 3
```

With a candidate BED:

```bash
uv run --no-editable XXEJ_scanner scan \
  --treated-bam bam/090-ETO-L189.sorted.bam \
  --control-bam bam/090-control.sorted.bam \
  --reference-fasta ref/genome.fa \
  --candidate-bed peaks/union_candidates.bed \
  --output-dir results/XXEJ_scanner
```

## Synthetic BAM Benchmark

For a first-pass correctness check, the `scripts/` directory contains a small
BAM-level simulator. It does not simulate FASTQ or run an aligner. Instead, it
writes a synthetic reference, candidate BED, treated/control BAMs, and
`truth.tsv` with 20 clean structural alleles: seven `LOCAL_INS`, seven
`LOCAL_DEL`, and six `BND_INTER` events. Four negative regions exercise
clip-only noise, a control-shared insertion, duplicate reads, and low-MAPQ reads.

Generate the benchmark:

```bash
uv run --no-editable python scripts/simulate_xxej_benchmark.py --force
```

Run the scanner and evaluator:

```bash
bash sim/xxej_benchmark_v1/run_scanner.sh
```

Or run the steps manually:

```bash
uv run --no-editable XXEJ_scanner scan \
  --treated-bam sim/xxej_benchmark_v1/treated.sorted.bam \
  --control-bam sim/xxej_benchmark_v1/control.sorted.bam \
  --reference-fasta sim/xxej_benchmark_v1/ref.fa \
  --output-dir sim/xxej_benchmark_v1/scanner_out \
  --min-alt-support 3 \
  --min-bnd-support 3 \
  --depth-count-method region

uv run --no-editable python scripts/evaluate_xxej_benchmark.py \
  --truth sim/xxej_benchmark_v1/truth.tsv \
  --events sim/xxej_benchmark_v1/scanner_out/events.tsv
```

## Outputs

The output directory contains:

```text
candidate_regions.bed
breakpoint_clusters.tsv
events.tsv
event_evidence.tsv
run_summary.json
raw_clip_sites.tsv
raw_discordant_pairs.tsv
raw_split_reads.tsv
igv_loci.bed
```

Coordinates are 0-based half-open for BED-like intervals. Breakpoint positions are reported as 0-based reference positions derived from alignment starts for left clips and alignment ends for right clips.

For `LOCAL_DEL`, `start`, `end`, and `deleted_length` describe the resolved junction. Additional `microhomology_*` columns annotate exact breakpoint sequence context, equivalent placements, and low-complexity status. Microhomology is annotation, not a repair-pathway assignment. `junction_evidence_support` and `junction_evidence_types` summarize read-level evidence such as CIGAR deletions, same-chromosome SA-tag split reads, or soft clips matching the opposite flank.

Depth columns in `breakpoint_clusters.tsv` and `events.tsv` report local depth around the clustered breakpoint window, not depth across the full candidate region. With `--depth-count-method pileup`, the scanner reports mean base-level pileup depth across that local window. With `--depth-count-method region`, it reports the number of unique read names overlapping the local window, which is similar to the original molecule/read support count.

Breakpoint clustering defaults to `--cluster-method window`, which groups soft-clipped read ends by genomic proximity only. The experimental `--cluster-method evidence-graph` mode builds a lightweight weighted graph over nearby clipped observations and uses deterministic community detection to link clips that are supported by local CIGAR indel, SA split-read, or discordant-pair evidence. It still emits the same `breakpoint_clusters.tsv` schema, so runs can be compared directly against the default window method.

By default, discovery filters use mapped, primary, non-secondary alignments that pass MAPQ and aligned-length thresholds. Duplicate reads are excluded unless `--allow-duplicates` is set. Supplementary alignments are excluded unless `--include-supplementary` is set; SA tags on primary alignments are still parsed for split-read evidence.

Candidate discovery always unions coverage/BED intervals with genome-wide strong structural evidence. Exact CIGAR indel alleles and bounded SA/discordant-pair clusters must meet their event support threshold; soft-clip-only evidence remains limited to coverage/BED intervals. This prevents low-coverage resolved junctions from being discarded while avoiding a genome-wide expansion driven only by clipping noise.

When no control BAM is supplied, coverage candidates use the treated sample's high-coverage bins. For broad no-control scans, consider stricter thresholds such as `--top-percentile 99`, `--min-treated-coverage 20`, and `--min-alt-support 5`. Avoid `--include-supplementary` unless supplementary records are specifically needed, because SA tags on primary alignments are already parsed.

## Event Types

`LOCAL_INS` reports an exact insertion position and sequence supported by CIGAR insertion reads. The label describes the observed allele and does not infer NHEJ, LIG3, or LIG4 causality.

`LOCAL_DEL` reports exact deletion endpoints supported by CIGAR deletion or compatible same-chromosome split reads. Any microhomology is reported separately as sequence context.

`BND_INTRA` and `BND_INTER` report distant same-chromosome or inter-chromosome breakends. `evidence_level` distinguishes resolved split-read junctions, multi-signal calls, and pair-only candidates.

## repair_evidence_fraction

`repair_evidence_fraction` is:

```text
ALT-like support / (ALT-like support + REF-like support)
```

ALT-like support is counted by unique read name across the evidence types assigned to the same structural allele. REF-like support is measured independently at each required breakend; the reported denominator uses the weaker covered endpoint for two-breakpoint events.

Because this is CUT&Tag-enriched data, this value is an enrichment-data-derived evidence fraction. It should not be interpreted as a true WGS allele fraction or VAF.

## Difference From WGS SV Callers

The scanner does not assume uniform whole-genome coverage and does not call repair events from coverage alone. Coverage, BED hints, and supported structural alignments define the search space; event calls still require allele-specific indel, split-read, or paired-end evidence.

## Difference From MEIGA-SR

MEIGA-SR uses REF/ALT logic for mobile element insertion genotyping. Here, `REF-like` means intact local spanning evidence around a candidate DSB breakpoint, while `ALT-like` means evidence supporting a repair-associated abnormal structure.
