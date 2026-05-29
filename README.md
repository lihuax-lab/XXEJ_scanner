# XXEJ_scanner

`XXEJ_scanner` is a first-version, rule-based scanner for LIG3/LIG4 CUT&Tag-enriched paired-end BAM files. It reports candidate DSB repair-associated events rather than definitive structural variants.

The pipeline keeps three concepts separate:

1. CUT&Tag-enriched candidate regions from coverage or user BED input.
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
  --coverage-bin-size 100 \
  --merge-distance 300 \
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

## Event Types

`NHEJ_INS` marks a candidate local end-joining event with small inserted or filler sequence evidence near a clipped breakpoint.

`MMEJ_DEL` marks a candidate local deletion event with paired local breakpoint evidence and reference microhomology context.

`NHEJ_BND_INS_INTRA` and `NHEJ_BND_INS_INTER` mark candidate wrong-end joining supported by local clipping plus distant same-chromosome or inter-chromosome split/paired evidence.

## repair_evidence_fraction

`repair_evidence_fraction` is:

```text
ALT-like support / (ALT-like support + REF-like support)
```

ALT-like support includes breakpoint-associated clipping, SA-tag split reads, discordant pairs, and local CIGAR insertion/deletion evidence. REF-like support includes reads spanning the candidate breakpoint or interval without strong clipping and with locally consistent alignment.

Because this is CUT&Tag-enriched data, this value is an enrichment-data-derived evidence fraction. It should not be interpreted as a true WGS allele fraction or VAF.

## Difference From WGS SV Callers

The scanner does not assume uniform whole-genome coverage and does not call repair events from coverage alone. Coverage only defines candidate enriched regions. Breakpoint and event calls require local clipped, indel, split-read, or paired-end evidence.

## Difference From MEIGA-SR

MEIGA-SR uses REF/ALT logic for mobile element insertion genotyping. Here, `REF-like` means intact local spanning evidence around a candidate DSB breakpoint, while `ALT-like` means evidence supporting a repair-associated abnormal structure.
