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
--depth-count-method pileup
--cluster-method window
--microhomology-search-window 5
--include-supplementary
--min-nhej-ins-indel-support 1
--allow-clip-only-nhej-ins
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

For `MMEJ_DEL`, the original `start`, `end`, and `deleted_length` columns keep the nominal paired-cluster interval. Additional `microhomology_*` columns report the breakpoint-adjusted microhomology placement, the implied deletion span when one microhomology copy is collapsed, equivalent-placement count, and low-complexity status. `junction_evidence_support` and `junction_evidence_types` summarize read-level evidence that is consistent with the local deletion junction, such as matched CIGAR deletions, same-chromosome SA-tag split reads, or soft clips that match the opposite flank.

Depth columns in `breakpoint_clusters.tsv` and `events.tsv` report local depth around the clustered breakpoint window, not depth across the full candidate region. With `--depth-count-method pileup`, the scanner reports mean base-level pileup depth across that local window. With `--depth-count-method region`, it reports the number of unique read names overlapping the local window, which is similar to the original molecule/read support count.

Breakpoint clustering defaults to `--cluster-method window`, which groups soft-clipped read ends by genomic proximity only. The experimental `--cluster-method evidence-graph` mode builds a lightweight weighted graph over nearby clipped observations and uses deterministic community detection to link clips that are supported by local CIGAR indel, SA split-read, or discordant-pair evidence. It still emits the same `breakpoint_clusters.tsv` schema, so runs can be compared directly against the default window method.

By default, discovery filters use mapped, primary, non-secondary alignments that pass MAPQ and aligned-length thresholds. Duplicate reads are excluded unless `--allow-duplicates` is set. Supplementary alignments are excluded unless `--include-supplementary` is set; SA tags on primary alignments are still parsed for split-read evidence.

When no control BAM is supplied, de novo candidate regions are selected from the treated sample's high-coverage bins. For broad no-control scans, consider stricter discovery thresholds such as `--top-percentile 99`, `--min-treated-coverage 20`, and `--min-alt-support 5`. Avoid `--include-supplementary` unless supplementary alignments are specifically needed, because SA tags on primary alignments are already parsed for split-read evidence.

## Event Types

`NHEJ_INS` marks a candidate local end-joining event with small inserted or filler sequence evidence near a clipped breakpoint. By default, it requires at least one nearby CIGAR insertion read (`--min-nhej-ins-indel-support 1`) so clip-only breakpoint clusters are not automatically reported as local insertions. Use `--allow-clip-only-nhej-ins` to restore clip-only filler-sequence calls.

`MMEJ_DEL` marks a candidate local deletion-like event with paired local breakpoint evidence and reference microhomology context. It is intended as a reviewable MMEJ-like candidate call, not a proof that every supporting read directly resolves the repaired junction. The left breakpoint must be supported by right-clipped reads, or by a mixed-side cluster, and the right breakpoint must be supported by left-clipped reads, or by a mixed-side cluster. This prevents a short doubly clipped alignment from being misreported as a deletion spanning its aligned seed.

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
