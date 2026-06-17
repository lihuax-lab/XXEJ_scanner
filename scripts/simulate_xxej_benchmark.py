#!/usr/bin/env python3
"""Create a BAM-level benchmark for XXEJ_scanner.

The benchmark intentionally skips FASTQ simulation and alignment. It writes
coordinate-sorted BAM records with the minimum evidence classes consumed by the
scanner: soft clips, local CIGAR indels, SA tags, and discordant mate fields.
"""

from __future__ import annotations

import argparse
import time
import csv
import random
from dataclasses import dataclass
from pathlib import Path

import pysam

CIGAR_MATCH = 0
CIGAR_INS = 1
CIGAR_DEL = 2
CIGAR_SOFT_CLIP = 4
MMEJ_MICROHOMOLOGY = "AGCTA"


@dataclass(frozen=True, slots=True)
class TruthEvent:
    truth_id: str
    event_type: str
    chrom: str
    bkp_A_pos: int
    bkp_B_chrom: str
    bkp_B_pos: str
    remote_chrom: str
    remote_pos: str
    inserted_sequence: str
    inserted_length: str
    deleted_length: str
    support_reads: int
    candidate_start: int
    candidate_end: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a 20-event synthetic BAM benchmark for XXEJ_scanner."
    )
    parser.add_argument(
        "--output-dir",
        default="sim/xxej_benchmark",
        help="Directory for ref.fa, BAMs, candidate BED, and truth.tsv.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed. Defaults to current time in nanoseconds if not provided.",
    )
    parser.add_argument(
        "--support",
        type=int,
        default=6,
        help="Support reads per synthetic evidence class.",
    )
    parser.add_argument(
        "--background-reads",
        type=int,
        default=24,
        help="Normal 150M reads to place around each candidate region per BAM.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite known benchmark files in the output directory.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    if args.seed is None:
        args.seed = time.time_ns() % (2**32)
    rng = random.Random(args.seed)
    _prepare_output_dir(output_dir, force=args.force)

    seqs = _build_reference(rng)
    ref_path = output_dir / "ref.fa"
    _write_reference(ref_path, seqs)

    events = _truth_events(args.support)
    _write_truth(output_dir / "truth.tsv", events)
    _write_candidates(output_dir / "candidates.bed", events)

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": chrom, "LN": len(seq)} for chrom, seq in seqs.items()],
    }
    alignment_header = pysam.AlignmentHeader.from_dict(header)
    tid_by_chrom = {chrom: idx for idx, chrom in enumerate(seqs)}

    treated_records: list[pysam.AlignedSegment] = []
    control_records: list[pysam.AlignedSegment] = []
    for event in events:
        treated_records.extend(
            _background_reads(
                seqs,
                alignment_header,
                tid_by_chrom,
                event.chrom,
                event.candidate_start,
                event.candidate_end,
                args.background_reads,
                prefix=f"treated_bg_{event.truth_id}",
            )
        )
        control_records.extend(
            _background_reads(
                seqs,
                alignment_header,
                tid_by_chrom,
                event.chrom,
                event.candidate_start,
                event.candidate_end,
                args.background_reads,
                prefix=f"control_bg_{event.truth_id}",
            )
        )

    for event in events:
        treated_records.extend(
            _event_records(seqs, alignment_header, tid_by_chrom, event, args.support)
        )

    _write_sorted_bam(
        output_dir / "treated.sorted.bam",
        treated_records,
        header,
        force=args.force,
    )
    _write_sorted_bam(
        output_dir / "control.sorted.bam",
        control_records,
        header,
        force=args.force,
    )
    _write_helper_commands(output_dir)

    print(f"Wrote synthetic benchmark to {output_dir}")
    print(f"  reference: {output_dir / 'ref.fa'}")
    print(f"  treated:   {output_dir / 'treated.sorted.bam'}")
    print(f"  control:   {output_dir / 'control.sorted.bam'}")
    print(f"  truth:     {output_dir / 'truth.tsv'}")
    print(f"  BED:       {output_dir / 'candidates.bed'}")
    print(f"  helper:    {output_dir / 'run_scanner.sh'}")
    return 0


def _prepare_output_dir(output_dir: Path, *, force: bool) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    known_files = [
        "ref.fa",
        "ref.fa.fai",
        "truth.tsv",
        "candidates.bed",
        "treated.sorted.bam",
        "treated.sorted.bam.bai",
        "control.sorted.bam",
        "control.sorted.bam.bai",
        "run_scanner.sh",
    ]
    existing = [
        output_dir / name for name in known_files if (output_dir / name).exists()
    ]
    if existing and not force:
        names = ", ".join(path.name for path in existing[:5])
        raise SystemExit(
            f"{output_dir} already contains benchmark files ({names}). "
            "Use --force to overwrite them."
        )
    for path in existing:
        path.unlink()


def _build_reference(rng: random.Random) -> dict[str, str]:
    seqs = {
        "chrSim1": list(_random_dna(rng, 220_000)),
        "chrSim2": list(_random_dna(rng, 80_000)),
    }
    # Make every MMEJ example carry an obvious 5 bp reference microhomology.
    for event in _truth_events(support=1):
        if event.event_type != "MMEJ_DEL":
            continue
        left_bkp = event.bkp_A_pos
        right_bkp = int(event.bkp_B_pos)
        seqs[event.chrom][
            left_bkp - len(MMEJ_MICROHOMOLOGY) : left_bkp
        ] = MMEJ_MICROHOMOLOGY
        seqs[event.chrom][
            right_bkp : right_bkp + len(MMEJ_MICROHOMOLOGY)
        ] = MMEJ_MICROHOMOLOGY
    return {chrom: "".join(seq) for chrom, seq in seqs.items()}


def _truth_events(support: int) -> list[TruthEvent]:
    nhej_specs = [
        (24_000, "GATTACA"),
        (50_000, "TACGGA"),
        (75_000, "CCGTTA"),
        (125_000, "ATGCCAA"),
        (150_000, "TTAGGC"),
        (175_000, "CGATAC"),
        (205_000, "GGAATTC"),
    ]
    mmej_specs = [
        (35_000, 35_240),
        (65_000, 65_210),
        (95_000, 95_275),
        (115_000, 115_180),
        (145_000, 145_260),
        (185_000, 185_220),
        (210_000, 210_240),
    ]
    bnd_specs = [
        (15_000, "chrSim2", 10_000),
        (85_000, "chrSim2", 20_000),
        (135_000, "chrSim2", 30_000),
        (165_000, "chrSim2", 40_000),
        (195_000, "chrSim2", 50_000),
        (215_000, "chrSim2", 60_000),
    ]

    events: list[TruthEvent] = []
    for idx, (breakpoint, inserted) in enumerate(nhej_specs, 1):
        events.append(
            TruthEvent(
                truth_id=f"SIM_NHEJ_INS_{idx:03d}",
                event_type="NHEJ_INS",
                chrom="chrSim1",
                bkp_A_pos=breakpoint,
                bkp_B_chrom="NA",
                bkp_B_pos="NA",
                remote_chrom="NA",
                remote_pos="NA",
                inserted_sequence=inserted,
                inserted_length=str(len(inserted)),
                deleted_length="NA",
                support_reads=support * 2,
                candidate_start=breakpoint - 300,
                candidate_end=breakpoint + 300,
            )
        )
    for idx, (left_bkp, right_bkp) in enumerate(mmej_specs, 1):
        events.append(
            TruthEvent(
                truth_id=f"SIM_MMEJ_DEL_{idx:03d}",
                event_type="MMEJ_DEL",
                chrom="chrSim1",
                bkp_A_pos=left_bkp,
                bkp_B_chrom="chrSim1",
                bkp_B_pos=str(right_bkp),
                remote_chrom="NA",
                remote_pos="NA",
                inserted_sequence="NA",
                inserted_length="NA",
                deleted_length=str(right_bkp - left_bkp),
                support_reads=support * 3,
                candidate_start=left_bkp - 300,
                candidate_end=right_bkp + 300,
            )
        )
    for idx, (breakpoint, remote_chrom, remote_pos) in enumerate(bnd_specs, 1):
        events.append(
            TruthEvent(
                truth_id=f"SIM_BND_INTER_{idx:03d}",
                event_type="NHEJ_BND_INS_INTER",
                chrom="chrSim1",
                bkp_A_pos=breakpoint,
                bkp_B_chrom=remote_chrom,
                bkp_B_pos=str(remote_pos),
                remote_chrom=remote_chrom,
                remote_pos=str(remote_pos),
                inserted_sequence="NA",
                inserted_length="NA",
                deleted_length="NA",
                support_reads=support,
                candidate_start=breakpoint - 300,
                candidate_end=breakpoint + 300,
            )
        )
    if len(events) != 20:
        raise AssertionError(f"Expected 20 benchmark events, found {len(events)}")
    return events


def _write_reference(path: Path, seqs: dict[str, str]) -> None:
    with path.open("w") as handle:
        for chrom, seq in seqs.items():
            handle.write(f">{chrom}\n")
            for start in range(0, len(seq), 80):
                handle.write(seq[start : start + 80] + "\n")
    pysam.faidx(str(path))


def _write_truth(path: Path, events: list[TruthEvent]) -> None:
    fields = list(TruthEvent.__dataclass_fields__)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for event in events:
            writer.writerow({field: getattr(event, field) for field in fields})


def _write_candidates(path: Path, events: list[TruthEvent]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for event in events:
            writer.writerow(
                [
                    event.chrom,
                    event.candidate_start,
                    event.candidate_end,
                    event.truth_id,
                    100,
                ]
            )


def _write_sorted_bam(
    path: Path,
    records: list[pysam.AlignedSegment],
    header: dict[str, object],
    *,
    force: bool,
) -> None:
    unsorted = path.with_name(f"{path.stem}.unsorted{path.suffix}")
    for candidate in (path, Path(f"{path}.bai"), unsorted, Path(f"{unsorted}.bai")):
        if candidate.exists():
            if not force:
                raise SystemExit(f"{candidate} exists. Use --force to overwrite it.")
            candidate.unlink()

    with pysam.AlignmentFile(unsorted, "wb", header=header) as bam:
        for record in records:
            bam.write(record)
    pysam.sort("-o", str(path), str(unsorted))
    pysam.index(str(path))
    unsorted.unlink()


def _background_reads(
    seqs: dict[str, str],
    header: pysam.AlignmentHeader,
    tid_by_chrom: dict[str, int],
    chrom: str,
    start: int,
    end: int,
    count: int,
    *,
    prefix: str,
    read_length: int = 150,
) -> list[pysam.AlignedSegment]:
    records: list[pysam.AlignedSegment] = []
    if count <= 0:
        return records
    available = max(1, end - start - read_length)
    for idx in range(count):
        offset = (
            int(idx * available / max(1, count - 1)) if count > 1 else available // 2
        )
        read_start = start + offset
        query = _ref(seqs, chrom, read_start, read_length)
        records.append(
            _segment(
                header,
                tid_by_chrom,
                query_name=f"{prefix}_{idx:03d}",
                chrom=chrom,
                start=read_start,
                cigar=[(CIGAR_MATCH, read_length)],
                query_sequence=query,
            )
        )
    return records


def _event_records(
    seqs: dict[str, str],
    header: pysam.AlignmentHeader,
    tid_by_chrom: dict[str, int],
    event: TruthEvent,
    support: int,
) -> list[pysam.AlignedSegment]:
    if event.event_type == "NHEJ_INS":
        return _nhej_ins_records(seqs, header, tid_by_chrom, event, support)
    if event.event_type == "MMEJ_DEL":
        return _mmej_del_records(seqs, header, tid_by_chrom, event, support)
    if event.event_type == "NHEJ_BND_INS_INTER":
        return _bnd_inter_records(seqs, header, tid_by_chrom, event, support)
    raise ValueError(f"Unsupported benchmark event type: {event.event_type}")


def _nhej_ins_records(
    seqs: dict[str, str],
    header: pysam.AlignmentHeader,
    tid_by_chrom: dict[str, int],
    event: TruthEvent,
    support: int,
) -> list[pysam.AlignedSegment]:
    chrom = event.chrom
    breakpoint = event.bkp_A_pos
    inserted = event.inserted_sequence
    prefix = event.truth_id.lower()
    records: list[pysam.AlignedSegment] = []

    for idx in range(support):
        clip_len = 24
        aligned = 96
        clip = _repeat_to_length("TTAACCGGATCG", clip_len)
        query = clip + _ref(seqs, chrom, breakpoint, aligned)
        records.append(
            _segment(
                header,
                tid_by_chrom,
                query_name=f"{prefix}_clip_{idx:03d}",
                chrom=chrom,
                start=breakpoint,
                cigar=[(CIGAR_SOFT_CLIP, clip_len), (CIGAR_MATCH, aligned)],
                query_sequence=query,
            )
        )

    for idx in range(support):
        left = 60
        right = 60
        start = breakpoint - left
        query = (
            _ref(seqs, chrom, start, left)
            + inserted
            + _ref(seqs, chrom, breakpoint, right)
        )
        records.append(
            _segment(
                header,
                tid_by_chrom,
                query_name=f"{prefix}_ins_{idx:03d}",
                chrom=chrom,
                start=start,
                cigar=[
                    (CIGAR_MATCH, left),
                    (CIGAR_INS, len(inserted)),
                    (CIGAR_MATCH, right),
                ],
                query_sequence=query,
            )
        )
    return records


def _mmej_del_records(
    seqs: dict[str, str],
    header: pysam.AlignmentHeader,
    tid_by_chrom: dict[str, int],
    event: TruthEvent,
    support: int,
) -> list[pysam.AlignedSegment]:
    chrom = event.chrom
    left_bkp = event.bkp_A_pos
    right_bkp = int(event.bkp_B_pos)
    deletion_length = right_bkp - left_bkp
    prefix = event.truth_id.lower()
    records: list[pysam.AlignedSegment] = []

    for idx in range(support):
        aligned = 90
        clip_len = 30
        start = left_bkp - aligned
        query = _ref(seqs, chrom, start, aligned) + _ref(
            seqs, chrom, right_bkp, clip_len
        )
        records.append(
            _segment(
                header,
                tid_by_chrom,
                query_name=f"{prefix}_leftclip_{idx:03d}",
                chrom=chrom,
                start=start,
                cigar=[(CIGAR_MATCH, aligned), (CIGAR_SOFT_CLIP, clip_len)],
                query_sequence=query,
            )
        )

    for idx in range(support):
        aligned = 90
        clip_len = 30
        query = _ref(seqs, chrom, left_bkp - clip_len, clip_len) + _ref(
            seqs, chrom, right_bkp, aligned
        )
        records.append(
            _segment(
                header,
                tid_by_chrom,
                query_name=f"{prefix}_rightclip_{idx:03d}",
                chrom=chrom,
                start=right_bkp,
                cigar=[(CIGAR_SOFT_CLIP, clip_len), (CIGAR_MATCH, aligned)],
                query_sequence=query,
            )
        )

    for idx in range(support):
        flank = 70
        start = left_bkp - flank
        query = _ref(seqs, chrom, start, flank) + _ref(seqs, chrom, right_bkp, flank)
        records.append(
            _segment(
                header,
                tid_by_chrom,
                query_name=f"{prefix}_del_{idx:03d}",
                chrom=chrom,
                start=start,
                cigar=[
                    (CIGAR_MATCH, flank),
                    (CIGAR_DEL, deletion_length),
                    (CIGAR_MATCH, flank),
                ],
                query_sequence=query,
            )
        )
    return records


def _bnd_inter_records(
    seqs: dict[str, str],
    header: pysam.AlignmentHeader,
    tid_by_chrom: dict[str, int],
    event: TruthEvent,
    support: int,
) -> list[pysam.AlignedSegment]:
    chrom = event.chrom
    breakpoint = event.bkp_A_pos
    remote_chrom = event.remote_chrom
    remote_pos = int(event.remote_pos)
    prefix = event.truth_id.lower()
    records: list[pysam.AlignedSegment] = []

    for idx in range(support):
        clip_len = 35
        aligned = 85
        mate_pos = remote_pos + (idx % 5)
        query = _ref(seqs, remote_chrom, mate_pos, clip_len) + _ref(
            seqs, chrom, breakpoint, aligned
        )
        record = _segment(
            header,
            tid_by_chrom,
            query_name=f"{prefix}_{idx:03d}",
            chrom=chrom,
            start=breakpoint,
            cigar=[(CIGAR_SOFT_CLIP, clip_len), (CIGAR_MATCH, aligned)],
            query_sequence=query,
            flag=0x1 | 0x20 | 0x40,
            next_chrom=remote_chrom,
            next_start=mate_pos,
        )
        if idx < max(1, support // 2):
            record.set_tag("SA", f"{remote_chrom},{mate_pos + 1},+,85M35S,60,0;")
        records.append(record)
    return records


def _segment(
    header: pysam.AlignmentHeader,
    tid_by_chrom: dict[str, int],
    *,
    query_name: str,
    chrom: str,
    start: int,
    cigar: list[tuple[int, int]],
    query_sequence: str,
    flag: int = 0,
    next_chrom: str | None = None,
    next_start: int = -1,
    mapq: int = 60,
) -> pysam.AlignedSegment:
    record = pysam.AlignedSegment(header)
    record.query_name = query_name
    record.query_sequence = query_sequence
    record.flag = flag
    record.reference_id = tid_by_chrom[chrom]
    record.reference_start = start
    record.mapping_quality = mapq
    record.cigartuples = cigar
    if next_chrom is None:
        record.next_reference_id = -1
        record.next_reference_start = -1
    else:
        record.next_reference_id = tid_by_chrom[next_chrom]
        record.next_reference_start = next_start
    record.template_length = 0
    record.query_qualities = pysam.qualitystring_to_array("I" * len(query_sequence))
    return record


def _write_helper_commands(output_dir: Path) -> None:
    script = output_dir / "run_scanner.sh"
    scanner_out = output_dir / "scanner_out"
    script.write_text(
        "\n".join(
            [
                "#!/usr/bin/env bash",
                "set -euo pipefail",
                'SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
                'PROJECT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"',
                'cd "${PROJECT_DIR}"',
                "",
                "uv run --no-editable XXEJ_scanner scan \\",
                f"  --treated-bam {output_dir / 'treated.sorted.bam'} \\",
                f"  --control-bam {output_dir / 'control.sorted.bam'} \\",
                f"  --reference-fasta {output_dir / 'ref.fa'} \\",
                f"  --candidate-bed {output_dir / 'candidates.bed'} \\",
                f"  --output-dir {scanner_out} \\",
                "  --min-alt-support 3 \\",
                "  --min-bnd-support 3 \\",
                "  --depth-count-method region",
                "",
                "uv run python scripts/evaluate_xxej_benchmark.py \\",
                f"  --truth {output_dir / 'truth.tsv'} \\",
                f"  --events {scanner_out / 'events.tsv'}",
                "",
            ]
        )
    )
    script.chmod(0o755)


def _ref(seqs: dict[str, str], chrom: str, start: int, length: int) -> str:
    if start < 0 or start + length > len(seqs[chrom]):
        raise ValueError(
            f"Reference slice out of bounds: {chrom}:{start}-{start + length}"
        )
    return seqs[chrom][start : start + length]


def _random_dna(rng: random.Random, length: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def _repeat_to_length(seed: str, length: int) -> str:
    return (seed * ((length // len(seed)) + 1))[:length]


if __name__ == "__main__":
    raise SystemExit(main())
