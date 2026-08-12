#!/usr/bin/env python3
"""Compare XXEJ_scanner events against the synthetic benchmark truth table."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate XXEJ_scanner synthetic benchmark recovery."
    )
    parser.add_argument("--truth", required=True, help="truth.tsv from the simulator.")
    parser.add_argument("--events", required=True, help="events.tsv from XXEJ_scanner.")
    parser.add_argument(
        "--breakpoint-window",
        type=int,
        default=25,
        help="Allowed absolute breakpoint error in bp.",
    )
    parser.add_argument(
        "--include-filtered",
        action="store_true",
        help="Count scanner events even if filter is not PASS.",
    )
    parser.add_argument(
        "--matches-out",
        default=None,
        help="Optional path for per-truth match details.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    truth_rows = _read_tsv(Path(args.truth))
    event_rows = _read_tsv(Path(args.events))
    if not args.include_filtered:
        event_rows = [row for row in event_rows if row.get("filter") == "PASS"]

    matches, unmatched_truth, unmatched_events = _match_events(
        truth_rows,
        event_rows,
        window=args.breakpoint_window,
    )
    _print_summary(truth_rows, event_rows, matches, unmatched_truth, unmatched_events)
    if args.matches_out:
        _write_matches(Path(args.matches_out), matches, unmatched_truth)
    return 0


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _match_events(
    truth_rows: list[dict[str, str]],
    event_rows: list[dict[str, str]],
    *,
    window: int,
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    remaining = list(event_rows)
    matches: list[dict[str, str]] = []
    unmatched_truth: list[dict[str, str]] = []

    for truth in truth_rows:
        candidates = [
            (idx, event, _match_distance(truth, event, window))
            for idx, event in enumerate(remaining)
            if _is_compatible(truth, event, window)
        ]
        if not candidates:
            unmatched_truth.append(truth)
            continue
        idx, event, distance = min(candidates, key=lambda item: item[2])
        remaining.pop(idx)
        matches.append(
            {
                "truth_id": truth["truth_id"],
                "truth_type": truth["event_type"],
                "event_id": event["event_id"],
                "event_type": event["event_type"],
                "filter": event.get("filter", "NA"),
                "truth_bkp_A": truth["bkp_A_pos"],
                "event_bkp_A": event["bkp_A_pos"],
                "truth_bkp_B": truth.get("bkp_B_pos", "NA"),
                "event_bkp_B": event.get("bkp_B_pos", "NA"),
                "truth_remote": truth.get("remote_pos", "NA"),
                "event_remote": event.get("remote_pos", "NA"),
                "distance": str(distance),
            }
        )
    return matches, unmatched_truth, remaining


def _is_compatible(
    truth: dict[str, str],
    event: dict[str, str],
    window: int,
) -> bool:
    if truth["event_type"] != event.get("event_type"):
        return False
    event_type = truth["event_type"]
    if event_type.startswith("BND_"):
        direct = (
            truth["chrom"] == event.get("bkp_A_chrom")
            and truth.get("remote_chrom") == event.get("bkp_B_chrom")
            and _near(truth.get("bkp_A_pos"), event.get("bkp_A_pos"), window)
            and _near(truth.get("remote_pos"), event.get("bkp_B_pos"), window)
        )
        reverse = (
            truth["chrom"] == event.get("bkp_B_chrom")
            and truth.get("remote_chrom") == event.get("bkp_A_chrom")
            and _near(truth.get("bkp_A_pos"), event.get("bkp_B_pos"), window)
            and _near(truth.get("remote_pos"), event.get("bkp_A_pos"), window)
        )
        return direct or reverse
    if truth["chrom"] != event.get("chrom"):
        return False
    if not _near(truth.get("bkp_A_pos"), event.get("bkp_A_pos"), window):
        return False

    if event_type == "LOCAL_INS":
        return truth.get("inserted_sequence") == event.get("inserted_sequence")
    if event_type == "LOCAL_DEL":
        return truth.get("bkp_B_chrom") == event.get("bkp_B_chrom") and _near(
            truth.get("bkp_B_pos"), event.get("bkp_B_pos"), window
        )
    return True


def _match_distance(
    truth: dict[str, str],
    event: dict[str, str],
    window: int,
) -> int:
    if truth["event_type"].startswith("BND_"):
        direct = (
            _distance(truth.get("bkp_A_pos"), event.get("bkp_A_pos"), window)
            + _distance(truth.get("remote_pos"), event.get("bkp_B_pos"), window)
            if truth["chrom"] == event.get("bkp_A_chrom")
            and truth.get("remote_chrom") == event.get("bkp_B_chrom")
            else 2 * (window + 1)
        )
        reverse = (
            _distance(truth.get("bkp_A_pos"), event.get("bkp_B_pos"), window)
            + _distance(truth.get("remote_pos"), event.get("bkp_A_pos"), window)
            if truth["chrom"] == event.get("bkp_B_chrom")
            and truth.get("remote_chrom") == event.get("bkp_A_chrom")
            else 2 * (window + 1)
        )
        return min(direct, reverse)

    distance = _distance(truth.get("bkp_A_pos"), event.get("bkp_A_pos"), window)
    if truth["event_type"] == "LOCAL_DEL":
        distance += _distance(truth.get("bkp_B_pos"), event.get("bkp_B_pos"), window)
    return distance


def _near(left: str | None, right: str | None, window: int) -> bool:
    return _distance(left, right, window + 1) <= window


def _distance(left: str | None, right: str | None, missing_value: int) -> int:
    try:
        return abs(int(str(left)) - int(str(right)))
    except (TypeError, ValueError):
        return missing_value


def _print_summary(
    truth_rows: list[dict[str, str]],
    event_rows: list[dict[str, str]],
    matches: list[dict[str, str]],
    unmatched_truth: list[dict[str, str]],
    unmatched_events: list[dict[str, str]],
) -> None:
    truth_by_type = Counter(row["event_type"] for row in truth_rows)
    event_by_type = Counter(row["event_type"] for row in event_rows)
    tp_by_type = Counter(row["truth_type"] for row in matches)
    fn_by_type = Counter(row["event_type"] for row in unmatched_truth)
    fp_by_type = Counter(row["event_type"] for row in unmatched_events)
    event_types = sorted(set(truth_by_type) | set(event_by_type))

    print("event_type\ttruth\tcalled\tTP\tFP\tFN\trecall\tprecision")
    for event_type in event_types:
        truth_count = truth_by_type[event_type]
        called_count = event_by_type[event_type]
        tp = tp_by_type[event_type]
        fp = fp_by_type[event_type]
        fn = fn_by_type[event_type]
        recall = tp / truth_count if truth_count else 0.0
        precision = tp / called_count if called_count else 0.0
        print(
            f"{event_type}\t{truth_count}\t{called_count}\t{tp}\t{fp}\t{fn}\t"
            f"{recall:.3f}\t{precision:.3f}"
        )

    print(
        f"TOTAL\t{len(truth_rows)}\t{len(event_rows)}\t{len(matches)}\t"
        f"{len(unmatched_events)}\t{len(unmatched_truth)}\t"
        f"{(len(matches) / len(truth_rows) if truth_rows else 0.0):.3f}\t"
        f"{(len(matches) / len(event_rows) if event_rows else 0.0):.3f}"
    )

    if unmatched_truth:
        missing = ", ".join(row["truth_id"] for row in unmatched_truth)
        print(f"Unmatched truth: {missing}")
    if unmatched_events:
        extras = ", ".join(row["event_id"] for row in unmatched_events[:10])
        print(f"Unmatched calls: {extras}")


def _write_matches(
    path: Path,
    matches: list[dict[str, str]],
    unmatched_truth: list[dict[str, str]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "truth_id",
        "truth_type",
        "event_id",
        "event_type",
        "filter",
        "truth_bkp_A",
        "event_bkp_A",
        "truth_bkp_B",
        "event_bkp_B",
        "truth_remote",
        "event_remote",
        "distance",
        "status",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for match in matches:
            writer.writerow({**match, "status": "TP"})
        for truth in unmatched_truth:
            writer.writerow(
                {
                    "truth_id": truth["truth_id"],
                    "truth_type": truth["event_type"],
                    "event_id": "NA",
                    "event_type": "NA",
                    "filter": "NA",
                    "truth_bkp_A": truth["bkp_A_pos"],
                    "event_bkp_A": "NA",
                    "truth_bkp_B": truth.get("bkp_B_pos", "NA"),
                    "event_bkp_B": "NA",
                    "truth_remote": truth.get("remote_pos", "NA"),
                    "event_remote": "NA",
                    "distance": "NA",
                    "status": "FN",
                }
            )


if __name__ == "__main__":
    raise SystemExit(main())
