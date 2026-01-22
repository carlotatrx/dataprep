#!/usr/bin/env python3
"""
Summarize SEF record lengths from filenames and merge consecutive/overlapping segments.

Assumes filenames contain a date block like YYYYMMDD-YYYYMMDD, e.g.:
stationname_start-end_var[_daily|_subdaily][_qc].tsv

We define a unique series as (station_source, variable), ignoring daily/subdaily and qc.

Outputs:
- overall start/end
- span_days / span_years
- covered_days (sum of merged blocks)
- n_files
- n_blocks (how many merged blocks remain)
- max_gap_days between blocks
- blocks string representation
- (optional) file_list

Usage:
    python sef_series_lengths.py /scratch3/PALAEO-RA/daily_data/final \
  -o sef_series_summary.csv --gap 1 --include-folder --write-file-list

"""

from __future__ import annotations

import argparse
import csv
import os
import re
from dataclasses import dataclass
from datetime import date, timedelta
from typing import Dict, Iterable, List, Optional, Tuple


DATE_BLOCK_RE = re.compile(r"(\d{8})-(\d{8})")


@dataclass(frozen=True)
class Segment:
    start: date
    end: date
    filepath: str


def parse_yyyymmdd(s: str) -> date:
    y = int(s[0:4])
    m = int(s[4:6])
    d = int(s[6:8])
    return date(y, m, d)


def parse_sef_filename(filename: str) -> Optional[Tuple[str, date, date, str]]:
    """
    Returns (station_source, start_date, end_date, variable) or None if not parseable.

    station_source = everything before the YYYYMMDD-YYYYMMDD block (stripped of trailing underscores)
    variable = first token after the date block, split on '_' (case preserved)
    """
    if not filename.lower().endswith(".tsv"):
        return None

    m = DATE_BLOCK_RE.search(filename)
    if not m:
        return None

    start_s, end_s = m.group(1), m.group(2)
    try:
        start_dt = parse_yyyymmdd(start_s)
        end_dt = parse_yyyymmdd(end_s)
    except ValueError:
        return None

    if end_dt < start_dt:
        return None

    left = filename[:m.start()].rstrip("_")
    right = filename[m.end():].lstrip("_")

    if right.lower().endswith(".tsv"):
        right = right[:-4]
    if not right:
        return None

    tokens = right.split("_")
    var = tokens[0]
    if not var or not left:
        return None

    return left, start_dt, end_dt, var


def merge_segments(segments: List[Segment], gap_tolerance_days: int = 1) -> List[Tuple[date, date]]:
    """
    Merge segments if they overlap or are within gap_tolerance_days.
    Returns merged blocks as list of (start, end), sorted.
    """
    if not segments:
        return []

    segs = sorted(segments, key=lambda s: (s.start, s.end))
    merged: List[Tuple[date, date]] = []
    gap = timedelta(days=gap_tolerance_days)

    cur_s, cur_e = segs[0].start, segs[0].end
    for s in segs[1:]:
        if s.start <= (cur_e + gap):
            if s.end > cur_e:
                cur_e = s.end
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = s.start, s.end

    merged.append((cur_s, cur_e))
    return merged


def max_gap_between_blocks(blocks: List[Tuple[date, date]]) -> int:
    """
    Computes the maximum number of gap days between consecutive blocks.
    If blocks overlap/are adjacent, gap is 0.
    If only one block, returns 0.
    """
    if len(blocks) <= 1:
        return 0

    blocks = sorted(blocks, key=lambda b: (b[0], b[1]))
    max_gap = 0
    for (s1, e1), (s2, e2) in zip(blocks[:-1], blocks[1:]):
        # gap days between e1 and s2 (exclusive). If s2 <= e1+1 => gap 0.
        gap = (s2 - e1).days - 1
        if gap > max_gap:
            max_gap = gap
    return max_gap


def blocks_to_string(blocks: List[Tuple[date, date]]) -> str:
    return "; ".join([f"{s.isoformat()}..{e.isoformat()}" for s, e in blocks])


def iter_tsv_files(root: str) -> Iterable[str]:
    for dirpath, _, filenames in os.walk(root):
        for fn in filenames:
            if fn.lower().endswith(".tsv"):
                yield os.path.join(dirpath, fn)


def main() -> None:
    ap = argparse.ArgumentParser(description="Summarize SEF record lengths by parsing filenames.")
    ap.add_argument("root", help="Root directory containing station folders and SEF .tsv files")
    ap.add_argument("-o", "--output", default="sef_series_summary.csv", help="Output CSV path")
    ap.add_argument("--gap", type=int, default=1, help="Gap tolerance in days for merging segments (default: 1)")
    ap.add_argument("--include-folder", action="store_true",
                    help="Include the immediate parent folder name as a column (useful if folder=station)")
    ap.add_argument("--write-file-list", action="store_true",
                    help="Include a 'file_list' column (can be very large).")
    args = ap.parse_args()

    groups: Dict[Tuple[str, str, Optional[str]], List[Segment]] = {}
    skipped = 0

    for path in iter_tsv_files(args.root):
        fn = os.path.basename(path)
        parsed = parse_sef_filename(fn)
        if not parsed:
            skipped += 1
            continue

        station_source, start_dt, end_dt, var = parsed
        folder = os.path.basename(os.path.dirname(path)) if args.include_folder else None

        key = (station_source, var, folder)
        groups.setdefault(key, []).append(Segment(start_dt, end_dt, path))

    fieldnames = [
        "station_source",
        "variable",
        "folder" if args.include_folder else None,
        "start",
        "end",
        "span_days",
        "span_years",
        "covered_days",
        "n_files",
        "n_blocks",
        "max_gap_days",
        "blocks",
        "file_list" if args.write_file_list else None,
    ]
    fieldnames = [f for f in fieldnames if f is not None]

    with open(args.output, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()

        for (station_source, var, folder), segs in sorted(
            groups.items(),
            key=lambda x: (x[0][0], x[0][1], x[0][2] or "")
        ):
            blocks = merge_segments(segs, gap_tolerance_days=args.gap)
            overall_start = min(b[0] for b in blocks)
            overall_end = max(b[1] for b in blocks)

            span_days = (overall_end - overall_start).days + 1
            span_years = span_days / 365.2425
            covered_days = sum((b[1] - b[0]).days + 1 for b in blocks)

            row = {
                "station_source": station_source,
                "variable": var,
                "start": overall_start.isoformat(),
                "end": overall_end.isoformat(),
                "span_days": span_days,
                "span_years": f"{span_years:.2f}",
                "covered_days": covered_days,
                "n_files": len(segs),
                "n_blocks": len(blocks),
                "max_gap_days": max_gap_between_blocks(blocks),
                "blocks": blocks_to_string(blocks),
            }
            if args.include_folder:
                row["folder"] = folder or ""
            if args.write_file_list:
                row["file_list"] = "|".join(sorted(s.filepath for s in segs))

            w.writerow(row)

    print(f"Wrote: {args.output}")
    if skipped:
        print(f"Skipped (unparseable) files: {skipped}")
    print(f"Parsed series: {len(groups)}")


if __name__ == "__main__":
    main()
