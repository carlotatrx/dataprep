#!/usr/bin/env python3
"""
Summarize SEF record lengths from filenames and split into separate series
when gaps exceed a threshold.

Unique identifier base: (station_source, variable), ignoring daily/subdaily and qc.

A "series" here = a merged cluster of segments where consecutive segments are
no more than --max-gap-days apart (or overlap).

Example:
  gap between end and next start > max_gap_days  => new series (new row)

Usage:
  python sef_series_lengths_v2.py /scratch3/PALAEO-RA/daily_data/final -o summary.csv --max-gap-days 30
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
    return date(int(s[0:4]), int(s[4:6]), int(s[6:8]))


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

    station_source = filename[:m.start()].rstrip("_")
    right = filename[m.end():].lstrip("_")
    if right.lower().endswith(".tsv"):
        right = right[:-4]

    if not station_source or not right:
        return None

    var = right.split("_", 1)[0]
    if not var:
        return None

    return station_source, start_dt, end_dt, var


def iter_tsv_files(root: str) -> Iterable[str]:
    for dirpath, _, filenames in os.walk(root):
        for fn in filenames:
            if fn.lower().endswith(".tsv"):
                yield os.path.join(dirpath, fn)


def cluster_into_series(segments: List[Segment], max_gap_days: int = 30) -> List[Tuple[date, date, List[Segment]]]:
    """
    Clusters segments into separate series based on max_gap_days.

    Two consecutive segments belong to the same series if:
      next.start <= current_end + max_gap_days

    Returns list of tuples:
      (series_start, series_end, contributing_segments)
    """
    if not segments:
        return []

    segs = sorted(segments, key=lambda s: (s.start, s.end))
    max_gap = timedelta(days=max_gap_days)

    out: List[Tuple[date, date, List[Segment]]] = []

    cur_start = segs[0].start
    cur_end = segs[0].end
    cur_list: List[Segment] = [segs[0]]

    for s in segs[1:]:
        if s.start <= (cur_end + max_gap):
            # same series; extend end if needed
            if s.end > cur_end:
                cur_end = s.end
            cur_list.append(s)
        else:
            # start new series
            out.append((cur_start, cur_end, cur_list))
            cur_start = s.start
            cur_end = s.end
            cur_list = [s]

    out.append((cur_start, cur_end, cur_list))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Summarize SEF record lengths by parsing filenames.")
    ap.add_argument("root", help="Root directory containing station folders and SEF .tsv files")
    ap.add_argument("-o", "--output", default="sef_series_summary.csv", help="Output CSV path")
    ap.add_argument("--max-gap-days", type=int, default=30,
                    help="If gap between segments is > this, start a new series (default: 30).")
    ap.add_argument("--include-folder", action="store_true",
                    help="Include the immediate parent folder name as a column (useful if folder=station).")
    ap.add_argument("--write-file-list", action="store_true",
                    help="Include a 'file_list' column (can be very large).")
    args = ap.parse_args()

    # Group raw segments by base key (station_source, variable) + optional folder
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
        "series_index",
        "start",
        "end",
        "span_days",
        "span_years",
        "n_files",
        "max_internal_gap_days",
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
            # split into separate series using max-gap-days
            series_list = cluster_into_series(segs, max_gap_days=args.max_gap_days)

            for idx, (s_start, s_end, contributing) in enumerate(series_list, start=1):
                # compute max gap *between segments* inside this series (for diagnostics)
                contributing_sorted = sorted(contributing, key=lambda s: (s.start, s.end))
                max_internal_gap = 0
                prev_end = contributing_sorted[0].end
                for seg in contributing_sorted[1:]:
                    gap_days = (seg.start - prev_end).days - 1
                    if gap_days > max_internal_gap:
                        max_internal_gap = gap_days
                    if seg.end > prev_end:
                        prev_end = seg.end

                span_days = (s_end - s_start).days + 1
                span_years = span_days / 365.2425

                blocks = f"{s_start.isoformat()}..{s_end.isoformat()}"

                row = {
                    "station_source": station_source,
                    "variable": var,
                    "series_index": idx,
                    "start": s_start.isoformat(),
                    "end": s_end.isoformat(),
                    "span_days": span_days,
                    "span_years": f"{span_years:.2f}",
                    "n_files": len(contributing),
                    "max_internal_gap_days": max_internal_gap,
                    "blocks": blocks,
                }
                if args.include_folder:
                    row["folder"] = folder or ""
                if args.write_file_list:
                    row["file_list"] = "|".join(sorted(s.filepath for s in contributing))

                w.writerow(row)

    print(f"Wrote: {args.output}")
    if skipped:
        print(f"Skipped (unparseable) files: {skipped}")
    print(f"Base groups (station_source+variable): {len(groups)}")
    print(f"Max gap days threshold: {args.max_gap_days}")


if __name__ == "__main__":
    main()
