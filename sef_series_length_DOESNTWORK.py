import pandas as pd

v1 = pd.read_csv("sef_series_summary_v1.csv")
v2 = pd.read_csv("sef_series_summary_v2.csv")

# Count series per (station_source, variable) in v2
v2_counts = v2.groupby(["station_source", "variable"]).size().reset_index(name="n_series")

# Keep only those split into more than 1 series
candidates = v2_counts[v2_counts["n_series"] > 1].sort_values("n_series", ascending=False)

# Enrich with the actual series details so you can inspect the gaps
candidates_full = v2.merge(candidates[["station_source", "variable"]], on=["station_source", "variable"])
candidates_full = candidates_full.sort_values(["station_source", "variable", "start"])

candidates_full.to_csv("grouping_candidates.csv", index=False)
print(f"Unique (station, variable) pairs to review: {len(candidates)}")
print(f"Total series rows involved: {len(candidates_full)}")

#!/usr/bin/env python3
"""
Unified SEF Series Processor
- Mode 1 (Summary): One row per station/variable. Shows total coverage and gaps.
- Mode 2 (Clustered): Splits station/variable into new rows if gaps exceed a threshold.

To get 1-row-per-station summary (equivalent to sef_series_length_v1.py):
python sef_series_length.py /scracth3/PALAEO-RA/daily_data/final/ -o summary.csv --mode summary

To get the continuous series, splitting by gaps (equivalent to sef_series_length_v2.py with 365-day threshold):
python sef_series_length.py /scracth3/PALAEO-RA/daily_data/final/ -o split.csv --mode split --threshold 365
"""

from __future__ import annotations
import argparse
import csv
import os
import re
from dataclasses import dataclass
from datetime import date, timedelta
from typing import Dict, Iterable, List, Optional, Tuple

DATE_BLOCK_RE = re.compile(r"(\d{8})[-_](\d{8})") # Updated to handle Durham's underscore

@dataclass(frozen=True)
class Segment:
    start: date
    end: date
    filepath: str

def parse_yyyymmdd(s: str) -> date:
    return date(int(s[0:4]), int(s[4:6]), int(s[6:8]))

def parse_sef_filename(filename: str) -> Optional[Tuple[str, date, date, str]]:
    if not filename.lower().endswith(".tsv"):
        return None
    m = DATE_BLOCK_RE.search(filename)
    if not m:
        return None
    try:
        start_dt, end_dt = parse_yyyymmdd(m.group(1)), parse_yyyymmdd(m.group(2))
    except ValueError:
        return None
    if end_dt < start_dt:
        return None

    left = filename[:m.start()].rstrip("_")
    right = filename[m.end():].lstrip("_").replace(".tsv", "").replace(".TSV", "")
    var = right.split("_")[0]
    return left, start_dt, end_dt, var

def get_metrics(blocks: List[Tuple[date, date]]):
    """Calculates gap and coverage metrics for a set of date blocks."""
    blocks = sorted(blocks, key=lambda b: b[0])
    overall_start, overall_end = blocks[0][0], blocks[-1][1]
    span_days = (overall_end - overall_start).days + 1
    covered_days = sum((b[1] - b[0]).days + 1 for b in blocks)
    
    max_gap = 0
    for i in range(len(blocks) - 1):
        gap = (blocks[i+1][0] - blocks[i][1]).days - 1
        max_gap = max(max_gap, gap)
    
    return overall_start, overall_end, span_days, covered_days, max(0, max_gap)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root", help="Directory with SEF files")
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--mode", choices=["summary", "split"], default="summary", 
                    help="'summary' = 1 row per station; 'split' = new row if gap > threshold")
    ap.add_argument("--threshold", type=int, default=365, help="Gap threshold in days")
    args = ap.parse_args()

    groups = {}
    for dirpath, _, filenames in os.walk(args.root):
        for fn in filenames:
            parsed = parse_sef_filename(fn)
            if parsed:
                key = (parsed[0], parsed[3]) # (station_source, variable)
                groups.setdefault(key, []).append(Segment(parsed[1], parsed[2], os.path.join(dirpath, fn)))

    fieldnames = ["station_source", "variable", "series_idx", "start", "end", 
                  "span_days", "span_years", "covered_days", "n_files", "max_gap", "blocks"]
    
    with open(args.output, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        
        for (stat, var), segs in sorted(groups.items()):
            segs = sorted(segs, key=lambda s: s.start)
            
            # Logic for splitting vs grouping
            clusters = []
            if args.mode == "summary":
                clusters.append(segs)
            else:
                current_cluster = [segs[0]]
                for i in range(1, len(segs)):
                    if (segs[i].start - current_cluster[-1].end).days - 1 > args.threshold:
                        clusters.append(current_cluster)
                        current_cluster = [segs[i]]
                    else:
                        current_cluster.append(segs[i])
                clusters.append(current_cluster)

            for idx, cluster in enumerate(clusters, 1):
                blocks = [(s.start, s.end) for s in cluster]
                s, e, d_span, d_cov, m_gap = get_metrics(blocks)
                w.writerow({
                    "station_source": stat, "variable": var, "series_idx": idx,
                    "start": s.isoformat(), "end": e.isoformat(),
                    "span_days": d_span, "span_years": round(d_span/365.25, 2),
                    "covered_days": d_cov, "n_files": len(cluster),
                    "max_gap": m_gap, "blocks": "; ".join([f"{b[0]}..{b[1]}" for b in blocks])
                })

if __name__ == "__main__":
    main()