#!/usr/bin/env python3
"""
Compute true per-year data completeness (years with at least one non-missing
value) for each of the 1803 series listed in metadata_summary_names.csv.

Unlike sef_active_years.csv (produced by sef_series_length_v2.py), which
dedupes active years to (station_source, variable) pairs, this keeps one
identity per SEF file/series so totals line up with the 1803-series count
used elsewhere (e.g. the pie charts in plots4presentation.ipynb).

Usage:
  python3 sef_active_years_1803.py
"""

import csv
import os

import pandas as pd

from sef_series_length_v2 import read_active_years

ROOT = "/scratch3/PALAEO-RA/daily_data/final"
METADATA_CSV = "metadata_summary_names.csv"
OUTPUT_CSV = "sef_active_years_1803.csv"


def main() -> None:
    df = pd.read_csv(METADATA_CSV)

    rows = []
    missing_files = 0
    for _, row in df.iterrows():
        path = os.path.join(ROOT, str(row["dirname"]), str(row["filename"]))
        if not os.path.exists(path):
            missing_files += 1
            continue
        for year in read_active_years(path):
            rows.append({
                "dirname": row["dirname"],
                "filename": row["filename"],
                "variable": row["Vbl"],
                "year": year,
            })

    with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["dirname", "filename", "variable", "year"])
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["dirname"], r["filename"], r["year"])))

    print(f"Processed {len(df)} series from {METADATA_CSV}")
    if missing_files:
        print(f"Missing files (skipped): {missing_files}")
    print(f"Wrote {len(rows)} active-year rows to {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
