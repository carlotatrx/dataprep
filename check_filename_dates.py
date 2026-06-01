"""
Checks that the start/end dates in each SEF .tsv filename match the
first and last observation dates in the file, and flags files where
the first or last data row has a NA value.

Usage:
    python check_filename_dates.py /scratch3/PALAEO-RA/daily_data/final/
"""

import os
import re
import sys
from datetime import date, datetime

DATE_PATTERN = re.compile(r"_(\d{8})-(\d{8})_")


def parse_filename_dates(filename):
    m = DATE_PATTERN.search(filename)
    if not m:
        return None, None
    return (datetime.strptime(m.group(1), "%Y%m%d").date(),
            datetime.strptime(m.group(2), "%Y%m%d").date())


def get_data_rows(filepath):
    rows = []
    in_data = False
    with open(filepath, encoding="utf-8", errors="replace") as f:
        for line in f:
            if not in_data:
                if line.startswith("Year\t"):
                    in_data = True
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 7:
                rows.append(parts)
    return rows


def row_date(parts):
    return date(int(parts[0]), int(parts[1]), int(parts[2]))


def is_na(parts):
    return parts[6].strip().upper() in ("NA", "NAN", "N/A", "", "NULL", "MISSING")


root = sys.argv[1]

for subdir in sorted(os.listdir(root)):
    subdir_path = os.path.join(root, subdir)
    if not os.path.isdir(subdir_path):
        continue
    for filename in sorted(os.listdir(subdir_path)):
        if not filename.endswith(".tsv"):
            continue
        filepath = os.path.join(subdir_path, filename)
        fn_start, fn_end = parse_filename_dates(filename)
        if fn_start is None:
            print(f"{filename}: cannot parse dates from filename")
            continue
        rows = get_data_rows(filepath)
        if not rows:
            print(f"{filename}: no data rows found")
            continue
        data_start = row_date(rows[0])
        data_end   = row_date(rows[-1])
        if data_start != fn_start or data_end != fn_end:
            print(f"{filename}: filename {fn_start}–{fn_end}, data {data_start}–{data_end}")
        if is_na(rows[0]):
            print(f"{filename}: first row has NA value")
        if is_na(rows[-1]):
            print(f"{filename}: last row has NA value")