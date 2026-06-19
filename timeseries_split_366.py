#!/usr/bin/env python3
"""
split_timeseries.py
===================================================================
Split SEF (.tsv) station files into separate continuous time series.

DEFINITION
-----------
A *series* is a run of observations at one station with no internal
gap longer than GAP_THRESHOLD days (366 by default). Wherever a file
contains a gap larger than that, it actually holds more than one
series and should be split into one file per series.

WHAT THIS SCRIPT DOES
---------------------
For every .tsv under DATA_DIR/<station>/ it:
  1. Reads the metaheader (everything up to and including the
     "Year  Month  Day ..." column-header line) and the data rows.
  2. Walks the data rows in file order and starts a new series each
     time the gap from the previous valid date exceeds GAP_THRESHOLD.
  3. POST-1900 RULE: once a gap would open a series starting in the
     20th century (year > LAST_19C_YEAR), it stops splitting. All
     remaining rows -- including any further post-1900 gaps -- are
     appended to that final 19th-century series instead of becoming
     their own files. This keeps every presented series rooted before
     1900 while preserving the original tail of the record.
  4. Writes one file per series, reusing the full metaheader verbatim.

NAMING
------
Each output reuses the original filename with its YYYYMMDD-YYYYMMDD
date span replaced by the span of that series, using the ACTUAL first
and last observation dates in the series (the final merged series
therefore carries its true 20th-century end date).
  e.g. WRDv2-2_Yarmouth_18610301-18750331_ta_barometer_subdaily.tsv
    -> WRDv2-2_Yarmouth_18610301-18631231_ta_barometer_subdaily.tsv
       WRDv2-2_Yarmouth_18720104-18750331_ta_barometer_subdaily.tsv

IN-PLACE REPLACEMENT
--------------------
Splits are written into the SAME station directory as the original,
and the original multi-series file is deleted. Files that contain a
single series are left untouched. Encoding is auto-detected (utf-8,
then latin-1) and reused on write so bytes round-trip losslessly.

DRY RUN
-------
With DRY_RUN = True nothing is written or deleted; the full plan is
only logged. A detailed report of every file and the series it would
become is always written to REPORT_FILE.

USAGE
-----
  1. Set DATA_DIR to the dataset root and REPORT_FILE if desired.
  2. Run with DRY_RUN = True:        python3 split_timeseries.py
  3. Review REPORT_FILE (timeseries_to_split.txt).
  4. Set DRY_RUN = False and rerun to apply the split in place.

The 19th/20th-century boundary lives in one place, LAST_19C_YEAR.
A gap whose far side lands in (LAST_19C_YEAR + 1) or later triggers
the tail-merge; set it to 1900 if you want 1900 itself to merge.
===================================================================
"""

import os
import re
import datetime
import pandas as pd

# ----------------------------- configuration -----------------------------
DATA_DIR = "/scratch3/PALAEO-RA/daily_data/final/"
REPORT_FILE = "timeseries_to_split.txt"   # detailed plan / log
GAP_THRESHOLD = 366                       # days; gaps longer than this break a series
LAST_19C_YEAR = 1899                      # a gap opening a series in a later year merges the tail
DRY_RUN = True                            # True = preview + report only; False = apply in place
# --------------------------------------------------------------------------

DATE_SPAN_RE = re.compile(r"\d{8}-\d{8}")


def read_lines(fpath):
    """Read a file as text, trying utf-8 then latin-1. Returns (lines, encoding)."""
    for enc in ("utf-8", "latin-1"):
        try:
            with open(fpath, "r", encoding=enc) as f:
                return f.readlines(), enc
        except UnicodeDecodeError:
            continue
    raise UnicodeDecodeError("codec", b"", 0, 1, f"could not decode {fpath}")


def plan_split(fpath, report):
    """Split one file's data rows into series.

    Returns (header, encoding, segments) or None if the file holds a
    single series (nothing to split) or has no column header.
    Each segment is a dict: {"lines": [...], "first": ts, "last": ts}.
    Raw data lines are preserved verbatim; rows with invalid dates
    (e.g. Feb 30) are kept in place and never trigger a split.
    """
    lines, enc = read_lines(fpath)

    # metaheader = everything up to and including the "Year  Month ..." row
    header_end = None
    for idx, line in enumerate(lines):
        if line.split("\t")[0].strip() == "Year":
            header_end = idx
            break
    if header_end is None:
        report.write(f"SKIP (no column header): {os.path.basename(fpath)}\n")
        return None

    header = lines[:header_end + 1]
    data = lines[header_end + 1:]

    segments = []
    current = {"lines": [], "first": None, "last": None}
    prev_date = None
    merging_tail = False   # once a 20th-century series would start, stop splitting

    for raw in data:
        parts = raw.strip().split("\t")
        ts = None
        if len(parts) >= 3 and parts[0].lstrip("-").isdigit():
            try:
                ts = pd.Timestamp(int(parts[0]), int(parts[1]), int(parts[2]))
            except ValueError:
                ts = None   # invalid calendar date: keep the line, no date

        if ts is not None and prev_date is not None and not merging_tail:
            if (ts - prev_date).days > GAP_THRESHOLD:
                if ts.year > LAST_19C_YEAR:
                    merging_tail = True          # append the rest, never split again
                else:
                    segments.append(current)
                    current = {"lines": [], "first": None, "last": None}

        current["lines"].append(raw)
        if ts is not None:
            if current["first"] is None:
                current["first"] = ts
            current["last"] = ts
            prev_date = ts

    if current["lines"]:
        segments.append(current)

    if len(segments) <= 1:
        return None
    return header, enc, segments


def process(fpath, report):
    """Plan and (unless DRY_RUN) apply the split for one file in place.

    Returns (files_split, series_made) counts for the run summary.
    """
    plan = plan_split(fpath, report)
    if plan is None:
        return 0, 0
    header, enc, segments = plan

    base = os.path.basename(fpath)
    if not DATE_SPAN_RE.search(base):
        report.write(f"SKIP (no date span in name): {base}\n")
        return 0, 0

    out_dir = os.path.dirname(fpath)
    tag = "[DRY] " if DRY_RUN else ""
    report.write(f"{tag}SPLIT {base} -> {len(segments)} series\n")

    for seg in segments:
        # naming uses the ACTUAL first/last observation dates of the series
        span = f"{seg['first'].strftime('%Y%m%d')}-{seg['last'].strftime('%Y%m%d')}"
        new_name = DATE_SPAN_RE.sub(span, base, count=1)
        marker = "  (+20th-c. tail)" if seg["last"].year > LAST_19C_YEAR else ""
        report.write(f"   {seg['first'].date()} -> {seg['last'].date()}  "
                     f"{len(seg['lines']):>6} rows  {new_name}{marker}\n")
        if not DRY_RUN:
            with open(os.path.join(out_dir, new_name), "w", encoding=enc) as f:
                f.writelines(header)
                f.writelines(seg["lines"])

    if DRY_RUN:
        report.write(f"   [DRY] would delete original: {base}\n")
    else:
        os.remove(fpath)
    return 1, len(segments)


def main():
    files_split = series_made = 0
    with open(REPORT_FILE, "w", encoding="utf-8") as report:
        report.write("# Time-series split plan\n")
        report.write(f"# generated: {datetime.datetime.now():%Y-%m-%d %H:%M:%S}\n")
        report.write(f"# mode: {'DRY RUN (no changes written)' if DRY_RUN else 'APPLIED IN PLACE'}\n")
        report.write(f"# data dir: {DATA_DIR}\n")
        report.write(f"# gap threshold: {GAP_THRESHOLD} days | last 19th-c. start year: {LAST_19C_YEAR}\n\n")

        for station in sorted(os.listdir(DATA_DIR)):
            sdir = os.path.join(DATA_DIR, station)
            if not os.path.isdir(sdir):
                continue
            for fname in sorted(os.listdir(sdir)):
                if fname.endswith(".tsv"):
                    fs, sm = process(os.path.join(sdir, fname), report)
                    files_split += fs
                    series_made += sm

        summary = (f"\n{'[DRY] ' if DRY_RUN else ''}"
                   f"{files_split} files split -> {series_made} series\n")
        report.write(summary)

    print(f"Wrote report to {REPORT_FILE}")
    print(summary.strip())


if __name__ == "__main__":
    main()
    