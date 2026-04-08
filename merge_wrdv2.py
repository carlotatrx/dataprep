# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
merge_wrdv2.py

For each station+variable combination, finds consecutive WRDv2-2 segments,
merges their data, writes the merged file to a final/ directory, copies all
originals to a tmp/ backup directory, and logs all merge operations to a
changes log file.

Directory layout (non-dry-run):
  TMP_DIR   = /scratch3/PALAEO-RA/daily_data/tmp/merge_WRDv2-2/
                All original input files are copied here under their original
                names, then removed from WORK_DIR.
  FINAL_DIR = /scratch3/PALAEO-RA/daily_data/final/<station>/
                Only the merged output file is written here.
  LOG_FILE  = /scratch2/ccorbella/code/dataprep/changes_log.txt
                Merge actions (not single-file skips) are appended here,
                with a station header block at the start of each run.

Overlap deduplication: when two consecutive files share dates, the
observations from the FIRST (older) file are kept for the overlap period,
and the duplicate rows from the second file are dropped.

Only the Meta line is updated in the output header; Lat/Lon/Alt are taken
from the file with the longest date coverage.

Usage:
    python3 merge_wrdv2.py /path/to/directory
    python3 merge_wrdv2.py /path/to/directory --dry-run

Options:
    --dry-run    Print what would happen without making any changes.
                 No files are written or copied, no log is updated.
"""

import os
import re
import sys
import shutil
from collections import defaultdict
from datetime import datetime
from io import StringIO

import pandas as pd

DRY_RUN = '--dry-run' in sys.argv

# ---------------------------------------------------------------------------
# fixed paths
# ---------------------------------------------------------------------------

TMP_DIR    = '/scratch3/PALAEO-RA/daily_data/tmp/merge_WRDv2-2'
FINAL_BASE = '/scratch3/PALAEO-RA/daily_data/final'
LOG_FILE   = '/scratch2/ccorbella/code/dataprep/changes_log.txt'

# ---------------------------------------------------------------------------
# parse directory argument
# ---------------------------------------------------------------------------

args = [a for a in sys.argv[1:] if not a.startswith('--')]

if len(args) == 0:
    WORK_DIR = '.'
elif len(args) == 1:
    WORK_DIR = args[0]
else:
    print("Usage: python3 merge_wrdv2.py /path/to/directory [--dry-run]")
    sys.exit(1)

if not os.path.isdir(WORK_DIR):
    print(f"Error: '{WORK_DIR}' is not a valid directory.")
    sys.exit(1)

print(f"Working directory: {os.path.abspath(WORK_DIR)}")
if DRY_RUN:
    print("[DRY RUN MODE — no files will be written or copied]\n")

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def full_path(fname):
    return os.path.join(WORK_DIR, fname)


def parse_header_and_data(fpath):
    """Return (header_lines, meta_dict, data_df) for a SEF file."""
    with open(fpath, encoding='utf-8') as f:
        lines = f.readlines()
    meta = {}
    data_start = None
    for i, line in enumerate(lines):
        parts = line.rstrip('\n').split('\t')
        if len(parts) >= 2:
            meta[parts[0]] = parts[1]
        if i > 0 and parts[0] == 'Year':
            data_start = i
            break
    header_lines = lines[:data_start]
    df = pd.read_csv(StringIO(''.join(lines[data_start:])), sep='\t')

    def safe_date(row):
        try:
            return pd.Timestamp(year=int(row['Year']), month=int(row['Month']),
                                day=int(row['Day']))
        except Exception:
            return pd.NaT

    df['_date'] = df.apply(safe_date, axis=1)
    df = df.dropna(subset=['_date'])
    return header_lines, meta, df


def is_daily(df):
    """True if the dataframe has at most one distinct Hour value (or no Hour col)."""
    if 'Hour' not in df.columns:
        return True
    return df['Hour'].dropna().nunique() <= 1


def date_str(ts):
    return ts.strftime('%Y%m%d')


def gap_days(end_date, start_date):
    return (start_date - end_date).days


def deduplicate_overlap(df_prev, df_curr):
    """
    Keep all rows from df_prev; drop rows from df_curr that fall within
    the overlap period. Returns (df_prev, df_curr_clean, info_str or None).
    """
    overlap_start = df_curr['_date'].min()
    overlap_end   = df_prev['_date'].max()

    if overlap_start > overlap_end:
        return df_prev, df_curr, None

    n_before = len(df_curr)
    df_curr_clean = df_curr[df_curr['_date'] > overlap_end].copy()
    n_dropped = n_before - len(df_curr_clean)

    info = (
        f"OVERLAP {overlap_start.date()} to {overlap_end.date()}: "
        f"{n_dropped} rows dropped from the later file "
        f"(kept earlier file's data for this period)"
    )
    return df_prev, df_curr_clean, info


def build_new_meta_line(group):
    """Build the new Meta field value for a merged group."""
    segments = []
    for fi in group:
        m = fi['meta']
        alias = None
        meta_parts = m.get('Meta', '')
        for part in meta_parts.split('|'):
            part = part.strip()
            if part.startswith('Alias='):
                alias = part[len('Alias='):]
        segments.append({
            'alias': alias,
            'lat':   m.get('Lat'),
            'lon':   m.get('Lon'),
            'alt':   m.get('Alt'),
            'start': date_str(fi['actual_start']),
            'end':   date_str(fi['actual_end']),
        })

    locs_change    = len(set((s['lat'], s['lon']) for s in segments)) > 1
    aliases_change = len(set(s['alias'] for s in segments)) > 1

    first_meta = group[0]['meta'].get('Meta', '')
    # preserve all existing fields except Alias and Coords (we rebuild those)
    base_fields = [p.strip() for p in first_meta.split('|')
                   if p.strip()
                   and not p.strip().startswith('Alias=')
                   and not p.strip().startswith('Coords=')]

    if aliases_change:
        alias_str = '; '.join(
            f"{s['alias']} ({s['start']}-{s['end']})" if s['alias']
            else f"unknown ({s['start']}-{s['end']})"
            for s in segments
        )
        base_fields.append(f"Alias={alias_str}")
    else:
        if segments[0]['alias']:
            base_fields.append(f"Alias={segments[0]['alias']}")

    if locs_change:
        coords_str = '; '.join(
            f"({s['lat']},{s['lon']},{s['alt']}m,{s['start']}-{s['end']})"
            for s in segments
        )
        base_fields.append(f"Coords={coords_str}")
    else:
        s = segments[0]
        base_fields.append(f"Coords=({s['lat']},{s['lon']},{s['alt']}m)")

    return ' | '.join(f for f in base_fields if f)


def append_to_log(station, log_lines):
    """Append a station block to the changes log file."""
    os.makedirs(os.path.dirname(LOG_FILE), exist_ok=True)
    with open(LOG_FILE, 'a', encoding='utf-8') as lf:
        lf.write('\n')
        lf.write('=' * 70 + '\n')
        lf.write(f'{station.upper()}\n')
        lf.write('=' * 70 + '\n')
        lf.write(f'Run: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n\n')
        for line in log_lines:
            lf.write(line + '\n')
        lf.write('\n')


# ---------------------------------------------------------------------------
# discover WRDv2-2 files
# ---------------------------------------------------------------------------

pattern = re.compile(
    r'^(WRDv2-2)_([^_]+)_(\d{8})-(\d{8})_(.+)\.tsv$'
)

groups = defaultdict(list)

for fname in sorted(os.listdir(WORK_DIR)):
    m = pattern.match(fname)
    if not m:
        continue
    project, station, date_start, date_end, var_suffix = m.groups()
    groups[(station, var_suffix)].append({
        'file': fname,
        'date_start_str': date_start,
        'date_end_str':   date_end,
    })

if not groups:
    print(f"No WRDv2-2 files found in '{WORK_DIR}'.")
    sys.exit(0)

stations_in_run = sorted(set(station for (station, _) in groups.keys()))

# ---------------------------------------------------------------------------
# process each group
# ---------------------------------------------------------------------------

GAP_THRESHOLD = 365  # days — gaps larger than this split into separate series

# accumulate log lines per station (only populated on actual runs)
station_log = defaultdict(list)

for (station, var_suffix), file_list in sorted(groups.items()):
    file_list.sort(key=lambda x: x['date_start_str'])

    if len(file_list) == 1:
        print(f"\n{'='*70}")
        print(f"SINGLE FILE -- no merge needed: {file_list[0]['file']}")
        print(f"  (no changes required)")
        continue

    # load all files
    loaded = []
    for fi in file_list:
        header_lines, meta, df = parse_header_and_data(full_path(fi['file']))
        loaded.append({
            'file':         fi['file'],
            'header_lines': header_lines,
            'meta':         meta,
            'df':           df,
            'actual_start': df['_date'].min(),
            'actual_end':   df['_date'].max(),
        })

    # split into merge-groups based on gap threshold and daily/subdaily mismatch
    merge_groups = [[loaded[0]]]
    for i in range(1, len(loaded)):
        prev = merge_groups[-1][-1]
        curr = loaded[i]
        gap = gap_days(prev['actual_end'], curr['actual_start'])
        period_mismatch = is_daily(prev['df']) != is_daily(curr['df'])
        if gap > GAP_THRESHOLD or period_mismatch:
            merge_groups.append([curr])
        else:
            merge_groups[-1].append(curr)

    # process each merge group
    for gi, group in enumerate(merge_groups):
        print(f"\n{'='*70}")

        if len(group) == 1:
            fi = group[0]
            print(f"SINGLE FILE (isolated by gap/period) -- no merge needed:")
            print(f"  {fi['file']}")
            continue

        # overlap detection
        overlap_warnings = []
        for i in range(len(group) - 1):
            gap = gap_days(group[i]['actual_end'], group[i+1]['actual_start'])
            if gap < 0:
                overlap_warnings.append((i, i+1, abs(gap)))

        # new filename and paths
        new_start  = date_str(group[0]['actual_start'])
        new_end    = date_str(group[-1]['actual_end'])
        new_name   = f"WRDv2-2_{station}_{new_start}-{new_end}_{var_suffix}.tsv"
        final_dir  = os.path.join(FINAL_BASE, station)
        final_path = os.path.join(final_dir, new_name)

        longest       = max(group, key=lambda x: x['df']['_date'].nunique())
        new_meta_line = build_new_meta_line(group)

        print(f"MERGE: {station} | {var_suffix}")
        print(f"  Files to merge ({len(group)}):")
        for fi in group:
            print(f"    {fi['file']}")
        print(f"  New filename : {new_name}")
        print(f"  Final path   : {final_path}")
        print(f"  Tmp backup   : {TMP_DIR}/")

        if overlap_warnings:
            print(f"\n  OVERLAP WARNINGS:")
            for i, j, days in overlap_warnings:
                print(f"    Files {i+1} and {i+2} overlap by {days} days "
                      f"({group[i]['actual_end'].date()} to "
                      f"{group[j]['actual_start'].date()})")
                print(f"    -> Rows from the later file in this period will be DROPPED")

        print(f"\n  File operations:")
        print(f"    copy originals (original names) -> {TMP_DIR}/")
        print(f"    remove originals from           -> {WORK_DIR}/")
        print(f"    write merged file            -> {final_path}")

        print(f"\n  New Meta line:")
        print(f"    Meta\t{new_meta_line}")

        if not DRY_RUN:
            # -- backup all originals to tmp under their original names, then remove --
            os.makedirs(TMP_DIR, exist_ok=True)
            for fi in group:
                src = full_path(fi['file'])
                dst = os.path.join(TMP_DIR, fi['file'])
                shutil.copy2(src, dst)
                print(f"  [BACKUP] {fi['file']} -> {dst}")
            for fi in group:
                os.remove(full_path(fi['file']))
                print(f"  [REMOVED] {fi['file']} from working directory")

            # -- concatenate with overlap deduplication --
            all_dfs   = [group[0]['df'].copy()]
            dedup_log = []
            for i in range(1, len(group)):
                df_prev = all_dfs[-1]
                df_curr = group[i]['df'].copy()
                df_prev, df_curr_clean, info = deduplicate_overlap(df_prev, df_curr)
                if info:
                    dedup_log.append(info)
                all_dfs[-1] = df_prev
                all_dfs.append(df_curr_clean)

            combined = pd.concat(all_dfs, ignore_index=True)
            combined = combined.drop(columns=['_date'])

            if dedup_log:
                print(f"\n  Deduplication applied:")
                for msg in dedup_log:
                    print(f"    {msg}")

            # -- build updated header: only Meta line is changed --
            new_header = []
            for line in group[0]['header_lines']:
                key = line.split('\t')[0]
                if key == 'Meta':
                    new_header.append(f"Meta\t{new_meta_line}\n")
                elif key == 'Lat':
                    new_header.append(f"Lat\t{longest['meta'].get('Lat')}\n")
                elif key == 'Lon':
                    new_header.append(f"Lon\t{longest['meta'].get('Lon')}\n")
                elif key == 'Alt':
                    new_header.append(f"Alt\t{longest['meta'].get('Alt')}\n")
                else:
                    new_header.append(line)

            # -- write merged file to final dir --
            os.makedirs(final_dir, exist_ok=True)
            with open(final_path, 'w', encoding='utf-8') as out:
                out.writelines(new_header)
                out.write(combined.to_csv(sep='\t', index=False))
            print(f"\n  [DONE] Merged file written to {final_path}")

            # -- verify filename dates match actual data start/end --
            _, _, df_check = parse_header_and_data(final_path)
            actual_start_check = date_str(df_check['_date'].min())
            actual_end_check   = date_str(df_check['_date'].max())
            if actual_start_check != new_start or actual_end_check != new_end:
                print(f"\n  [ERROR] Filename date mismatch in {new_name}!")
                print(f"    Filename says : {new_start} to {new_end}")
                print(f"    Data contains : {actual_start_check} to {actual_end_check}")
                print(f"  Aborting. Originals are safe in {TMP_DIR}")
                sys.exit(1)
            print(f"  [CHECK] Filename dates verified: {new_start} to {new_end} ✓")

            # -- accumulate log entry for this merge --
            entry = [
                f"MERGE: {station} | {var_suffix}",
                f"  Merged {len(group)} files:",
            ]
            for fi in group:
                entry.append(f"    {fi['file']}")
            entry.append(f"  Output: {final_path}")
            entry.append(f"  New Meta: {new_meta_line}")
            if dedup_log:
                for msg in dedup_log:
                    entry.append(f"  {msg}")
            entry.append("")
            station_log[station].extend(entry)

        else:
            print(f"\n  [DRY RUN] No files were written or copied.")

# ---------------------------------------------------------------------------
# write log (non-dry-run only, once per station at the end)
# ---------------------------------------------------------------------------

if not DRY_RUN:
    for station in stations_in_run:
        if station_log[station]:
            append_to_log(station, station_log[station])
            print(f"\n[LOG] Changes for {station} appended to {LOG_FILE}")

print(f"\n{'='*70}")
print("Done.")
if not DRY_RUN:
    print(f"  Originals backed up in : {TMP_DIR}")
    print(f"  Merged files written to: {FINAL_BASE}/<station>/")
    print(f"  Changes logged in      : {LOG_FILE}")
    print("  Remember to manually verify the Meta line in each merged file.")
