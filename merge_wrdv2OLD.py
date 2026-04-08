# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
merge_wrdv2.py

For each station+variable combination, finds consecutive WRDv2-2 segments,
merges their data into the oldest file (renamed with the full date span),
renames the remaining files to 2.tsv, 3.tsv etc. for easy deletion,
and prints the new metahead for manual review.

Overlap deduplication: when two consecutive files share dates, the
observations from the FIRST (older) file are kept for the overlap period,
and the duplicate rows from the second file are dropped. A warning is
printed with the overlap details.

Usage:
    python3 merge_wrdv2.py /path/to/directory
    python3 merge_wrdv2.py /path/to/directory --dry-run

Options:
    --dry-run    Print what would happen without making any changes.
"""

import os
import re
import sys
from collections import defaultdict
from io import StringIO

import pandas as pd

DRY_RUN = '--dry-run' in sys.argv

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

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def full_path(fname):
    return os.path.join(WORK_DIR, fname)


def parse_header_and_data(fname):
    """Return (header_lines, meta_dict, data_df) for a SEF file."""
    with open(full_path(fname), encoding='utf-8') as f:
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


def date_str(ts):
    return ts.strftime('%Y%m%d')


def gap_days(end_date, start_date):
    return (start_date - end_date).days


def deduplicate_overlap(df_prev, df_curr):
    """
    When df_prev and df_curr have overlapping dates, keep all rows from
    df_prev and drop rows from df_curr that fall within the overlap period.
    Returns (df_prev, df_curr_clean, overlap_info_str or None).
    """
    overlap_start = df_curr['_date'].min()
    overlap_end = df_prev['_date'].max()

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
        'date_end_str': date_end,
    })

if not groups:
    print(f"No WRDv2-2 files found in '{WORK_DIR}'.")
    sys.exit(0)

# ---------------------------------------------------------------------------
# process each group
# ---------------------------------------------------------------------------

GAP_THRESHOLD = 365  # days -- groups with a gap larger than this are split

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
        header_lines, meta, df = parse_header_and_data(fi['file'])
        loaded.append({
            'file': fi['file'],
            'header_lines': header_lines,
            'meta': meta,
            'df': df,
            'actual_start': df['_date'].min(),
            'actual_end': df['_date'].max(),
        })

    # obs per day to detect daily vs subdaily
    # obs_per_day = [round(len(x['df']) / x['df']['_date'].nunique(), 1) for x in loaded]

    # split into merge-groups based on gap threshold and period consistency
    merge_groups = [[loaded[0]]]
    for i in range(1, len(loaded)):
        prev = merge_groups[-1][-1]
        curr = loaded[i]
        gap = gap_days(prev['actual_end'], curr['actual_start'])
        prev_is_daily = prev['df']['Hour'].dropna().nunique() <= 1 if 'Hour' in prev['df'].columns else True
        curr_is_daily = curr['df']['Hour'].dropna().nunique() <= 1 if 'Hour' in curr['df'].columns else True
        # prev_is_daily = obs_per_day[loaded.index(prev)] <= 1.1
        # curr_is_daily = obs_per_day[loaded.index(curr)] <= 1.1
        period_mismatch = prev_is_daily != curr_is_daily
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

        # check and report overlaps between consecutive files in the group
        overlap_warnings = []
        for i in range(len(group) - 1):
            gap = gap_days(group[i]['actual_end'], group[i+1]['actual_start'])
            if gap < 0:
                overlap_warnings.append((i, i+1, abs(gap)))

        # build new filename
        new_start = date_str(group[0]['actual_start'])
        new_end = date_str(group[-1]['actual_end'])
        new_name = f"WRDv2-2_{station}_{new_start}-{new_end}_{var_suffix}.tsv"
        new_path = full_path(new_name)

        print(f"MERGE: {station} | {var_suffix}")
        print(f"  Files to merge ({len(group)}):")
        for fi in group:
            print(f"    {fi['file']}")
        print(f"  New filename: {new_name}")

        if overlap_warnings:
            print(f"\n  OVERLAP WARNINGS:")
            for i, j, days in overlap_warnings:
                print(f"    Files {i+1} and {i+2} overlap by {days} days "
                      f"({group[i]['actual_end'].date()} to "
                      f"{group[j]['actual_start'].date()})")
                print(f"    -> Rows from the later file in this period will be DROPPED")

        # mv commands (shown with full paths for clarity)
        print(f"\n  mv commands:")
        print(f"    mv {full_path(group[0]['file'])} {new_path}")
        for idx, fi in enumerate(group[1:], start=2):
            # print(f"    mv {full_path(fi['file'])} {full_path(str(idx) + '.tsv')}")
            print(f"    mv {full_path(fi['file'])} {full_path(str(idx) + '-' + var_suffix + '.tsv')}")

        # build new metahead
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
                'lat': m.get('Lat'),
                'lon': m.get('Lon'),
                'alt': m.get('Alt'),
                'start': date_str(fi['actual_start']),
                'end': date_str(fi['actual_end']),
                'meta_raw': meta_parts,
            })

        locs_change = len(set((s['lat'], s['lon']) for s in segments)) > 1
        aliases_change = len(set(s['alias'] for s in segments)) > 1

        first_meta = group[0]['meta'].get('Meta', '')
        base_fields = [p.strip() for p in first_meta.split('|')
                       if not p.strip().startswith('Alias=')]

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

        new_meta_line = ' | '.join(f for f in base_fields if f)

        longest = max(group, key=lambda x: x['df']['_date'].nunique())

        print(f"\n  New metahead (for manual update):")
        print(f"  ---------------------------------")
        for line in group[0]['header_lines']:
            key = line.split('\t')[0]
            if key == 'Lat':
                print(f"    Lat\t{longest['meta'].get('Lat')}")
            elif key == 'Lon':
                print(f"    Lon\t{longest['meta'].get('Lon')}")
            elif key == 'Alt':
                print(f"    Alt\t{longest['meta'].get('Alt')}")
            elif key == 'Meta':
                print(f"    Meta\t{new_meta_line}")
            else:
                print(f"    {line.rstrip()}")

        # perform file operations
        if not DRY_RUN:
            # concatenate with overlap deduplication
            all_dfs = [group[0]['df'].copy()]
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

            # build updated header
            new_header = []
            for line in group[0]['header_lines']:
                key = line.split('\t')[0]
                if key == 'Lat':
                    new_header.append(f"Lat\t{longest['meta'].get('Lat')}\n")
                elif key == 'Lon':
                    new_header.append(f"Lon\t{longest['meta'].get('Lon')}\n")
                elif key == 'Alt':
                    new_header.append(f"Alt\t{longest['meta'].get('Alt')}\n")
                else:
                    new_header.append(line)

            with open(full_path(group[0]['file']), 'w', encoding='utf-8') as out:
                out.writelines(new_header)
                out.write(combined.to_csv(sep='\t', index=False))

            os.rename(full_path(group[0]['file']), new_path)
            print(f"\n  [DONE] Renamed and data merged into {new_name}")

            for idx, fi in enumerate(group[1:], start=2):
                # os.rename(full_path(fi['file']), full_path(f"{idx}.tsv"))
                os.rename(full_path(fi['file']), full_path(f"{idx}-{var_suffix}.tsv"))
                print(f"  [DONE] Renamed {fi['file']} -> {idx}.tsv")
        else:
            print(f"\n  [DRY RUN] No files were changed.")

print(f"\n{'='*70}")
print("Done. Remember to manually update the Meta line in each merged file.")
