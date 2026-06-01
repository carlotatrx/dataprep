#!/usr/bin/env python3
"""
check_dd_directions.py

Audit wind-direction (dd) SEF files: verify that the numeric Value (degrees)
matches the original textual direction recorded in the Meta column.

- Recursively finds every file whose name contains "_dd_" under a root dir.
- For each data row, reads the original direction from Meta. The original may
  be written as `orig=...` or `orig.dd=...` (but NOT `orig.time=...`).
- If a row has no original direction, it is skipped (not checked).
- Recognises both English (N, NW, W, WNW, ...) and French (N, NO, O, ONO, ...)
  16-point compass abbreviations; they map to the same degrees by position and
  do not collide (English uses W, French uses O).
- Tokens that are not standard compass abbreviations (e.g. named local winds
  such as "Solano", "Gallego") cannot be verified automatically and are
  reported separately as "unrecognized" for manual review, NOT as errors.

Read-only: this script never modifies any file. It prints a summary and
writes a CSV of every issue found.
"""

import argparse
import csv
import os
import re
import sys

# --- Compass abbreviation -> degrees (English + French 16-point rose) -------
# English and French agree on N/E/S letters; they differ only on West (W vs O),
# and those never produce conflicting degree values.
DIR_TO_DEG = {
    # cardinal / shared
    'N': 0, 'NNE': 22.5, 'NE': 45, 'ENE': 67.5,
    'E': 90, 'ESE': 112.5, 'SE': 135, 'SSE': 157.5,
    'S': 180,
    # English west half
    'SSW': 202.5, 'SW': 225, 'WSW': 247.5,
    'W': 270, 'WNW': 292.5, 'NW': 315, 'NNW': 337.5,
    # French west half (O = Ouest)
    'SSO': 202.5, 'SO': 225, 'OSO': 247.5,
    'O': 270, 'ONO': 292.5, 'NO': 315, 'NNO': 337.5,
}

# Capture orig= or orig.dd= but never orig.time= (the '=' must follow orig or orig.dd)
ORIG_RE = re.compile(r'orig(?:\.dd)?=([^|]+)')


def find_dd_files(root):
    out = []
    for dirpath, _dirs, files in os.walk(root):
        for name in files:
            if '_dd_' in name and name.lower().endswith(('.tsv', '.txt', '.sef')):
                out.append(os.path.join(dirpath, name))
    return sorted(out)


def parse_sef(path):
    """Return (header_index, value_col, meta_col, lines) or raise ValueError."""
    with open(path, encoding='utf-8') as fh:
        lines = fh.read().splitlines()
    hdr = None
    for i, line in enumerate(lines):
        if line.startswith('Year\t'):
            hdr = i
            break
    if hdr is None:
        raise ValueError("no 'Year\\t...' header row found")
    cols = lines[hdr].split('\t')
    try:
        value_col = cols.index('Value')
        meta_col = cols.index('Meta')
    except ValueError:
        raise ValueError("header missing Value or Meta column")
    return hdr, value_col, meta_col, lines


def degrees_match(stored, expected, tol):
    """True if stored value equals expected within tol (handles 0/360 + rounding)."""
    try:
        v = float(stored)
    except (TypeError, ValueError):
        return False  # NA / non-numeric never matches a real direction
    candidates = [expected]
    if expected == 0:
        candidates.append(360.0)
    if expected == 360:
        candidates.append(0.0)
    return any(abs(v - c) <= tol for c in candidates)


def check_file(path, tol):
    """Return dict with per-file stats and a list of issue rows."""
    hdr, vcol, mcol, lines = parse_sef(path)
    issues = []
    stats = dict(rows=0, checked=0, ok=0, mismatch=0, missing_value=0,
                 unrecognized=0, no_orig=0)
    for line in lines[hdr + 1:]:
        if not line.strip():
            continue
        stats['rows'] += 1
        c = line.split('\t')
        meta = c[mcol] if len(c) > mcol else ''
        value = c[vcol] if len(c) > vcol else ''
        m = ORIG_RE.search(meta)
        if not m:
            stats['no_orig'] += 1
            continue                      # no original -> don't check
        orig_raw = m.group(1).strip()
        token = orig_raw.upper()
        ymd = (c[0], c[1], c[2], c[3] if len(c) > 3 else '',
               c[4] if len(c) > 4 else '')
        if token not in DIR_TO_DEG:
            stats['unrecognized'] += 1
            issues.append((path, *ymd, orig_raw, value, '', 'unrecognized_orig'))
            continue
        stats['checked'] += 1
        expected = DIR_TO_DEG[token]
        exp_str = str(int(expected) if float(expected).is_integer() else expected)
        if value.strip().upper() in ('', 'NA'):
            stats['missing_value'] += 1
            issues.append((path, *ymd, orig_raw, value, exp_str, 'orig_present_value_missing'))
        elif degrees_match(value, expected, tol):
            stats['ok'] += 1
        else:
            stats['mismatch'] += 1
            issues.append((path, *ymd, orig_raw, value, exp_str, 'value_mismatch'))
    return stats, issues


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('root', help='Root directory to scan recursively')
    ap.add_argument('-o', '--out', default='dd_direction_check_issues.csv',
                    help='CSV file for issues (default: %(default)s)')
    ap.add_argument('--tol', type=float, default=1.0,
                    help='Degree tolerance to absorb rounding of .5 values '
                         '(default: %(default)s)')
    args = ap.parse_args()

    files = find_dd_files(args.root)
    if not files:
        print(f"No '_dd_' files found under {args.root!r}")
        return

    print(f"Scanning {len(files)} dd file(s) under {args.root}\n")
    all_issues = []
    totals = dict(rows=0, checked=0, ok=0, mismatch=0, missing_value=0,
                  unrecognized=0, no_orig=0)
    files_with_problems = 0

    for path in files:
        try:
            stats, issues = check_file(path, args.tol)
        except ValueError as e:
            print(f"  [SKIP] {path}: {e}")
            continue
        all_issues.extend(issues)
        for k in totals:
            totals[k] += stats[k]
        problem = stats['mismatch'] + stats['missing_value']
        flag = ''
        if problem:
            files_with_problems += 1
            flag = f"  <-- {problem} problem row(s)"
        elif stats['unrecognized']:
            flag = f"  ({stats['unrecognized']} unrecognized orig token(s) to review)"
        rel = os.path.relpath(path, args.root)
        print(f"  {rel}: checked={stats['checked']} ok={stats['ok']} "
              f"mismatch={stats['mismatch']} missing={stats['missing_value']} "
              f"unrecognized={stats['unrecognized']}{flag}")

    # write issues CSV
    if all_issues:
        with open(args.out, 'w', newline='', encoding='utf-8') as fh:
            w = csv.writer(fh)
            w.writerow(['file', 'Year', 'Month', 'Day', 'Hour', 'Minute',
                        'orig', 'stored_value', 'expected_value', 'issue'])
            w.writerows(all_issues)

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  files scanned ............ {len(files)}")
    print(f"  files with problems ...... {files_with_problems}")
    print(f"  data rows ................ {totals['rows']}")
    print(f"  rows checked ............. {totals['checked']}")
    print(f"    matching ............... {totals['ok']}")
    print(f"    value mismatch ......... {totals['mismatch']}")
    print(f"    orig present, value NA . {totals['missing_value']}")
    print(f"  unrecognized orig token .. {totals['unrecognized']} (manual review)")
    print(f"  rows without an original . {totals['no_orig']} (skipped)")
    if all_issues:
        print(f"\n  -> {len(all_issues)} issue row(s) written to {args.out}")
    else:
        print("\n  -> No issues found.")


if __name__ == '__main__':
    main()