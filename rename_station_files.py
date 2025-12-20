#!/usr/bin/env python3

# Usage: python rename_sef.py /path/to/root_folder

import os
import re
import sys

ALLOWED_VARS = {"ta", "p", "eee", "fs", "w", "dd", "tb", "rr", "rh"}
DATE_RE = re.compile(r"^\d{8}-\d{8}$")

def get_start_end_dates(filepath):
    """Read SEF file and return (start_yyyymmdd, end_yyyymmdd)."""
    with open(filepath, "r", encoding="utf-8") as f:
        lines = f.readlines()

    # find data header line: "Year\tMonth\tDay..."
    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith("Year\tMonth\tDay"):
            header_idx = i
            break
    if header_idx is None:
        raise ValueError(f"No data header found in {filepath}")

    header_cols = lines[header_idx].rstrip("\n").split("\t")
    try:
        y_i = header_cols.index("Year")
        m_i = header_cols.index("Month")
        d_i = header_cols.index("Day")
    except ValueError:
        raise ValueError(f"Year/Month/Day not found in header of {filepath}")

    data_lines = lines[header_idx + 1 :]

    # drop empty lines at the end
    while data_lines and not data_lines[-1].strip():
        data_lines.pop()
    if not data_lines:
        raise ValueError(f"No data lines in {filepath}")

    first = data_lines[0].rstrip("\n").split("\t")
    last = data_lines[-1].rstrip("\n").split("\t")

    def fmt_date(fields):
        y = int(fields[y_i])
        m = int(fields[m_i])
        d = int(fields[d_i])
        return f"{y:04d}{m:02d}{d:02d}"

    return fmt_date(first), fmt_date(last)


def process_root(root):
    for dirpath, _, filenames in os.walk(root):
        for fname in filenames:
            if not fname.endswith(".tsv"):
                continue
            if not fname.endswith("_qc.tsv"):
                continue

            base = fname[:-4]  # drop .tsv
            parts = base.split("_")

            # find variable position
            var_idx = None
            for i, p in enumerate(parts):
                if p in ALLOWED_VARS:
                    var_idx = i
                    break
            if var_idx is None:
                continue  # not a SEF we care about

            # check if there is already a date range before the var
            has_date = any(DATE_RE.match(p) for p in parts[:var_idx])
            if has_date:
                continue  # already has yyyymmdd-yyyymmdd

            old_path = os.path.join(dirpath, fname)

            try:
                start, end = get_start_end_dates(old_path)
            except Exception as e:
                print(f"SKIP (error reading dates) {old_path}: {e}")
                continue

            date_segment = f"{start}-{end}"
            new_parts = parts[:var_idx] + [date_segment] + parts[var_idx:]
            new_base = "_".join(new_parts)
            new_fname = new_base + ".tsv"
            new_path = os.path.join(dirpath, new_fname)

            if old_path == new_path:
                continue

            if os.path.exists(new_path):
                print(f"SKIP (target exists) {new_path}")
                continue

            print(f"RENAMING:\n  {old_path}\n  -> {new_path}")
            os.rename(old_path, new_path)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} ROOT_DIR")
        sys.exit(1)
    process_root(sys.argv[1])
    
