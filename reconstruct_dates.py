"""
reconstruct_dates.py

Reconstructs missing (NA) dates in SEF .tsv files, and cleans up junk rows
that carry an impossible calendar date.

NA-date reconstruction
----------------------
Some series contain blocks of consecutive data rows where Year/Month/Day are
NA but a Value is present. When such a block is bounded by real dates and
contains exactly the right number of rows to fill the gap with no missing
days, the dates are reconstructed by counting forward one day at a time from
the last real date before the block.

SAFETY: for every NA-date block the script checks that filling sequentially
lands exactly on the day before the first real date AFTER the block. If it
does not (row count implies hidden internal gaps), date filling for that file
is REFUSED and the NA dates are left as-is.

Impossible dates
----------------
A row whose date is numerically impossible (e.g. 1712-9-31) is handled as
follows:
  * value is NA AND meta column empty  -> the row is DELETED (junk).
  * otherwise                          -> FLAGGED for manual review, never
                                          deleted (it carries real data).

When a file contains ANY unparseable/impossible date, the "hidden gap" block
messages for that file are suppressed (they are unreliable once the date
sequence is broken). The unparseable dates themselves are always reported.

Usage:
    python reconstruct_dates.py <root_directory> [--apply]

Without --apply it is a DRY RUN: reports what would change, writes nothing.
"""

import sys
from datetime import date, timedelta
from pathlib import Path


def is_na(s):
    return s.strip().upper() in ("NA", "NAN", "N/A", "", "NULL", "MISSING")


def parse_file(path):
    """Return (header_lines, data_rows, trailing_newline)."""
    with open(path, encoding="utf-8") as f:
        raw = f.read()
    trailing_newline = raw.endswith("\n")
    lines = raw.split("\n")
    if trailing_newline:
        lines = lines[:-1]

    header, data, in_data = [], [], False
    for line in lines:
        if not in_data:
            header.append(line)
            cols = line.split("\t")
            if cols and cols[0] == "Year" and len(cols) > 1 and cols[1] == "Month":
                in_data = True
            continue
        data.append(line.split("\t"))
    return header, data, trailing_newline


def numeric_impossible(y, m, d):
    """True if y/m/d are integers but do not form a valid calendar date."""
    try:
        yi, mi, di = int(y), int(m), int(d)
    except ValueError:
        return False
    try:
        date(yi, mi, di)
        return False
    except ValueError:
        return True


def reconstruct(path):
    header, data, trailing_newline = parse_file(path)

    r = dict(
        filled=0, blocks=[], recon_ok=True,
        block_errors=[],     # "hidden gap" messages (suppressed if date problems)
        deletions=[],        # (idx, datestr) impossible date + NA value + empty meta
        flags=[],            # (idx, reason) unparseable, NOT safe to delete
        other_errors=[],
        header=header, trailing_newline=trailing_newline,
        delete_idx=set(), new_data=None,
    )

    last_real = None
    expected_after = None
    in_block = False
    blk_idx = None
    blk_count = 0
    blk_start = None

    filled = [list(row) for row in data]   # working copy with fills

    for i, row in enumerate(data):
        if len(row) < 7:
            if in_block:
                r["other_errors"].append(f"row {i}: malformed/blank line inside an NA-date block")
                r["recon_ok"] = False
            continue

        y, m, d = row[0], row[1], row[2]

        # --- NA-date row: candidate for reconstruction ---
        if is_na(y) or is_na(m) or is_na(d):
            if last_real is None:
                r["block_errors"].append(
                    f"row {i}: NA-date block at start of series (no preceding real date) "
                    f"-> cannot reconstruct.")
                r["recon_ok"] = False
                continue
            if not in_block:
                in_block, blk_idx, blk_count = True, i, 0
                blk_start = last_real + timedelta(days=1)
            last_real = last_real + timedelta(days=1)
            blk_count += 1
            expected_after = last_real + timedelta(days=1)
            filled[i][0] = str(last_real.year)
            filled[i][1] = str(last_real.month)
            filled[i][2] = str(last_real.day)
            r["filled"] += 1
            if r["blocks"] and r["blocks"][-1][0] == blk_idx:
                s = r["blocks"][-1]
                r["blocks"][-1] = (s[0], blk_count, blk_start, last_real)
            else:
                r["blocks"].append((blk_idx, blk_count, blk_start, last_real))
            continue

        # --- real-date row: try to parse ---
        try:
            cur = date(int(y), int(m), int(d))
            parsed = True
        except ValueError:
            parsed = False

        if not parsed:
            val = row[6] if len(row) > 6 else ""
            meta = row[7] if len(row) > 7 else ""
            if numeric_impossible(y, m, d) and is_na(val) and meta.strip() == "":
                r["deletions"].append((i, f"{y}\t{m}\t{d}"))
                r["delete_idx"].add(i)
            else:
                reason = f"unparseable date {y}\\t{m}\\t{d}"
                if not is_na(val) or meta.strip() != "":
                    reason += f" (carries data: value={val!r} meta={meta!r}) -> manual review"
                r["flags"].append((i, reason))
            # a bad date cannot anchor or validate a block; end any open block
            # without emitting a (spurious) block error, and refuse filling.
            if in_block:
                in_block = False
                r["recon_ok"] = False
            continue

        # parsed OK
        if in_block:
            if expected_after is not None and cur != expected_after:
                r["block_errors"].append(
                    f"block ending before row {i}: after filling, next real date should be "
                    f"{expected_after} but file has {cur}. Row count implies hidden gaps "
                    f"-> refusing to fill this block.")
                r["recon_ok"] = False
            in_block = False
        last_real = cur

    # --- build output ---
    has_dates_problem = bool(r["deletions"]) or bool(r["flags"])
    do_fill = r["recon_ok"] and r["filled"] > 0
    base = filled if do_fill else [list(row) for row in data]  # fills only if safe

    if r["deletions"] or do_fill:
        r["new_data"] = [row for idx, row in enumerate(base) if idx not in r["delete_idx"]]

    r["has_dates_problem"] = has_dates_problem
    r["do_fill"] = do_fill
    return r


def write_file(path, header, data, trailing_newline):
    text = "\n".join(list(header) + ["\t".join(row) for row in data])
    if trailing_newline:
        text += "\n"
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


def main():
    if len(sys.argv) < 2:
        print("Usage: python reconstruct_dates.py <root_directory> [--apply]")
        return
    root = sys.argv[1]
    apply = len(sys.argv) > 2 and sys.argv[2] == "--apply"

    print("=== APPLYING — files will be modified ===" if apply
          else "=== DRY RUN — no files will be changed ===")
    print()

    n_recon = n_delfiles = n_review = n_clean = 0
    n_delrows = 0
    skipped = []

    for path in sorted(Path(root).rglob("*.tsv")):
        fname = path.name
        try:
            r = reconstruct(str(path))
        except UnicodeDecodeError as e:
            skipped.append((fname, f"not UTF-8 (byte 0x{e.object[e.start]:02x} at position {e.start})"))
            continue
        except Exception as e:
            skipped.append((fname, f"{type(e).__name__}: {e}"))
            continue

        interesting = (r["do_fill"] or r["deletions"] or r["flags"]
                       or r["block_errors"] or r["other_errors"])
        if not interesting:
            n_clean += 1
            continue

        print(fname)

        # impossible-date rows to delete
        if r["deletions"]:
            n_delfiles += 1
            n_delrows += len(r["deletions"])
            verb = "deleted" if apply else "would delete"
            print(f"  {len(r['deletions'])} impossible-date row(s) — {verb} (NA value, empty meta):")
            for idx, ds in r["deletions"]:
                print(f"    - row {idx}: {ds}")

        # unparseable dates needing manual review
        if r["flags"]:
            n_review += 1
            print(f"  {len(r['flags'])} unparseable date(s) — MANUAL REVIEW (not deleted):")
            for idx, reason in r["flags"]:
                print(f"    ! row {idx}: {reason}")

        # block / hidden-gap errors: only when there are NO date problems
        if r["block_errors"] and not r["has_dates_problem"]:
            n_review += 1
            print("  date reconstruction REFUSED (inconsistent blocks):")
            for e in r["block_errors"]:
                print(f"    ! {e}")

        for e in r["other_errors"]:
            print(f"    ! {e}")

        # successful reconstruction
        if r["do_fill"]:
            n_recon += 1
            for (_, count, d0, d1) in r["blocks"]:
                print(f"  block: {count} rows  ->  {d0}  ..  {d1}")
            print(f"  total rows reconstructed: {r['filled']}")

        if apply and r["new_data"] is not None:
            write_file(str(path), r["header"], r["new_data"], r["trailing_newline"])
            print("  [written]")
        print()

    if skipped:
        print(f"--- SKIPPED {len(skipped)} file(s) (could not read — check manually) ---")
        for fn, reason in skipped:
            print(f"  {fn}")
            print(f"      {reason}")
        print()

    print(f"Summary: {n_recon} reconstructed, "
          f"{n_delrows} junk row(s) removed across {n_delfiles} file(s), "
          f"{n_review} file(s) need review, "
          f"{n_clean} clean/no-op, {len(skipped)} skipped")
    if not apply and (n_recon or n_delrows):
        print("(dry run — re-run with --apply to write changes)")


if __name__ == "__main__":
    main()