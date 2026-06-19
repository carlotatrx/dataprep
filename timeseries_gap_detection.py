import os
import pandas as pd

DATA_DIR = "/scratch3/PALAEO-RA/daily_data/final/"
GAP_THRESHOLD = 366  # days

for station in os.listdir(DATA_DIR):
    station_dir = os.path.join(DATA_DIR, station)
    if not os.path.isdir(station_dir):
        continue

    for fname in os.listdir(station_dir):
        if not fname.endswith(".tsv"):
            continue
        fpath = os.path.join(station_dir, fname)

        # Read with utf-8, fall back to latin-1
        try:
            with open(fpath, "r", encoding="utf-8") as f:
                lines = f.readlines()
        except UnicodeDecodeError:
            try:
                with open(fpath, "r", encoding="latin-1") as f:
                    lines = f.readlines()
            except Exception as e:
                print(f"SKIP (unreadable): {fpath} — {e}")
                continue

        # Parse data rows
        rows = []
        for line in lines:
            parts = line.strip().split("\t")
            if len(parts) >= 3 and parts[0].lstrip("-").isdigit():
                try:
                    y, m, d = int(parts[0]), int(parts[1]), int(parts[2])
                    rows.append((y, m, d))
                except ValueError:
                    continue

        if len(rows) < 2:
            continue

        # Build sorted unique dates, skipping invalid ones
        dates_set = set()
        for y, m, d in rows:
            try:
                dates_set.add(pd.Timestamp(y, m, d))
            except ValueError:
                continue
        dates = sorted(dates_set)

        if len(dates) < 2:
            continue

        # Check consecutive gaps
        for i in range(1, len(dates)):
            if dates[i - 1].year > 1900: # don't flag any gaps after 1900
                continue
            delta = (dates[i] - dates[i - 1]).days
            if delta > GAP_THRESHOLD:
                print(f"GAP {delta}d | {dates[i-1].date()} → {dates[i].date()} | {fname}")
