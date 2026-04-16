#!/usr/bin/env python3
"""
SEF Spike Detective
===================
Scans a database of SEF files organized as:
    database_root/
        STATION_NAME/
            file1.tsv
            file2.tsv
            ...

For each SEF file it extracts:
  - Station ID, Name, Lat, Lon, Variable (Vbl), Source, Link
  - The year range of the actual data rows

Then it:
  1. Counts how many active series exist per year, per variable
  2. Computes year-over-year changes to detect spikes
  3. For user-specified spike windows (e.g. 1855-1870, 1875-1890),
     lists exactly which stations START and STOP in those windows
  4. Groups results by Source to identify bulk database ingestions
  5. Exports a CSV of all series metadata for further analysis

Usage:
    python filename.py /scratch3/PALAEO-RA/daily_data/final/

    Or edit the CONFIG section below.
"""

import os
import sys
import csv
from collections import defaultdict, Counter
from pathlib import Path

class Tee:
    """Duplicate all stdout to both console and a log file."""
    def __init__(self, filepath, mode='w'):
        self.file = open(filepath, mode, encoding='utf-8')
        self.stdout = sys.stdout

    def write(self, data):
        self.stdout.write(data)
        self.file.write(data)

    def flush(self):
        self.stdout.flush()
        self.file.flush()

    def close(self):
        self.file.close()
        sys.stdout = self.stdout
        
# ============================================================
# CONFIG - Edit these to match your needs
# ============================================================

# Path to the root of your SEF database
DATABASE_ROOT = sys.argv[1] if len(sys.argv) > 1 else "/scratch3/PALAEO-RA/daily_data/final/"

# Variables to investigate (SEF Vbl codes for temperature)
# Common SEF variable codes:
#   ta = air temperature (mean), tx = max temp, tn = min temp
#   p = pressure, rr = precipitation, dd = wind direction, ff = wind speed
# Set to None to analyze ALL variables
TARGET_VARIABLES = None  # e.g. ['ta', 'tx', 'tn'] or None for all

# Spike windows to investigate (start_year, end_year) inclusive
SPIKE_WINDOWS = [
    (1855, 1870, "1860s spike"),
    (1875, 1895, "1880s spike"),
    # Add more windows as needed:
    # (1770, 1785, "Palatinate Society spike"),
]

# Minimum year-over-year change to flag as a spike
SPIKE_THRESHOLD = 5  # number of series appearing/disappearing in one year

# Output files
OUTPUT_DIR = "excel_files"

# ============================================================
# SEF PARSER
# ============================================================

def parse_sef_file(filepath):
    """
    Parse a SEF file and extract header metadata + year range.
    Returns a dict with metadata, or None if parsing fails.
    """
    metadata = {
        'filepath': str(filepath),
        'folder': filepath.parent.name,
    }
    
    header_keys = {
        'SEF', 'ID', 'Name', 'Lat', 'Lon', 'Alt',
        'Source', 'Link', 'Vbl', 'Stat', 'Units', 'Meta'
    }
    
    try:
        with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
            lines = f.readlines()
    except Exception as e:
        metadata['error'] = str(e)
        return metadata
    
    data_start = None
    for i, line in enumerate(lines):
        line = line.strip()
        if not line:
            continue
        
        # Check if this is a header line (Key\tValue)
        parts = line.split('\t')
        if len(parts) >= 2 and parts[0] in header_keys:
            metadata[parts[0]] = parts[1].strip()
        
        # Detect the data header row
        if line.startswith('Year') and 'Month' in line:
            data_start = i + 1
            break
    
    if data_start is None:
        metadata['error'] = 'No data header found'
        return metadata
    
    # Extract year range from data rows
    years = []
    for line in lines[data_start:]:
        line = line.strip()
        if not line:
            continue
        parts = line.split('\t')
        try:
            year = int(parts[0])
            years.append(year)
        except (ValueError, IndexError):
            continue
    
    if years:
        metadata['year_start'] = min(years)
        metadata['year_end'] = max(years)
        metadata['n_records'] = len(years)
        
        # Also get the set of unique years for gap detection
        unique_years = sorted(set(years))
        metadata['unique_years'] = unique_years
        metadata['n_unique_years'] = len(unique_years)
        
        # Detect gaps > 1 year
        gaps = []
        for j in range(1, len(unique_years)):
            gap = unique_years[j] - unique_years[j-1]
            if gap > 1:
                gaps.append((unique_years[j-1], unique_years[j], gap))
        metadata['gaps'] = gaps
    else:
        metadata['error'] = 'No valid data rows'
    
    return metadata


# ============================================================
# MAIN ANALYSIS
# ============================================================

def main():
    root = Path(DATABASE_ROOT)
    
    if not root.exists():
        print(f"ERROR: Database root not found: {root}")
        print(f"\nUsage: python {sys.argv[0]}  /scratch3/PALAEO-RA/daily_data/final/")
        sys.exit(1)
    
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    log_path = os.path.join(OUTPUT_DIR, "spike_detective_report.txt")
    tee = Tee(log_path)
    sys.stdout = tee
        
    print(f"SEF Spike Detective")
    print(f"{'='*60}")
    print(f"Database: {root}")
    print(f"Scanning for SEF files...\n")
    
    # ----------------------------------------------------------
    # 1. Scan all SEF files
    # ----------------------------------------------------------
    all_series = []
    n_files = 0
    n_errors = 0
    
    station_dirs = sorted([d for d in root.iterdir() if d.is_dir()])
    print(f"Found {len(station_dirs)} station folders")
    
    for i, station_dir in enumerate(station_dirs):
        if (i + 1) % 100 == 0:
            print(f"  Scanning station {i+1}/{len(station_dirs)}: {station_dir.name}")
        
        for sef_file in sorted(station_dir.iterdir()):
            if sef_file.is_file():
                n_files += 1
                meta = parse_sef_file(sef_file)
                if 'error' in meta:
                    n_errors += 1
                all_series.append(meta)
    
    print(f"\nScanned {n_files} files ({n_errors} with errors)")
    
    # Filter to valid series
    valid_series = [s for s in all_series if 'year_start' in s]
    print(f"Valid series with data: {len(valid_series)}")
    
    # Filter by target variables if specified
    if TARGET_VARIABLES:
        filtered = [s for s in valid_series if s.get('Vbl') in TARGET_VARIABLES]
        print(f"Filtered to variables {TARGET_VARIABLES}: {len(filtered)} series")
    else:
        filtered = valid_series
    
    # ----------------------------------------------------------
    # 2. Count active series per year, per variable
    # ----------------------------------------------------------
    print(f"\n{'='*60}")
    print("ACTIVE SERIES COUNT PER YEAR")
    print(f"{'='*60}\n")
    
    # Determine year range
    all_starts = [s['year_start'] for s in filtered]
    all_ends = [s['year_end'] for s in filtered]
    if not all_starts:
        print("No valid series found!")
        sys.exit(1)
    
    year_min = min(all_starts)
    year_max = max(all_ends)
    
    # Group by variable
    by_variable = defaultdict(list)
    for s in filtered:
        vbl = s.get('Vbl', 'unknown')
        by_variable[vbl].append(s)
    
    print(f"Variables found: {dict((k, len(v)) for k, v in sorted(by_variable.items()))}\n")
    
    # For each variable, count active series per year
    # A series is "active" in a year if year_start <= year <= year_end
    # But also track series that HAVE DATA in that specific year (using unique_years)
    
    yearly_counts = {}  # {vbl: {year: count}}
    yearly_series = {}  # {vbl: {year: [list of series]}}
    
    for vbl, series_list in by_variable.items():
        yearly_counts[vbl] = defaultdict(int)
        yearly_series[vbl] = defaultdict(list)
        
        for s in series_list:
            # Use unique_years for precise counting (not just start-end range)
            for y in s.get('unique_years', []):
                yearly_counts[vbl][y] += 1
                yearly_series[vbl][y].append(s)
    
    # ----------------------------------------------------------
    # 3. Detect year-over-year spikes
    # ----------------------------------------------------------
    print(f"{'='*60}")
    print(f"YEAR-OVER-YEAR SPIKES (threshold: {SPIKE_THRESHOLD}+ series change)")
    print(f"{'='*60}\n")
    
    spike_report = []
    
    for vbl in sorted(yearly_counts.keys()):
        counts = yearly_counts[vbl]
        if not counts:
            continue
        
        years_sorted = sorted(counts.keys())
        
        print(f"\n--- Variable: {vbl} ({len(by_variable[vbl])} total series) ---")
        
        found_spikes = False
        for j in range(1, len(years_sorted)):
            y_prev = years_sorted[j-1]
            y_curr = years_sorted[j]
            
            # Only look at consecutive years
            if y_curr - y_prev > 1:
                continue
            
            change = counts[y_curr] - counts[y_prev]
            
            if abs(change) >= SPIKE_THRESHOLD:
                found_spikes = True
                direction = "UP" if change > 0 else "DOWN"
                print(f"  {y_prev} -> {y_curr}: {counts[y_prev]:4d} -> {counts[y_curr]:4d}  "
                      f"(change: {change:+d} {direction})")
                
                spike_report.append({
                    'variable': vbl,
                    'year_from': y_prev,
                    'year_to': y_curr,
                    'count_from': counts[y_prev],
                    'count_to': counts[y_curr],
                    'change': change,
                })
        
        if not found_spikes:
            print(f"  No spikes above threshold found.")
    
    # ----------------------------------------------------------
    # 4. Investigate specific spike windows
    # ----------------------------------------------------------
    print(f"\n\n{'='*60}")
    print("SPIKE WINDOW INVESTIGATION")
    print(f"{'='*60}")
    
    for win_start, win_end, win_label in SPIKE_WINDOWS:
        print(f"\n{'─'*60}")
        print(f"Window: {win_label} ({win_start}-{win_end})")
        print(f"{'─'*60}")
        
        for vbl in sorted(yearly_counts.keys()):
            series_list = by_variable[vbl]
            
            # Series that START within this window
            starting = [s for s in series_list
                        if win_start <= s['year_start'] <= win_end]
            
            # Series that END within this window
            ending = [s for s in series_list
                      if win_start <= s['year_end'] <= win_end]
            
            if not starting and not ending:
                continue
            
            print(f"\n  Variable: {vbl}")
            
            if starting:
                print(f"    Series STARTING in {win_start}-{win_end}: {len(starting)}")
                
                # Group by Source
                by_source = defaultdict(list)
                for s in starting:
                    by_source[s.get('Source', 'unknown')].append(s)
                
                for source, slist in sorted(by_source.items(), 
                                            key=lambda x: -len(x[1])):
                    print(f"      Source: {source} ({len(slist)} series)")
                    # Show first few station names
                    names = sorted(set(s.get('Name', s.get('ID', '?')) 
                                      for s in slist))
                    for name in names[:10]:
                        # Find this series to show details
                        example = next(s for s in slist 
                                      if s.get('Name', s.get('ID', '?')) == name)
                        print(f"        - {name} "
                              f"({example['year_start']}-{example['year_end']}) "
                              f"[{example.get('folder', '?')}]")
                    if len(names) > 10:
                        print(f"        ... and {len(names)-10} more")
            
            if ending:
                print(f"    Series ENDING in {win_start}-{win_end}: {len(ending)}")
                
                by_source = defaultdict(list)
                for s in ending:
                    by_source[s.get('Source', 'unknown')].append(s)
                
                for source, slist in sorted(by_source.items(), 
                                            key=lambda x: -len(x[1])):
                    print(f"      Source: {source} ({len(slist)} series)")
                    names = sorted(set(s.get('Name', s.get('ID', '?')) 
                                      for s in slist))
                    for name in names[:10]:
                        example = next(s for s in slist 
                                      if s.get('Name', s.get('ID', '?')) == name)
                        print(f"        - {name} "
                              f"({example['year_start']}-{example['year_end']}) "
                              f"[{example.get('folder', '?')}]")
                    if len(names) > 10:
                        print(f"        ... and {len(names)-10} more")
    
    # ----------------------------------------------------------
    # 5. Source summary across ALL series
    # ----------------------------------------------------------
    print(f"\n\n{'='*60}")
    print("SOURCE SUMMARY (all variables)")
    print(f"{'='*60}\n")
    
    source_stats = defaultdict(lambda: {
        'count': 0, 'variables': set(), 
        'year_min': 9999, 'year_max': 0,
        'stations': set()
    })
    
    for s in valid_series:
        src = s.get('Source', 'unknown')
        source_stats[src]['count'] += 1
        source_stats[src]['variables'].add(s.get('Vbl', '?'))
        source_stats[src]['year_min'] = min(source_stats[src]['year_min'], 
                                            s['year_start'])
        source_stats[src]['year_max'] = max(source_stats[src]['year_max'], 
                                            s['year_end'])
        source_stats[src]['stations'].add(s.get('Name', s.get('ID', '?')))
    
    print(f"{'Source':<30s} {'#Series':>8s} {'#Stations':>10s} "
          f"{'Years':>15s} {'Variables'}")
    print(f"{'─'*90}")
    
    for src, stats in sorted(source_stats.items(), 
                             key=lambda x: -x[1]['count']):
        vbls = ','.join(sorted(stats['variables']))
        print(f"{src:<30s} {stats['count']:>8d} {len(stats['stations']):>10d} "
              f"{stats['year_min']:>6d}-{stats['year_max']:<6d}  {vbls}")
    
    # ----------------------------------------------------------
    # 6. Export full metadata CSV
    # ----------------------------------------------------------
    csv_path = os.path.join(OUTPUT_DIR, "all_series_metadata.csv")
    csv_fields = [
        'folder', 'ID', 'Name', 'Lat', 'Lon', 'Alt', 'Vbl', 'Stat',
        'Units', 'Source', 'Link', 'year_start', 'year_end', 
        'n_records', 'n_unique_years', 'filepath'
    ]
    
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=csv_fields, extrasaction='ignore')
        writer.writeheader()
        for s in sorted(valid_series, key=lambda x: (x.get('Vbl',''), 
                                                      x.get('year_start', 0))):
            writer.writerow(s)
    
    print(f"\n\nFull metadata exported to: {csv_path}")
    
    # ----------------------------------------------------------
    # 7. Export yearly counts CSV (for re-plotting)
    # ----------------------------------------------------------
    counts_path = os.path.join(OUTPUT_DIR, "yearly_series_counts.csv")
    all_years = sorted(set(
        y for vbl_counts in yearly_counts.values() 
        for y in vbl_counts.keys()
    ))
    all_vbls = sorted(yearly_counts.keys())
    
    with open(counts_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['Year'] + all_vbls + ['Total'])
        for year in all_years:
            row = [year]
            total = 0
            for vbl in all_vbls:
                c = yearly_counts[vbl].get(year, 0)
                row.append(c)
                total += c
            row.append(total)
            writer.writerow(row)
    
    print(f"Yearly counts exported to: {counts_path}")
    
    # ----------------------------------------------------------
    # 8. Export spike window details CSV
    # ----------------------------------------------------------
    detail_path = os.path.join(OUTPUT_DIR, "spike_window_series.csv")
    detail_fields = [
        'window', 'event', 'folder', 'ID', 'Name', 'Vbl', 'Source', 
        'Link', 'year_start', 'year_end', 'Lat', 'Lon'
    ]
    
    with open(detail_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=detail_fields, extrasaction='ignore')
        writer.writeheader()
        
        for win_start, win_end, win_label in SPIKE_WINDOWS:
            for s in valid_series:
                if win_start <= s['year_start'] <= win_end:
                    row = dict(s)
                    row['window'] = win_label
                    row['event'] = 'START'
                    writer.writerow(row)
                if win_start <= s['year_end'] <= win_end:
                    row = dict(s)
                    row['window'] = win_label
                    row['event'] = 'END'
                    writer.writerow(row)
    
    print(f"Spike window details exported to: {detail_path}")
    
    print(f"\n{'='*60}")
    print(f"\n{'='*60}")
    print(f"DONE! Check the {OUTPUT_DIR}/ folder for CSV exports.")
    print(f"Full report saved to: {log_path}")
    print(f"{'='*60}")

    tee.close()

if __name__ == '__main__':
    main()