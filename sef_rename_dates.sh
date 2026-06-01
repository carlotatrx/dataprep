#!/usr/bin/env bash
# Usage: bash rename_files.sh <root_directory> [--dry-run]
#
# Renames SEF .tsv files whose filename dates don't match the actual data dates.
# With --dry-run (or no second argument), just prints what would happen.

ROOT="${1:?Usage: $0 <root_directory> [--dry-run]}"
DRY_RUN=true
[[ "$2" == "--apply" ]] && DRY_RUN=false

$DRY_RUN && echo "=== DRY RUN — no files will be changed ===" || echo "=== APPLYING RENAMES ==="
echo ""

COUNT=0
for filepath in "$ROOT"/*/*.tsv; do
    filename=$(basename "$filepath")
    dir=$(dirname "$filepath")

    # Extract dates from filename (format: _YYYYMMDD-YYYYMMDD_)
    if [[ ! "$filename" =~ _([0-9]{8})-([0-9]{8})_ ]]; then continue; fi
    fn_start="${BASH_REMATCH[1]}"
    fn_end="${BASH_REMATCH[2]}"

    # Find the data header line and read first/last data rows
    data_line=$(grep -n "^Year	" "$filepath" | cut -d: -f1)
    [[ -z "$data_line" ]] && continue

    first_row=$(awk -v n="$((data_line+1))" 'NR==n' "$filepath")
    last_row=$(tail -n 1 "$filepath")

    [[ -z "$first_row" || -z "$last_row" ]] && continue

    # Parse year/month/day from first and last rows (tab-separated)
    IFS=$'\t' read -r yr mo dy _ <<< "$first_row"
    data_start=$(printf "%04d%02d%02d" "$yr" "$mo" "$dy")

    IFS=$'\t' read -r yr mo dy _ <<< "$last_row"
    data_end=$(printf "%04d%02d%02d" "$yr" "$mo" "$dy")

    # Compare
    if [[ "$data_start" == "$fn_start" && "$data_end" == "$fn_end" ]]; then continue; fi

    new_filename="${filename/_${fn_start}-${fn_end}_/_${data_start}-${data_end}_}"
    COUNT=$((COUNT + 1))

    if $DRY_RUN; then
        echo "  $filename"
        echo "→ $new_filename"
        echo ""
    else
        mv "$dir/$filename" "$dir/$new_filename"
        echo "Renamed: $filename → $new_filename"
    fi
done

echo "Total: $COUNT file(s) to rename"