#!/usr/bin/env bash
# Usage: bash strip_leading_na.sh <root_directory> [--apply]
#
# For each SEF .tsv file, removes consecutive data rows at the START of the
# series whose Value column is NA, stopping at the first real value.
# Rows with NA in the MIDDLE of the series are left untouched.
#
# Without --apply it is a DRY RUN: it only reports what would be removed.

ROOT="${1:?Usage: $0 <root_directory> [--apply]}"
APPLY=false
[[ "$2" == "--apply" ]] && APPLY=true

$APPLY && echo "=== APPLYING — files will be modified ===" || echo "=== DRY RUN — no files will be changed ==="
echo ""

COUNT=0
for filepath in "$ROOT"/*/*.tsv; do
    # awk does the work. It prints, on stderr, a summary of leading-NA rows.
    # On stdout it writes the cleaned file content.
    summary=$(awk -F'\t' '
        function is_na(v) {
            gsub(/^[ \t]+|[ \t]+$/, "", v)
            v = toupper(v)
            return (v=="NA" || v=="NAN" || v=="N/A" || v=="" || v=="NULL" || v=="MISSING")
        }
        BEGIN { in_data=0; started=0; removed=0; first_old=""; first_new="" }
        {
            if (!in_data) {
                print > "/dev/stdout"
                if ($1=="Year" && $2=="Month") in_data=1
                next
            }
            # in data section
            if (!started && NF>=7) {
                if (is_na($7)) {
                    if (first_old=="") first_old=$1"-"$2"-"$3
                    last_removed=$1"-"$2"-"$3
                    removed++
                    next            # skip this leading-NA row
                } else {
                    started=1
                    first_new=$1"-"$2"-"$3
                }
            }
            print > "/dev/stdout"
        }
        END {
            if (removed>0)
                printf("REMOVED %d  old_start=%s  last_removed=%s  new_start=%s",
                       removed, first_old, last_removed, first_new) > "/dev/stderr"
        }
    ' "$filepath" 2>/tmp/_na_summary >/tmp/_na_cleaned)

    info=$(cat /tmp/_na_summary)
    [[ -z "$info" ]] && continue        # nothing removed -> skip

    fname=$(basename "$filepath")
    COUNT=$((COUNT+1))

    if $APPLY; then
        mv /tmp/_na_cleaned "$filepath"
        echo "$fname"
        echo "  $info  [applied]"
    else
        echo "$fname"
        echo "  $info"
    fi
    echo ""
done

echo "Total: $COUNT file(s) with leading NA rows"
$APPLY || echo "(dry run — re-run with --apply to actually remove them)"