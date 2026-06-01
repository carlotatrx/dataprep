#!/usr/bin/env bash
# Usage: bash strip_trailing_na.sh <root_directory> [--apply]
#
# For each SEF .tsv file, removes consecutive data rows at the END of the
# series whose Value column is NA, stopping at the last real value.
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
    awk -F'\t' '
        function is_na(v) {
            gsub(/^[ \t]+|[ \t]+$/, "", v)
            v = toupper(v)
            return (v=="NA" || v=="NAN" || v=="N/A" || v=="" || v=="NULL" || v=="MISSING")
        }
        BEGIN { in_data=0; n=0; last_good=0 }
        {
            if (!in_data) {
                print > "/dev/stdout"
                if ($1=="Year" && $2=="Month") in_data=1
                next
            }
            # buffer every data line
            n++
            lines[n]=$0
            # remember last line index that has a real value
            if (NF>=7 && !is_na($7)) {
                last_good=n
                good_date=$1"-"$2"-"$3
            }
        }
        END {
            # print kept data rows
            for (i=1; i<=last_good; i++) print lines[i] > "/dev/stdout"

            removed = n - last_good
            if (removed>0) {
                # date of the first removed row and the last removed row
                split(lines[last_good+1], a, "\t"); first_rm=a[1]"-"a[2]"-"a[3]
                split(lines[n],          b, "\t"); last_rm =b[1]"-"b[2]"-"b[3]
                printf("REMOVED %d  new_end=%s  first_removed=%s  old_end=%s",
                       removed, good_date, first_rm, last_rm) > "/dev/stderr"
            }
        }
    ' "$filepath" 2>/tmp/_na_summary >/tmp/_na_cleaned

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

echo "Total: $COUNT file(s) with trailing NA rows"
$APPLY || echo "(dry run — re-run with --apply to actually remove them)"