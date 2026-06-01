## using find
SEARCH_TERM="reclido"
find . -type f -iname "*${SEARCH_TERM}*" | awk -F'[_/]' '{
    # find the field that matches the date pattern YYYYMMDD-YYYYMMDD
    for (i=1; i<=NF; i++) {
        if ($i ~ /^[0-9]{8}-[0-9]{8}$/) {
            split($i, dates, "-")
            start = dates[1]
            end   = dates[2]
            var   = $(i+1)
            if (!(var in min_s) || start < min_s[var]) min_s[var] = start
            if (!(var in max_e) || end   > max_e[var]) max_e[var] = end
        }
    }
}
END {
    for (v in min_s) printf "%-6s  %s  %s\n", v, min_s[v], max_e[v]
}' | sort

## using grep
SEARCH_TERM="CHIMES"
grep -ri $SEARCH_TERM */* | grep -oP '^[^:]+' | awk -F'[_/]' '{
    for (i=1; i<=NF; i++) {
        if ($i ~ /^[0-9]{8}-[0-9]{8}$/) {
            split($i, dates, "-")
            start = dates[1]
            end   = dates[2]
            var   = $(i+1)
            if (!(var in min_s) || start < min_s[var]) min_s[var] = start
            if (!(var in max_e) || end   > max_e[var]) max_e[var] = end
        }
    }
}
END {
    for (v in min_s) printf "%-6s  %s  %s\n", v, min_s[v], max_e[v]
}' | sort

#########################################################################################
######## ispd_in_histdaily.txt ##########################################################
# First, extract unique station names from the last column
stations=$(awk '{print $NF}' ../tmp/ispd_in_histdaily.txt | sort -u)

echo "=== PRESSURE SERIES FOUND ==="
all_results=""

for station in $stations; do
    # find matching directories (case-insensitive, prefix match)
    matches=$(find . -maxdepth 1 -type d -iname "${station}*")
    
    if [ -z "$matches" ]; then
        echo "NO DIRECTORY: $station" >&2
        continue
    fi
    
    # find pressure files in matching dirs
    results=$(find $matches -type f -iname "*_p_*" -o -type f -iname "*_p.*" | \
        grep -iP '_p_(daily|subdaily)' | \
        awk -F'[_/]' '{
            for (i=1; i<=NF; i++) {
                if ($i ~ /^[0-9]{8}-[0-9]{8}$/) {
                    split($i, dates, "-")
                    print dates[1], dates[2], $0
                }
            }
        }')
    
    all_results="${all_results}${results}"$'\n'
done

echo "$all_results" | grep -v '^$' | sort -k1,1

echo ""
echo "=== MISSING DIRECTORIES ==="
for station in $stations; do
    matches=$(find . -maxdepth 1 -type d -iname "${station}*")
    if [ -z "$matches" ]; then
        echo "$station"
    fi
done