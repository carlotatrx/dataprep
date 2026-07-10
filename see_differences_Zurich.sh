#!/bin/bash

# Define file names
SMALL_FILE="PALAEO-RA_Zurich_Weiss_18231201-18240630_ta_subdaily_qc.tsv"
BIG_FILE="GCOS_Zurich_Weiss_18231201-18350818_ta_subdaily_qc.tsv"

awk -F'\t' -v tolerance=5 '
# 1. Process the Big File (FNR == NR means first file)
NR == FNR {
    if (FNR > 12) {
        date = $1 "-" $2 "-" $3
        # Store values for each date in a string separated by "|"
        big_data[date] = (big_data[date] == "" ? $7 : big_data[date] "|" $7)
    }
    next
}

# 2. Process the Small File
FNR > 12 {
    date = $1 "-" $2 "-" $3
    val = $7
    found_date = 0
    match_found = 0
    
    if (date in big_data) {
        found_date = 1
        split(big_data[date], big_vals, "|")
        for (i in big_vals) {
            diff = val - big_vals[i]
            if (diff < 0) diff = -diff
            if (diff <= tolerance) {
                match_found = 1
                break
            }
        }
    }

    if (!found_date) {
        missing_dates[date] = 1
    } else if (match_found) {
        same_count++
    } else {
        diff_count++
        # Save the mismatch details to show later
        diff_report = diff_report "Date: " date " | Small: " val " | Big Vals: [" big_data[date] "]\n"
    }
}

# 3. Print Results
END {
    printf "\n--- Statistics ---\n"
    printf "Records with matching values (+/- %d): %d\n", tolerance, same_count + 0
    printf "Records with different values: %d\n", diff_count + 0
    
    if (diff_count > 0) {
        print "\n--- Details of Different Entries ---"
        printf "%s", diff_report
    } else if (same_count > 0) {
        print "(All overlapping records match within the tolerance.)"
    }

    print "\n--- Dates in Small File missing from Big File ---"
    # Sort and print missing dates
    for (d in missing_dates) print d | "sort -t- -k1,1n -k2,2n -k3,3n"
}' "$BIG_FILE" "$SMALL_FILE"
