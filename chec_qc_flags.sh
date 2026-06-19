# Shared awk logic: classifies a qc token as a real flag or not.
# Non-flags: numeric codes (external QC), Y, NA, manually_accepted, parallel_record, indoor.
AWK_FLAG_LOGIC='
  function is_real_flag(tok) {
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", tok)
    if (tok ~ /^[0-9]+$/) return 0
    if (tok == "Y" || tok == "NA") return 0
    if (tok == "manually_accepted" || tok == "parallel_record" || tok == "indoor") return 0
    return 1
  }
  function row_is_flagged(meta,    n, parts, i, val, toks, nt, j) {
    n = split(meta, parts, /[[:space:]]\|[[:space:]]/)
    for (i=1; i<=n; i++) {
      if (parts[i] ~ /^qc=/) {
        val = substr(parts[i], 4)
        gsub(/\r$/, "", val)
        nt = split(val, toks, ";")
        for (j=1; j<=nt; j++) if (is_real_flag(toks[j])) return 1
      }
    }
    return 0
  }
'

# Total non-NA values and flagged values across the whole dataset
find /scratch3/PALAEO-RA/daily_data/final -type f -name "*.tsv" | while read f; do
  awk -F'\t' "$AWK_FLAG_LOGIC"'
    FNR>1 && /^[0-9]/ && $7!="NA" {
      total++
      if (row_is_flagged($NF)) flagged++
    }
    END { print total+0, flagged+0 }
  ' "$f"
done | awk '{t+=$1; f+=$2} END {
  printf "Total non-NA values:   %d\nFlagged values:         %d\nPercentage:             %.2f%%\n", t, f, (f/t)*100
}'

# Total non-NA values:   28077006
# Flagged values:         82771
# Percentage:             0.29%

# Per-variable totals and flagged counts
find /scratch3/PALAEO-RA/daily_data/final -type f -name "*.tsv" | while read f; do
  awk -F'\t' "$AWK_FLAG_LOGIC"'
    /^Vbl\t/ { vbl = $2; gsub(/\r$/, "", vbl) }
    FNR>1 && /^[0-9]/ && $7!="NA" {
      total++
      if (row_is_flagged($NF)) flagged++
    }
    END { print vbl, total+0, flagged+0 }
  ' "$f"
done | awk '{
  t[$1]+=$2; f[$1]+=$3
} END {
  printf "%-6s %14s %14s %10s\n", "Var", "Total", "Flagged", "Pct"
  for (v in t) printf "%-6s %14d %14d %9.2f%%\n", v, t[v], f[v], (f[v]/t[v])*100
}'

# Var             Total        Flagged        Pct
# w              254066           2745      1.08%
# eee             16236            111      0.68%
# dd            1191400          38148      3.20%
# ta            7403947          19816      0.27%
# rrt              7448              0      0.00%
# rr            7998148           2389      0.03%
# p             8639984          13400      0.16%
# fs              13951             28      0.20%
# Tx             983528           3271      0.33%
# rh             561889            472      0.08%
# Tn            1006409           2391      0.24%


# Flagged and own-digitized (Source contains WeaR or PALAEO) breakdown
find /scratch3/PALAEO-RA/daily_data/final -type f -name "*.tsv" | while read f; do
  awk -F'\t' "$AWK_FLAG_LOGIC"'
    /^Source\t/ {
      src = tolower($2)
      if (src ~ /wear/ || src ~ /palaeo/) is_own = 1
    }
    FNR>1 && /^[0-9]/ && $7!="NA" {
      total++
      if (is_own) own_total++
      if (row_is_flagged($NF)) {
        flagged++
        if (is_own) own_flagged++
      }
    }
    END { print total+0, own_total+0, flagged+0, own_flagged+0 }
  ' "$f"
done | awk '{t+=$1; ot+=$2; f+=$3; of+=$4} END {
  printf "Total non-NA values:                %d\n", t
  printf "Own-digitized non-NA values:        %d (%.2f%%)\n", ot, (ot/t)*100
  printf "Flagged values:                     %d (%.2f%%)\n", f, (f/t)*100
  printf "Flagged from own-digitized:         %d (%.2f%% of flagged)\n", of, (of/f)*100
}'

# Total non-NA values:                28077006
# Own-digitized non-NA values:        1724444 (6.14%)
# Flagged values:                     82771 (0.29%)
# Flagged from own-digitized:         11832 (14.29% of flagged)

# Frequency table of real flag types across the whole dataset
find /scratch3/PALAEO-RA/daily_data/final -type f -name "*.tsv" | while read f; do
  awk -F'\t' "$AWK_FLAG_LOGIC"'
    FNR>1 && /^[0-9]/ && $7!="NA" {
      meta=$NF
      n = split(meta, parts, /[[:space:]]\|[[:space:]]/)
      for (i=1; i<=n; i++) {
        if (parts[i] ~ /^qc=/) {
          val = substr(parts[i], 4)
          gsub(/\r$/, "", val)
          nt = split(val, toks, ";")
          for (j=1; j<=nt; j++) if (is_real_flag(toks[j])) print toks[j]
        }
      }
    }
  ' "$f"
done | sort | uniq -c | sort -rn


  # 37379 duplicate_columns
  # 34898 subdaily_repetition
  # 10287 wmo_time_consistency
  # 10001 duplicate_times
  #  8452 daily_repetition
  #  3493 climatic_outliers
  #  1229 wmo_gross_errors
  #   405 duplicate_dates
  #   278 plot_subdaily
  #   258 subdaily_out_of_range
  #   144 Tn>Tx
  #    65 impossible_values
  #    56 outlier
  #    50 in_sunlight
  #    39 internal_consistency
  #    22 inconsistent_with_2nd_thermometer
  #    18 date-time_doubtful
  #    17 temporal_coherence
  #     2 time_uncertain
  #     2 spatially_inconsistent
