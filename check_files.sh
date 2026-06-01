# count number of flagged values

find /scratch3/PALAEO-RA/daily_data/final -type f -name "*.tsv" | while read f; do
  awk 'NR>1 && /^[0-9]/ && $7!="NA" {
    total++
    meta=$NF
    # extract qc value
    if (match(meta, /qc=([^[:space:]]+)/, arr)) {
      qc=arr[1]
      # not flagged: qc=Y, qc=0, qc=1, qc=5, no qc at all
      if (qc != "Y" && qc != "0" && qc != "1" && qc != "5" && qc != "NA") {
        flagged++
      }
    }
  }
  END {print total, flagged}' "$f"
done | awk '{t+=$1; f+=$2} END {
  printf "Total values:   %d\nFlagged values: %d\nPercentage:     %.2f%%\n", t, f, (f/t)*100
}'

# show wrongly flagged values:
find /scratch3/PALAEO-RA/daily_data/final -type f -name "*.tsv" | while read f; do
  result=$(awk 'NR>1 && /^[0-9]/ {print $NF}' "$f" | \
    grep -oP 'qc=[^[:space:]]+' | \
    sed 's/qc=//' | \
    tr ';' '\n' | \
    grep -P '^(\d+|[0-9]\|PTC|manually_acc[ep]+ted|.*PTC=Y.*)$' | \
    sort | uniq -c)
  [ -n "$result" ] && echo "=== $f ===" && echo "$result"
done

# see what different types of qc there are:
find /scratch3/PALAEO-RA/daily_data/final -type f -name "*.tsv" | while read f; do
  awk 'NR>1 && /^[0-9]/ {print $NF}' "$f"
done | \
grep -oP 'qc=[^[:space:]]+' | \
sed 's/qc=//' | \
tr ';' '\n' | \
sort | uniq -c | sort -rn

# find total number of flagged values
find /scratch3/PALAEO-RA/daily_data/final -type f -name "*.tsv" | while read f; do
  awk 'NR>1 && /^[0-9]/ && $7!="NA" {
    total++
    meta=$NF
    # extract qc value
    if (match(meta, /qc=([^[:space:]]+)/, arr)) {
      qc=arr[1]
      # not flagged: qc=Y, qc=0, qc=1, qc=5, no qc at all
      if (qc != "Y" && qc != "0" && qc != "1" && qc != "5" && qc != "NA") {
        flagged++
      }
    }
  }
  END {print total, flagged}' "$f"
done | awk '{t+=$1; f+=$2} END {
  printf "Total values:   %d\nFlagged values: %d\nPercentage:     %.2f%%\n", t, f, (f/t)*100
}'

# find total number of NA values
find /scratch3/PALAEO-RA/daily_data/final -type f -name "*.tsv" | while read f; do
  awk 'NR>1 && /^[0-9]/ && $7=="NA" {count++} END {print count+0}' "$f"
done | awk '{sum+=$1} END {print "NA values: " sum}'

# find total number of values
find /scratch3/PALAEO-RA/daily_data/final -type f -name "*.tsv" | while read f; do
  awk 'NR>1 && /^[0-9]/ {print $7}' "$f"
done | awk '{total++; if($1=="NA") na++} END {
  printf "Total values:    %d\nNA values:       %d\nNon-NA values:   %d\n", total, na, total-na
}'