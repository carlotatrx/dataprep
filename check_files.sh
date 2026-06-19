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
# This is wrong because it only captures the first QC and doesn't filter NA values
find /scratch3/PALAEO-RA/daily_data/final -type f -name "*.tsv" | while read f; do
  awk 'NR>1 && /^[0-9]/ {print $NF}' "$f"
done | \
grep -oP 'qc=[^[:space:]]+' | \
sed 's/qc=//' | \
tr ';' '\n' | \
sort | uniq -c | sort -rn

# find total number of flagged values
find /scratch3/PALAEO-RA/daily_data/final -type f -name "*.tsv" | while read f; do
  awk -F'\t' '
    /^Source\t/ {
      src = tolower($2)
      if (src ~ /wear/ || src ~ /palaeo/) is_own = 1
    }
    FNR>1 && /^[0-9]/ && $7!="NA" {
      total++
      meta=$NF
      if (match(meta, /qc=([^[:space:]]+)/, arr)) {
        qc=arr[1]
        if (qc != "Y" && qc != "0" && qc != "1" && qc != "5" && qc != "NA") {
          flagged++
          if (is_own) own_flagged++
        }
      }
    }
    END {print total+0, flagged+0, own_flagged+0}
  ' "$f"
done | awk '{t+=$1; f+=$2; o+=$3} END {
  printf "Total non-NA values:                 %d\n", t
  printf "Flagged values:                      %d\n", f
  printf "Flagged values from own-digitized:   %d\n", o
  printf "Own-digitized share of flagged:      %.2f%%\n", (o/f)*100
  printf "Own-digitized share of all values:   %.2f%%\n", (o/t)*100
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

# find number of flagged values per variable
find /scratch3/PALAEO-RA/daily_data/final -type f -name "*.tsv" | while read f; do
  awk -F'\t' '
    /^Vbl\t/ { vbl = $2; gsub(/\r$/, "", vbl) }
    FNR>1 && /^[0-9]/ && $7!="NA" {
      total++
      meta=$NF
      if (match(meta, /qc=([^[:space:]]+)/, arr)) {
        qc=arr[1]
        if (qc != "Y" && qc != "0" && qc != "1" && qc != "5" && qc != "NA") {
          flagged++
        }
      }
    }
    END { print vbl, total+0, flagged+0 }
  ' "$f"
done | awk '{
  t[$1] += $2
  f[$1] += $3
} END {
  printf "%-6s %14s %14s %10s\n", "Var", "Total", "Flagged", "Pct"
  for (v in t) {
    printf "%-6s %14d %14d %9.2f%%\n", v, t[v], f[v], (f[v]/t[v])*100
  }
}'

# Per variable, shows what share of all non-NA values (not just flagged ones) come from own-digitized series

find /scratch3/PALAEO-RA/daily_data/final -type f -name "*.tsv" | while read f; do
  awk -F'\t' '
    /^Vbl\t/ { vbl = $2; gsub(/\r$/, "", vbl) }
    /^Source\t/ {
      src = tolower($2)
      if (src ~ /wear/ || src ~ /palaeo/) is_own = 1
    }
    FNR>1 && /^[0-9]/ && $7!="NA" {
      total++
      if (is_own) own_total++
    }
    END { print vbl, total+0, own_total+0 }
  ' "$f"
done | awk '{
  t[$1] += $2
  o[$1] += $3
} END {
  printf "%-6s %14s %14s %10s\n", "Var", "Total", "Own-digitized", "Pct"
  for (v in t) {
    printf "%-6s %14d %14d %9.2f%%\n", v, t[v], o[v], (o[v]/t[v])*100
  }
}'

# Counts occurrences of each qc flag type across the whole dataset

find /scratch3/PALAEO-RA/daily_data/final -type f -name "*.tsv" | while read f; do
  awk -F'\t' '
    FNR>1 && /^[0-9]/ && $7!="NA" {
      meta=$NF
      if (match(meta, /qc=([^[:space:]]+)/, arr)) {
        qc=arr[1]
        gsub(/\r$/, "", qc)
        if (qc != "Y" && qc != "0" && qc != "1" && qc != "5" && qc != "NA") {
          n = split(qc, tests, ";")
          for (i=1; i<=n; i++) print tests[i]
        }
      }
    }
  ' "$f"
done | sort | uniq -c | sort -rn

# Finds qc-test/variable combinations that fall outside each test's documented scope
# (climatic_outliers, wmo_gross_errors, wmo_time_consistency, internal_consistency, temporal_coherence)

find /scratch3/PALAEO-RA/daily_data/final -type f -name "*.tsv" | while read f; do
  awk -F'\t' -v fname="$f" '
    BEGIN {
      valid_for["climatic_outliers"]    = "ta,Tx,Tn,rr,snow,p,eee"
      valid_for["wmo_gross_errors"]     = "ta,Tx,Tn,td,p,mslp,dd,w"
      valid_for["wmo_time_consistency"] = "ta,td,p"
      valid_for["internal_consistency"] = "Tx,Tn,w,dd,sc,sd,fs"
      valid_for["temporal_coherence"]   = "Tx,Tn,w,sc,sd,fs"
      for (test in valid_for) {
        n = split(valid_for[test], arr, ",")
        for (i=1; i<=n; i++) is_valid[test SUBSEP arr[i]] = 1
      }
    }
    /^Vbl\t/ { vbl = $2; gsub(/\r$/, "", vbl) }
    FNR>1 && /^[0-9]/ && $7!="NA" {
      meta=$NF
      if (match(meta, /qc=([^[:space:]]+)/, m)) {
        qcstr = m[1]
        gsub(/\r$/, "", qcstr)
        n2 = split(qcstr, tests, ";")
        for (j=1; j<=n2; j++) {
          test = tests[j]
          if ((test in valid_for) && !((test SUBSEP vbl) in is_valid)) {
            bad[test]++
          }
        }
      }
    }
    END {
      for (test in bad) print fname, vbl, test, bad[test]
    }
  ' "$f"
done

# Finds files where climatic_outliers is flagged on a variable outside its intended scope (ta/rr/fs/p/Tx/Tn)
find /scratch3/PALAEO-RA/daily_data/final -type f -name "*.tsv" | while read f; do
  awk -F'\t' -v fname="$f" '
    /^Vbl\t/ { vbl = $2; gsub(/\r$/, "", vbl) }
    FNR>1 && /^[0-9]/ && $7!="NA" {
      meta=$NF
      if (match(meta, /qc=([^[:space:]]+)/, arr)) {
        qc=arr[1]
        gsub(/\r$/, "", qc)
        if (qc ~ /climatic_outliers/) {
          count++
        }
      }
    }
    END {
      if (count > 0) {
        valid = (vbl=="ta" || vbl=="Tx" || vbl=="Tn" || vbl=="rr" || vbl=="fs" || vbl=="p")
        if (!valid) {
          print fname, vbl, count+0
        }
      }
    }
  ' "$f"
done


# Extract Lat, Lon, Name from each file's SEF header and store alongside the filepath
find "$DATA_DIR" -type f -name "*.tsv" -print0 | while IFS= read -r -d '' filepath; do
    lat=$(awk -F'\t' '$1=="Lat"{print $2; exit}' "$filepath")
    lon=$(awk -F'\t' '$1=="Lon"{print $2; exit}' "$filepath")
    name=$(awk -F'\t' '$1=="Name"{print $2; exit}' "$filepath")

    # skip files missing coordinates (e.g. malformed headers)
    [[ -z "$lat" || -z "$lon" ]] && continue

    printf '%s\t%s\t%s\t%s\n' "$lat" "$lon" "$name" "$filepath"
done | awk -F'\t' '
{
    lat=$1; lon=$2; name=$3; file=$4
    if (NR==1 || lat>maxlat) { maxlat=lat; n_lat=lat; n_lon=lon; n_name=name; n_file=file }
    if (NR==1 || lat<minlat) { minlat=lat; s_lat=lat; s_lon=lon; s_name=name; s_file=file }
    if (NR==1 || lon>maxlon) { maxlon=lon; e_lat=lat; e_lon=lon; e_name=name; e_file=file }
    if (NR==1 || lon<minlon) { minlon=lon; w_lat=lat; w_lon=lon; w_name=name; w_file=file }
}
END {
    printf "Northernmost: lat=%s lon=%s  name=%s  file=%s\n", n_lat, n_lon, n_name, n_file
    printf "Southernmost: lat=%s lon=%s  name=%s  file=%s\n", s_lat, s_lon, s_name, s_file
    printf "Easternmost:  lat=%s lon=%s  name=%s  file=%s\n", e_lat, e_lon, e_name, e_file
    printf "Westernmost:  lat=%s lon=%s  name=%s  file=%s\n", w_lat, w_lon, w_name, w_file
}'


# find daily and subdaily files, and count how many of each there are, and how many non-conforming files there are. Non-conforming files are those that do not have "_daily_" or "_subdaily_" in their filename, or do not end with ".tsv".
DATA_DIR="/scratch3/PALAEO-RA/daily_data/final"

daily_count=0
subdaily_count=0
nonconforming_count=0

while IFS= read -r -d '' filepath; do
    fname=$(basename "$filepath")

    if [[ "$fname" == *_subdaily_* || "$fname" == *_subdaily.tsv ]]; then
        ((subdaily_count++))
    elif [[ "$fname" == *_daily_* || "$fname" == *_daily.tsv ]]; then
        ((daily_count++))
    else
        ((nonconforming_count++))
        echo "Non-conforming: $filepath"
    fi
done < <(find "$DATA_DIR" -type f -name "*.tsv" -print0)

echo ""
echo "Daily series:         $daily_count"
echo "Subdaily series:      $subdaily_count"
echo "Non-conforming files: $nonconforming_count"


# replace some lines
grep -rl  'qc= Tmax<=Tmin' . --include='*.tsv' | xargs sed -i 's/qc= Tmax<=Tmin/qc=Tn>Tx/g'
