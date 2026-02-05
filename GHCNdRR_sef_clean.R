rm(list=ls())
library(openxlsx)

# directory with SEF files
indir <- "/scratch3/PALAEO-RA/daily_data/tmp/GHCNdRR"

files <- list.files(indir, pattern = "^GHCN-D_.*\\.tsv$", full.names = TRUE)

get_meta <- function(f) {
  h <- readLines(f, n = 20)
  
  data.frame(
    ID   = sub("^ID\\t",   "", h[grep("^ID\\t", h)]),
    Name = sub("^Name\\t", "", h[grep("^Name\\t", h)]),
    Lat  = as.numeric(sub("^Lat\\t", "", h[grep("^Lat\\t", h)])),
    Lon  = as.numeric(sub("^Lon\\t", "", h[grep("^Lon\\t", h)])),
    stringsAsFactors = FALSE
  )
}

out <- do.call(rbind, lapply(files, get_meta))

write.xlsx(out, file = "/scratch3/PALAEO-RA/daily_data/tmp/GHCNdRR/GHCN_SEF_station_metadata_Europe.xlsx", overwrite = TRUE)
write.table(
  out,
  file = "/scratch3/PALAEO-RA/daily_data/tmp/GHCNdRR/GHCN_SEF_station_metadata.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# remove files not in the domain -----------------------------------------

awk -F'\t' '
NR==1{next}  # skip header
($3 < 22) || ($4 < -51) || ($4 > 81) {
  printf("rm -v -- *%s*\n", $1)
}
' GHCN_SEF_station_metadata.tsv | bash


# rename the files -------------------------------------------------------
for f in GHCN-D_*_rr_daily.tsv; do
[ -f "$f" ] || continue

id=$(awk -F'\t' '$1=="ID"{print $2; exit}' "$f")
name=$(awk -F'\t' '$1=="Name"{print $2; exit}' "$f")

# skip if not a SEF file / missing header fields
[ -n "$id" ] || { echo "SKIP (no ID): $f"; continue; }
[ -n "$name" ] || { echo "SKIP (no Name): $f"; continue; }

station=$(printf "%s" "$name" |
            sed 's/[()]/ /g; s/_/ /g; s/[^[:alnum:] ]/ /g; s/[[:space:]]\+/ /g; s/^ //; s/ $//' |
            perl -pe 's/\b([[:alnum:]])([[:alnum:]]*)/\U$1\L$2/g' |
            sed 's/ /-/g; s/-\+/-/g; s/^-//; s/-$//'
)

daterest=$(echo "$f" | sed -E 's/^GHCN-D_[^_]+_//')   # everything after ID_
daterange=$(echo "$daterest" | cut -d_ -f1)          # YYYYMMDD-YYYYMMDD
suffix=$(echo "$daterest" | cut -d_ -f2-)            # rr_daily.tsv

new="GHCN-D_${id}_${station}_${daterange}_${suffix}"

if [ "$f" = "$new" ]; then
  echo "OK (already): $f"
else
  mv -n -- "$f" "$new"
fi
done


# move GHCN-D_RR files to their corresponding station ----------
#!/bin/bash
set -euo pipefail

TMP="/scratch3/PALAEO-RA/daily_data/tmp/GHCNdRR"
ORIGTMP="/scratch3/PALAEO-RA/daily_data/tmp/GHCNdRR/original_GHCNdRR"

FINAL_BASE="/scratch3/PALAEO-RA/daily_data/final"
ORIGINAL_BASE="/scratch3/PALAEO-RA/daily_data/original"

while IFS= read -r station; do
  [ -n "$station" ] || continue

  final_dir="$FINAL_BASE/$station"
  orig_dir="$ORIGINAL_BASE/$station"
  
  if [ ! -d "$final_dir" ]; then
    echo "SKIP (missing final dir): $final_dir"
    continue
  fi
  if [ ! -d "$orig_dir" ]; then
    echo "SKIP (missing original dir): $orig_dir"
    continue
  fi

  # 1) warn if final already has rr
  if ls "$final_dir"/*_rr* >/dev/null 2>&1; then
    echo "WARNING: $station already has *_rr* in $final_dir"
  fi

  # 2) find renamed rr_daily files for this station
  shopt -s nullglob
  files=( "$TMP"/GHCN-D_*_"$station"_*_rr_daily.tsv )
  shopt -u nullglob

  if [ ${#files[@]} -eq 0 ]; then
    echo "WARNING: no rr_daily file found for $station in $TMP"
    continue
  fi
  
  for f in "${files[@]}"; do
    base=$(basename "$f")
    code=$(echo "$base" | cut -d_ -f2)   # GHCN-D_<CODE>_<Station>_...
    
    mv -v -- "$f" "$final_dir/"
    
    # 4) find original file(s) by code
    shopt -s nullglob
    orig_matches=( "$ORIGTMP"/GHCN-D_"$code"_*_rr.tsv )
    shopt -u nullglob
    
    if [ ${#orig_matches[@]} -eq 0 ]; then
      echo "WARNING: no original file found for code=$code (station=$station) in $ORIGTMP"
      continue
    fi
    
    for of in "${orig_matches[@]}"; do
      mv -v -- "$of" "$orig_dir/"
    done
  done
      
done < already_files.txt


# move seeing if filename coincides --------------------------------
#!/bin/bash
set -euo pipefail

TMP="/scratch3/PALAEO-RA/daily_data/tmp/GHCNdRR"
ORIGTMP="/scratch3/PALAEO-RA/daily_data/tmp/GHCNdRR/original_GHCNdRR"

FINAL_BASE="/scratch3/PALAEO-RA/daily_data/final"
ORIGINAL_BASE="/scratch3/PALAEO-RA/daily_data/original"

resolve_target_station() {
  local station="$1"
  local cand first
  
  # A) remove hyphens (Port-Said -> PortSaid)
  cand="${station//-/}"
  if [ -d "$FINAL_BASE/$cand" ]; then
    echo "$cand"; return 0
  fi
  
  # B) fallback to first token ONLY if first token is NOT in the exception list
  first="${station%%-*}"
  case "$first" in
    Port|Isle|Island|Pte)
      return 1 ;;   # do NOT try first-token fallback
    *)
    cand="$first"
    if [ -d "$FINAL_BASE/$cand" ]; then
      echo "$cand"; return 0
    fi
    ;;
    esac
    
    return 1
  }

shopt -s nullglob
for f in "$TMP"/GHCN-D_*_*_rr_daily.tsv; do
  base=$(basename "$f")
  station=$(echo "$base" | cut -d_ -f3)
  code=$(echo "$base" | cut -d_ -f2)
  
  if target=$(resolve_target_station "$station"); then
    final_dir="$FINAL_BASE/$target"
    orig_dir="$ORIGINAL_BASE/$target"
  
    if [ ! -d "$orig_dir" ]; then
      echo "SKIP (missing original dir): $orig_dir  | file=$base"
      continue
    fi
    
    # 1) warn if final already has rr
    if ls "$final_dir"/*_rr* >/dev/null 2>&1; then
      echo "WARNING: $target already has *_rr* in $final_dir"
    fi
    
    # 2) move final file
    mv -v -- "$f" "$final_dir/"
    
    # 3) move original(s) by code
    orig_matches=( "$ORIGTMP"/GHCN-D_"$code"_*_rr.tsv )
    if [ ${#orig_matches[@]} -eq 0 ]; then
      echo "WARNING: no original file found for code=$code (from $base) in $ORIGTMP"
      continue
    fi
      
    for of in "${orig_matches[@]}"; do
      mv -v -- "$of" "$orig_dir/"
    done
  else
    echo "NO MATCH: station=$station  | file=$base"
  fi
done
shopt -u nullglob


# copy making folders -----------------------------------------------------


#!/bin/bash
set -euo pipefail

TMP="/scratch3/PALAEO-RA/daily_data/tmp/GHCNdRR"
ORIGTMP="/scratch3/PALAEO-RA/daily_data/tmp/GHCNdRR/original_GHCNdRR"

FINAL_BASE="/scratch3/PALAEO-RA/daily_data/final"
ORIGINAL_BASE="/scratch3/PALAEO-RA/daily_data/original"

STATIONS_OUT="all_stations.txt"
: > "$STATIONS_OUT"

shopt -s nullglob
for f in "$TMP"/GHCN-D_*_*_rr_daily.tsv; do
base=$(basename "$f")
code=$(echo "$base" | cut -d_ -f2)
station=$(echo "$base" | cut -d_ -f3)

# save station name
echo "$station" >> "$STATIONS_OUT"

final_dir="$FINAL_BASE/$station"
orig_dir="$ORIGINAL_BASE/$station"

mkdir -p -- "$final_dir"
mkdir -p -- "$orig_dir"

if [ -d "$final_dir" ] && ls "$final_dir"/*_rr* >/dev/null 2>&1; then
echo "WARNING: $station already has *_rr* in $final_dir"
fi

mv -v -- "$f" "$final_dir/"
orig_matches=( "$ORIGTMP"/GHCN-D_"$code"_*_rr.tsv )
if [ ${#orig_matches[@]} -eq 0 ]; then
  echo "WARNING: no original file found for code=$code (from $base) in $ORIGTMP"
  continue
  fi
  
  for of in "${orig_matches[@]}"; do
  mv -v -- "$of" "$orig_dir/"
  done
  done
  shopt -u nullglob
  
  # unique + sorted station list
  sort -u "$STATIONS_OUT" -o "$STATIONS_OUT"
  
  
    