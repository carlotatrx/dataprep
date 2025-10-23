#!/bin/bash
# ------------------------------------------------------------
# Fix missing Hour and Minute in SEF files.
# Replaces:
#   Hour (NA)   → 24
#   Minute (NA) → 0
# Adds IMPROVE Link to metadata
#
# Usage: bash unix_commands.sh input.sef output.sef
# ------------------------------------------------------------


for f in *.tsv; do
  f2="${f%.tsv}_new.tsv"
  awk -F'\t' '                 # read file line-by-line and split by tabs
  BEGIN { OFS = "\t" }         # Output Field Separator is tab
  {
    if ($1=="Source" && $2=="GCOS") sou=1
    if ($1=="Link" && sou==1) $2="https://doi.pangaea.de/10.1594/PANGAEA.948258"
    print
    next
    # add link to metadata
  }     
# {
#   if ($4 == "NA") $4 = 24;   # Hour column
#   if ($5 == "NA") $5 = 0;    # Minute column
#   print
# }
  ' "$f" > "$f2"
done