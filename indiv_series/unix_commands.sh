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


awk -F'\t' '                 # read file line-by-line and split by tabs
BEGIN { OFS = "\t" }         # Output Field Separator is tab
NR <= 13 {
  if ($1=="Source" && $2=="IMPROVE") sou=1
  if ($1=="Link" && sou==1) $2="https://link.springer.com/book/10.1007/978-94-010-0371-1"
  print; next
  # add link to metadata
}     
{
  if ($4 == "NA") $4 = 24;   # Hour column
  if ($5 == "NA") $5 = 0;    # Minute column
  print
}
' "$1" > "$2"