#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./rename_chimes.sh /path/to/root
# If no path is given, it runs in the current directory.
ROOT="${1:-.}"

LOG="rename_chimes_$(date +%Y%m%d_%H%M%S).log"
touch "$LOG"

echo "Root: $ROOT" | tee -a "$LOG"
echo "Log:  $LOG"  | tee -a "$LOG"
echo "" | tee -a "$LOG"

# Find matching files (recursively), NUL-delimited for safety
# Pattern: CHIMES_ + [A-Z][A-Z] + [0-9][0-9] + _ + anything + .tsv
find "$ROOT" -type f -regextype posix-extended \
  -regex '.*/CHIMES_[A-Z]{2}[0-9]{2}_.+\.tsv' -print0 |
while IFS= read -r -d '' SRC; do
  DIR="$(dirname "$SRC")"
  BASE="$(basename "$SRC")"

  # Remove the AAXX_ part after CHIMES_
  # Example: CHIMES_AB12_filename.tsv -> CHIMES_filename.tsv
  DST_BASE="$(printf '%s' "$BASE" | sed -E 's/^CHIMES_[A-Z]{2}[0-9]{2}_/CHIMES_/')"
  DST="$DIR/$DST_BASE"

  if [[ "$SRC" == "$DST" ]]; then
    echo "SKIP (no change): $SRC" | tee -a "$LOG"
    continue
  fi

  if [[ -e "$DST" ]]; then
    echo "CONFLICT (exists, not renamed):" | tee -a "$LOG"
    echo "  src: $SRC" | tee -a "$LOG"
    echo "  dst: $DST" | tee -a "$LOG"
    continue
  fi

  mv -- "$SRC" "$DST"
  echo "RENAMED:" | tee -a "$LOG"
  echo "  $SRC" | tee -a "$LOG"
  echo "  -> $DST" | tee -a "$LOG"
done

echo "" | tee -a "$LOG"
echo "Done. Review: $LOG"
