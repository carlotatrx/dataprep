# Usage: bash scale_rr.sh input.tsv output.tsv
# Multiplies the precipitation Value column (7th column) by 10 in SEF files.
# to use in some KNMI files where I didn't convert properly.

# Read input and output file names
infile="$1"
outfile="$2"

awk '
BEGIN {
    FS=OFS="\t"   # Use tab as field separator for input and output
}

# Copy all lines until the table header (Year Month ...)
# These lines form the SEF metadata block
/^Year[[:space:]]+Month/ {header=1}
header == 0 {print; next}

# After the header: process data rows
header == 1 {
    if ($7 == "Value") { print; next }  # Print header line as is
    if (NF < 7) {print; next}   # Skip malformed or empty rows
    $7 = ($7 * 10)              # Multiply the Value column by 10
    print                       # Output modified row
}
' "$infile" > "$outfile"

# example usage:
# bash /scratch2/ccorbella/code/dataprep/scale_rr_KNMI.sh /scratch3/PALAEO-RA/daily_data/final/Bergen/Bergen_17410904-17460518_rr_subdaily_qc.tsv /scratch3/PALAEO-RA/daily_data/final/Bergen/Bergen_17410904-17460518_rr_subdaily.tsv