## procedure to make GCOS files from Yuri compliant
## with the format I'm using

# add link if it's NA
for f in *.tsv; do
  f2="${f%.tsv}_subdaily.tsv"
  awk -F'\t' '                 # read file line-by-line and split by tabs
  BEGIN { OFS = "\t" }         # Output Field Separator is tab
  {
    if ($1=="Source" && $2=="GCOS") sou=1
    if ($1=="Link" && sou==1) $2="https://doi.pangaea.de/10.1594/PANGAEA.948258"
    print
    next
    # add link to metadata
  }' "$f" > "$f2"    
done

# clean ugly NAs
sed -i 's/NANA/NA/g' *
sed -i 's/orig.time=NA:NA/orig.time=NA/g' *
sed -i 's/NAR/NA/g' *
sed -i 's/NAC/NA/g' *
sed -i 's/NA.NAPin/NA/g' *
sed -i 's/NAmmHg/NA/g' *

sed -i 's/|/ | /g' *           # add spaces like I like
sed -i 's/  |  / | /g' *       # clean metaHead double spaces
sed -i 's/minimum/min/g' *     # shorten minimum to min
sed -i 's/maximum/max/g' *     # shorten maximum to max

mkdir old

# move all files that are not subdaily to a subdir
find . -maxdepth 1 -type f ! -name "*_subdaily.tsv" -exec mv {} subdir/ \;

# rename subdaily files to daily for Tn, Tx, rr
rename 's/_subdaily\.tsv/_daily.tsv/' *{Tn,Tx,rr}_subdaily.tsv

# other useful commands
rename 's/\.tsv$/_subdaily.tsv/' *