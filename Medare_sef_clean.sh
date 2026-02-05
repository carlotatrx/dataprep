# output: ID<TAB>Name
: > output.txt
for f in Medare_*_rr_daily.tsv; do
  id=$(awk -F'\t' '$1=="ID"{print $2; exit}' "$f")
  name=$(awk -F'\t' '$1=="Name"{print $2; exit}' "$f")
  printf "%s\t%s\n" "$id" "$name" >> output.txt
done

# --------------------------
while IFS=$'\t' read -r id name; do
  clean_name=$(echo "$name" | tr ' _' '--')
  for f in "Medare_${id}_"*"_rr_daily.tsv"; do
    [ -e "$f" ] || continue
    new="${f/Medare_${id}_/Medare_${clean_name}_}"
    mv -n -- "$f" "$new"
  done
done

