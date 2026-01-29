ELE=46 # elevation for the files
SITE=Nottingham # to rename the files

for f in *Tn*.tsv *Tx*.tsv *rr*.tsv; do 

  [ -f "$f" ] || continue 

  awk -F'\t' 'BEGIN{OFS="\t"; indata=0} 

  /^Year\tMonth\tDay\tHour\tMinute\tPeriod\tValue\tMeta/ {indata=1; print; next} 

  !indata {print; next} 

  { 

    $4="NA"; $5="NA" 

    print 

  }' "$f" > "$f.tmp" && mv "$f.tmp" "$f" 

done
 

# Add the source to all the files in the current dir. 

for f in *.tsv; do 

  [ -f "$f" ] || continue 

  awk -F'\t' 'BEGIN{OFS="\t"; done=0} 

    /^Link(\t|$)/ && !done { 

      done=1 

      if (NF==1 || $2=="" ) $2="DOI 10.5281/zenodo.5940391" 

      print 

      next 

    } 

    {print} 

  ' "$f" > "${f}.tmp" && mv "${f}.tmp" "$f" 

done 

# Add elevation to all the files in cd 

for f in *.tsv; do 

  [ -f "$f" ] || continue 

  awk -F'\t' -v ele="$ELE" 'BEGIN{OFS="\t"; done=0}

    /^Alt(\t|$)/ && !done { 

      done=1 

      if (NF==1 || $2=="" ) $2=ele 

      print 

      next 

    } 

    {print} 

  ' "$f" > "${f}.tmp" && mv "${f}.tmp" "$f" 

done 

 

 

# make that all the Values with more than 2 decimal places only have two decimal places. If it has one decimal or is an integer, don't change it 

for f in *.tsv; do 

  awk -F'\t' 'BEGIN{OFS="\t"; indata=0} 

  /^Year\tMonth\tDay\tHour\tMinute\tPeriod\tValue\tMeta/ {indata=1; print; next} 

  !indata {print; next} 

  { 

    if ($7 ~ /^-?[0-9]+(\.[0-9]{3,})([eE][+-]?[0-9]+)?$/) $7=sprintf("%.2f",$7) 

    print 

  }' "$f" > "$f.tmp" && mv "$f.tmp" "$f" 

done 

 

rename "s/DWR_UKMO_DWRUK_${SITE^^}/WRDv2.2_${SITE}/" *

rename 's/mslp/p/' * 

rename 's/tb/ta_barometer/' WRD*.tsv 

rename 's/.tsv/_subdaily.tsv/' WRD*.tsv

sed -i 's/p1day/day/g' * 

sed -i 's/mslp/p/' *_p*.tsv 

sed -i 's/|/ | /g' WRD*.tsv 

sed -i 's/tb/ta/' *_ta_bar*.tsv 

rename 's/rr_subdaily/rr_daily/' WRD*.tsv 

rename 's/Tn_subdaily/Tn_daily/' WRD*.tsv 

rename 's/Tx_subdaily/Tx_daily/' WRD*.tsv 