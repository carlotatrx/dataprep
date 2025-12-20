# usage: bash Gurzelen.sh
awk -F'\t' 'BEGIN{OFS="\t"}
NR==1 {print; next}
{
  if ($8 ~ /orig.time=Abend/)  {$4=19; $5=0}
  if ($8 ~ /orig.time=Morgen/) {$4=8;  $5=0}
  if ($8 ~ /orig.time=Mittag/) {$4=12; $5=0}
  print
}' /scratch3/PALAEO-RA/daily_data/final/Gurzelen/Gurzelen_17710101-17841231_dd_subdaily.tsv > /scratch3/PALAEO-RA/daily_data/final/Gurzelen/Gurzelen_17710101-17841231_dd_subdaily_qc2.tsv
