## Only temperature formatted

library(dataresqc)
library(dplyr)
library(stringr)

source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')
sef_test_path <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/sef_tests/'
outpath_preprocessed <- "/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed"


# Delemont ----------------------------------------------------------------

ta <- read.delim('/scratch3/PALAEO-RA/Instrumental/SEF/Daily/QCed/Chimes/TA/CHIMES_JU01_Delemont_18011222-18321230_ta.tsv',
                 header=T, sep="\t", skip=12)
meta <- read_meta('/scratch3/PALAEO-RA/Instrumental/SEF/Daily/QCed/Chimes/TA/CHIMES_JU01_Delemont_18011222-18321230_ta.tsv')

ta <- ta %>%
  mutate(
    Hour = case_when(
      str_detect(Meta, "matin") ~ 8,
      str_detect(Meta, "apresmidi") ~ 14,
      str_detect(Meta, "soir") ~ 20,
      TRUE ~ NA_real_
    ),
    
    Meta = str_replace_all(Meta, "\\|", " | "),
    Meta = str_replace_all(Meta, "orig=", "orig_ta="),
    Meta = str_replace_all(Meta, "orig.time", "orig_time")
  )

ta$Minute <- 0
ta$Period <- NULL
write_sef_f(Data=ta, outfile="Delemont_ta_subdaily.tsv",
            outpath=outpath_preprocessed,
            cod=meta[["id"]],
            metaHead=meta[["meta"]],
            variable="ta",
            nam=meta[["name"]], link=meta[['link']],
            lat=meta[["lat"]],
            lon=meta[["lon"]], alt=meta[["alt"]], sou=meta[["source"]],
            units="C", stat="point",keep_na = F
)

p <- read.delim('/scratch3/PALAEO-RA/Instrumental/SEF/Daily/QCed/Chimes/P/CHIMES_JU01_Delemont_18011222-18321230_p.tsv',
                 header=T, sep="\t", skip=12)
meta <- read_meta('/scratch3/PALAEO-RA/Instrumental/SEF/Daily/QCed/Chimes/P/CHIMES_JU01_Delemont_18011222-18321230_p.tsv')

p <- p %>%
  mutate(
    Hour = case_when(
      str_detect(Meta, "matin") ~ 8,
      str_detect(Meta, "apresmidi") ~ 14,
      str_detect(Meta, "soir") ~ 20,
      TRUE ~ NA_real_
    ),
    
    Meta = str_replace_all(Meta, "\\|", " | "),
    Meta = str_replace_all(Meta, "orig=", "orig_p="),
    Meta = str_replace_all(Meta, "orig.time", "orig_time")
  )

p$Minute <- 0
p$Period <- NULL

write_sef_f(Data=ta, outfile="Delemont_p_subdaily.tsv",
            outpath=outpath_preprocessed,
            cod=meta[["id"]],
            metaHead="Observer=François-Joseph Helg | PTC=N | PGC=Y",
            variable="p",
            nam=meta[["name"]], link=meta[['link']],
            lat=meta[["lat"]],
            lon=meta[["lon"]], alt=meta[["alt"]], sou=meta[["source"]],
            units="hPa", stat="point",keep_na = F
)


# Delsberg ----------------------------------------------------------------

lat <- 47.366
lon <- 7.35
ele <- 435

## Read data
template <- read.table("/scratch3/PALAEO-RA/DataRescue/Incoming/DigiHom/Delsberg_Suppl_DC/Delsberg_1802_1832_Suppl_dc.txt", header=TRUE, 
                       sep="\t", stringsAsFactors=FALSE)

## Initialize data frame
out <- data.frame(year = template[,1], 
                  month = template[,2], 
                  day = template[,3], 
                  hour = 24, 
                  minute = 0)

## Temperature
out$ta <- as.numeric(template$Temperatur)
out$p <- as.numeric(template$Barometer)
out$meta.p <- paste0("orig_p=",out$p,"mmHg | atb=",out$ta,"C")
out$p.corrected <- round(convert_pressure(out$p, lat=lat, alt=ele, atb=out$ta))

## Write SEF files
write_sef_f(
  Data=out[,c("year","month","day","hour","minute","ta")], outfile="Delsberg_ta_daily.tsv",
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod="Delsberg",
  variable="ta",
  nam="Delsberg",
  lat=lat, lon=lon, alt=ele, period="day",
  units="C", sou="DigiHom", stat="mean",keep_na = F
)

write_sef_f(
  Data=out[,c("year","month","day","hour","minute","p.corrected")], outfile="Delsberg_p_daily.tsv",
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod="Delsberg",
  variable="p",
  nam="Delsberg",
  lat=lat, lon=lon, alt=ele, period='day',
  meta=out$meta.p, metaHead = "PGC=Y | PTC=N",
  units="hPa", sou="DigiHom", stat="mean",keep_na = F
)
