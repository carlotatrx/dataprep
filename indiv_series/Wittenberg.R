## Wittenberg.R

rm(list=ls())
library(dataresqc)
library(XLConnect)
library(lubridate)

source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')

name <- "Wittenberg"
lat <-	51.86666667
lon	<- 12.8
alt	<- 72
code <- "Europe_Wittenberg_1"
metahead <- "Observer=Weidler"

## Define function to fill missing values with the value on the previous row
fill_variable <- function(timeseries, strings) {
  for (i in 2:length(timeseries)) {
    if (timeseries[i] %in% strings) {
      if(!timeseries[i-1] %in% strings) {
        timeseries[i] <- timeseries[i-1]
      }
    }
  }
  return(timeseries)
}


file <- "/scratch3/PALAEO-RA/daily_data/original/Europe_T3_DE_Wittenberg_1728-1729_subdaily.xls"

## Read data
raw <- readWorksheetFromFile(file, sheet=1, startRow=7, header=T, colTypes="character")
raw[,1] <- fill_variable(raw[,1], NA)
raw[,2] <- fill_variable(raw[,2], NA)
raw[,3] <- fill_variable(raw[,3], NA)

times <- strptime(raw[,4], format="%Y-%m-%d %H:%M")

meta_time <- case_when(
  is.na(times) ~ "orig.time=NA",
  nchar(raw[,4]) < 10 ~ paste0("orig.time", raw[,4]),
  TRUE ~ paste0("orig.time=", substr(raw[,4],12,13))
)

hours <- hour(times)
hours[which(is.na(hours))] <- 6
minutes <- minute(times)
minutes[which(is.na(minutes))] <- 0

out <- data.frame(Years = raw[,1],
                  Month = raw[,2],
                  Day   = raw[,3],
                  Hour  = hours,
                  Minute = minutes,
                  Value = raw[,7],
                  meta  = paste0(meta_time, " | orig.ta=", raw[,7]),
                  stringsAsFactors = F)

## Write SEF files
write_sef_f(out,
            outfile=paste0(name, "_ta_subdaily.tsv"),
            outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
            variable = "ta",
            cod = code,
            nam = name,
            lat = lat,
            lon = lon,
            alt = alt,
            sou = "PALAEO-RA",
            units = "unknown",
            stat = "point",
            meta = out$meta,
            metaHead = metahead,
            period = 0,
            keep_na = F)