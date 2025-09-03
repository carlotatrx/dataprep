rm(list=ls())
library(dataresqc)
library(readxl)
library(dplyr)
library(tidyr)
source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')

file <- '/scratch3/PALAEO-RA/daily_data/original/Europe_T2_RU_Kaliningrad_1702_subdaily.xls'

lat	<- 54.71666667
lon	<- 20.5
alt	<- 23

## NOTES:
## For the barometer Rheinländer Zoll is assumed, as it best fits the data, but
## there is no metadata supporting that.
## Thermometer and hygrometer were not formatted.

raw <- read_excel(file, sheet=1, skip=5)

# pressure conversion
conversions <- list(p  = function(x) round(convert_pressure(as.numeric(x),f = 26.154,lat = lat,alt = alt), 1))

## Wind direction
directions <- c("N", "NtO", "NNO", "NOtN", "NO", "NOtO", "ONO", "OtN", "O", "OtS", "OSO", "SOtO",
                "SO", "SOtS", "SSO", "StO", "S", "StW", "SSW", "SWtS", "SW", "SWtW", "WSW", "WtS",
                "W", "WtN", "WNW", "NWtW", "NW", "NWtN", "NNW", "NtW")

# Keep only the first 10 columns and rename
df <- raw %>%
  select(1:9) %>%
  rename(
    Year = 1, Month = 2, Day = 3, Hour.orig = 4, dd.orig = 5,
    p.zoll = 6, p.in = 7, ta.A = 8, ta.B = 9
  ) %>%
  # Fill down Year/Month/Day
  mutate(across(where(is.character), ~ na_if(trimws(.x), ""))) %>%
  fill(Year, Month, Day, .direction = "down") %>%
  mutate(
    Year = as.integer(Year),
    Month = as.integer(Month),
    Day = as.integer(Day),
    ta.A = as.numeric(ta.A),
    ta.B = suppressWarnings(as.numeric(ta.B)),
    Hour = case_when(
      Hour.orig=="M"~8,
      Hour.orig=="N"~12,
      Hour.orig=="A"~16,
      TRUE ~ NA_real_
      ),
    Minute = 0,
    Hour = suppressWarnings(as.integer(Hour)),
    p = as.numeric(p.zoll) + as.numeric(p.in)/12,
    p = conversions[['p']](p),
    dd = round(11.25 * (match(dd.orig, directions) -1), 0),
    ta = if_else(ta.A != 0, ta.A, ta.B),
    meta.time = paste0("orig.time=", Hour.orig)
  )
head(df)


# save

# ta
meta.col.ta <- paste0(df$meta.time," | orig=", df$ta)
write_sef_f(
  Data=as.data.frame(df[, c("Year", "Month", "Day", "Hour", "Minute", "ta")]),
  outfile="Kaliningrad_ta_subdaily.tsv",
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod="Kaliningrad",
  lat=lat, lon=lon, alt=alt,
  variable="ta",
  nam="Kaliningrad",
  meta=meta.col.ta,
  metaHead="Observer=Gottsched",
  units="unknown", sou="PALAEO-RA", stat="point",keep_na = F
)

# p
meta.col.p <- paste0(df$meta.time," | orig=",df$p.zoll,".",df$p.in,"Rh.in")
write_sef_f(
  Data=as.data.frame(df[, c("Year", "Month", "Day", "Hour", "Minute", "p")]),
  outfile="Kaliningrad_p_subdaily.tsv",
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod="Kaliningrad",
  lat=lat, lon=lon, alt=alt,
  variable="p",
  nam="Kaliningrad",
  meta=meta.col.p,
  metaHead="Observer=Gottsched | PGC=Y | PTC=N",
  units="hPa", sou="PALAEO-RA", stat="point",keep_na = F
)

# dd
meta.col.dd <- paste0(df$meta.time," | orig=",df$dd.orig)
write_sef_f(
  Data=as.data.frame(df[, c("Year", "Month", "Day", "Hour", "Minute", "dd")]),
  outfile="Kaliningrad_dd_subdaily.tsv",
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod="Kaliningrad",
  lat=lat, lon=lon, alt=alt,
  variable="p",
  nam="Kaliningrad",
  meta=meta.col.dd,
  metaHead="Observer=Gottsched",
  units="deg", sou="PALAEO-RA", stat="point",keep_na = F
)
