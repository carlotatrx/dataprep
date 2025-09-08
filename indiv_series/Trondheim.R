rm(list=ls())
library(dataresqc)
library(readxl)
library(dplyr)
library(tidyr)
source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')

file <- '/scratch3/PALAEO-RA/daily_data/original/Europe_T3_NO_Trondheim_1762-1764_subdaily.xls'
raw <- read_excel(file, sheet = 1, skip  = 6)

lat	<- 63.4305149
lon	<-10.3950528
alt	<- NA
name <- "Trondheim"


# pressure conversion
conversions <- list(p  = function(x, ta) round(convert_pressure(as.numeric(x), f = 27.07, lat = lat, alt = alt, atb=ta), 1),
                    dd = function(x) 22.5 * (match(toupper(x), directions) - 1))
time.offset <- as.numeric(lon)*12/180

###### Wind direction ########
directions <- c("N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", 
                "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW")

# Keep only the first 10 columns and rename
df <- raw %>%
  select(1:7) %>%
  rename(
    Year = 1, Month = 2, Day = 3, p.zoll = 4, p.in = 5, ta.R = 6, dd.orig = 7
  ) %>%
  # Fill down Year/Month/Day
  mutate(across(where(is.character), ~ na_if(trimws(.x), ""))) %>%
  fill(p.zoll, .direction = "down") %>%
  mutate(
    Year = as.integer(Year),
    Month = as.integer(Month),
    Day = as.integer(Day),
    Hour = NA,
    Minute = NA,
    
    dd = conversions[['dd']](dd.orig),
    p = suppressWarnings(as.numeric(p.zoll)) + suppressWarnings(as.numeric(p.in))/12,
    ta = round(suppressWarnings(as.numeric(ta.R)) * 1.25,2),
    p = conversions[['p']](p, ta),
    meta.p = paste0("orig.p=", p.zoll,".",p.in,"Pin"),
    meta.ta = paste0("orig.ta=",ta.R,"R"),
    meta.dd = paste0("orig.dd=",dd.orig)
  )
head(df)

# save

# ta
write_sef_f(
  Data=as.data.frame(df[, c("Year", "Month", "Day", "Hour", "Minute", "ta")]),
  outfile=paste0(name,"_ta_subdaily.tsv"),
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod=name,
  lat=lat, lon=lon, alt=alt,
  variable="ta",
  nam=name,
  meta=df$meta.ta,
  units="C", sou="PALAEO-RA", stat="point",keep_na = F
)

# p
write_sef_f(
  Data=as.data.frame(df[, c("Year", "Month", "Day", "Hour", "Minute", "p")]),
  outfile=paste0(name,"_p_subdaily.tsv"),
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod=name,
  lat=lat, lon=lon, alt=alt,
  variable="p",
  nam=name,
  meta=df$meta.p,
  metaHead="PGC=Y | PTC=Y",
  units="hPa", sou="PALAEO-RA", stat="point",keep_na = F
)

# dd
write_sef_f(
  Data=as.data.frame(df[, c("Year", "Month", "Day", "Hour", "Minute", "dd")]),
  outfile=paste0(name,"_dd_subdaily.tsv"),
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod=name,
  lat=lat, lon=lon, alt=alt,
  variable="dd",
  nam=name,
  meta=df$meta.dd,
  units="deg", sou="PALAEO-RA", stat="point",keep_na = F
)
