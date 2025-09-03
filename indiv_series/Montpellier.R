rm(list=ls())
library(dataresqc)
library(readxl)
library(dplyr)
library(tidyr)
source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')

file <- '/scratch3/PALAEO-RA/daily_data/original/Europe_T3_FR_Montpellier_1705-1748_subdaily.xls'
raw <- read_excel(file, sheet = 1, skip  = 6)

lat	<- 43.610769
lon	<- 3.876716
alt	<- NA
metaHead <- "Observer=Bon"

## NOTES:
## For the barometer Rheinländer Zoll is assumed, as it best fits the data, but
## there is no metadata supporting that.
## Thermometer and hygrometer were not formatted.


# pressure conversion
conversions <- list(p  = function(x) round(convert_pressure(as.numeric(x), f = 27.07, lat = lat, alt = alt), 1),
                    dd = function(x) 22.5 * (match(toupper(x), directions) - 1))
time.offset <- as.numeric(lon)*12/180

###### Wind direction ########
directions <- c("N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", 
                "SSO", "SO", "OSO", "O", "ONO", "NO", "NNO")

# Keep only the first 10 columns and rename
df <- raw %>%
  select(1:11) %>%
  rename(
    Year = 1, Month = 2, Day = 3, Hour.orig = 4, ta.p = 5, ta.l = 6, ta.R = 7,
    p.zoll = 8, p.in = 9, dd.orig = 10, notes = 11
  ) %>%
  # Fill down Year/Month/Day
  mutate(across(where(is.character), ~ na_if(trimws(.x), ""))) %>%
  mutate(
    Year = as.integer(Year),
    Month = as.integer(Month),
    Day = as.integer(Day),
    
    Hour.clean = case_when(
      grepl("(?i)^morning$", Hour.orig) ~ "08:00",
      grepl("(?i)^afternoon$", Hour.orig) ~ "15:00",
      TRUE ~ Hour.orig
    ),
    flag.hours = if_else(
      !is.na(Hour.orig) & is.na(suppressWarnings(as.numeric(Hour.orig))) &
        !grepl("(?i)^(morning|afternoon)$", Hour.orig),
      TRUE, FALSE
    ),
    
    
    Hour.decimal = suppressWarnings(as.numeric(Hour.orig)) * 24,
    Hour = floor(Hour.decimal),
    Minute = round((Hour.decimal-Hour)*60),
    
    dd = conversions[['dd']](dd.orig),
    p = suppressWarnings(as.numeric(p.zoll)) + suppressWarnings(as.numeric(p.in))/12,
    p = conversions[['p']](p),
    ta = round(suppressWarnings(as.numeric(ta.R)) * 1.25,2),
    
    meta.time = if_else(
      is.na(Hour),
      NA,
      paste0("orig.time=", sprintf("%02d", Hour),":", sprintf("%02d",Minute))
    )
  )
df %>% filter(flag.hours)
head(df)

# save

# ta
meta.col.ta <- ifelse(
  is.na(df$Hour),
  paste0("orig=", df$ta.R, "R"),
  paste0(df$meta.time," | orig=",df$ta.R,"R")
)

write_sef_f(
  Data=as.data.frame(df[, c("Year", "Month", "Day", "Hour", "Minute", "ta")]),
  outfile="Montpellier_ta_subdaily.tsv",
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod="Montpellier",
  lat=lat, lon=lon, alt=alt,
  variable="ta",
  nam="Montpellier",
  meta=meta.col.ta,
  metaHead=metaHead,
  time_offset = time.offset,
  units="C", sou="PALAEO-RA", stat="point",keep_na = F
)

# p
meta.col.p <- ifelse(
  is.na(df$Hour),
  paste0("orig=", df$p.zoll, ".", df$p.in, "P.in"),
  paste0(df$meta.time," | orig=",df$p.zoll,".",df$p.in,"P.in")
)

write_sef_f(
  Data=as.data.frame(df[, c("Year", "Month", "Day", "Hour", "Minute", "p")]),
  outfile="Montpellier_p_subdaily.tsv",
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod="Montpellier",
  lat=lat, lon=lon, alt=alt,
  variable="p",
  nam="Montpellier",
  meta=meta.col.p,
  metaHead=paste0(metaHead, " | PGC=Y | PTC=N"),
  time_offset=time.offset,
  units="hPa", sou="PALAEO-RA", stat="point",keep_na = F
)

# dd
meta.col.dd <- ifelse(
  is.na(df$Hour),
  paste0("orig=", df$dd.orig),
  paste0(df$meta.time," | orig=",df$dd.orig)
)
  
write_sef_f(
  Data=as.data.frame(df[, c("Year", "Month", "Day", "Hour", "Minute", "dd")]),
  outfile="Montpellier_dd_subdaily.tsv",
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod="Montpellier",
  lat=lat, lon=lon, alt=alt,
  variable="dd",
  nam="Montpellier",
  meta=meta.col.dd,
  metaHead=metaHead,
  time_offset=time.offset,
  units="deg", sou="PALAEO-RA", stat="point",keep_na = F
)
