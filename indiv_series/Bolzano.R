rm(list=ls())
library(dataresqc)
library(readxl)
library(tidyr)
library(dplyr)
library(lubridate)
library(tibble)
library(stringr)
source('/scratch2/ccorbella/code/dataprep/helpfun.R')

outdir <- "/scratch3/PALAEO-RA/daily_data/final/Bolzano"

# -------------------------------------------------------------------------
lat <- 46.493
lon <- 11.341
alt <- 310
name <- "PALAEO-RA_Bolzano"
code <- "Bolzano"
source <- "PALAEO-RA"
link <- ""
metaHead <- "coords=approx."
raw <- read.csv("/scratch3/PALAEO-RA/daily_data/original/Bolzano/Bolzano_1842-1849_raw.csv", na.strings = c("NA","?"))

df <- raw %>%
  mutate(
    dt = ymd_hms(datetime),
    Year   = year(dt),
    Month  = month(dt),
    Day    = day(dt),
    Hour   = hour(dt),
    Minute = minute(dt),
    
    ta = round(Thermometer_Reaumur*1.25,2),
    p = round(convert_pressure((Barometer_Zoll+Barometer_Linien/12), f=27.07, lat=lat, alt=alt, atb=ta),2),
    meta.time = meta_time(Hour, Minute),
    meta.p = paste0(meta.time, " | orig=", Barometer_Zoll, "in", Barometer_Linien, "in | atb=", ta, "C"),
    meta.ta = paste0(meta.time, " | orig=", Thermometer_Reaumur, "R")
  ) %>% select(Year, Month, Day, Hour, Minute, p, ta, meta.p, meta.ta)

head(df)


# save
var <-"ta"
write_sef_f(
  as.data.frame(df[,c("Year","Month", "Day", "Hour", "Minute", "ta"),]),
  outfile = outfile.name(name, var, df, TRUE),
  outpath = outdir,
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "point",
  units   = units(var),
  metaHead = metaHead,
  meta    = df$meta.ta,
  keep_na = TRUE
)

# save
var <-"p"
write_sef_f(
  as.data.frame(df[,c("Year","Month", "Day", "Hour", "Minute", "p"),]),
  outfile = outfile.name(name, var, df, TRUE),
  outpath = outdir,
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "point",
  units   = units(var),
  metaHead = metaHead,
  meta    = df$meta.p,
  keep_na = TRUE
)


