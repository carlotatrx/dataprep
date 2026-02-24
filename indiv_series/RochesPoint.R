rm(list=ls())
library(dataresqc)
library(readxl)
library(tidyr)
library(dplyr)

source('/scratch2/ccorbella/code/dataprep/helpfun.R')

indir <- "/scratch3/PALAEO-RA/daily_data/original/"
name <- "RochesPoint"

outdir <- paste0("/scratch3/PALAEO-RA/daily_data/final/",name)

source <- "Mateus, C., Potito, A., & Curley, M. (2020). Reconstruction of a long‐term historical daily maximum and minimum air temperature network dataset for Ireland (1831‐1968). Geoscience Data Journal, 7(2), 102-115."
link   <- "https://www.edepositireland.ie/entities/publication/7271837d-dbe6-463f-b103-422f42a15446"
lat <- 51.795
lon <- -8.252
alt <- 13.1

metaHead <- paste0(
  "Observer=William Kennedy (1873-1902), Station staff (1902-1920) | Instrument=Stevenson screen (verified 1872-1920),",
  "Max/Min thermometer"
)

raw <- read.csv(paste0(indir, name, "/Roches Point_1872-1920.csv"), skip=1, header=FALSE,
                col.names=c("year", "month", "day", "maxF", "minF", "maxC", "minC", "notes1", "notes2"))
head(raw)
tail(raw)

df <- raw %>%
  mutate(
    hour=NA,
    minute=NA,
    meta.Tx = case_when(
      # 1. 1872 to June 1908: 8:00 a.m.
      (year >= 1872 & year < 1908) | (year == 1908 & month <= 6) ~ 
        paste0("orig=", maxF, "F | obs.time=8a.m."),
      
      # 2. July 1908 to February 1912: 7:00 a.m.
      (year == 1908 & month >= 7) | (year > 1908 & year < 1912) | (year == 1912 & month <= 2) ~ 
        paste0("orig=", maxF, "F | obs.time=7a.m."),
      
      # 3. March 1912 to 1920: 7:00 a.m. and 6:00 p.m.
      (year == 1912 & month >= 3) | (year > 1912 & year <= 1920) ~ 
        paste0("orig=", maxF, "F | obs.time=7a.m./6p.m."),
      
      # Catch-all for missing data
      is.na(maxF) ~ "obs.time=9a.m.",
      
      TRUE ~ NA_character_
    ),
    meta.Tn = case_when(
      (year >= 1872 & year < 1908) | (year == 1908 & month <= 6) ~ 
        paste0("orig=", minF, "F | obs.time=8a.m."),
      (year == 1908 & month >= 7) | (year > 1908 & year < 1912) | (year == 1912 & month <= 2) ~ 
        paste0("orig=", minF, "F | obs.time=7a.m."),
      (year == 1912 & month >= 3) | (year > 1912 & year <= 1920) ~ 
        paste0("orig=", minF, "F | obs.time=7a.m./6p.m."),
      is.na(minF) ~ "obs.time=9a.m.",
      TRUE ~ NA_character_
    ),
  ) %>% select(year, month, day, hour, minute, Tx=maxC, Tn=minC, meta.Tx, meta.Tn)

head(df)
tail(df)

var <- "Tn"
write_sef_f(
  as.data.frame(df[c("year","month", "day", "hour", "minute","Tn")]),
  outfile = outfile.name(paste0("ILMMT_",name), var, df, FALSE),
  outpath = outdir,
  cod     = paste0("ILMMT-",name),
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "minimum",
  period    = "day",
  units   = units(var),
  meta = df$meta.Tn,
  metaHead=metaHead,
  keep_na=FALSE
)


var <- "Tx"
write_sef_f(
  as.data.frame(df[c("year","month", "day", "hour", "minute",var)]),
  outfile = outfile.name(paste0("ILMMT_",name), var, df, FALSE),
  outpath = outdir,
  cod     = paste0("ILMMT-",name),
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "maximum",
  period    = "day",
  units   = units(var),
  metaHead=metaHead,
  meta=df$meta.Tx,
  keep_na=FALSE
)

