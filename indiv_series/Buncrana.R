rm(list=ls())
library(dataresqc)
library(readxl)
library(tidyr)
library(dplyr)

source('/scratch2/ccorbella/code/dataprep/helpfun.R')

indir <- "/scratch3/PALAEO-RA/daily_data/original/"
name <- "Buncrana"

outdir <- paste0("/scratch3/PALAEO-RA/daily_data/final/",name)

source <- "Mateus, C., Potito, A., & Curley, M. (2020). Reconstruction of a long‐term historical daily maximum and minimum air temperature network dataset for Ireland (1831‐1968). Geoscience Data Journal, 7(2), 102-115."
link   <- "https://www.edepositireland.ie/entities/publication/7271837d-dbe6-463f-b103-422f42a15446"
lat <- 55.133 # 55° 8' N
lon <- -7.450 # 7° 27' W
alt <- 15     # 48 feet converted to meters
metaHead <- "Observer=Coast-guard Boatmen | Instrument=Standardized RIA supply (George Yeates) | Location=Buncrana Coast-guard Station, Co. Donegal | DataSource=RIA Manuscripts 12M16"


raw <- read.csv(paste0(indir, name, "/Buncrana_1850-1851.csv"),
                col.names=c("year", "month", "day", "maxF", "minF", "maxC", "minC"))
head(raw)
tail(raw)

df <- raw %>%
  mutate(
    hour = NA,
    minute = NA,
    obs_time_label = "9a.m.",
    # Meta logic including the standardized instrument source
    meta.Tx = case_when(
      !is.na(maxF) ~ paste0("orig=", maxF, "F | obs.time=", obs_time_label),
      is.na(maxF)  ~ ""
    ),
    meta.Tn = case_when(
      !is.na(minF) ~ paste0("orig=", minF, "F | obs.time=", obs_time_label),
      is.na(minF)  ~ ""
    )
  ) %>% 
  select(year, month, day, hour, minute, Tx = maxC, Tn = minC, meta.Tx, meta.Tn)

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
  keep_na=TRUE
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
  keep_na=TRUE
)

