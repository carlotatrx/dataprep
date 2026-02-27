rm(list=ls())
library(dataresqc)
library(readxl)
library(tidyr)
library(dplyr)

source('/scratch2/ccorbella/code/dataprep/helpfun.R')

indir <- "/scratch3/PALAEO-RA/daily_data/original/"
name <- "Scattery"

outdir <- paste0("/scratch3/PALAEO-RA/daily_data/final/",name)

source <- "Mateus, C., Potito, A., & Curley, M. (2020). Reconstruction of a long‐term historical daily maximum and minimum air temperature network dataset for Ireland (1831‐1968). Geoscience Data Journal, 7(2), 102-115."
link   <- "https://www.edepositireland.ie/entities/publication/7271837d-dbe6-463f-b103-422f42a15446"
lat <- 52.617 # (Approximate, based on Scattery Island lighthouse)
lon <- -9.517 # (Approximate, based on Scattery Island lighthouse)
alt <- 10 # (Estimated elevation for coastal station)
metaHead <- paste0(
  "Observer=Coast-guard boatmen | Instrument=Self-registering thermometers | ",
  "DataSource=Royal Irish Academy Manuscript 12 M 16 | ",
  "Location=Scattery Island",
  "Notes=primarily whole degree precision; additional note in metadata regarding high Min bias in early 1852."
)
metaHead


raw <- read.csv(paste0(indir, name, "/Scattery_1850-1852.csv"), header=FALSE,
                skip=1,fileEncoding = "latin1",
                col.names=c("year", "month", "day", "maxF", "minF", "maxC", "minC", "extra"))

head(raw)
tail(raw)

df <- raw %>%
  mutate(
    hour = NA,
    minute = NA,
    obs_time_label = "9a.m.",
    meta.Tx = case_when(
      !is.na(maxF) ~ paste0("orig=", maxF, "F | obs.time=", obs_time_label),
      is.na(maxF)  ~ ""
    ),
    meta.Tn = case_when(
      !is.na(minF) ~ paste0("orig=", minF, "F | obs.time=", obs_time_label),
      is.na(minF)  ~ ""
    )
  ) %>% 
  select(year, month, day, hour, minute, Tx = maxC, Tn = minC, meta.Tx, meta.Tn) %>%
  filter(!is.na(Tx) | !is.na(Tn))

head(df)
tail(df)

var <- "Tx"
df.Tx <- df[c("year","month", "day", "hour", "minute","Tx", "meta.Tx")] %>% filter(!is.na(Tx))
write_sef_f(
  as.data.frame(df.Tx),
  outfile = outfile.name(paste0("ILMMT_",name), var, df.Tx, FALSE),
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
  meta = df.Tx$meta.Tx,
  metaHead=metaHead,
  keep_na=TRUE
)


var <- "Tn"
df.Tn <- df[c("year","month", "day", "hour", "minute","Tn", "meta.Tn")] %>% filter(!is.na(Tn))
write_sef_f(
  as.data.frame(df.Tn),
  outfile = outfile.name(paste0("ILMMT_",name), var, df.Tn, FALSE),
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
  meta = df.Tn$meta.Tn,
  metaHead=metaHead,
  keep_na=TRUE
)


