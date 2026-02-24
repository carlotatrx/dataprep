rm(list=ls())
library(dataresqc)
library(readxl)
library(tidyr)
library(dplyr)

source('/scratch2/ccorbella/code/dataprep/helpfun.R')

indir <- "/scratch3/PALAEO-RA/daily_data/original/"
name <- "Killarney"

outdir <- paste0("/scratch3/PALAEO-RA/daily_data/final/",name)

source <- "Mateus, C., Potito, A., & Curley, M. (2020). Reconstruction of a long‐term historical daily maximum and minimum air temperature network dataset for Ireland (1831‐1968). Geoscience Data Journal, 7(2), 102-115."
link   <- "https://www.edepositireland.ie/entities/publication/7271837d-dbe6-463f-b103-422f42a15446"
lat <- 52.062
lon <- -9.507
alt <- 20.7
metaHead <- paste0(
  "Observer=V. A. G. R. Wynne (1881-1898), G. J. Wynne (1920) | ",
  "Instrument=Stevenson screen, Max/Min thermometers | ",
  "Location=The Park, Killarney; ",
  "0.1F precision used in primary manuscript series."
)

raw <- read.csv(paste0(indir, name, "/Killarney_1881-1933.csv"),
                col.names=c("year", "month", "day", "maxF", "minF", "maxC", "minC"))
head(raw)
tail(raw)

df <- raw %>%
  mutate(
    hour=NA,
    minute=NA,
    meta.Tx=ifelse(year<=1911, paste0("orig=",maxF, "F | obs.time=9a.m."),paste0("orig=",maxF, "F")),
    meta.Tn=ifelse(year<=1911, paste0("orig=",minF, "F | obs.time=9a.m."),paste0("orig=",minF, "F")),
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
  metaHead = metaHead,
  meta = df$meta.Tn
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
  metaHead = metaHead,
  meta=df$meta.Tx
)
