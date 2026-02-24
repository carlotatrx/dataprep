rm(list=ls())
library(dataresqc)
library(readxl)
library(tidyr)
library(dplyr)

source('/scratch2/ccorbella/code/dataprep/helpfun.R')

indir <- "/scratch3/PALAEO-RA/daily_data/original/"
name <- "Markree"

outdir <- paste0("/scratch3/PALAEO-RA/daily_data/final/",name)

source <- "Mateus, C., Potito, A., & Curley, M. (2020). Reconstruction of a long‐term historical daily maximum and minimum air temperature network dataset for Ireland (1831‐1968). Geoscience Data Journal, 7(2), 102-115."
link   <- "https://www.edepositireland.ie/entities/publication/7271837d-dbe6-463f-b103-422f42a15446"
lat <- 54.173
lon <- -8.457
alt <- 38

metaHead <- paste0(
  "Observer=Edward J. Cooper (1833-1863), Albert Marth & W. Doberck (1874-1883), ",
  "F.W. Henkel (1898-1902) | Instrument=Stevenson screen, Max/Min thermometers"
)

raw <- read.csv(paste0(indir, name, "/Markree_1850-1968.csv"),
                col.names=c("year", "month", "day", "maxF", "minF", "maxC", "minC"))
head(raw)
tail(raw)

df <- raw %>%
  mutate(
    hour=NA,
    minute=NA,
    meta.Tx=case_when(
      year<=1863 ~ paste0("orig=",maxF, "F | obs.time=8:30a.m./3:30p.m."),
      TRUE~ paste0("obs.time=9a.m.")
    ),
    meta.Tn=case_when(
      year<=1863 ~ paste0("orig=",minF, "F | obs.time=8:30a.m./3:30p.m."),
      TRUE~ paste0("obs.time=9a.m.")
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

