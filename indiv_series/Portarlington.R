rm(list=ls())
library(dataresqc)
library(readxl)
library(tidyr)
library(dplyr)

source('/scratch2/ccorbella/code/dataprep/helpfun.R')

indir <- "/scratch3/PALAEO-RA/daily_data/original/"
name <- "Portarlington"

outdir <- paste0("/scratch3/PALAEO-RA/daily_data/final/",name)

source <- "Mateus, C., Potito, A., & Curley, M. (2020). Reconstruction of a long‐term historical daily maximum and minimum air temperature network dataset for Ireland (1831‐1968). Geoscience Data Journal, 7(2), 102-115."
link   <- "https://www.edepositireland.ie/entities/publication/7271837d-dbe6-463f-b103-422f42a15446"
lat <- 53.150 # (53° 9' N)
lon <- -7.200 # (7° 12' W)
alt <- 70.1 # (230 ft)

metaHead <- paste0(
  "Observer=Dr. Hanlon (1845-1864) | Instrument=Max/Min thermometers | ",
  "DataSource=Dublin Medical Press newspaper; high-quality manuscripts from the Royal Irish Academy network (1851-1852, overlap) | ",
  "Notes=",
  "whole degree precision"
)


metaHead


raw <- read.csv(paste0(indir, name, "/Portarlington_1845-1864.csv"), header=FALSE,
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
      !is.na(maxF) ~ paste0("orig=", maxF),
      is.na(maxF)  ~ ""
    ),
    meta.Tn = case_when(
      !is.na(minF) ~ paste0("orig=", minF),
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

indir <- '/scratch3/PALAEO-RA/daily_data/final/'
files <- list.files(paste0(indir, name), pattern = "daily\\.tsv$", full.names = TRUE)
files

for (f in files) {
  qc(f, outpath = paste0(indir,name))
  
  qcfile <- file.path(
    indir,
    paste0("qc_", tools::file_path_sans_ext(basename(f)), ".txt")
  )
}
qc_files  <- list.files(paste0(indir, name), pattern="^qc_.*\\.txt$", full.names=TRUE)
qc_files


