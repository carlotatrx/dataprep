rm(list=ls())
library(dataresqc)
library(readxl)
library(tidyr)
library(dplyr)

source('/scratch2/ccorbella/code/dataprep/helpfun.R')

indir <- "/scratch3/PALAEO-RA/daily_data/original/"
name <- "Glendooen"

outdir <- paste0("/scratch3/PALAEO-RA/daily_data/final/",name)

source <- "Mateus, C., Potito, A., & Curley, M. (2020). Reconstruction of a long‐term historical daily maximum and minimum air temperature network dataset for Ireland (1831‐1968). Geoscience Data Journal, 7(2), 102-115."
link   <- "https://www.edepositireland.ie/entities/publication/7271837d-dbe6-463f-b103-422f42a15446"

lat <- 54.950 # Estimated for Glendooen/Letterkenny area
lon <- -7.733 
alt <- 50     # Rough estimate for the Letterkenny area (needs verification if specific site known)
metaHead <- "Observer=Rev. H. Kingsmill | DataSource=Armagh Observatory MS M32"

raw <- read.csv(paste0(indir, name, "/Glendooen_1865-1870.csv"),
                col.names=c("year", "month", "day", "minF","minC"))
head(raw)
tail(raw)

df <- raw %>%
  mutate(
    hour = NA,
    minute = NA,
    obs_time_label = "9a.m.",
    meta.Tn = case_when(
      !is.na(minF) ~ paste0("orig=", minF, "F | obs.time=", obs_time_label),
      is.na(minF)  ~ ""
    )
  ) %>% 
  select(year, month, day, hour, minute, Tn = minC, meta.Tn)

head(df)
tail(df)

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