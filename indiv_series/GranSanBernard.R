rm(list=ls())
library(dplyr)
library(readxl)
library(purrr)
library(tidyr)
source('/scratch2/ccorbella/code/dataprep/helpfun.R')

outdir <- "/scratch3/PALAEO-RA/daily_data/final/GrandStBernard"

raw <- read.csv("/scratch3/PALAEO-RA/daily_data/tmp/GCOS_GrandStBernard_18170914-18251231_p_subdailyNO.tsv",
                skip=12, sep="\t")
meta <- read_meta("/scratch3/PALAEO-RA/daily_data/tmp/GCOS_GrandStBernard_18170914-18251231_p_subdailyNO.tsv")

head(raw)

df <- raw %>%
  mutate(
    p=round(convert_pressure(Value),2)
  ) %>% select(Year, Month, Day, Hour, Minute, p, Meta)
head(df)
meta

var <-"p"

write_sef_f(
  df,
  outpath=outdir,
  outfile=outfile.name("GCOS_GrandStBernard",var, df, TRUE),
  cod     = meta[["id"]],
  lat     = meta[["lat"]],
  lon     = meta[["lon"]],
  alt     = meta[["alt"]],
  sou     = meta[["source"]],
  link    = "https://doi.pangaea.de/10.1594/PANGAEA.948258",
  nam     = meta[["name"]],
  var     = var,
  stat    = "point",
  units   = units(var),
  metaHead = meta[["meta"]],
  meta=df$Meta,
  keep_na = TRUE
)




raw <- read_sef("/scratch3/PALAEO-RA/daily_data/tmp/DIGIHOM_GrandStBernard_18620101-18631130_p_daily.tsv")
meta <- read_meta("/scratch3/PALAEO-RA/daily_data/tmp/DIGIHOM_GrandStBernard_18620101-18631130_p_daily.tsv")
head(raw)

df <- raw %>%
  mutate(
    Hour=NA,
    Minute=NA, 
    p=round(convert_pressure(Value),2)
  ) %>% select(Year, Month, Day, Hour, Minute, p)
head(df)

var <-"p"

write_sef_f(
  df,
  outpath=outdir,
  outfile=outfile.name("DIGIHOM_GrandStBernard",var, df, FALSE),
  cod     = meta[["id"]],
  lat     = meta[["lat"]],
  lon     = meta[["lon"]],
  alt     = meta[["alt"]],
  sou     = meta[["source"]],
  link    = "",
  nam     = meta[["name"]],
  var     = var,
  stat    = "mean",
  period = "day",
  units   = units(var),
  metaHead = meta[["meta"]],
  keep_na = TRUE
)
