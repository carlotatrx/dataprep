rm(list=ls())
library(dataresqc)
library(readr)
library(dplyr)
library(lubridate)
source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')

file <- '/scratch3/PALAEO-RA/daily_data/original/Trieste/Trieste_PPPP_daily.tab'

raw <- suppressWarnings(read_delim(
  file = file,
  delim = '\t',
  skip = 15,
  col_types = cols(.default = col_guess()),
  trim_ws = T,
  guess_max = 100000
))%>% rename(
  Date = 1,
  p = 2,
  n_miss = 3,
  corr = 4
)

df <- raw %>%
  transmute(
    Year = year(Date),
    Month = month(Date),
    Day = day(Date),
    Hour = NA,
    Minute = NA,
    Value = p,
    meta = paste0("hom.adj=", corr)
  )

write_sef_f(as.data.frame(df),
            outfile=paste0('Trieste_', get_date_range(df), "_p_daily.tsv"),
            outpath="/scratch3/PALAEO-RA/daily_data/final/Trieste/",
            variable = "p",
            cod = "Trieste",
            nam = "Trieste",
            lat = 45.646400,
            lon = 13.764000,
            alt = NA,
            sou = "Raicich, Fabio; Colucci, Renato R (2021): Daily mean-sea-level atmospheric pressure from 1841 to 2018 at Trieste, Italy [dataset]",
            link = "https://doi.org/10.1594/PANGAEA.926896",
            units = "C",
            stat = "mean",
            period = "day",
            metaHead = "PTC=Y | PGC=Y",
            meta = df$meta,
            keep_na = F)
