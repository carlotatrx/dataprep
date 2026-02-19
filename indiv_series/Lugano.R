rm(list=ls())
library(dplyr)
library(readxl)
library(purrr)
library(tidyr)
library(lubridate)
source('/scratch2/ccorbella/code/dataprep/helpfun.R')

name <- "Lugano"
code <-"LUG"

outdir <- paste0("/scratch3/PALAEO-RA/daily_data/final/", name)

raw <- read.csv(paste0("/scratch3/PALAEO-RA/daily_data/original/", name, "/", code,"_SFP.csv"), na=c("","NA"))

head(raw)

lat <- 46.000
lon <- 8.970
alt <- 273

df <- raw %>%
  mutate(
    year=year(dates),
    month=month(dates),
    day=day(dates),
    Hour=NA,
    minute=NA,
  )

df <- df %>%
  filter(!is.na(SFP))
head(df)


var<-"p"
write_sef_f(
  as.data.frame(df[c("year","month", "day", "Hour", "minute","SFP")]),
  outfile = outfile.name(paste0("Pfister_",name), var, df, FALSE),
  outpath = outdir,
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = "Pfister, L., Wilhelm, L., Brugnara, Y., Imfeld, N., & Brönnimann, S. (2024). Weather type reconstruction using machine learning approaches. EGUsphere, 2024, 1-33.",
  link    = "https://doi.org/10.5194/wcd-6-571-2025",
  nam     = name,
  var     = var,
  stat    = "mean",
  period    = "day",
  units   = units(var),
  metaHead = "PTC=Y | PGC=N | homogenized",
)
