rm(list=ls())
library(dplyr)
library(readxl)
library(purrr)
library(tidyr)
library(lubridate)
source('/scratch2/ccorbella/code/dataprep/helpfun.R')

outdir <- "/scratch3/PALAEO-RA/daily_data/SEF/DeBilt"

raw <- read.csv("/scratch3/PALAEO-RA/daily_data/original/DeBilt/DBL_SFP.csv", na=c("","NA"))

head(raw)

name <-	"Pfister_DeBilt"
code<-"DeBilt"
lat <- 52.100
lon <- 5.180
alt <- 1

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
  outfile = outfile.name(name, var, df, FALSE),
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
  metaHead = "PTC=Y | PGC=N | merged from Zwanenburg, Haarlem, Den Helder, Delft | homogenized | QC=Y",
)
