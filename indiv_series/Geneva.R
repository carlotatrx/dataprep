rm(list=ls())
library(dataresqc)
library(lubridate)
library(dplyr)

source('/home/ccorbella/scratch2_symboliclink/code/dataprep/helpfun.R')
outdir <- "/scratch3/PALAEO-RA/daily_data/final/Geneva"


# Cadiz from Lucas --------------------------------------------------------


raw <- read.csv("/scratch3/PALAEO-RA/daily_data/original/Geneva/GVA_SFP.csv", na=c("","NA"))

head(raw)

name <-	"Pfister_Geneva"
code<-"Geneva"
lat <- 46.248
lon <- 6.128
alt <- 410

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
  metaHead = "PTC=Y | PGC=N | merged from various observers/sources in Geneva | homogenized | QC=Y",
)

