rm(list=ls())
library(dplyr)
library(readxl)
library(purrr)
library(tidyr)
source('/scratch2/ccorbella/code/dataprep/helpfun.R')

outdir <- "/scratch3/PALAEO-RA/daily_data/final/Upminster"

raw <- read_excel("/scratch3/PALAEO-RA/daily_data/original/Upminster/Upminster_1697.xlsx", skip=2,
                  na=c("","NA"))

head(raw)

name <-	"Upminster"
code<-"Upminster in Essex"
lat <- 51.55591
lon <-	0.248894
alt <-	19

df <- raw %>%
  fill(day, .direction="down") %>%
  mutate(
    minute=0,
    p=round(convert_pressure(pressure,f=25.04, lat=lat, alt=alt),2),
    meta.p=paste0(meta_time(time, minute), " | orig=", pressure, " English inch")
  )

head(df)

var<-"p"
write_sef_f(
  as.data.frame(df[c("year","month", "day", "time", "minute","p")]),
  outfile = outfile.name(name, var, df, TRUE),
  outpath = outdir,
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = "",
  link    = "",
  nam     = name,
  var     = var,
  stat    = "point",
  units   = units(var),
  metaHead = "PTC=N | PGC=Y",
  meta    = df$meta.p,
  time_offset = time.offset(lon),
  keep_na = TRUE
)

var<-"rr"
write_sef_f(
  as.data.frame(df[c("year","month", "day", "time", "minute","rain")]),
  outfile = outfile.name(name, var, df, TRUE),
  outpath = outdir,
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = "",
  link    = "",
  nam     = name,
  var     = var,
  stat    = "point",
  units   = "unknown",
  metaHead = "",
  keep_na = TRUE
)
