rm(list=ls())
library(dataresqc)
library(readxl)
library(tidyr)
library(dplyr)
library(lubridate)
library(tibble)
library(stringr)
source('/scratch2/ccorbella/code/dataprep/helpfun.R')

outdir <- "/scratch3/PALAEO-RA/daily_data/final/Waldenburg"
name <- "Waldenburg"
code <- "GCOS_Waldenburg"
lat	<- 47.385911
lon	<- 7.748753
alt	<- 500
source <- "PALAEO-RA"
link   <- ""
metaHead <- "approximate coordinates | obs.times=unknown"


# 1717-1726  --------------------------------------------------------------

raw <- read_excel("/scratch3/PALAEO-RA/daily_data/original/Waldenburg/BL05_Waldenburg_1776-1790.xlsx",
                  skip=3,
                  col_names=c("Year", "Month", "Day", "pzoll", "plinien", "ta0", "ta_plus", "ta_minus", "notes")
)
raw <- raw %>%
  fill(Year, Month, Day, pzoll, .direction="down")
head(raw)

# add observation numbers
df <- raw %>%
  group_by(Year, Month, Day) %>%
  mutate(obs.num = row_number()) %>%
  ungroup()

df <- df %>%
  mutate(
    Hour=NA_integer_,
    Minute=NA_integer_,
    
    ta.orig = case_when(
      !is.na(ta0) ~ 0,
      !is.na(ta_plus) ~ as.numeric(ta_plus),
      !is.na(ta_minus) ~ -as.numeric(ta_minus),
      TRUE ~ NA
    ),
    ta = round((ta.orig-12)/4,1),
    p = round(convert_pressure(p=(pzoll+plinien/12), f=27.07, lat = lat, alt = alt), 1),
    
    meta.time = paste0("obs.num=", obs.num),
    meta.ta = paste0(meta.time, " | orig.ta=",ta.orig,"MdC"),
    meta.p = paste0(meta.time, " | orig.p=", pzoll, "in", plinien, "l")
  ) %>% select(Year, Month, Day, Hour, Minute, ta, p,meta.ta, meta.p)

head(df)

base_cols <- df[c("Year","Month","Day","Hour","Minute")]

# ta
var <-"ta"
dat <- cbind(base_cols, setNames(df[var], var))
write_sef_f(
  dat,
  outfile = outfile.name(name, var, dat, TRUE),
  outpath = outdir,
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "point",
  units   = "C",
  meta    = df$meta.ta,
  metaHead = paste0(metaHead, " | Instrument=Micheli du Crest thermometer | conversion=(orig.ta-12)/4"),
  keep_na = TRUE
)

# p

var <-"p"
dat <- cbind(base_cols, setNames(df[var], var))
write_sef_f(
  dat,
  outfile = outfile.name(name, var, dat, TRUE),
  outpath = outdir,
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "point",
  units   = units(var),
  meta    = df$meta.p,
  metaHead = paste0(metaHead, " | PGC=Y | PTC=N"),
  keep_na = TRUE
)



