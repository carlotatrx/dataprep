rm(list=ls())
library(dataresqc)
library(readxl)
library(tidyr)
library(dplyr)
library(lubridate)
source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')

outdir <- "/scratch3/PALAEO-RA/daily_data/final/Granada"
name <- "Granada"
code <- "Gr1728-1730"
lat	<- 37.166519
lon	<- -3.600014
alt	<- 600
source <- "Rodrigo, Fernando S.. “The climate of Granada (southern Spain) during the first third of the 18th century (1706–1730) according to documentary sources.” Climate of the Past (2019)."
link   <- "http://hdl.handle.net/10835/6248"
metaHead <- "Observer=Beintema van Peima | Instrument=florentine thermometer"

file <- "/scratch3/PALAEO-RA/daily_data/original/Granada/NavarreteData.xlsx"
raw <- read_excel(file, sheet=3, skip=8)

df <- raw %>%
  rename(
    ta.orig = `T`,
    ta = `T(ºC)`,
    er.ta = `e(T)(°C)`
  ) %>%
  mutate(
    Date = as.Date(Date, format="%d-%b-%Y"),
    Year = year(Date),
    Month = month(Date),
    Day = day(Date),
    Hour = NA_integer_,
    Minute = NA_integer_,
    
    ta = round(ta, 2),
    
    meta = paste0("orig.ta=", ta.orig, " | error.ta=", er.ta, "C")
  ) %>%
  select(Year, Month, Day, Hour, Minute, ta, meta)

head(df)

var <- "ta"

write_sef_f(
  as.data.frame(df),
  outfile  = outfile.name(code, var, df, FALSE),
  outpath  = outdir,
  cod      = code,
  lat      = lat,
  lon      = lon,
  alt      = alt,
  sou      = source,
  link     = link,
  nam      = name,
  var      = var,
  stat     = "point",
  period   = "day",
  units    = "C",
  metaHead = metaHead,
  meta     = df$meta,
  keep_na  = FALSE
)