rm(list = ls())
library(lubridate)
library(dplyr)
library(readr)
library(tidyr)
library(dataresqc)
library(readxl)
library(XLConnect)

source("/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R")

code <-	"Küsnacht"
name <-	"Küsnacht"
lat <-	47.31701
lon	<- 8.58323
alt	<- 420
source <-	"CHIMES"
link <-	"https://doi.pangaea.de/10.1594/PANGAEA.948258"
metaHead <- "Observer=Zollinger,others"


file <- "/scratch3/PALAEO-RA/daily_data/original/Kuesnacht/ZH08_Kuesnacht_1852-1856.xlsx"
raw <- readWorksheetFromFile(file, sheet=1, startRow=3, header=T) %>%
  select(Year, Month, Day, `time.1`:`time.5`, D1:D5, S1:S5, H1:H5)

head(raw)

df <- raw %>%
  mutate(
    across(`time.1`:`time.5`, as.numeric)
  ) %>%
  pivot_longer(
    cols = matches("^(time\\.|[DSH])[0-9]"),
    names_to = c(".value", "tod"),
    names_pattern = "^([A-Za-z]+)[\\.]?([0-9])$"
  ) %>% 
  mutate(
    Hour = time,
    Minute = 0L,
    
    ddnorm = dd_normalize(D),
    ddfinal = dd2deg(ddnorm),
    
    rh = ifelse(H=="NA", NA_integer_, H),
    
    meta.dd = paste0("orig.time=", sprintf('%02d', Hour), ":", sprintf('%02d', Minute),
                     " | orig.dd=", D, " | force=", S),
    meta.time = paste0("orig.time=", sprintf('%02d', Hour), ":", sprintf('%02d', Minute))
  ) %>% select(Year, Month, Day, Hour, Minute, ddfinal, meta.dd, rh, meta.time)
head(df)

dat <- as.data.frame(df[c("Year", "Month", "Day", "Hour", "Minute", "ddfinal")])
var <- "dd"
write_sef_f(
  dat,
  outfile = outfile.name(name, var, dat, TRUE),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Kuesnacht',
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
  meta    = df$meta.dd,
  metaHead = metaHead,
  time_offset = time.offset(lon),
  keep_na = FALSE
)

dat <- as.data.frame(df[c("Year", "Month", "Day", "Hour", "Minute", "rh")])
var <- "rh"
write_sef_f(
  dat,
  outfile = outfile.name(name, var, dat, TRUE),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Kuesnacht',
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
  meta    = df$meta.time,
  metaHead = metaHead,
  time_offset = time.offset(lon),
  keep_na = FALSE
)
