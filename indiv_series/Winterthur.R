rm(list = ls())
library(lubridate)
library(dplyr)
library(readr)
library(tidyr)
library(dataresqc)
library(readxl)
library(XLConnect)

source("/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R")

code <-	"Winterthur"
name <-	"Winterthur"
lat <-	47.49922
lon	<- 8.72897
alt	<- 443
source <-	"GCOS"
link <-	"https://doi.pangaea.de/10.1594/PANGAEA.948258"
metaHead <- "Observer=Furrer"


file <- "/scratch3/PALAEO-RA/daily_data/original/Winterthur/ZH03_Winterthur_1857-1867.xlsx"
raw <- readWorksheetFromFile(file, sheet=1, startRow=3, header=T) %>%
  select(Year, Month, Day, c(21:23)) %>%
  rename(
    dd1 = 4,
    dd2 = 5,
    dd3 = 6
  )

head(raw)

df <- raw %>%
  pivot_longer(
    cols = matches("dd[1-3]"),
  ) %>% 
  mutate(
    Hour = case_when(
      name == "dd1" ~ 9,
      name == "dd2" ~ 12,
      name == "dd3" ~ 16
    ),
    Minute = 0L,
    
    ddnorm = dd_normalize(value),
    ddfinal = dd2deg(ddnorm),
    

    meta.dd = paste0("orig.time=", sprintf('%02d', Hour), ":", sprintf('%02d', Minute),
                     " | orig.dd=", value),
  ) %>% select(Year, Month, Day, Hour, Minute, ddfinal, meta.dd)
head(df)

dat <- as.data.frame(df[c("Year", "Month", "Day", "Hour", "Minute", "ddfinal")])
var <- "dd"
write_sef_f(
  dat,
  outfile = outfile.name(name, var, dat, TRUE),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Winterthur',
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
