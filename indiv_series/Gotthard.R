rm(list = ls())
library(lubridate)
library(dplyr)
library(readr)
library(tidyr)
library(dataresqc)
library(readxl)
library(XLConnect)

source("/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R")

files <- list.files("/scratch3/PALAEO-RA/daily_data/original/Gotthard",
                    full.names = TRUE)

code <-	"TI01_Gotthard"
name <-	"Gotthardpass"
lat <-	46.55535
lon	<- 8.56791
alt	<- 2093
source <-	"CHIMES"
link <-	"https://doi.org/10.5194/cp-15-1345-2019"
metaHead <- "Observer=Pater Onuphrius, Pater Laurentius Mediolanensis, Jos. Belmas de Caladray"


file <- files[1]
raw <- readWorksheetFromFile(file, sheet=1, startRow=12, header=T) %>%
  rename(
    dd1 = 19,
    dd2 = 20,
    dd3 = 21,
    w1  = 22,
    w2  = 23,
    w3  = 24
  ) %>% select(Year, Month, Day, dd1, dd2, dd3, w1, w2, w3)

head(raw)

df <- raw %>%
  pivot_longer(
    cols = c(dd1, dd2, dd3, w1, w2, w3),
    names_to = c(".value", "tod"),
    names_pattern = "^([a-z]+)([123])$",
  ) %>% 
  mutate(
    Hour = case_when(
      tod=="1" ~ 7L,
      tod=="2" ~ 14L,
      tod=="3" ~ 21L
    ),
    Minute = 0L,
    
    ddnorm = dd_normalize(dd),
    ddfinal = dd2deg(ddnorm),
    meta.dd = paste0("orig.time=", sprintf('%02d', Hour), ":", sprintf('%02d', Minute),
                     " | orig.dd=", dd, " | force=", w)
  ) %>% select(Year, Month, Day, Hour, Minute, Value=ddfinal, meta.dd)
head(df)

var <- "dd"
write_sef_f(
  as.data.frame(df),
  outfile = outfile.name(name, var, df, TRUE),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Gotthard',
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

