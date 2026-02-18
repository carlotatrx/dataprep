rm(list=ls())
library(dplyr)
library(readxl)
library(purrr)
library(tidyr)
library(dataresqc)
source('/scratch2/ccorbella/code/dataprep/helpfun.R')

outdir <- "/scratch3/PALAEO-RA/daily_data/final/Armagh"


lat <- round(54+21.2/60, 4)
lon <- round(6+38.9/60,4)
alt <- 64
name <- "Yuri_Armagh"
code <- "Armagh"
source <- "Butler, C. J., García Suárez, A. M., Coughlin, A. D. S., & Morrell, C. (2005). Air temperatures at Armagh observatory, northern Ireland, from 1796 to 2002. International Journal of Climatology: A Journal of the Royal Meteorological Society, 25(8), 1055-1079."
link   <- "DOI: 10.1002/joc.1148"
raw <- read_excel("/scratch3/PALAEO-RA/daily_data/original/Armagh/Armagh.xls",
                  skip=2,
                  col_names = c("Year", "Month", "Day", "porig8", "porig12", "porig20", "taborig8", "taborig12", "taborig20", "taorig8", "taorig12", "taorig20", "notes"))
df <- raw %>%
  pivot_longer(
    cols = porig8:taorig20, 
    names_to = c(".value", "Hour"), 
    names_pattern = "(porig|taborig|taorig)(\\d+)"
  ) %>% mutate(
    Hour=ifelse(is.na(notes), as.integer(Hour), 14L),
    Minute = 0L,

  p = round(convert_pressure(porig, f=27.07, lat=lat, alt=alt, atb=(taborig-32)/1.8),2),
  ta = round((taorig- 32) / 1.8,2),
  
  meta.p = paste0(meta_time(Hour, Minute), " | orig=",porig,"in | atb=", taborig, "F"),
  meta.ta = paste0(meta_time(Hour, Minute), " | orig=",taorig, "F"),
  ) %>% select(Year, Month, Day, Hour, Minute, p, ta, meta.ta, meta.p)

head(df)

# save
var <-"p"
write_sef_f(
  as.data.frame(df[,c("Year","Month", "Day", "Hour", "Minute", "p")]),
  outfile = outfile.name(name, var, df, TRUE),
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
  metaHead = "PTC=Y | PGC=Y",
  meta    = df$meta.p,
  keep_na = TRUE
)

var <-"ta"
write_sef_f(
  as.data.frame(df[,c("Year","Month", "Day", "Hour", "Minute", "ta")]),
  outfile = outfile.name(name, var, df, TRUE),
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
  meta    = df$meta.ta,
  keep_na = TRUE
)

