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


# Granada_Dalmau ----------------------------------------------------------

lat <- round(37 + 10/60,4)
lon <- round(-(3+36/60),4)
alt <- 660
code <- "Gr1796-1797"
metaHead <- "Observer=Francisco Dalmau"
source   <- "El Mensagero Económico y Erudito de Granada. Universidad de Granada FLA F-9-9-3 http://www.bibliotecavirtualdeandalucia.es"
file <- list.files("/scratch3/PALAEO-RA/daily_data/original/Southern_Spain", full.names=TRUE)
outdir <- "/scratch3/PALAEO-RA/daily_data/final/Granada"

raw <- read_excel(file, sheet=code, skip=8)

df <- raw %>%
  rename(
    ta.orig = `T(°R)`,
    pin = `P(p)`,
    pl  = `p(l)`,
  ) %>%
  mutate(
    Hour = 12L,
    Minute = 0L,
    
    ta = round(ta.orig*1.25, 1),
    p  = round(convert_pressure(pin+pl/12, f=27.07, lat=lat, alt=alt, atb=ta),1),
    dd.norm = dd_normalize(W),
    dd = dd2deg(dd.norm),
    
    meta.ta = paste0("orig.time=noon | orig.ta=", ta.orig, "R"),
    meta.p  = paste0("orig.time=noon | orig.p=", pin,"in", pl,"l | atb=", ta, "C"),
    meta.dd = paste0("orig.time=noon | orig.dd=", W)
  ) %>%
  select(Year, Month, Day, Hour, Minute, ta, p, dd, meta.ta, meta.p, meta.dd)

head(df)

meta_map  <- list(ta = df$meta.ta, p = df$meta.p, dd = df$meta.dd)
base_cols <- df[c("Year","Month","Day","Hour","Minute")]
vars <- c("ta", "p", "dd")

for (var in vars) {
  meta  <- meta_map[[var]]
  dat   <- cbind(base_cols, setNames(df[var], var))
  write_sef_f(
    dat,
    outfile = outfile.name(name, var, df, subdaily="noon", obs_name="_Dalmau_"),
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
    metaHead= ifelse(var=="p", paste0(metaHead, " | PTC=Y | PGC=Y"), metaHead),
    meta    = meta,
    time_offset = time.offset(lon),
    keep_na = TRUE
  )
}



# Granada_diario_1813 -----------------------------------------------------


source <- "Diario de Granada: El Publicista. 1812-1813. Biblioteca de Andalucia http://www.bibliotecavirtualdeandalucia.es"
code   <- "Gr1813"

raw <- read_excel(file, sheet=code, skip=8)

df <- raw %>%
  rename(
    ta.orig = `T(ºR)`,
  ) %>%
  mutate(
    Hour = 12L,
    Minute = 0L,
    
    ta = round(ta.orig*1.25, 1),
    
    meta.ta = paste0("orig.time=noon | orig.ta=", ta.orig, "R"),
  ) %>%
  select(Year, Month, Day, Hour, Minute, ta, meta.ta)

head(df)

var <- "ta"

write_sef_f(
  as.data.frame(df),
  outfile  = outfile.name(code, var, df, subdaily="noon"),
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
  period   = "point",
  units    = "C",
  meta     = df$meta.ta,
  time_offset = time.offset(lon),
  keep_na  = FALSE
)



# Granada_diario_1820 -----------------------------------------------------

code   <- "Gr1820"
source <- "Diario Constitucional de Granada. 1820. Biblioteca de Andalucia http://www.bibliotecavirtualdeandalucia.es"

raw <- read_excel(file, sheet=code, skip=8)

df <- raw %>%
  rename(
    ta.orig = `T(ºR)`,
    pin = `P(FI)`,
    pl  = `P(Fl)`,
  ) %>%
  mutate(
    Minute = 0L,
    
    ta = round(ta.orig*1.25, 1),
    p  = round(convert_pressure(pin+pl/12, f=27.07, lat=lat, alt=alt, atb=ta),1),
    dd.norm = dd_normalize(W),
    dd = dd2deg(dd.norm),
    
    meta.ta = paste0("orig.time=",Hour," | orig.ta=", ta.orig, "R"),
    meta.p  = paste0("orig.time=",Hour," | orig.p=", pin,"French in", pl,"l | atb=", ta, "C"),
    meta.dd = paste0("orig.time=",Hour," | orig.dd=", W)
  ) %>%
  select(Year, Month, Day, Hour, Minute, ta, p, dd, meta.ta, meta.p, meta.dd)

head(df)

meta_map  <- list(ta = df$meta.ta, p = df$meta.p, dd = df$meta.dd)
base_cols <- df[c("Year","Month","Day","Hour","Minute")]
vars <- c("ta", "p", "dd")

for (var in vars) {
  meta  <- meta_map[[var]]
  dat   <- cbind(base_cols, setNames(df[var], var))
  write_sef_f(
    dat,
    outfile = outfile.name(name, var, df, subdaily=TRUE),
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
    metaHead= ifelse(var=="p", "PTC=Y | PGC=Y", ""),
    meta    = meta,
    time_offset = time.offset(lon),
    keep_na = TRUE
  )
}




