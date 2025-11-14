rm(list=ls())
library(dataresqc)
library(readxl)
library(tidyr)
library(dplyr)
library(lubridate)
source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')

outdir <- "/scratch3/PALAEO-RA/daily_data/final/Southern_Spain"
source <- "Rodrigo, Fernando S.. “Early Meteorological Observations in Southern Spain Version 2.” (2019)."
link   <- "http://hdl.handle.net/10835/6806"

file <- list.files("/scratch3/PALAEO-RA/daily_data/original/Southern_Spain", full.names=TRUE)


# Sevilla_Matute ----------------------------------------------------------


raw <- read_excel(file, sheet=15, skip=8)
name <- "Sevilla"
lat <- 36 + 23/60
lon <- -(5 + 59/60)
alt <- 11
code <- "Se1803-1830"
metaHead <- "Obsever=Justino Matute"
source <- "Mature, J. (ed.). 1803-1807. El Correo Literario y Económico de Sevilla. Biblioteca de la Universidad de Sevilla, A059/042-051 (http://fondosdigitales.us.es)"

df <- raw %>%
  rename(
    pl  = `p(lin)`,
    pin = `p(EI)`,
    Tx.orig = `Tmax(°R)`,
    Tn.orig = `Tmin(°R)`,
  ) %>%
  mutate(
    Hour = 12,
    Minute = 0,
    
    Tx = round(Tx.orig*1.25),
    Tn = round(Tn.orig*1.25),
    
    p.orig = pin + pl / 12,
    p = round(convert_pressure(p.orig, f=25.4, lat=lat, alt=alt),1),
    
    dd.norm = dd_normalize(W),
    dd = dd2deg(dd.norm),
    
    meta.Tx = paste0("obs.time=noon | orig.Tx=", Tx.orig, "R"),
    meta.Tn = paste0("obs.time=noon | orig.Tn=", Tn.orig, "R"),
    
    meta.p  = paste0("obs.time=noon | orig.p=", pin,"English in", pl,"l"),
    meta.dd = paste0("obs.time=noon | orig.dd=", W, ifelse(is.na(W2), "", paste0(" | alternative.dd=", W2)))
  )

head(df)


meta_map  <- list(Tx = df$meta.Tx, p = df$meta.p, Tn=df$meta.Tn, dd = df$meta.dd)
base_cols <- df[c("Year","Month","Day","Hour","Minute")]
vars <- c("Tn", "Tx", "p", "dd")

for (var in vars) {
  meta  <- meta_map[[var]]
  dat   <- cbind(base_cols, setNames(df[var], var))
  write_sef_f(
    dat,
    outfile = outfile.name(name, var, df, FALSE, obs_name="_JMChacon_"),
    outpath = '/scratch3/PALAEO-RA/daily_data/final/Cadiz',
    cod     = code,
    lat     = lat,
    lon     = lon,
    alt     = alt,
    sou     = source,
    link    = link,
    nam     = name,
    var     = var,
    stat    = case_when(
      var=="Tx" ~ "maximum",
      var=="Tn" ~ "minimum",
      TRUE ~ "point"
    ),
    units   = units(var),
    metaHead = ifelse(var=="p", paste0(metaHead, " | PTC=Y | PGC=Y"), metaHead),
    meta    = meta,
    time_offset= time.offset(lon),
    keep_na = F
  )
}
