rm(list=ls())
library(dataresqc)
library(readxl)
library(tidyr)
library(dplyr)
library(lubridate)
source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')

outdir <- "/scratch3/PALAEO-RA/daily_data/final/Sevilla"
source <- "Rodrigo, Fernando S.. “Early Meteorological Observations in Southern Spain Version 2.” (2019)."
link   <- "http://hdl.handle.net/10835/6806"

file <- list.files("/scratch3/PALAEO-RA/daily_data/original/Southern_Spain", full.names=TRUE)


# Sevilla_Matute ----------------------------------------------------------


raw <- read_excel(file, sheet=15, skip=8)
name <- "Sevilla"
lat <- round(36 + 23/60, 4)
lon <- round(-(5 + 59/60), 4)
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
    outfile = outfile.name(name, var, df, FALSE, obs_name="_Matute_"),
    outpath = outdir,
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


# Sevilla_diario ---------------------------------------------------------

code <- "Se1812-1830"
source <- " Diario del Gobierno de Sevilla/Diario de Sevilla/Diario de Sevilla de Comercio, Artes y Literatura. 1812, 1813, 1829, 1830. Biblioteca Nacional de España, R60042(11), (http://www.hemerotecadigital.bne.es)."

raw <- read_excel(file, sheet=17, skip=8)

df <- raw %>%
  rename(
    p.orig = `p(EI)`,
    ta.orig = `T(ºR)`,
    hh.orig = Hour,
  ) %>%
  mutate(
    Date = as.Date(sprintf("%04d-%02d-%02d", Year, Month, Day)), # will need later
    ta = round(ta.orig*1.25, 1),
    p  = round(convert_pressure(p.orig, f=25.4, lat=lat, alt=alt, atb=ta), 1),
    dd.norm = dd_normalize(W),
    dd = dd2deg(dd.norm),
    
    orig_time_fmt = ifelse(
      grepl("^[0-9.]+$", hh.orig),
      ifelse(as.numeric(hh.orig) <= 1,
             sprintf("%02d:%02d",
                     floor(as.numeric(hh.orig) * 24),
                     round((as.numeric(hh.orig) * 24 - floor(as.numeric(hh.orig) * 24)) * 60)),
             hh.orig),
      hh.orig
    ),
    meta.time = paste0("orig.time=", orig_time_fmt),
    
    meta.p  = paste0(meta.time, " | orig.p=", p.orig, "English in | atb=", ta, "C"),
    meta.ta = paste0(meta.time, " | orig.ta=", ta.orig, "R"),
    meta.dd = paste0(meta.time, " | orig.dd=", W)
  ) %>% select(Date, Year, Month, Day, hh.orig, ta, p, dd, meta.p, meta.ta, meta.dd)

df.hours <- df %>%
  distinct(Date) %>%
  rowwise() %>%
  mutate(sun = list(suncalc::getSunlightTimes(date = Date, lat = lat, lon = lon,
                                              keep = c("sunrise","sunset")))) %>%
  unnest_wider(sun) %>%
  select(Date, sunrise, sunset)


df.final <- df %>%
  left_join(df.hours, by = "Date") %>%
  mutate(
    Hour_min = dplyr::case_when(
      hh.orig == "Sunrise" ~ format(sunrise, "%H:%M"),
      hh.orig == "Sunset"  ~ format(sunset,  "%H:%M"),
      hh.orig == "Midday"  ~ "12:00",
      grepl(":", hh.orig)  ~ hh.orig,                       # already "H:MM"
      grepl("^[0-9.]+$", hh.orig) ~ {                    # numeric hours or fraction of day
        v <- as.numeric(hh.orig)
        v <- ifelse(v <= 1, v * 24, v)               # treat <=1 as fraction of day
        sprintf("%02d:%02d", floor(v), round((v - floor(v)) * 60))
      },
      TRUE ~ NA_character_
    )
  ) %>%
  tidyr::separate(Hour_min, into = c("Hour","Minute"), sep = ":", convert = TRUE, remove = FALSE) %>%
  select(-sunrise, -sunset, -hh.orig, Hour_min)

head(df.final)


meta_map  <- list(ta = df.final$meta.ta, p = df.final$meta.p, dd=df.final$meta.dd)
base_cols <- df.final[c("Year","Month","Day","Hour","Minute")]
vars <- c("ta", "p", "dd")

for (var in vars) {
  meta  <- meta_map[[var]]
  dat   <- cbind(base_cols, setNames(df.final[var], var))
  write_sef_f(
    dat,
    outfile = outfile.name(name, var, df.final, TRUE, obs_name="_Matute_"),
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
    metaHead = ifelse(var=="p", paste0(metaHead, " | PTC=Y | PGC=Y"), metaHead),
    meta    = meta,
    time_offset= time.offset(lon),
    keep_na = TRUE
  )
}





