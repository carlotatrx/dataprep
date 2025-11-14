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


# Cadiz_Chacón -------------------------------------------------------------------

raw <- read_excel(file, sheet=2, skip=8)
name <- "Cadiz"
lat <- 36 + 31/60
lon <- -(6 + 17/60)
alt <- 13
code <- "Ca1799-1800"
metaHead <- "Obsever=José María Chacón | obs.time=noon"
source <- "Aréjula, J.M. 1806. Breve descripción de la fiebre amarilla padecida en Cádiz. Biblioteca de Andalucía ANT-XIX-614 http://www.biblitecavirtualdeandalucia.es"

df <- raw %>%
  rename(
    pin = `p(FI)`,
    pl = `p(Fl)`,
    ta.orig = `T(ºF)`
  ) %>%
  mutate(
    Hour = 12,
    Minute = 0,
    
    ta = round((ta.orig-32)/1.8, 1),
    
    p.orig = pin + pl / 12,
    p = round(convert_pressure(p.orig, f=27.07, lat=lat, alt=alt, atb=ta),1),
    
    dd.norm = dd_normalize(W),
    dd = dd2deg(dd.norm),
    
    meta.ta = paste0("orig.ta=", ta.orig, "F"),
    meta.p  = paste0("orig.p=", pin,"in", pl,"l | atb=", ta, "C"),
    meta.dd = paste0("orig.dd=", W, ifelse(is.na(W2), "", paste0(" | alternative.dd=", W2)))
  )

head(df)

meta_map  <- list(ta = df$meta.ta, p = df$meta.p, dd = df$meta.dd)
base_cols <- df[c("Year","Month","Day","Hour","Minute")]
vars <- c("ta", "p", "dd")

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
    stat    = "point",
    period  = "day",
    units   = units(var),
    metaHead = ifelse(var=="p", paste0(metaHead, " | PTC=Y | PGC=Y | orig.unit=French inches"), metaHead),
    meta    = meta,
    keep_na = F
  )
}



# Cadiz_González ----------------------------------------------------------


raw <- read_excel(file, sheet=3, skip=8)
name <- "Cadiz"
lat <- 36 + 31/60
lon <- -(6 + 17/60)
alt <- 13
code <- "Ca1800"
metaHead <- "Obsever=Pedro María González | obs.time=noon"
source <- "González, PM. 1801. Disertación médica sobre la calentura maligna contagiosa que reynó en Cádiz en el año de 1800. Biblioteca Universidad Complutense BH FG 1272"

df <- raw %>%
  rename(
    w1 = `W(sunrise)`,
    w2 = `W(sunset)`,
    ta.orig = `T(ºF)`
  ) %>%
  mutate(
    Hour = 12,
    Minute = 0,
    
    ta = round((ta.orig-32)/1.8, 1),
    
    meta.ta = paste0("orig.ta=", ta.orig, "F"),

  )

head(df)

base_cols <- df[c("Year","Month","Day","Hour","Minute")]
vars <- c("ta")

for (var in vars) {
  meta  <- meta_map[[var]]
  dat   <- cbind(base_cols, setNames(df[var], var))
  write_sef_f(
    dat,
    outfile = outfile.name(name, var, df, TRUE, obs_name="_PMGonzalez_"),
    outpath = '/scratch3/PALAEO-RA/daily_data/final/Cadiz',
    cod     = code,
    lat     = lat,
    lon     = lon,
    alt     = alt,
    sou     = source,
    link    = link,
    nam     = name,
    var     = var,
    stat    = "point",
    period  = "day",
    units   = units(var),
    metaHead = metaHead,
    meta    = df$meta.ta,
    keep_na = F
  )
}


df.dd <- df %>%
  pivot_longer(
    cols = c(w1,w2),
    names_to = "tod",
    names_prefix = "w",
    values_to = "dd.orig"
  ) %>% mutate(
    Date = as.Date(paste(Year, Month, Day, sep = "-")),
    meta.time = case_when(
      tod=="1" ~ "sunrise",
      tod=="2" ~ "sunset",
    ),
    
    dd.norm = dd_normalize(dd.orig),
    dd = dd2deg(dd.norm),
    
    meta   = paste0("obs.time=", meta.time, " | orig.dd=", dd.orig)
  ) %>% select(Date, Year, Month, Day, Hour, Minute, tod, Value=dd, meta)

df.hours <- df %>%
  mutate(Date = as.Date(paste(Year, Month, Day, sep = "-"))) %>%
  rowwise() %>%
  mutate(
    # get sunrise/sunset for that day
    sun = list(getSunlightTimes(date = Date, lat = lat, lon = lon, keep = c("sunrise", "sunset")))
  ) %>%
  unnest_wider(sun) %>%
  select(Date, sunrise, sunset)

df.final <- df.dd %>%
  left_join(df.hours, by = "Date") %>%
  mutate(
    # assign actual datetime of observation based on sunrise/sunset
    obs.time = case_when(
      tod == "1" ~ sunrise,
      tod == "2" ~ sunset
    ),
    Hour = hour(obs.time),
    Minute = minute(obs.time)
  ) %>%
  select(Year, Month, Day, Hour, Minute, Value, meta)

head(df.final)

var <- "dd"

write_sef_f(
  as.data.frame(df.final),
  outfile = outfile.name(name, var, df.final, TRUE, obs_name="_PMGonzalez_"),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Cadiz',
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
  metaHead = metaHead,
  meta    = df.final$meta,
  keep_na = F
)


# Cadiz_diario ------------------------------------------------------------

raw <- read_excel(file, sheet=4, skip=8)
name <- "Cadiz"
lat <- 36 + 31/60
lon <- -(6 + 17/60)
alt <- 13
code <- "Ca1802-1803"
source <- "Diario Comercial de Cádiz (1802-1830), Cádiz Library, sgn: FL-PP-Est.59. http://www.bibliotecavirtualdeandalucia.es"

df <- raw %>%
  rename(
    pin = `P(FI)`,
    pl = `P(Fl)`,
    ta.orig = `T(ºF)`
  ) %>%
  mutate(
    Minute = 0,
    ta = round((ta.orig-32)/1.8, 1),
    
    p.orig = pin + pl / 12,
    p = round(convert_pressure(p.orig, f=27.07, lat=lat, alt=alt, atb=ta),1),
    
    dd.norm = dd_normalize(W),
    dd = dd2deg(dd.norm),
    
    meta.time = paste0("orig.time=", Hour),
    meta.ta = paste0(meta.time, " | orig.ta=", ta.orig, "F"),
    meta.p  = paste0(meta.time, " | orig.p=", pin,"in", pl,"l | atb=", ta, "C"),
    meta.dd = paste0(meta.time, " | orig.dd=", W)
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
    outfile = outfile.name(name, var, df, TRUE),
    outpath = '/scratch3/PALAEO-RA/daily_data/final/Cadiz',
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
    meta    = meta,
    metaHead= ifelse(var=="p", "PTC=Y | PGC=Y | orig.unit=French inches", ""),
    time_offset = time.offset(lon),
    keep_na = F
  )
}


# Cadiz_Aréjula -----------------------------------------------------------

raw <- read_excel(file, sheet=5, skip=8)
name <- "Cadiz"
lat <- 36 + 31/60
lon <- -(6 + 17/60)
alt <- 13
code <- "CaA1803"
source <- "Aréjula, J.M. 1806. Breve descripción de la fiebre amarilla padecida en Cádiz. Biblioteca de Andalucía ANT-XIX-614 http://www.biblitecavirtualdeandalucia.es"

df <- raw %>%
  rename(
    pin = `p(EI)`,
    pl = `p(EL)`,
    pt = `p(ELt)`,
    ta.orig = `T(ºR)`
  ) %>%
  mutate(
    Year = 1803,
    Minute = 0,
    
    ta = ta.orig*1.25,
    
    p.orig = pin + pl / 12 + pt/120,
    p = round(convert_pressure(p.orig, f=25.4, lat=lat, alt=alt, atb=ta),1),
    
    dd.norm = dd_normalize(WD),
    dd = dd2deg(dd.norm),
    
    meta.time = paste0("orig.time=", Hour),
    meta.ta = paste0(meta.time, " | orig.ta=", ta.orig, "R"),
    meta.p  = paste0(meta.time, " | orig.p=", pin,"in", pl+pt/10,"l | atb=", ta, "C"),
    meta.dd = paste0(meta.time, " | orig.dd=", WD, " | force=", WF)
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
    outfile = outfile.name(name, var, df, TRUE, "_JMdeArejula_"),
    outpath = '/scratch3/PALAEO-RA/daily_data/final/Cadiz',
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
    metaHead= ifelse(var=="p", "Observer=Juan Manuel de Aréjula | orig.unit=English inches | PTC=Y | PGC=Y", "Observer=Juan Manuel de Aréjula"),
    meta    = meta,
    time_offset = time.offset(lon),
    keep_na = F
  )
}



# Cadiz_deMolina ----------------------------------------------------------

raw <- read_excel(file, sheet=6, skip=8)
name <- "Cadiz"
lat <- 36 + 31/60
lon <- -(6 + 17/60)
alt <- 13
code <- "CaU1803"
source <- "Ureña, M. 1804. Observaciones meteorológicas hechas en la Isla de León en 1803. Anales de Ciencias Naturales. Biblioteca Universidad Complutense BH Rev 55 (5-11) http://www.books.google.es"

df <- raw %>%
  mutate(
    Year = 1803,
    Minute = minute(Hour),
    Hour = hour(Hour),
    
    ta = round(Tout*1.25,1),
    
    meta.time = paste0("orig.time=", meta_time(Hour, Minute)),
    meta.ta = paste0(meta.time, " | orig.ta=", Tout, "R")
  ) %>%
  select(Year, Month, Day, Hour, Minute, ta, meta.ta)

head(df)

var <- "ta"
write_sef_f(
  as.data.frame(df),
  outfile = outfile.name(name, var, df, TRUE, "_deMolina-deUrena_"),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Cadiz',
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
  metaHead= "Observer=Gaspar de Molina, Marqués de Ureña",
  meta    = df$meta.ta,
  time_offset = time.offset(lon),
  keep_na = F
)


# Cadiz_PerezLasso --------------------------------------------------------

raw <- read_excel(file, sheet=7, skip=8)
name <- "Cadiz"
lat <- 36 + 31/60
lon <- -(6 + 17/60)
alt <- 13
code <- "Ca1820-1822"
metaHead <- "Obsever=Francisco Javier Pérez Lasso | obs.time=midday"
source <- "Pérez-Laso, FJ. (ed.), 1820-1822, 1824. Periódico de la Sociedad Médico-Quirúrgica de Cádiz. Biblioteca Universidad Complutense BH Rev 63-1, 63-2, 63-3, 63-4 (http://www.books.google.com)"

df <- raw %>%
  rename(
    pin = `p(FI/EI)`,
    pl  = `p(l)`,
    rrin = `R(inch)`,
    rrl  = `R(line)`,
    ta.orig = `T(ºR)`
  ) %>%
  mutate(
    Hour = 12,
    Minute = 0,
    
    ta = round(ta.orig*1.25, 1),
    rr = round((as.numeric(rrin) + as.numeric(rrl)/12)*.22, 1),
    
    p_in_total = case_when(
      Year < 1821 ~ as.numeric(pin) + as.numeric(pl)/12,
      Year == 1821 & Month < 10 ~ as.numeric(pin) + as.numeric(pl)/12,
      Year > 1821 | (Year == 1821 & Month >= 10) ~ as.numeric(pin),
      TRUE ~ NA_real_
    ),
    p_mmHg = case_when(
      Year < 1821 ~ p_in_total * 27.07,   # French inches
      Year == 1821 & Month < 10 ~ p_in_total * 27.07,
      Year > 1821 | (Year == 1821 & Month >= 10) ~ p_in_total * 25.4,  # English inches
      TRUE ~ NA_real_
    ),
    p_mmHg = case_when(
      Year < 1821 ~ p_in_total * 27.07,   # French inches
      Year == 1821 & Month < 10 ~ p_in_total * 27.07,
      Year > 1821 | (Year == 1821 & Month >= 10) ~ p_in_total * 25.4,  # English inches
      TRUE ~ NA_real_
    ),
    
    p = round(convert_pressure(p_mmHg, lat=lat, alt=alt, atb=ta),1),
    
    meta.ta = paste0("orig.ta=", ta.orig, "R"),
    meta.rr = paste0("orig.rr=", ifelse(is.na(rr), paste0(as.numeric(rrl), "l"), paste0(rrin, "in",rrl,"l"))),
    meta.p  = case_when(
      Year < 1821 ~ paste0("orig.p=", pin,"French inch", pl, "l | atb=", ta, "C"),   # French inches
      Year == 1821 & Month < 10 ~ paste0("orig.p=", pin," French inch", pl, "l | atb=", ta, "C"),
      Year > 1821 | (Year == 1821 & Month >= 10) ~ paste0("orig.p=", pin,"English in | atb=", ta, "C"),  # English inches
      TRUE ~ NA_character_
    ),
  )

head(df)

base_cols <- df[c("Year","Month","Day","Hour","Minute")]
meta_map  <- list(ta = df$meta.ta, p = df$meta.p, rr = df$meta.rr)

vars <- c("ta", "rr", "p")

for (var in vars) {
  meta  <- meta_map[[var]]
  dat   <- cbind(base_cols, setNames(df[var], var))
  write_sef_f(
    dat,
    outfile = outfile.name(name, var, dat, TRUE, obs_name="_PerezLasso_"),
    outpath = '/scratch3/PALAEO-RA/daily_data/final/Cadiz',
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
    keep_na = F
  )
}


df.dd <- df %>%
  pivot_longer(
    cols = c(Wmorning,Wafternoon, Wnight),
    names_to = "tod",
    values_to = "dd.orig"
  ) %>% mutate(
    meta.time = case_when(
      tod=="Wmorning" ~ "morning",
      tod=="Wafternoon" ~ "afternoon",
      tod=="Wnight" ~ "night"
    ),
    
    Hour = case_when(
      tod=="Wmorning" ~ 8L,
      tod=="Wafternoon" ~ 15L,
      tod=="Wnight" ~ 21L
    ),
    Minute = 0L,
    
    
    dd.norm = dd_normalize(dd.orig),
    dd = dd2deg(dd.norm),
    
    meta   = paste0("obs.time=", meta.time, " | orig.dd=", dd.orig)
  ) %>% select(Year, Month, Day, Hour, Minute, Value=dd, meta)


head(df.dd)

var <- "dd"

write_sef_f(
  as.data.frame(df.dd),
  outfile = outfile.name(name, var, df.dd, TRUE, obs_name="_PerezLasso_"),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Cadiz',
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
  metaHead = metaHead,
  meta    = df.dd$meta,
  keep_na = F
)


# Cadiz_Diario-1816-1830 --------------------------------------------------

raw <- read_excel(file, sheet=9, skip=8)
name <- "Cadiz"
lat <- 36 + 31/60
lon <- -(6 + 17/60)
alt <- 13
code <- "Ca1816-1830"
source <- "Diario Comercial de Cádiz (1802-1830), Cádiz Library, sgn: FL-PP-Est.59. http://www.bibliotecavirtualdeandalucia.es"

df <- raw %>%
  rename(
    ta.orig = `T(ºF)`
  ) %>%
  mutate(
    Minute = 0,
    ta = round((ta.orig-32)/1.8, 1),
    
    dd.norm = dd_normalize(W),
    dd = dd2deg(dd.norm),
    
    meta.time = paste0("orig.time=", Hour),
    meta.ta = paste0(meta.time, " | orig.ta=", ta.orig, "F"),
    meta.dd = paste0(meta.time, " | orig.dd=", W)
  ) %>%
  select(Year, Month, Day, Hour, Minute, ta, dd, meta.ta, meta.dd)

head(df)

meta_map  <- list(ta = df$meta.ta, dd = df$meta.dd)
base_cols <- df[c("Year","Month","Day","Hour","Minute")]
vars <- c("ta", "dd")

for (var in vars) {
  meta  <- meta_map[[var]]
  dat   <- cbind(base_cols, setNames(df[var], var))
  write_sef_f(
    dat,
    outfile = outfile.name(name, var, df, TRUE),
    outpath = '/scratch3/PALAEO-RA/daily_data/final/Cadiz',
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
    meta    = meta,
    time_offset = time.offset(lon),
    keep_na = F
  )
}

