rm(list = ls())
library(lubridate)
library(dplyr)
library(readr)
library(tidyr)
library(dataresqc)
source("/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R")

files <- list.files("/scratch3/PALAEO-RA/daily_data/original/Utrecht", full.names = TRUE)

source <- "KNMI"
link   <- "https://www.knmi.nl/nederland-nu/klimatologie/daggegevens/antieke-waarnemingen"

time.offset <- function(lon) {as.numeric(lon)*12/180}

outfile.name <- function(name, var, df, subdaily=TRUE) {
  subdaily.str <- ifelse(subdaily, "subdaily", "daily")
  df.short <- df[!is.na(df[[var]]), , drop=FALSE]
  paste0(name,"_",
         get_date_range(if (nrow(df.short)) df.short else df),
         "_", var,"_", subdaily.str)
}

units <- function(var) {
  case_when(
    var == "ta" ~ "C",
    var == "p"  ~ "hPa",
    var == "dd" ~ "deg",
    var == "rr" ~ "mm",
    var == "rh" ~ "perc",
    TRUE ~ "unknown"
  )
}

# wind conversion
directions <- c("N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
                "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW")

dd_normalize <- function(x) {
  y <- toupper(trimws(x))
  y <- recode(y,
              "NWE"  = "WNW",
              "NWW"  = "WNW",
              "SWW"  = "WSW",
              "WENW" = "WNW",
              "ENW"  = "ENE",
              "NEE"  = "ENE",
              "NN0"  = "NNW",   # OCR '0' -> 'W'
              "NNO"  = "NNW",   # 'O' (Oeste) -> W
              "SWS"  = "SSW",
              "ES"   = "SE",
              "SES"  = "SE",
              "NSW"  = "SSW",
              "C"    = "calm",
              .default = y,
              .missing = NA_character_
  )
  ifelse(y %in% c(directions, "calm"), y, NA_character_)
}
dd2deg <- function(x) {
  deg <- rep(NA_character_, length(x))
  i_dir <- !is.na(x) & x %in% directions
  i_calm <- !is.na(x) & x == "calm"
  deg[i_dir] <- as.character(22.5 * (match(x[i_dir], directions) - 1))
  deg[i_calm] <- "calm"
  deg
}

dir_deg <- setNames(seq(0, 337.5, by=22.5), directions)
# base Dutch -> English mapping (no 'T')
dutch_map <- c(
  "N"="N", "NO"="NE", "NNO"="NNE",
  "ONO"="ENE", "O"="E", "OZO"="ESE",
  "ZO"="SE", "ZZO"="SSE", "Z"="S",
  "ZZW"="SSW", "ZW"="SW", "WZW"="WSW",
  "W"="W", "WNW"="WNW", "NW"="NW", "NNW"="NNW",
  # Iberian variants
  "O"="E", "NO"="NE", "SO"="SE", "Z0"="S", "N0"="N",
  # common typos
  "WZW"="WSW", "WZW"="WSW"
)

# helper: English 16-pt to degrees and back
to_deg <- function(code) dir_deg[[code]]
nearest_dir <- function(deg) names(dir_deg)[which.min(abs(((deg - unname(dir_deg) + 180) %% 360) - 180))]

# Normalize one token (Dutch → English 16-pt or "calm"/NA), handling "T"
norm_dutch_token <- function(tok) {
  if (is.na(tok) || !nzchar(tok)) return(NA_character_)
  y <- toupper(trimws(tok))
  if (y %in% c("C","CALM","V")) return("calm")   # C=Calm (Calmo/Calma), V=Variable → keep as calm label
  # direct map first
  if (y %in% names(dutch_map)) return(dutch_map[[y]])
  # if contains 'T' (range like "ZWTZ" = SW to S, "NOTN" = NE to N)
  if (grepl("T", y)) {
    parts <- strsplit(y, "T", fixed = TRUE)[[1]]
    if (length(parts) != 2) return(NA_character_)
    left  <- dutch_map[toupper(parts[1])]
    right <- dutch_map[toupper(parts[2])]
    if (any(is.na(c(left, right)))) return(NA_character_)
    la <- to_deg(left); ra <- to_deg(right)
    # angular midpoint (shortest arc)
    diff <- ((ra - la + 540) %% 360) - 180
    mid  <- (la + diff/2) %% 360
    return(nearest_dir(mid))
  }
  # other single Dutch combos like "ONO","OZO","ZZO","ZZW","NNO" covered above;
  # handle simple two-letter Dutch not in map by translating letters and retry
  y2 <- chartr("OZ", "ES", y) # O->E, Z->S
  if (y2 %in% directions) return(y2)
  NA_character_
}

# Vectorized normalizer
dd_normalize_nl <- function(x) vapply(x, norm_dutch_token, character(1))




# Haarlem -------------------------------------------------------------------
files <- list.files("/scratch3/PALAEO-RA/daily_data/original/Haarlem", full.names = TRUE)
files <- files[c(4, 1, 2, 3)] # rearrange so first one is year 1735
code <- "KNMI-Haarlem"
lat	 <- 52.38
lon	 <- 4.64
alt	 <- 5
name <- "Haarlem"


raw_list <- lapply(files, function(f)
  read.csv(
    file      = f,
    skip      = 54,
    header = TRUE,
    fill = TRUE,
    strip.white = TRUE,
    stringsAsFactors = FALSE
  ) %>% rename(
    station = 1,
    Date = 2,
    obsnum = 3,
    porig  = 4,
    taorig = 5,
    ddorig = 6,
    fh     = 7,
    ww     = 8,
    rrorig = 9,
    w2     = 10
  )
)

raw <- bind_rows(raw_list)
head(raw)

df <- raw %>% mutate(
  Date= ymd(Date),
  Year = year(Date),
  Month = month(Date),
  Day = day(Date),
  Hour = case_when(
    obsnum==1 ~ 8L,
    obsnum==2 ~ 13L,
    obsnum==3 ~ 22L
  ),
  Minute = 0L,
  meta.time = paste0("obs.num=", obsnum),
  
  # wind
  # check that all wind directions are okay
  dd_norm = dd_normalize_nl(ddorig),
  stopifnot(all(dd_norm[!is.na(dd_norm)] %in% c(directions, "calm"))),
  dd  = dd2deg(dd_norm),
  
  s   = sprintf("%05d", as.integer(porig)),
  ii  = as.numeric(substr(s, 1, 2)),
  ll  = as.numeric(substr(s, 3, 4)),
  q   = as.numeric(substr(s, 5, 5)),
  pinch = ii + ll/12 + q/48,
  ta  = round((taorig/10-32)*(5/9),1),
  p   = round(convert_pressure(pinch, f=26.2, lat=lat, alt=alt, atb=ta),2),
  rr  = rrorig/10,
  
  meta.ta = paste0("orig.ta=", taorig, "F | ", meta.time),
  meta.dd = paste0("orig.dd=", ddorig, " | ", meta.time),
  meta.p  = paste0("orig.p=", ii, ".", ll, ".", q, "Rijnlandse inch.line.quarter | atb=", ta, "C | ", meta.time)
  
) %>% select(Year, Month, Day, Hour, Minute, p, ta, dd, rr, meta.time, meta.p, meta.dd, meta.ta)

head(df)

base_cols <- df[c("Year","Month","Day","Hour","Minute")]
meta.list <- list(dd=df$meta.dd,
                  rr=df$meta.time, ta=df$meta.ta,
                  p=df$meta.p)
vars <- c("dd", "p", "ta", "rr")

for (var in vars) {
  dat <- cbind(base_cols, setNames(df[var], var))
  write_sef_f(
    dat,
    outfile = outfile.name(name, var, dat, TRUE),
    outpath = '/scratch3/PALAEO-RA/daily_data/final/Haarlem',
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
    meta    = meta.list[[var]],
    metaHead = ifelse(var=="p", "PTC=Y | PGC=Y | alt.1735-1742=2m", "alt.1735-1742=2m"),
    keep_na = FALSE
  )
}



# Utrecht -------------------------------------------------------------------
files <- list.files("/scratch3/PALAEO-RA/daily_data/original/Utrecht", full.names = TRUE)

f <- files[5]
code <- "KNMI-43_Utrecht"
lat	<- 52.09
lon	<- 5.12
alt	<- 5
name <- "Utrecht"

raw <- read.csv(
  file      = f,
  skip      = 51,
  header = TRUE,
  fill = TRUE,
  strip.white = TRUE,
  stringsAsFactors = FALSE
) %>% rename(
  station = 1,
  Date = 2,
  obsnum = 3,
  porig = 4,
  taorig = 5,
  ddorig = 6,
  rrorig = 9
)


df <- raw %>% mutate(
  Date= ymd(Date),
  Year = year(Date),
  Month = month(Date),
  Day = day(Date),
  Hour = case_when(
    obsnum==1 ~ 8L,
    obsnum==2 ~ 13L,
    obsnum==3 ~ 22L
  ),
  Minute = 0L,
  
  # wind
  # check that all wind directions are okay
  dd_norm = dd_normalize_nl(ddorig),
  stopifnot(all(dd_norm[!is.na(dd_norm)] %in% c(directions, "calm"))),
  dd  = dd2deg(dd_norm),
  meta.dd = paste0("orig.dd=", ddorig),
  
  # precip
  rr = round(as.numeric(rrorig)/12,1),
  meta.rr = paste0("orig.rr=", rrorig, "l")
) %>% select(Year, Month, Day, Hour, Minute, dd,  meta.dd, rr, meta.rr)

base_cols <- df[c("Year","Month","Day","Hour","Minute")]
meta.list <- list(dd=df$meta.dd, rr=df$meta.rr)

for (var in c("dd", "rr")){
  
  dat   <- cbind(base_cols, setNames(df[var], var))
  write_sef_f(
    dat,
    outfile = outfile.name(name, var, dat, TRUE),
    outpath = '/scratch3/PALAEO-RA/daily_data/final/Utrecht',
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
    meta    = meta.list[[var]],
    keep_na = FALSE
  )
}

f <- files[6]
raw <- read.csv(
  file      = f,
  skip      = 18,
  header = TRUE,
  fill = TRUE,
  strip.white = TRUE,
  stringsAsFactors = FALSE
) %>% rename(
  Date = 1,
  hh = 2,
  porig = 3,
  taorig = 4,
  ddorig = 5,
)

df <- raw %>% mutate(
  Date= ymd(Date),
  Year = year(Date),
  Month = month(Date),
  Day = day(Date),
  Hour = as.numeric(hh)/100,
  Minute = 0L,
  
  # wind
  # check that all wind directions are okay
  dd_norm = dd_normalize_nl(ddorig),
  stopifnot(all(dd_norm[!is.na(dd_norm)] %in% c(directions, "calm"))),
  dd  = dd2deg(dd_norm),
  meta.dd = paste0("orig.dd=", ddorig),
  
) %>% select(Year, Month, Day, Hour, Minute, dd,  meta.dd)

var <- "dd"
write_sef_f(
  df,
  outfile = outfile.name(name, var, df, TRUE),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Utrecht',
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
  keep_na = FALSE
)

code <- "KNMI-155_Utrecht"

raw_list <- lapply(files[2:4], function(f)
  read.csv(
    file      = f,
    skip      = 25,
    header = TRUE,
    fill = TRUE,
    strip.white = TRUE,
    stringsAsFactors = FALSE
  ) %>% rename(
    station = 1,
    Date = 2,
    hhmm = 3,
    ddorig = 4,
    taorig = 8,
    tmin   = 9,
    tmax   = 10,
    rrorig = 11,
    porig  = 12,
    rh     = 14
  )
)

raw <- bind_rows(raw_list)
head(raw)

df <- raw %>% mutate(
  Date= ymd(Date),
  Year = year(Date),
  Month = month(Date),
  Day = day(Date),
  Hour = as.numeric(hhmm)/100,
  Minute = 0L,
  
  # wind
  # check that all wind directions are okay
  dd_norm = dd_normalize_nl(ddorig),
  stopifnot(all(dd_norm[!is.na(dd_norm)] %in% c(directions, "calm"))),
  dd  = dd2deg(dd_norm),
  
  ta  = taorig/10,
  p   = round(convert_pressure(porig/100, lat=lat, alt=alt, atb=ta),2),
  rr  = rrorig/10,
  meta.time = paste0("orig.time=", sprintf("%02d", Hour), ":00"),
  
  meta.dd = paste0("orig.dd=", ddorig, " | ", meta.time),
  meta.p  = paste0("orig.p=", porig/100, "mmHg | atb=", ta, "C | ", meta.time)
  
) %>% select(Year, Month, Day, Hour, Minute,p, ta, tmin, tmax, dd, rr, rh, meta.time, meta.p, meta.dd)

head(df)

base_cols <- df[c("Year","Month","Day","Hour","Minute")]
meta.list <- list(dd=df$meta.dd, rh=df$meta.time,
                  rr=df$meta.time, ta=df$meta.time, # it's the same meta
                  p=df$meta.p, rr=df$meta.time)
vars <- c("dd", "rh", "p", "ta", "rr")

for (var in vars) {
  dat <- cbind(base_cols, setNames(df[var], var))
  write_sef_f(
    dat,
    outfile = outfile.name(name, var, dat, TRUE),
    outpath = '/scratch3/PALAEO-RA/daily_data/final/Utrecht',
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
    meta    = meta.list[[var]],
    metaHead = ifelse(var=="p", "PTC=Y | PGC=Y", ""),
    time_offset = time.offset(lon),
    keep_na = FALSE
  )
}

df.daily <- raw %>% mutate(
  Date= ymd(Date),
  Year = year(Date),
  Month = month(Date),
  Day = day(Date),
  Hour = 24L,
  Minute = 0L,
  
  tmax = tmax/10,
  tmin = tmin/10
  
) %>% select(Year, Month, Day, Hour, Minute, tmax, tmin)

for (var in c("tmax", "tmin")) {
  dat <- df.daily[c("Year","Month","Day","Hour","Minute", var)]
  write_sef_f(
    dat,
    outfile = outfile.name(name, var, dat, FALSE),
    outpath = '/scratch3/PALAEO-RA/daily_data/final/Utrecht',
    cod     = code,
    lat     = lat,
    lon     = lon,
    alt     = alt,
    sou     = source,
    link    = link,
    nam    