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

outfile.name <- function(name, var, df, subdaily=T) {
  subdaily.str <- ifelse(subdaily, "subdaily", "daily")
  paste0(name,"_",get_date_range(df), "_", var,"_", subdaily.str)
}

units <- function(var) {
  case_when(
    var == "ta" ~ "C",
    var == "p"  ~ "hPa",
    var == "dd" ~ "deg",
    var == "rr" ~ "mm",
    T ~ "unknown"
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
  # common typos seen
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

# Breda -------------------------------------------------------------------

f <- "/scratch3/PALAEO-RA/daily_data/original/Breda/Breda_1726-1740.txt"
lat <- 51.57
lon <- 4.77
alt <- 3
code <- "KNMI-56_Breda"
name <- "Breda"

raw <- read.csv(
  file      = f,
  skip      = 55,
  header = FALSE,
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
  rhorig = 9
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
  ta = round((taorig/10-32)/1.8,2),
  p = round(convert_pressure(as.numeric(porig)/12, f=26.2, lat=lat, alt=alt, atb=ta),2),
  ddnorm = dd_normalize_nl(ddorig),
  dd  = dd2deg(ddnorm),
  meta.p = paste0("orig.p=",porig/10, "inch(duim)", porig %%10, "l | atb=", ta,"C"),
  meta.dd = paste0("orig.dd=", ddorig)
) %>% select(Year, Month, Day, Hour, Minute, dd, ta, p, meta.p, meta.dd)


base_cols <- df[c("Year","Month","Day","Hour","Minute")]
meta.list <- list(dd=df$meta.dd, p=df$meta.p)

for (var in c("dd", "p")){
  
  dat   <- cbind(base_cols, setNames(df[var], var))
  write_sef_f(
    dat,
    outfile = outfile.name(name, var, dat, TRUE),
    outpath = '/scratch3/PALAEO-RA/daily_data/final/Breda',
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
    keep_na = FALSE
  )
}

# Breda2 -------------------------------------------------------------------

f <- "/scratch3/PALAEO-RA/daily_data/original/Breda/Breda2_1778-1781.txt"

raw <- read.csv(
  file      = f,
  skip      = 55,
  header = FALSE,
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
  ddnorm = dd_normalize_nl(ddorig),
  stopifnot(all(ddnorm[!is.na(ddnorm)] %in% c(directions, "calm"))), # safety check
  dd  = dd2deg(ddnorm),
  meta.dd = paste0("orig.dd=", ddorig),
  
  # precip
  rr = round(rrorig/12,1),
  meta.rr = paste0("orig.rr=", rrorig, "l")
) %>% select(Year, Month, Day, Hour, Minute, dd,  meta.dd, rr, meta.rr)

base_cols <- df[c("Year","Month","Day","Hour","Minute")]
meta.list <- list(dd=df$meta.dd, rr=df$meta.rr)

vars <- c("dd", "rr")
for (var in vars) {
  write_sef_f(
    dat <- cbind(base_cols, setNames(df[var], var)),
    outfile = outfile.name(name, var, dat, TRUE),
    outpath = '/scratch3/PALAEO-RA/daily_data/final/Breda',
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



