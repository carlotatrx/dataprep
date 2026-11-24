rm(list = ls())
library(lubridate)
library(dplyr)
library(readr)
library(tidyr)
library(dataresqc)
library(purrr)
source("/scratch2/ccorbella/code/dataprep/helpfun.R")

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

# check that array is sorted by date
date_sorted_check <- function(x, nm = "Date") stopifnot(all(x[[nm]] == sort(x[[nm]])))



#############################################################################
# Senguerdius (Leiden) ------------------------------------------------------
#############################################################################

name <- "Leiden"
observer <- "Observer=Wolferd Senguerd"
lat <- 52.1601144
lon <- 4.49700970
code <- paste0("KNMI-", name)

files <- "/scratch3/PALAEO-RA/daily_data/original/Leiden/senguerdius.dat"

raw <- read.csv(
  file      = files,
  skip      = 33,
  header    = FALSE,
  fill = TRUE,
  strip.white = TRUE,
  stringsAsFactors = FALSE
) %>% rename(
  station = 1,
  Date = 2,
  porigA  = 3,
  porigB  = 4,
  taorigA = 5,
  taorigB = 6,
  taorigC = 7,
  taorigD = 8,
  ddorig  = 9,
  worig   = 10,
  # pA = 12,
  # pB = 13,
  taC = 14,
  taB = 15,
  pA  = 16,
  pB  = 17,
  p   = 20
)
  
  
head(raw)

raw <- raw[!duplicated(raw), ]
raw <- raw[!is.na(raw$Date), ]
raw <- raw[order(raw$Date), ]
stopifnot(all(raw$Date==sort(raw$Date)))

df <- raw %>% mutate(
  # fix non-existing 29th Feb 1697
  feb29 = Date==16970229,
  Date = if_else(feb29, 16970228, Date),
  Date = ymd(Date),
  Year = year(Date),
  Month = month(Date),
  Day = day(Date),
  Hour = NA_integer_,
  Minute = NA_integer_,
  
  # wind
  # check that all wind directions are okay
  dd_norm = dd_normalize_nl(ddorig),
  stopifnot(all(dd_norm[!is.na(dd_norm)] %in% c(directions, "calm"))),
  dd  = dd2deg(dd_norm),
  
  ta  = taC/10,
  p = p/10,

  meta.dd = paste0("orig.dd=", ddorig, " | force=", worig),
  meta.p  = paste0("orig.p=", porigA/100, "in", porigA%%100, "l Rijnlandse | atb=", taB/10, "C"),
  
  meta.ta = if_else(feb29, "orig.date=1697-02-29", ""),
  meta.dd = if_else(feb29, paste0(meta.dd, " | orig.date=1697-02-29"), meta.dd),
  meta.p  = if_else(feb29, paste0(meta.p,  " | orig.date=1697-02-29"), meta.p)
) %>% select(Year, Month, Day, Hour, Minute, ta, dd, p, meta.p, meta.dd, meta.ta)

head(df)


base_cols <- df[c("Year","Month","Day","Hour","Minute")]
meta.list <- list(dd=df$meta.dd,
                  p=df$meta.p, ta=df$meta.ta)
metaHead <- list(dd=observer,
                 p=paste0(observer, " | PTC=Y | PGC=Y"),
                 ta=paste0(observer, " | Instrument=liquid thermometer"))
vars <- c("dd", "ta", "p")

for (var in vars) {
  dat <- cbind(base_cols, setNames(df[var], var))
  write_sef_f(
    dat,
    outfile = outfile.name(name, var, dat, FALSE),
    outpath = paste0('/scratch3/PALAEO-RA/daily_data/SEF/', name),
    cod     = code,
    lat     = lat,
    lon     = lon,
    alt     = NA,
    sou     = source,
    link    = link,
    nam     = name,
    var     = var,
    stat    = "point",
    period  = "day",
    units   = units(var),
    meta    = meta.list[[var]],
    metaHead = metaHead[[var]],
    keep_na = FALSE
  )
}


#############################################################################
# Bergen -------------------------------------------------------------------
#############################################################################

name <- "Bergen"
code <- paste0("KNMI-", name)
lat	 <- 52.67
lon	 <- 4.72
alt	 <- 1
files <- list.files(paste0("/scratch3/PALAEO-RA/daily_data/original/", name), full.names = TRUE)

raw <- read.csv(
  file      = files,
  skip      = 51,
  header = TRUE,
  fill = TRUE,
  strip.white = TRUE,
  stringsAsFactors = FALSE
) %>% rename(
  station = 1,
  Date    = 2,
  obsnum = 3,
  porig  = 4,
  taorig = 5,
  ddorig = 6,
  fh     = 7,
  ww     = 8,
  rrorig = 9,
  w2     = 10
)

stopifnot(all(raw$Date==sort(raw$Date)))

df <- raw %>% mutate(
  Date = ymd(Date),
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
  
  ta  = round((taorig/10-32)*(5/9),1),

  rr  = round(rrorig*2.2, 1),
  
  meta.ta = paste0(meta.time, " | orig.ta=", taorig/10, "F"),
  meta.dd = paste0(meta.time, " | orig.dd=", ddorig),
  meta.rr = paste0(meta.time, " | orig.rr=", rrorig, "l"),

) %>% select(Year, Month, Day, Hour, Minute, ta, dd, rr, meta.rr, meta.dd, meta.ta)

head(df)

base_cols <- df[c("Year","Month","Day","Hour","Minute")]
meta.list <- list(dd=df$meta.dd,
                  rr=df$meta.rr, ta=df$meta.ta)
vars <- c("dd", "ta", "rr")

for (var in vars) {
  dat <- cbind(base_cols, setNames(df[var], var))
  write_sef_f(
    dat,
    outfile = outfile.name(name, var, dat, TRUE),
    outpath = paste0('/scratch3/PALAEO-RA/daily_data/SEF/', name),
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

#############################################################################
# Alkmaar -------------------------------------------------------------------
#############################################################################

name <- "Alkmaar"
code <- paste0("KNMI-", name)
lat	 <- 52.63
lon	 <- 4.75
alt	 <- 1
files <- list.files(paste0("/scratch3/PALAEO-RA/daily_data/original/", name), full.names = TRUE)

raw <- read.csv(
  file      = files,
  skip      = 51,
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

head(raw)

stopifnot(all(raw$Date==sort(raw$Date)))

df <- raw %>% mutate(
  Date = ymd(Date),
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
  
  # no pressure available
  # wind
  # check that all wind directions are okay
  dd_norm = dd_normalize_nl(ddorig),
  stopifnot(all(dd_norm[!is.na(dd_norm)] %in% c(directions, "calm"))),
  dd  = dd2deg(dd_norm),
  
  ta  = round((taorig/10-32)*(5/9),1),
  
  rr  = round(rrorig*2.2,1),
  
  meta.ta = paste0(meta.time, " | orig.ta=", taorig/10, "F"),
  meta.dd = paste0(meta.time, " | orig.dd=", ddorig),
  meta.rr = paste0(meta.time, " | orig.rr=", rrorig, "l"),

) %>% select(Year, Month, Day, Hour, Minute, ta, dd, rr, meta.rr, meta.dd, meta.ta)

head(df)

base_cols <- df[c("Year","Month","Day","Hour","Minute")]
meta.list <- list(dd=df$meta.dd,
                  rr=df$meta.rr, ta=df$meta.ta)
vars <- c("dd", "ta", "rr")

for (var in vars) {
  dat <- cbind(base_cols, setNames(df[var], var))
  write_sef_f(
    dat,
    outfile = outfile.name(name, var, dat, TRUE),
    outpath = paste0('/scratch3/PALAEO-RA/daily_data/SEF/', name),
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


#############################################################################
# Leiden -------------------------------------------------------------------
#############################################################################

name <- "Leiden"
code <- paste0("KNMI-", name)
lat	 <- 52.16
lon	 <- 4.5
alt	 <- 2
files <- list.files(paste0("/scratch3/PALAEO-RA/daily_data/original/", name), full.names = TRUE)

raw_list <- lapply(files, function(f)
  read.csv(
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

stopifnot(all(raw$Date==sort(raw$Date)))

df <- raw %>% mutate(
  Date = ymd(Date),
  Year = year(Date),
  Month = month(Date),
  Day = day(Date),
  Hour = case_when(
    obsnum==1 ~ 8L,
    obsnum==2 ~ 13L,
    obsnum==3 ~ 22L
  ),
  Minute = 0L,
  meta.time = paste0("obs.num=", obsnum)
)

# test before not passed
df <- df[order(df$Date), ]

# check again
stopifnot(all(df$Date==sort(df$Date)))

df <- df %>% mutate(
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
  
  p   = round(convert_pressure(pinch, f=26.2, lat=lat, alt=alt, atb=ta),1),
  rr  = round(rrorig*2.2,1),
  
  meta.ta = paste0(meta.time, " | orig.ta=", taorig/10, "F"),
  meta.dd = paste0(meta.time, " | orig.dd=", ddorig),
  meta.rr = paste0(meta.time, " | orig.rr=", rrorig, "l"),
  meta.p  = paste0(meta.time, " | orig.p=", ii, ".", ll, ".", q, " Rijnlandse_inch.line.quarter | atb=", ta, "C")
  
) %>% select(Year, Month, Day, Hour, Minute, p, ta, dd, rr, meta.rr, meta.p, meta.dd, meta.ta)

head(df)

base_cols <- df[c("Year","Month","Day","Hour","Minute")]
meta.list <- list(dd=df$meta.dd,
                  rr=df$meta.rr, ta=df$meta.ta,
                  p=df$meta.p)
vars <- c("dd", "p", "ta", "rr")

for (var in vars) {
  dat <- cbind(base_cols, setNames(df[var], var))
  write_sef_f(
    dat,
    outfile = outfile.name(name, var, dat, TRUE),
    outpath = paste0('/scratch3/PALAEO-RA/daily_data/SEF/', name),
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

#############################################################################
# Zwanenburg -------------------------------------------------------------------
#############################################################################

name <- "Zwanenburg"
code <- paste0("KNMI-", name)
lat	 <- 52.32
lon	 <- 4.74
alt	 <- 1
files <- list.files(paste0("/scratch3/PALAEO-RA/daily_data/original/", name),
                    pattern = "^his_42-.*\\.dat$", full.names = TRUE)

raw_list <- lapply(files, function(f)
  read.csv(
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
  Date = ymd(Date),
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
  
  inch_mm = if_else(Date < ymd("1809-05-01"), 25.4, 26.2),
  p   = round(convert_pressure(pinch, f=inch_mm, lat=lat, alt=alt, atb=ta),1),
  rr  = rrorig/10,
  
  meta.ta = paste0(meta.time, " | orig.ta=", taorig/10, "F"),
  meta.dd = paste0(meta.time, " | orig.dd=", ddorig),
  meta.rr = paste0(meta.time, " | orig.rr=", rrorig, "l"),
  meta.p = paste0(meta.time, " | orig.p=", ii, ".", ll, ".", q, " ",
                  if_else(Date < ymd("1809-05-01"), "English_inch.line.quarter",
                          "Rijnlandse_inch.line.quarter"),
                  " | atb=", ta, "C")
  
) %>% select(Year, Month, Day, Hour, Minute, p, ta, dd, rr, meta.rr, meta.p, meta.dd, meta.ta)

head(df)

stopifnot(all(df$Date==sort(df$Date)))

base_cols <- df[c("Year","Month","Day","Hour","Minute")]
meta.list <- list(dd=df$meta.dd,
                  rr=df$meta.rr, ta=df$meta.ta,
                  p=df$meta.p)
vars <- c("dd", "p", "ta", "rr")

for (var in vars) {
  dat <- cbind(base_cols, setNames(df[var], var))
  write_sef_f(
    dat,
    outfile = outfile.name(name, var, dat, TRUE),
    outpath = paste0('/scratch3/PALAEO-RA/daily_data/SEF/', name),
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
    metaHead = ifelse(var == "p",
                      paste0(
                        "PTC=Y | PGC=Y | MSLP=N | ",
                        "units=English inch (1738-01-01--1809-04-30); ",
                        "Rhineland inch (1809-05-01--1860-12-31) | ",
                        "Instrument=Prins barometer with unboiled mercury, English scale (1735--1809); ",
                        "Solari & Butti barometer, Rhineland scale (from 1809-04)"
                      ),
                      ""),
    keep_na = FALSE
  )
}

#############################################################################
# Haarlem -------------------------------------------------------------------
#############################################################################

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

stopifnot(all(raw$Date==sort(raw$Date)))

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

  # scale in use: station 48 (1735-42) Rhineland 12 lines; 1788-11 to 1793-10-22 English12 lines;
  # from 1793-10-23 English with decimal lines (10 per inch)
  scl = case_when(
    Date <  ymd("1788-11-01") ~ "rhine12",
    Date <  ymd("1793-10-23") ~ "eng12",
    TRUE                      ~ "eng10"),
  pinch = case_when(          
    scl == "eng10" & ll >= 10 ~ NA_real_,          # impossible with decimal lines: digitization error
    scl == "eng10"            ~ ii + ll/10 + q/40,
    TRUE                      ~ ii + ll/12 + q/48),
  inch_mm = if_else(scl == "rhine12", 26.2, 25.4),
  
  ta  = round((taorig/10-32)*(5/9),1),
  p   = round(convert_pressure(pinch, f=inch_mm, lat=lat, alt=alt, atb=ta),1),
  rr  = round(rrorig*2.2, 1),
  
  meta.ta = paste0(meta.time, " | orig.ta=", taorig/10, "F"),
  meta.dd = paste0(meta.time, " | orig.dd=", ddorig),
  meta.rr = paste0(meta.time, " | orig.rr=", rrorig, "l"),
  meta.p = paste0(meta.time, " | orig.p=", ii, ".", ll, ".", q, " ",
                  case_when(scl == "rhine12" ~ "Rijnlandse_inch.line.quarter",
                            scl == "eng12"   ~ "English_inch.line.quarter",
                            TRUE             ~ "English_inch.decimalline.quarter"),
                  " | atb=", ta, "C")   
  
) %>% select(Year, Month, Day, Hour, Minute, p, ta, dd, rr, meta.rr, meta.p, meta.dd, meta.ta)

head(df)

base_cols <- df[c("Year","Month","Day","Hour","Minute")]
meta.list <- list(dd=df$meta.dd, p=df$meta.p,
                  rr=df$meta.rr, ta=df$meta.ta)

# check: should be roughly 950-1050hPa
range(df$p, na.rm=TRUE)

vars <- c("dd", "p", "ta", "rr")

for (var in vars) {
  dat <- cbind(base_cols, setNames(df[var], var))
  write_sef_f(
    dat,
    outfile = outfile.name(name, var, dat, TRUE),
    outpath = '/scratch3/PALAEO-RA/daily_data/SEF/Haarlem',
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
    metaHead = ifelse(var=="p",
                      paste0("PTC=Y | PGC=Y | MSLP=N | alt.=2m (1735-1742) | units=Rhineland inch with 12 lines ",
                             "(1735-01-01--1742-12-31); English inch with 12 lines (1788-11-01--1793-10-22); English ",
                             "inch with decimal lines (1793-10-23--1841-12-31) | Instrument=Barometer with English ",
                             "scale (1788-11--1793-10); barometer by Van Marum, ivory English scale with decimal lines (1794-04--1841)"),
                      ""),
    keep_na = FALSE
  )
}



for (var in vars[2]) {
  dat <- cbind(base_cols, setNames(df[var], var))
  write_sef_f(
    dat,
    outfile = outfile.name(name, var, dat, TRUE),
    outpath = '/scratch3/PALAEO-RA/daily_data/HIST-DAILY_v2/fixes/',
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
    metaHead = ifelse(var=="p",
                      paste0("PTC=Y | PGC=Y | MSLP=N | alt.=2m (1735-1742) | units=Rhineland inch with 12 lines ",
                      "(1735-01-01--1742-12-31); English inch with 12 lines (1788-11-01--1793-10-22); English ",
                      "inch with decimal lines (1793-10-23--1841-12-31) | Instrument=Barometer with English ",
                      "scale (1788-11--1793-10); barometer by Van Marum, ivory English scale with decimal lines (1794-04--1841)"),
                      ""),
    keep_na = FALSE
  )
}
#############################################################################
# Utrecht -------------------------------------------------------------------
#############################################################################

files <- list.files("/scratch3/PALAEO-RA/daily_data/original/Utrecht", full.names = TRUE)

f <- "/scratch3/PALAEO-RA/daily_data/original/Utrecht/his_43.dat"
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
  Minute = 0L
)
  
df <- df[order(df$Date), ]

df <- df %>% mutate(
  # wind
  # check that all wind directions are okay
  dd_norm = dd_normalize_nl(ddorig),
  stopifnot(all(dd_norm[!is.na(dd_norm)] %in% c(directions, "calm"))),
  dd  = dd2deg(dd_norm),
  
  # ta
  ta = round((taorig/10 - 32) * 5/9, 1),
  
  # p
  s   = sprintf("%05d", as.integer(porig)),
  ii  = as.numeric(substr(s, 1, 2)),
  ll  = as.numeric(substr(s, 3, 4)),
  q   = as.numeric(substr(s, 5, 5)),
  pinch = ii + ll/12 + q/48,
  p   = round(convert_pressure(pinch, f=26.2, lat=lat, alt=alt, atb=ta),1),
  
  # precip
  rr = round(as.numeric(rrorig)*2.2, 1),
  
  meta.time = paste0("obs.num=", obsnum),
  
  meta.p  = paste0(meta.time, " | orig.p=", ii, ".", ll, ".", q, "Rijnlandse inch.line.quarter | atb=", ta, "C"),
  meta.ta = paste0(meta.time, " | orig.ta=", taorig/10, "F"),
  meta.dd = paste0(meta.time, " | orig.dd=", ddorig),
  meta.rr = paste0(meta.time, " | orig.rr=", rrorig, "l")
  
) %>% select(Year, Month, Day, Hour, Minute, dd, p, ta, rr, meta.dd, meta.rr, meta.p, meta.ta)

base_cols <- df[c("Year","Month","Day","Hour","Minute")]
meta.list <- list(dd=df$meta.dd, rr=df$meta.rr,
                  ta=df$meta.ta, p=df$meta.p)
vars <- c("dd", "p", "ta", "rr")

for (var in vars){
  
  dat   <- cbind(base_cols, setNames(df[var], var))
  write_sef_f(
    dat,
    outfile = outfile.name(name, var, dat, TRUE),
    outpath = paste0("/scratch3/PALAEO-RA/daily_data/SEF/", name),
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

#############################################################################

f <- "/scratch3/PALAEO-RA/daily_data/original/Utrecht/Utrecht_1836-1846.txt"
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
  Minute = 0L
)

df <- df[order(df$Date), ]

df <- df %>% mutate(
  
  ta = case_when(
    Year == 1836                ~ round((taorig/10 - 32) * 5/9, 1),
    Year >= 1837 & Year <= 1846 ~ taorig/10,
    TRUE                        ~ NA_real_   # outside documented range
  ),
  
  porig = as.numeric(porig),
  p   = round(convert_pressure(porig/10, lat=lat, alt=alt, atb=ta),1),
  
  # wind
  # check that all wind directions are okay
  dd_norm = dd_normalize_nl(ddorig),
  stopifnot(all(dd_norm[!is.na(dd_norm)] %in% c(directions, "calm"))),
  dd  = dd2deg(dd_norm),
  
  meta.time = paste0("orig.time=", sprintf("%02d", Hour), ":", sprintf("%02d", Minute)),
  meta.ta  = case_when(
    Year == 1836                ~ paste0(meta.time, " | orig.ta=", taorig/10, "F"),
    Year >= 1837 & Year <= 1846 ~ paste0(meta.time, " | orig.ta=", taorig/10, "C"),
    TRUE                        ~ "orig.ta=unknown_scale"
  ),
  meta.dd = paste0(meta.time, " | orig.dd=", ddorig),
  meta.p  = paste0(meta.time, " | orig.p=", porig/10, "mmHg | atb=", ta, "C")
  
) %>% select(Year, Month, Day, Hour, Minute, dd, p, ta, meta.dd, meta.p, meta.ta)

head(df)

vars <- c("dd", "p", "ta")
base_cols <- df[c("Year","Month","Day","Hour","Minute")]
meta.list <- list(dd=df$meta.dd, p=df$meta.p, ta=df$meta.ta)

for (var in vars){
  
  dat   <- cbind(base_cols, setNames(df[var], var))
  write_sef_f(
    dat,
    outfile = outfile.name(name, var, dat, TRUE),
    outpath = paste0('/scratch3/PALAEO-RA/daily_data/SEF/', name),
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
    time_offset = time.offset(lon),
    metaHead = ifelse(var=="p", "PTC=Y | PGC=Y | Observer=Prof. PJI de Fremery", "Observer=Prof. PJI de Fremery"),
    keep_na = FALSE
  )
}

#############################################################################

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
  Minute = 0L
)

df <- df[order(df$Date), ]
  
df <- df %>% mutate(
  # wind
  # check that all wind directions are okay
  dd_norm = dd_normalize_nl(ddorig),
  stopifnot(all(dd_norm[!is.na(dd_norm)] %in% c(directions, "calm"))),
  dd  = dd2deg(dd_norm),
  
  ta  = taorig/10,
  p   = round(convert_pressure(porig/100, lat=lat, alt=alt, atb=ta),1),
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
    outpath = paste0('/scratch3/PALAEO-RA/daily_data/SEF/', name),
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
    outpath = '/scratch3/PALAEO-RA/daily_data/SEF/Utrecht',
    cod     = code,
    lat     = lat,
    lon     = lon,
    alt     = alt,
    sou     = source,
    link    = link,
    nam     = name,
    var     = var,
    period  = "day",
    stat    = substr(var, 2, 4),
    units   = "C",
    keep_na = FALSE
  )
}

#############################################################################
# Breda -------------------------------------------------------------------
#############################################################################

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
  ta = round((taorig/10-32)/1.8,1),
  p = round(convert_pressure(as.numeric(porig)/12, f=26.2, lat=lat, alt=alt, atb=ta),1),
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
    outpath = '/scratch3/PALAEO-RA/daily_data/SEF/Breda',
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

#############################################################################
# Breda2 -------------------------------------------------------------------
#############################################################################

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

# order check
stopifnot(all(raw$Date==sort(raw$Date)))

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
  
  # precip
  rr = round(rrorig*2.2,1),
  
  meta.time = paste0("obs.num=", obsnum),
  meta.dd = paste0(meta.time, " | orig.dd=", ddorig),
  meta.rr = paste0(meta.time, " | orig.rr=", rrorig, "l")
) %>% select(Year, Month, Day, Hour, Minute, dd,  meta.dd, rr, meta.rr)

head(df)

base_cols <- df[c("Year","Month","Day","Hour","Minute")]
meta.list <- list(dd=df$meta.dd, rr=df$meta.rr)

vars <- c("dd", "rr")
for (var in vars) {
  write_sef_f(
    dat <- cbind(base_cols, setNames(df[var], var)),
    outfile = outfile.name(name, var, dat, TRUE),
    outpath = '/scratch3/PALAEO-RA/daily_data/SEF/Breda',
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



