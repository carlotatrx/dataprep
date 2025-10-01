rm(list=ls())
library(lubridate)
library(dplyr)
library(readr)
library(tidyr)
library(dataresqc)
source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')

files <- list.files('/scratch3/PALAEO-RA/daily_data/original/Extremadura',
                    pattern = ".tab$", full.names = TRUE)

link <- 'https://doi.pangaea.de/10.1594/PANGAEA.925805'
source <- 'Vaquero, José Manuel; Bravo-Paredes, Nieves; Obregón, María Angeles; Sánchez-Carrasco, Víctor Manuel; Valente, Maria Antónia; Trigo, Ricardo M; Domínguez-Castro, Fernando; Montero-Martín, Javier; Vaquero-Martínez, Javier; Antón, Manuel; García, José Agustín; Gallego, María Cruz (2020): Early meteorological records from BAALAN1 (1866). PANGAEA.'

time.offset <- function(lon) {as.numeric(lon)*12/180}

outfile.name <-function(name, var, df, subdaily=T) {
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
directions <- c("N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", 
                "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW")
dd2deg <- function(x) 22.5 * (match(toupper(x), directions) - 1)
map2dd <- function(x) {
  x2 <- x %>% trimws()
  x2[x2==""] <- NA_character_
  x2.low <- tolower(x2)  # per si de cas
  out <- x2
  out[!is.na(x2.low) & x2.low=="gallego"] <- "NW"   # gallego: Dicho del viento: Procedente del noroeste, de la parte de Galicia.(RAE)
  out[!is.na(x2.low) & x2.low=="solano"]  <- "W"    # solano: Viento que sopla de donde nace el sol.(RAE)
  toupper(out)
}

# BABADA2 -----------------------------------------------------------------
code <- "BABADA2"
lat  <- 38.900000
lon  <- -6.962500
name <- "Extremadura-Badajoz"

f <- files[grepl(code, files)] # select file from list

raw <- suppressWarnings(read_tsv(
  file      = f,
  skip      = 15,
  col_types = cols(.default = col_guess()),
  trim_ws   = T,
  guess_max = 100000
)) %>%
  rename(
    Date         = 1,
    p.orig       = 2,
    ta.orig      = 3,
    dd.orig      = 4,
    dd.descr     = 5
  ) 

df <- raw %>%
  mutate(
    Year   = year(Date),
    Month  = month(Date),
    Day    = day(Date),
    Hour   = 24L,
    Minute = 0L,
    p.orig  = suppressWarnings(as.numeric(p.orig)),
    ta.orig = suppressWarnings(as.numeric(ta.orig)),
    p  = round(convert_pressure(p.orig, lat=lat, atb=ta.orig),2),
    dd = dd2deg(dd.orig),
    meta.p  = paste0("orig.p=", p.orig, "mmHg | atb=", ta.orig, "C"),
    meta.dd = paste0("orig.dd=", dd.orig, " | force=", dd.descr)
  ) %>% select(Year, Month, Day, Hour, Minute, p, dd, ta.orig, meta.p, meta.dd)

head(df)

value_map <- c(ta = "ta.orig", p = "p", dd = "dd")
meta_map  <- list(ta = "", p = df$meta.p, dd = df$meta.dd)
base_cols <- df[c("Year","Month","Day","Hour","Minute")]


for (var in names(value_map)) {
  value <- value_map[[var]]
  meta  <- meta_map[[var]]
  dat   <- cbind(base_cols, setNames(df[value], value))
  write_sef_f(
    dat,
    outfile = outfile.name(name, var, df, F),
    outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
    cod     = code,
    lat     = lat,
    lon     = lon,
    alt     = NA,
    sou     = source,
    link    = link,
    nam     = name,
    var     = var,
    stat    = "mean",
    period  = "day",
    units   = units(var),
    metaHead = ifelse(var=="p", "PTC=Y | PGC=N", ""),
    meta    = meta,
    keep_na = F
  )
}

# BABADA6 -----------------------------------------------------------------
code <- "BABADA6"
f <- files[grepl(code, files)] # select file from list

raw <- read.delim(
  file          = f,
  skip          = 25,
  sep           = "\t",
  fileEncoding  = "UTF-8",
  check.names   = FALSE,    # keep original column names
  stringsAsFactors = FALSE,
  quote         = ""        # columns contain parentheses/commas but no quotes
) %>%
  rename(
    Date = 1,
    ta.1 = 2,
    ta.2 = 3,
    ta.3 = 4,
    p.zoll.1 = 5,
    p.line.1 = 6,
    p.zoll.2 = 7,
    p.line.2 = 8,
    p.zoll.3 = 9,
    p.zoll.4 = 10,
    dd.orig.1 = 11,
    dd.orig.2 = 13,
    dd.orig.3 = 15
  )

head(raw)


head(raw)

# BAALAN1 & 2 -----------------------------------------------------------------
code <- "BAALAN2" # "BAALAN1" or "BALAAN2"
lat  <- 38.784444
lon  <- -6.243611
name <- "Extremadura-Alange"

f <- files[grepl(code, files)] # select file from list

raw <- suppressWarnings(read_tsv(
  file      = f,
  skip      = 13,
  col_types = cols(.default = col_guess()),
  trim_ws   = T,
  guess_max = 100000
  )) %>%
  rename(
    Date         = 1,
    p.orig       = 2,
    ta.orig      = 3,
    rh           = 4
  )

df <- raw %>%
  mutate(
    Year   = year(Date),
    Month  = month(Date),
    Day    = day(Date),
    Hour   = hour(Date),
    Minute = minute(Date),
    p      = round(convert_pressure(p.orig, lat=lat, atb=ta.orig),2),
    meta.time = paste0("orig.time=", sprintf("%02d", Hour),":", sprintf("%02d",Minute)),
    meta.p = paste0("orig.p=",p.orig,"mmHg | ", meta.time)
  )
head(df)

# save_ta
var <- "ta"
write_sef_f(
  as.data.frame(df[, c("Year", "Month", "Day", "Hour", "Minute", "ta.orig")]),
  outfile = outfile.name(name, var, df, T),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = NA,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "point",
  units   = units(var),
  meta    = df$meta.time,
  metaHead = "Instrument=Thermometer",
  time_offset = time.offset(lon),
  keep_na = F
)

# save p
var <- "p"
write_sef_f(
  as.data.frame(df[, c("Year", "Month", "Day", "Hour", "Minute", "p")]),
  outfile = outfile.name(name, var, df, T),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = NA,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "point",
  units   = units(var),
  meta    = df$meta.p,
  metaHead = "PGC=N | PTC=Y",
  time_offset = time.offset(lon),
  keep_na = F
)

# save rh
var <- "rh"
units.rh <- ifelse(code=="BAALAN1", "relative units", "%")
write_sef_f(
  as.data.frame(df[, c("Year", "Month", "Day", "Hour", "Minute", "rh")]),
  outfile = outfile.name(name, var, df, T),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = NA,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "point",
  units   = units.rh,
  meta    = df$meta.time,
  metaHead = "Instrument=Saussure hygrometer",
  time_offset = time.offset(lon),
  keep_na = F
)


# BAALAN3 -----------------------------------------------------------------
code <- "BAALAN3"

f <- files[grepl(code, files)] # select file from list

raw <- suppressWarnings(read_tsv(
  file      = f,
  skip      = 13,
  col_types = cols(.default = col_guess()),
  trim_ws   = T,
  guess_max = 100000
)) %>%
  rename(
    Date         = 1,
    ta.orig      = 2,
    p.orig       = 3
  )


head(raw)
df <- raw %>%
  mutate(
    Year   = year(Date),
    Month  = month(Date),
    Day    = day(Date),
    Hour   = NA,
    Minute = NA,
    ta     = ta.orig*1.25,
    p.num  = suppressWarnings(as.numeric(p.orig)),
    p      = round(convert_pressure(p.num, lat=lat, atb=ta),2),
    meta.p = paste0("orig.p=",p.orig,"mmHg"),
    meta.ta = paste0("orig.ta=", ta.orig,"R")
  )
head(df)

# save_ta
var <- "ta"
write_sef_f(
  as.data.frame(df[, c("Year", "Month", "Day", "Hour", "Minute", "ta")]),
  outfile = outfile.name(name, var, df, F),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = NA,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "point",
  units   = units(var),
  meta    = df$meta.ta,
  metaHead = "Instrument=Reaumur thermometer",
  time_offset = time.offset(lon),
  keep_na = F
)

# save_p
var <- "p"
write_sef_f(
  as.data.frame(df[, c("Year", "Month", "Day", "Hour", "Minute", "p")]),
  outfile = outfile.name(name, var, df, F),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = NA,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "point",
  units   = units(var),
  meta    = df$meta.p,
  metaHead = "PTC=Y | PGC=N",
  time_offset = time.offset(lon),
  keep_na = F
)


# BAALAN4 -----------------------------------------------------------------
code <- "BAALAN4"

f <- files[grepl(code, files)] # select file from list

raw <- read.delim(
  file          = f,
  skip          = 20,
  sep           = "\t",
  fileEncoding  = "UTF-8",
  check.names   = FALSE,    # keep original column names
  stringsAsFactors = FALSE,
  quote         = ""        # columns contain parentheses/commas but no quotes
  ) %>%
  rename(
    Date         = 1,
    ta.09        = 2,
    ta.14        = 3,
    dd.descr.09  = 4,
    dd.descr.14  = 5,
    dd.descr.mean = 6,
    dd.09        = 7,
    dd.14        = 8,
    dd.mean      = 9,
    p.orig.09    = 10,
    p.orig.14    = 11
  ) %>%
  mutate(
    Year   = year(Date),
    Month  = month(Date),
    Day    = day(Date)
  )

head(raw)

############ wind #############

df.dd <- raw %>%
  pivot_longer(
    cols = c(dd.descr.09, dd.descr.14, dd.descr.mean),
    names_to = "tod",
    names_prefix = "dd\\.descr\\.",
    values_to = "dd.orig"
  ) %>% mutate(
    Hour = case_when(
      tod=="09" ~ 9L,
      tod=="14" ~ 14L,
      tod=="mean" ~ 24L
    ),
    Minute = 0,
    dd.clean = map2dd(dd.orig),
    Value  = dd2deg(dd.clean),
    meta   = ifelse(is.na(Value),"", paste0("orig.dd=", dd.orig))
  ) %>% select(Year, Month, Day, Hour, Minute, Value, meta, dd.orig)

head(df.dd)

df.dd.daily <- df.dd %>%
  filter(Hour==24L, !is.na(Value)) 

df.dd.subdaily <- df.dd %>%
  filter(Hour != 24L, !is.na(Value)) %>%
  mutate(
    meta = paste0(meta, " | orig.time", sprintf("%02d", Hour),":", sprintf("%02d",Minute))
  )

# save_dd_daily
var <- "dd"
write_sef_f(
  as.data.frame(df.dd.daily[, c("Year", "Month", "Day", "Hour", "Minute", "Value")]),
  outfile = outfile.name(name, var, df.dd.daily, F),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = NA,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "mean",
  period  = "day",
  units   = "deg",
  meta    = df.dd.daily$meta,
  keep_na = F
)

write_sef_f(
  as.data.frame(df.dd.subdaily[, c("Year", "Month", "Day", "Hour", "Minute", "Value")]),
  outfile = outfile.name(name, var, df.dd.subdaily, T),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = NA,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "point",
  units   = "deg",
  meta    = df.dd.subdaily$meta,
  time_offset = time.offset(lon),
  keep_na = F
)

############# ta #######################

df.ta <- raw %>%
  pivot_longer(
    cols = c(ta.09, ta.14),
    names_to = "tod",
    names_prefix = "ta\\.",
    values_to = "ta.orig"
  ) %>% mutate(
    Hour = case_when(
      tod=="09" ~ 9L,
      tod=="14" ~ 14L
    ),
    Minute = 0,
    Value  = ta.orig*1.25,
    meta   = paste0("orig.ta=", ta.orig,"R | orig.time=", sprintf("%02d", Hour),":", sprintf("%02d",Minute))
  ) %>% select(Year, Month, Day, Hour, Minute, Value, meta)
head(df.ta)

var <- "ta"
write_sef_f(
  as.data.frame(df.ta[, c("Year", "Month", "Day", "Hour", "Minute", "Value")]),
  outfile = outfile.name(name, var, df.ta, T),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = NA,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "point",
  units   = "C",
  meta    = df.ta$meta,
  time_offset = time.offset(lon),
  metaHead = "Instrument=Reaumur thermometer",
  keep_na = F
)


############# p #######################

df.p <- raw %>%
  pivot_longer(
    cols = c(p.orig.09, p.orig.14),
    names_to = "tod",
    names_prefix = "p\\.orig\\.",
    values_to = "p.orig"
  ) %>% mutate(
    Hour = case_when(
      tod=="09" ~ 9L,
      tod=="14" ~ 14L
    ),
    Minute = 0,
    Value  = round(convert_pressure(p.orig, lat=lat, atb=df.ta$Value),2),
    meta   = paste0("orig.p=", p.orig,"mmHg | orig.time=", sprintf("%02d", Hour),":", sprintf("%02d",Minute), " | atb=", df.ta$Value, "C")
  ) %>% select(Year, Month, Day, Hour, Minute, Value, p.orig, meta)
head(df.p)

var <- "p"
write_sef_f(
  as.data.frame(df.p[, c("Year", "Month", "Day", "Hour", "Minute", "Value")]),
  outfile = outfile.name(name, var, df.p, T),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = NA,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "point",
  units   = "hPa",
  meta    = df.p$meta,
  time_offset = time.offset(lon),
  metaHead = "PTC=Y | PGC=N",
  keep_na = F
)

