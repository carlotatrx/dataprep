rm(list=ls())
library(dataresqc)
library(readxl)
library(dplyr)
library(stringr)
source('/scratch2/ccorbella/code/dataprep/helpfun.R')


# Warsaw-from-Rajmund -----------------------------------------------------

outdir <- "/scratch3/PALAEO-RA/daily_data/final/Warsaw"

name <- "Przybylak_Warsaw"
code <- "Warsaw"
lat <- 52.23
lon <- 21
alt	<- 124
source <- "Wahlen E., 1887, Wahre Tagesmittel und Tagliche Variation der Temperatur, an 18 Stationen des Russischen Reiches, Dritter Supplementband zum Repertorium für Meteorologie, Kaiserliche Academie der Wissenschaften, St. Petersburg."
link <- ""
metaHead <- "coords=approx"

file <- "/scratch3/PALAEO-RA/daily_data/original/Warsaw/Warsaw_daily_air_temperature_1761-1882.xls"


# Read all sheets
sheets <- excel_sheets(file)

read_one_month_raw <- function(sheet) {
  raw <- read_excel(
    file, sheet = sheet,
    skip = 1,                 # skip title line above "Days"
    col_names = TRUE,
    col_types = "text",       # IMPORTANT: consistent types across sheets
    na = "-",                 # "-" -> NA
    .name_repair = "unique",
    trim_ws = TRUE
  )
  
  # Standardize first column name
  names(raw)[1] <- "Year"
  
  # DROP rows like "Mean"
  raw <- raw %>%
    filter(str_detect(Year, "^\\d{4}$"))
  
  # Keep only Year + day columns 1..31 (some sheets can have extra junk columns)
  wanted_days <- as.character(1:31)
  keep <- c("Year", intersect(names(raw), wanted_days))
  raw <- raw[, keep, drop = FALSE]
  
  # Ensure missing day columns exist (e.g., if a sheet only goes to 30)
  missing_days <- setdiff(wanted_days, names(raw))
  if (length(missing_days) > 0) {
    raw[missing_days] <- NA_character_
  }
  
  # Order columns Year, 1..31 and add Month label
  raw <- raw %>%
    select(Year, all_of(wanted_days)) %>%
    mutate(Month = sheet, .before = Year)
  
  raw
}

raw_all <- bind_rows(lapply(sheets, read_one_month_raw))
raw <- raw_all %>%
  mutate(
    Year = as.integer(Year),
    across(all_of(as.character(1:31)), ~ suppressWarnings(as.numeric(.x)))
  )

# sanity checks
# all years should be valid
stopifnot(nrow(raw %>% filter(is.na(Year)))==0)

# Month name -> month number (handles "January", "May", etc.)
month_to_num <- function(x) match(tolower(x), tolower(month.name))


df <- raw %>%
  pivot_longer(
    cols = all_of(as.character(1:31)),
    names_to = "Day",
    values_to = "Value"
  ) %>%
  mutate(
    Day   = as.integer(Day),
    Month = month_to_num(Month),
    Date  = make_date(Year, Month, Day)
  ) %>%
  mutate(
    Value = round(Value, 1)
  ) %>%
  select(Date, Year, Month, Day, Value) %>%
  arrange(Year, Month, Day)

# check false dates
false_dates <- df %>%
  filter(is.na(Date)) %>%
  select(Year, Month, Day) %>%
  distinct() %>%
  arrange(Year, Month, Day)
false_dates

false_dates_with_values <- df %>%
  filter(is.na(Date) & !is.na(Value)) %>%
  select(Year, Month, Day, Value) %>%
  arrange(Year, Month, Day)

# sanity check
# there shouldn't be any false date with a non-NA value
stopifnot(nrow(false_dates_with_values)==0)

# once all is good, drop non-valid Dates
df2 <- df %>%
  filter(!is.na(Date)) %>%
  mutate(Hour=NA_integer_, Minute=NA_integer_) %>%
  select(!Date)
View(df2)


# save
var <-"ta"
write_sef_f(
  as.data.frame(df2[,c("Year","Month", "Day", "Hour", "Minute", "Value"),]),
  outfile = outfile.name(name, var, df2, FALSE),
  outpath = outdir,
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "unknown",
  period  = "day",
  units   = units(var),
  metaHead = metaHead,
  keep_na = TRUE
)


# Warsaw-from-PALAEO-RA ---------------------------------------------------

file <- '/scratch3/PALAEO-RA/daily_data/original/Europe_T3_PL_Warsaw_1725-1728_subdaily.xlsx'
raw <- read_excel(file, sheet = 1, skip  = 6)

lat	<- 52.2296756
lon	<- 21.0122287
alt	<- NA
name <- "Warsaw"
metaHead <- "Observer=Erndtel"


time.offset <- as.numeric(lon)*12/180

###### Wind direction ########
# Convert Warsaw-style wind codes (O=East) to daily degrees
# Strategy: split on '.'/' ' separators, map each token to degrees, circular mean per row.

# Base 16-point (with O=East) + 32-point extensions + a few forgiving shortcuts
dir_deg <- c(
  "N"= 0, "NNO"= 22.5, "NO"= 45, "ONO"= 67.5, "O"= 90, "OSO"=112.5,
  "SO"=135, "SSO"=157.5, "S"=180, "SSW"=202.5, "SW"=225, "WSW"=247.5,
  "W"=270, "WNW"=292.5, "NW"=315, "NNW"=337.5,
  # Forgiving shorthands occasionally seen in manuscripts:
  # 'WN' ~ WNW (292.5), 'WS' ~ WSW (247.5), 'NO' and 'SO' already in table.
  "WN"=292.5, "WS"=247.5
)

# helper: map a single token to degrees
map_token <- function(tok){
  if (is.na(tok)) return(NA_real_)
  t0 <- toupper(tok)
  t0 <- gsub("[^A-Z]", "", t0, perl=T)  # keep letters only
  if (!nzchar(t0)) return(NA_real_)     # if t0 is empty
  
  # Direct hit
  hit <- unname(dir_deg[t0])
  if (!is.na(hit)) return(hit)
  
  # Try collapsing repeated letters like 'NNOO' -> 'NNO' (rare OCR oddity)
  t2 <- gsub("(.)\\1+", "\\1", t0, perl=T)
  hit <- unname(dir_deg[t2])
  if (!is.na(hit)) return(hit)
  
  # else return NA  
  NA_real_
}

# circular mean in degrees (0-360)
circ_mean_deg <- function(deg){
  deg <- deg[!is.na(deg)]
  if (!length(deg)) return(NA_real_)
  rad <- deg * pi/180
  mX <- mean(cos(rad)); mY <- mean(sin(rad))
  out <- (atan2(mY, mX) * 180/pi) %% 360
  out
}

# Main function: vectorized over a character vector of codes
wind_to_dd <- function(x){
  unknown <- character(0)
  out <- vapply(x, function(s){
    if (is.na(s)) return(NA_real_)
    s2 <- toupper(trimws(as.character(s)))
    if (!nzchar(s2)) return(NA_real_)
    # standardize separators '.' and spaces
    parts <- unlist(strsplit(s2, "[\\.[:space:]]+", perl=T))
    parts <- parts[nzchar(parts)]
    if (!length(parts)) return(NA_real_)
    vals <- vapply(parts, function(p){
      d <- map_token(p)
      if (is.na(d) && nzchar(p)) unknown <<- unique(c(unknown, p))
      d
    }, numeric(1))
    circ_mean_deg(vals)
  }, numeric(1))
  attr(out, "unknown_tokens") <- sort(unique(unknown))
  out
}


# Keep only the first 10 columns and rename
df <- raw %>%
  select(1:8) %>%
  rename(
    Year = 1, Month = 2, Day = 3, p.morning.zoll = 4, p.morning.in = 5, p.afternoon.zoll = 6, p.afternoon.in = 7, dd.orig = 8
  ) %>%
  # Fill down Year/Month/Day
  mutate(across(where(is.character), ~ na_if(trimws(.x), ""))) %>%
  mutate(
    Year = as.integer(Year),
    Month = as.integer(Month),
    Day = as.integer(Day),
    Hour = NA,
    Minute = NA,
    
    dd = round(wind_to_dd(dd.orig),2),
    meta.dd = paste0("orig.dd=",dd.orig)
  )
head(df)

# save
# dd
write_sef_f(
  Data=as.data.frame(df[, c("Year", "Month", "Day", "Hour", "Minute", "dd")]),
  outfile=paste0(name,"_dd_subdaily.tsv"),
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod=name,
  lat=lat, lon=lon, alt=alt,
  variable="dd",
  nam=name,
  meta=df$meta.dd,
  units="deg", sou="PALAEO-RA", stat="point",keep_na = F
)
