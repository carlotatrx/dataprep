rm(list=ls())
library(dataresqc)
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(lubridate)
source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')

file <- '/scratch3/PALAEO-RA/daily_data/original/Europe_T3_DE_Halle_1700_subdaily.xlsx'
raw <- read_excel(file, sheet = 1, skip  = 5)

lat <- NA
lon <- NA
ele <- NA

# Keep only the first 10 columns and rename
df <- raw %>%
  select(1:8) %>%
  rename(
    Year = 1, Month = 2, Day = 3, Hour.orig = 4, dd.orig = 5,
    ta.orig = 6, p.orig = 7, notes = 8
  ) %>%
  # Fill down Year/Month/Day
  mutate(across(where(is.character), ~ na_if(trimws(.x), ""))) %>%
  fill(Year, Month, Day, .direction = "down") %>%
  mutate(
    Year = as.integer(Year),
    Month = as.integer(Month),
    Day = as.integer(Day),
    BaseDate = as.Date(sprintf("%04d-%02d-%02d", Year, Month, Day)),
    Hour_raw = Hour.orig,
    Minute = 0,
    Hour.orig = suppressWarnings(as.integer(Hour.orig)),
    p.orig = suppressWarnings(as.numeric(p.orig)),
    ta.orig = suppressWarnings(as.numeric(ta.orig))
  )

adjust_hours <- function(df) {
  n <- nrow(df)
  stopifnot(n!=0)
  
  date_corr <- df$BaseDate
  Hour24 <- df$Hour.orig
  
  prev_date <- date_corr[1]
  prev_hour <- Hour24[1]
  
  # if (!is.na(prev_hour)) prev_hour <- prev_hour %% 24 # in case
  
  for (i in seq_len(n)) {
    h  <- df$Hour.orig[i]
    bd <- df$BaseDate[i]
    
    # if hour is missing, move forward base date (e.g. next day's header started)
    # and advance prev_date
    if (is.na(h)) {
      if (!is.na(bd) && bd > prev_date) {
        prev_date <- bd
        prev_hour <- -1L
      }
      date_corr[i] <- prev_date
      Hour24[i] <- NA_integer_ # set Hour to NA for this row
      next  # skip the rest of the loop
    }
    
    # if the base date jumped to later day, reset our day context
    if (!is.na(bd) && bd > prev_date) {
      prev_date <- bd
      prev_hour <- -1L # reset so we can have small hours at the start of the day
    }
    
    h1 <- h  # corrected hour
    
    # if new hour is smaller than last one
    if (prev_hour >= 0L && h1 < prev_hour) {
      
      # convert 12a.m. to 24:00h or switch to p.m. time
      pm_candidate <- if(h==12L) 24L else h + 12L
      
      if (pm_candidate >= prev_hour && pm_candidate <= 24L) {
        # accept PM interpretation on the SAME day
        h1 <- pm_candidate
      } else {
        # or advance one day if it's larger than 25
        prev_date <- prev_date + 1L
        h1 <- h
      }
    }
    
    date_corr[i] <- prev_date
    Hour24[i]    <- h1
    prev_hour    <- ifelse(is.na(h1), prev_hour, h1)
    
  }
  
  minute_vector <- integer(n)
  
  # if there's two consecutive hours in the same day, the second is the half
  for (i in 2:n) {
    same_day  <- !is.na(date_corr[i]) && !is.na(date_corr[i-1]) && date_corr[i] == date_corr[i-1]
    same_hour <- !is.na(Hour24[i]) && !is.na(Hour24[i-1]) && Hour24[i] == Hour24[i-1]
    if (same_day && same_hour) {
      cat(date_corr[i], "has", Hour24[i], "h repeated.")
      minute_vector[i] <- ifelse(minute_vector[i-1]==30L, 0, 30L)
    }
    if (!is.na(Hour24[i]) && Hour24[i] == 24L) minute_vector[i] <- 0L
  }
  
  out <- df %>%
    mutate(date_corrected = date_corr,
           Hour24 = Hour24,
           Minute = minute_vector,
           Year = year(date_corr),
           Month = month(date_corr),
           Day = day(date_corr)
           )
  out
}

df.adj <- adjust_hours(df)

head(df.adj)

## Define wind directions
dirs <- c("N", "NNO", "NO", "ONO", "O", "OSO", "SO", "SSO", "S", 
          "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW")

library(stringr)
# Build regex: longest-first, match as a standalone token (no letters around)
dir.ordered <- dirs[order(nchar(dirs), decreasing = TRUE)]
dir.pattern <- paste0("(?<![A-Z])(", paste(dir.ordered, collapse="|"), ")(?![A-Z])")

df.all <- df.adj %>%
  mutate(
    dd.orig = na_if(dd.orig, "--"),
    dd.clean = str_extract(toupper(trimws(dd.orig)),
                            regex(dir.pattern, ignore_case=F)),
    dd = case_when(
      is.na(dd.clean) ~ NA_real_,
      TRUE ~ round(22.5 * (match(dd.clean, dirs) - 1), 0)
    ),
    
    meta.dd = if_else(!is.na(dd.orig), paste0("orig.dd=", dd.orig), NA_character_),
    meta.p  = if_else(!is.na(p.orig), paste0("orig.p=", p.orig), NA_character_),
    meta.ta = if_else(!is.na(ta.orig), paste0("orig.ta=", ta.orig), NA_character_),
    
    # --- CHANGES ONLY when hour/date actually changed --- #
    meta.hour = if_else(
      suppressWarnings(as.integer(Hour.orig)) != suppressWarnings(as.integer(Hour24)),
      paste0("orig.hour=",Hour.orig),
      NA_character_
    ),
    
    meta.date = if_else(
      !is.na(date_corrected) & !is.na(BaseDate) & date_corrected != BaseDate,
      paste0("orig.date=", BaseDate),
      NA_character_
    )
    
  )

# put the metas properly
# Helper to paste only non-empty, non-NA parts with " | "
combine_meta <- function(...) {
  x <- c(...)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x)) paste(x, collapse = " | ") else NA_character_
}

df.all <- df.all %>%
  mutate(
    meta.ta.full = pmap_chr(list(meta.ta, meta.date, meta.hour), combine_meta),
    meta.p.full  = pmap_chr(list(meta.p,  meta.date, meta.hour), combine_meta),
    meta.dd.full = pmap_chr(list(meta.dd, meta.date, meta.hour), combine_meta)
  )

head(df.all)

# save

# ta
write_sef_f(
  Data=as.data.frame(df.all[, c("Year", "Month", "Day", "Hour24", "Minute", "ta.orig")]),
  outfile="Halle_ta_subdaily.tsv",
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod="Halle",
  lat=lat, lon=lon, alt=ele,
  variable="ta",
  nam="Halle",
  meta=df.all$meta.ta.full,
  units="unknown", sou="Ephemeris Barometrico Meteorologica", stat="point",keep_na = F
)

# p
write_sef_f(
  Data=as.data.frame(df.all[, c("Year", "Month", "Day", "Hour24", "Minute", "p.orig")]),
  outfile="Halle_p_subdaily.tsv",
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod="Halle",
  variable="p",
  nam="Halle",
  lat=lat, lon=lon, alt=ele,
  meta=df.all$meta.p.full,
  metaHead = "PGC=N | PTC=N",
  units="unknown", sou="Ephemeris Barometrico Meteorologica", stat="point", keep_na = F
)

# dd
write_sef_f(
  Data=as.data.frame(df.all[, c("Year", "Month", "Day", "Hour24", "Minute", "dd")]),
  outfile="Halle_dd_subdaily.tsv",
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod="Halle",
  variable="dd",
  lat=lat, lon=lon, alt=ele,
  nam="Halle",
  meta=df.all$meta.dd.full,
  units="deg", sou="Ephemeris Barometrico Meteorologica", stat="point", keep_na = F
)
