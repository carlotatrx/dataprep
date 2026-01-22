rm(list=ls())
library(dataresqc)
library(readxl)
library(tidyr)
library(lubridate)
library(dplyr)
library(stringr)
source('/scratch2/ccorbella/code/dataprep/helpfun.R')

outdir <- "/scratch3/PALAEO-RA/daily_data/final/Vilnius"

name <- "Przybylak_Vilnius"
code <- "Vilnius"
lat <- 54.69
lon <- 25.28
alt	<- 104
source <- "Wahlen E., 1887, Wahre Tagesmittel und Tagliche Variation der Temperatur, an 18 Stationen des Russischen Reiches, Dritter Supplementband zum Repertorium für Meteorologie, Kaiserliche Academie der Wissenschaften, St. Petersburg."
link <- ""
metaHead <- "coords=approx"

file <- "/scratch3/PALAEO-RA/daily_data/original/Vilnius/Vilniu_daily_air_temperature_1777-1882.xls"


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


