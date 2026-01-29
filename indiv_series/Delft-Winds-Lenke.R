rm(list=ls())
# Delft wind directions (Lenke Table 4) -> SEF dd (degrees)
# ---------------------------------------------------------
# Packages
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
source('/scratch2/ccorbella/code/dataprep/helpfun.R')

outdir <- "/scratch3/PALAEO-RA/daily_data/final/Delft"
infile <- "/scratch3/PALAEO-RA/daily_data/original/Berlin/Lenke_Tab4_1708_1709.xlsx"

name <- "Lenke_Delft"
code <- "Delft"
lat	<- 43.2965
lon	<- 5.36978
alt	<- 44
source <- "Untersuchung der ältesten Temperaturmessungen mit Hilfe des strengen Winters 1708-1709, Walter Lenke"

# --- helpers ---
month_de_to_int <- function(x) {
  m <- c("Januar","Februar","März","April","Mai","Juni","Juli","August",
         "September","Oktober","November","Dezember")
  match(x, m)
}

# --- read sheet without assuming fixed header rows ---
raw <- read_excel(infile, sheet = 1, skip=4, col_names = TRUE)



# let's do this dd2deg bc AI suggested it and it will be more reliable
# than the gazillion names of my function
dd_to_deg <- function(dd) {
  dd <- toupper(dd)
  map <- c(
    "N"=0, "NNE"=22.5, "NE"=45, "ENE"=67.5, "E"=90, "ESE"=112.5, "SE"=135,
    "SSE"=157.5, "S"=180, "SSW"=202.5, "SW"=225, "WSW"=247.5, "W"=270,
    "WNW"=292.5, "NW"=315, "NNW"=337.5
  )
  unname(map[dd])
}


df <- raw %>%
  filter(!is.na(...1), !is.na(...2), ...2 == "Delft") %>%
  mutate(
    ...1 = stringr::str_trim(...1),
    Month = month_de_to_int(...1),
    Year  = if_else(Month >= 10, 1708L, 1709L)
  ) %>%
  tidyr::pivot_longer(
    cols = `1`:`31`,
    names_to = "Day",
    values_to = "orig_dd"
  ) %>%
  mutate(
    Day = as.integer(Day),
    # orig_dd = stringr::str_trim(as.character(orig_dd)),
    # orig_dd = dplyr::na_if(orig_dd, "NA"),
    # orig_dd = dplyr::na_if(orig_dd, ""),
    # dd_txt  = stringr::str_extract(toupper(orig_dd),
    #                                "NNW|WNW|WSW|SSW|NNE|ENE|ESE|SSE|NW|NE|SW|SE|N|E|S|W|C"),
    Value = dd_to_deg(orig_dd),
    Meta  = dplyr::case_when(
      is.na(orig_dd) ~ NA_character_,
      TRUE           ~ paste0("orig=", orig_dd)
    ),
    Hour = NA_integer_,
    Minute = NA_integer_,
  ) %>%
  select(Year, Month, Day, Hour, Minute, Value, Meta)

# save
var <-"dd"
write_sef_f(
  as.data.frame(df[,c("Year","Month", "Day", "Hour", "Minute", "Value"),]),
  outfile = outfile.name(name, var, df, FALSE),
  outpath = outdir,
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = "",
  nam     = name,
  var     = var,
  stat    = "point",
  period  = "day",
  units   = units(var),
  keep_na = TRUE,
  meta    = df$Meta
)

