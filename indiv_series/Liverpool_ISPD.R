rm(list=ls())
library(dplyr)
source('/scratch2/ccorbella/code/dataprep/helpfun.R')

outdir <- "/scratch3/PALAEO-RA/daily_data/SEF/Liverpool"

lat  <- 53.4
lon  <- -2.99
alt  <- 20
name <- "Liverpool"
code <- "004002"
source_name <- "ISPD"
link <- "DOI: 10.1002/joc.1335"

raw <- read.csv(
  "/scratch3/PALAEO-RA/daily_data/original/ISPD/ISPD_004002_1768_1793.csv",
  colClasses = "character",
  na.strings = c("", "NA"),
  stringsAsFactors = FALSE
)

df <- raw %>%
  mutate(
    Year   = as.integer(year),
    Month  = as.integer(month),
    Day    = as.integer(day),
    Hour   = NA_integer_,
    Minute = NA_integer_,
    p      = round(as.numeric(surface_pressure_hpa), 2),
    orig_units = recode(trimws(original_surface_pressure_units), "inchHg" = "inHg"),
    meta.p = paste0("orig.p=", trimws(original_surface_pressure), orig_units)
  ) %>%
  filter(!is.na(p)) %>%
  arrange(Year, Month, Day) %>%
  select(Year, Month, Day, Hour, Minute, p, meta.p)

var <- "p"
write_sef_f(
  as.data.frame(df[, c("Year","Month","Day","Hour","Minute","p")]),
  outfile  = paste0("ISPD_Liverpool_", get_date_range(df), "_p_daily"),
  outpath  = outdir,
  cod      = code,
  lat      = lat,
  lon      = lon,
  alt      = alt,
  sou      = source_name,
  link     = link,
  nam      = name,
  variable = var,
  stat     = "mean",
  units    = units(var),
  period   = "day",
  metaHead = "PGC=N | PTC=N | QC=Y",
  meta     = df$meta.p,
  keep_na  = FALSE
)
