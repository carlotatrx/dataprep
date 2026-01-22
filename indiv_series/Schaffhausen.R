
rm(list=ls())
suppressWarnings(suppressMessages({
  library(readr)
  library(dplyr)
  library(stringr)
}))
source('/scratch2/ccorbella/code/dataprep/helpfun.R')


# --------- user inputs ----------
infile  <- "/scratch3/PALAEO-RA/daily_data/original/Schaffhausen/SH01_Schaffhausen.tsv"
outdir <- "/scratch3/PALAEO-RA/daily_data/final/Schaffhausen/"

code <- "Schaffhausen"
name  <- "CHIMES_Schaffhausen"
lat <- 47.696
lon <- 8.639
alt <- 400
source <- "CHIMES"
link   <- "https://doi.org/10.5194/cp-15-1345-2019"

# --------------------------------

# 1) read tsv
raw <- readr::read_tsv(infile, show_col_types = FALSE, na = c("NA", "", "NaN"))

# -------------------------------------------------------------------------


# 3) keep only the requested columns
keep_cols <- c("year","month","day","hour","minutes","wind_from_direction")


df <- raw %>%
  select(all_of(keep_cols)) %>%
  mutate(
    year    = as.integer(year),
    month   = as.integer(month),
    day     = as.integer(day),
    hour    = as.integer(hour),
    minutes = as.integer(minutes),
    
    dd.norm = dd_normalize(wind_from_direction),
    dd = dd2deg(dd.norm),
    meta = paste0("orig.time=", hour, ":00 | orig=", wind_from_direction)
  )

df

# save
var <-"dd"
write_sef_f(
  as.data.frame(df[,c("year","month", "day", "hour", "minutes", "dd"),]),
  outfile = outfile.name(name, var, df, TRUE),
  outpath = outdir,
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "point",
  meta    = df$meta,
  time_offset = time.offset(lon),
  units   = units(var),
  metaHead = "Observer=Johann Christoph Schalch",
  keep_na = TRUE
)
