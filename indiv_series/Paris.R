rm(list=ls())
library(dataresqc)
library(readr)
library(dplyr)
library(lubridate)
library(tidyr)
library(purrr)
library(stringr)
source('/home/ccorbella/scratch2_symboliclink/code/dataprep/helpfun.R')


# Pressure_Lucas ----------------------------------------------------------
name <- "Paris"
code <-"PAR"

outdir <- paste0("/scratch3/PALAEO-RA/daily_data/final/", name)

raw <- read.csv(paste0("/scratch3/PALAEO-RA/daily_data/original/", name, "/", code,"_SLP.csv"), na=c("","NA"))

head(raw)

lat <- 48.817
lon <- 2.322
alt <- 77

df <- raw %>%
  mutate(
    year=year(dates),
    month=month(dates),
    day=day(dates),
    Hour=NA,
    minute=NA,
  )

df <- df %>%
  filter(!is.na(SLP))
head(df)


var<-"p"
write_sef_f(
  as.data.frame(df[c("year","month", "day", "Hour", "minute","SLP")]),
  outfile = outfile.name(paste0("Pfister_",name), var, df, FALSE),
  outpath = outdir,
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = "Pfister, L., Wilhelm, L., Brugnara, Y., Imfeld, N., & Brönnimann, S. (2024). Weather type reconstruction using machine learning approaches. EGUsphere, 2024, 1-33.",
  link    = "https://doi.org/10.5194/wcd-6-571-2025",
  nam     = name,
  var     = var,
  stat    = "mean",
  period    = "day",
  units   = units(var),
  metaHead = "PTC=Y | PGC=Y | homogenized=Y | QC=Y",
)



# Pressure_Cornes ---------------------------------------------------------

file <- '/scratch3/PALAEO-RA/daily_data/original/Paris/paris_MSLP_1670_2007.txt'

raw <- suppressWarnings(read.table(
  file = file,
  skip = 16,
  #col_types = cols(.default = col_guess()),
  col.names = c("Year", "Month", "Day", "p.orig", "p.hom", "diff", "code"),
  header = FALSE,
  na.strings = NA
))

head(raw)

df <- raw %>%
  mutate(
    Hour = NA_integer_,
    Minute = NA_integer_,
  ) %>% select(Year, Month, Day, Hour, Minute, p.orig)

write_sef_f(as.data.frame(df),
            outfile=paste0('Paris_Cornes_', get_date_range(df), "_p_daily.tsv"),
            outpath="/scratch3/PALAEO-RA/daily_data/final/Paris/",
            variable = "p",
            cod = "Paris_Morin",
            nam = "Paris",
            lat = 48.86,
            lon = 2.337,
            alt = 67,
            sou = "Cornes, R. C., Jones, P. D., Briffa, K. R. and Osborn, T. J. (2012) A daily series of mean sea-level pressure for Paris, 1670-2007. International Journal of Climatology 32, 1135-1150. doi: 10.1002/joc.2349",
            link = "https://crudata.uea.ac.uk/cru/data/parislondon/",
            units = "hPa",
            period = "day",
            stat = "mean",
            metaHead = "PTC=Y | PGC=Y | a homogenized version exists, see `raw` file",
            keep_na = FALSE)

# Precipitation_Pliemon ---------------------------------------------------

file <- '/scratch3/PALAEO-RA/daily_data/original/Paris/Supplement_dataset_Pliemon_etal_CP_2023.txt'

raw <- suppressWarnings(read_delim(
  file = file,
  delim = '\t',
  skip = 47,
  col_types = cols(.default = col_guess()),
  trim_ws = T,
  guess_max = 100000
))

mk_meta <- function(v1, v2, v3, v4, v5, v6){
  vals <- list(v1, v2, v3, v4, v5, v6)
  notes <- purrr::map_chr(seq_along(vals), function(i){
    x <- vals[[i]]
    if (is.na(x) || (is.character(x) && trimws(x) == "")) return(NA_character_)
    sx <- as.character(x)
    
    # Special codes:
    if (sx %in% c("0", 0))   return(paste0("obs", i, ":RI=0(light_rainfall)"))
    if (sx %in% c("-1", -1)) return(paste0("obs", i, ":RD=-1(snow)"))
    
    # Default: two digits = intensity(d1), duration(d2)
    digits <- stringr::str_extract_all(sx, "\\d")[[1]]
    if (length(digits) < 2) return(NA_character_)
    paste0("obs", i, ":RI=", digits[1], ",RD=", digits[2])
  })
  notes <- notes[!is.na(notes)]
  paste(notes, collapse = " | ")
}

raw <- raw %>%
  mutate(
    meta = pmap_chr(
      list(Prec_raw_1, Prec_raw_2, Prec_raw_3, Prec_raw_4, Prec_raw_5, Prec_raw_6),
      mk_meta
    )
  )

df <- raw %>%
  transmute(
    Year = year(Date),
    Month = month(Date),
    Day = day(Date),
    Hour = 24L,
    Minute = 0L,
    Value = Prec,
    meta = meta
  ) %>% drop_na(Value)

head(df)
write_sef_f(as.data.frame(df),
            outfile=paste0('Paris_', get_date_range(df), "_rr_daily.tsv"),
            outpath="/scratch3/PALAEO-RA/daily_data/final/Paris/",
            variable = "rr",
            cod = "Paris_Morin",
            nam = "Paris",
            lat = 48.861347,
            lon = 2.350561,
            alt = NA,
            sou = "Pliemon, T., Foelsche, U., Rohr, C., & Pfister, C. (2023). Precipitation reconstructions for Paris based on the observations by Louis Morin, 1665–1713 CE [Data set]. Zenodo.",
            link = "https://doi.org/10.5281/zenodo.7404635",
            units = "mm",
            stat = "sum",
            period = "day",
            metaHead = "loc: Rue Quinquempoix (1665-1685), Hôtel de Rohan-Soubise (1685-1688), Abbey of Saint-Victor (1688-1715) | legend: RI=rain intensity, RD=rain duration",
            meta = df$meta,
            keep_na = F)
