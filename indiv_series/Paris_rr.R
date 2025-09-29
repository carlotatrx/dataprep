rm(list=ls())
library(dataresqc)
library(readr)
library(dplyr)
library(lubridate)
library(tidyr)
library(purrr)
library(stringr)
source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')

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
