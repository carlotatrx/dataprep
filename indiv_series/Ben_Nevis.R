rm(list=ls())
library(dataresqc)
library(readxl)
library(dplyr)
source('/scratch2/ccorbella/code/dataprep/helpfun.R')

outdir <- "/scratch3/PALAEO-RA/daily_data/final/Ben_Nevis"

name <- "Weather-Rescue-Data-v2-2_Ben-Nevis"
code <- "BEN NEVIS SUMMIT OBSERVATORY"
lat <- 56.8000
lon <- -5.0000
alt	<- 1345
source <- "Ed Hawkins. (2024). ed-hawkins/weather-rescue-data: Weather Rescue Data v2.2 (v2.2). Zenodo."
link <- "https://doi.org/10.5281/zenodo.11001125"

p_to_slp <- function(p_hpa, z_m, T_C, lapse_K_per_m = 0.0065) {
  g0 <- 9.80665
  Rd <- 287.05
  T0 <- T_C + 273.15
  # mean temperature of the layer (simple lapse-rate approximation)
  T_layer <- T0 - 0.5 * lapse_K_per_m * z_m
  p_hpa * exp((g0 * z_m) / (Rd * T_layer))
}


# for the qcs

indir <- '/scratch3/PALAEO-RA/daily_data/final/'
dirname <- "Ben_Nevis/"

# hourly-file -------------------------------------------------------------

file <- "/scratch3/PALAEO-RA/daily_data/original/Ben_Nevis/ben_nevis_summit_hourly_v2.csv"

raw <- read.delim(file,
                  skip =4,
                  sep=",",
                  na="-9999") %>%
  rename(
    p.orig=5,
    rr=6,
    ta=7,
    tawet=8,
    w=9,
    rh=10
  ) %>%
  mutate(
    Minute=0,
    p = round(p_to_slp(p.orig, z_m = alt, T_C = ta), 2),
    meta_time=paste0("orig.time=",Hour),
    meta_p = paste0(meta_time, " | orig=", p.orig, "mbar | atb=",ta, "C")
  )

head(raw)


for (var in c('rh', 'w')){#r', 'ta', 'p', 'w', 'rh')) {
  meta_col <- if (var == "p") "meta_p" else "meta_time"
  
  df <- raw %>%
    select(Year, Month, Day, Hour, Minute, all_of(var), all_of(meta_col)) %>%
    rename(Value = all_of(var), Meta = all_of(meta_col)) %>%
    drop_na(Value)
  
  filename <- outfile.name(name, var, df, TRUE)
  filename_full <- paste0(filename, ".tsv")
  
  write_sef_f(
    df,
    outfile = filename,
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
    units   = ifelse(var=="w", "Ben Nevis scale", units(var)),
    meta    = df$Meta,
    metaHead = ifelse(var=="p", "PTC=Y | PGC=Y", ""),
    keep_na = FALSE
  )
  
  
  qc(glue(indir,dirname,filename_full), outpath=glue(indir, dirname))
  qcfilename <- paste0("qc_BEN NEVIS SUMMIT OBSERVATORY_",var,"_subdaily.txt")
  
  write_flags_f(infile=glue(indir,dirname,filename_full),
                qcfile=glue(indir,dirname, qcfilename),
                outpath=glue(indir,dirname),
                match=FALSE)
}

# daily-file --------------------------------------------------------------

file <- "/scratch3/PALAEO-RA/daily_data/original/Ben_Nevis/ben_nevis_summit_daily_v2.csv"

raw <- read.delim(file,
                  skip =4,
                  sep=",",
                  na="-9999") %>%
  rename(
    rr=4,
    Tn=5,
    Tx=6
  ) %>%
  mutate(
    Hour=NA_integer_,
    Minute=NA_integer_
  )


for (var in c('rr', 'Tn', 'Tx')) {
  df <- as.data.frame(raw[,c("Year","Month", "Day", "Hour", "Minute", var)]) %>% drop_na(any_of(var))
  
  filename <- outfile.name(name, var, df, FALSE)
  filename_full <- paste0(filename, ".tsv")
  
  write_sef_f(
    df,
    outfile = filename,
    outpath = outdir,
    cod     = code,
    lat     = lat,
    lon     = lon,
    alt     = alt,
    sou     = source,
    link    = link,
    nam     = name,
    var     = var,
    stat    = case_when(
      var=="rr" ~ "sum",
      var=="Tn" ~ "minimum",
      var=="Tx" ~ "maximum"
    ),
    period  = "day",
    units   = units(var),
    keep_na = FALSE
  )
  
  qc(glue(indir,dirname,filename_full), outpath=glue(indir, dirname))
  extradaily <- ifelse(var=='rr', "subdaily", "daily")
  qcfilename <- paste0("qc_BEN NEVIS SUMMIT OBSERVATORY_",var,"_", extradaily,".txt")
  
  write_flags_f(infile=glue(indir,dirname,filename_full),
                qcfile=glue(indir,dirname, qcfilename),
                outpath=glue(indir,dirname),
                match=FALSE)
}







