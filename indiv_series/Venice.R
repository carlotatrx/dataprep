rm(list=ls())
library(dataresqc)
library(readxl)
library(tidyr)
library(dplyr)
library(lubridate)
library(tibble)
source('/scratch2/ccorbella/code/dataprep/helpfun.R')

outdir <- "/scratch3/PALAEO-RA/daily_data/final/Venice"
name <- "Venice"
code <- "Venezia"
lat	<- 45.44
lon	<- 12.33
alt	<- NA
source <- " Camuffo, D., della Valle, A., & Becherini, F. (2025). Temperature and Pressure Observations by Tommaso Temanza from 1751 to 1769 in Venice, Italy. Climate, 13(10), 217."
link   <- " https://doi.org/10.3390/cli13100217 "
metaHead.p <- "Observer=Tommaso Temanza | Instrumet=Amontons thermometer (1751-1759), Réaumur thermometer (1761-1769) | exact location unknown"
metaHead.ta <- "Observer=Tommaso Temanza | PGC=N | exact location unknown"


# pressure ----------------------------------------------------------------


file <- "/scratch3/PALAEO-RA/daily_data/original/Venice/Venice_pressure_by_Temanza.csv"
raw <- read.csv(file)

colnames(raw)[4] <- "p"

head(raw)

df <- raw %>%
  add_column(Hour=12L, Minute=0L, .before="p") %>%
  mutate(

    meta=paste0(
    "obs.time=solar_noon",
    ifelse(outlier==0,""," | qc=outlier")
    )
  )

head(df)

var <- "p"

write_sef_f(
  as.data.frame(df),
  outfile  = outfile.name(code, var, df, subdaily="noon"),
  outpath  = outdir,
  cod      = code,
  lat      = lat,
  lon      = lon,
  alt      = alt,
  sou      = source,
  link     = link,
  nam      = name,
  var      = var,
  stat     = "point",
  units    = units(var),
  metaHead = metaHead.p,
  meta     = df$meta,
  time_offset = time.offset(lon),
  keep_na  = FALSE
)


# temperature -------------------------------------------------------------

file <- "/scratch3/PALAEO-RA/daily_data/original/Venice/Venice_temperature_by_Temanza.csv"
raw <- read.csv(file)

colnames(raw)[4] <- "ta"

head(raw)

df <- raw %>%
  add_column(Hour=12, Minute=0, .before="ta") %>%
  mutate(meta="obs.time=solar_noon")

head(df)

var <- "ta"

write_sef_f(
  as.data.frame(df),
  outfile  = outfile.name(code, var, df, subdaily="noon"),
  outpath  = outdir,
  cod      = code,
  lat      = lat,
  lon      = lon,
  alt      = alt,
  sou      = source,
  link     = link,
  nam      = name,
  var      = var,
  stat     = "point",
  units    = "C",
  metaHead = metaHead.ta,
  meta     = df$meta,
  keep_na  = FALSE
)
