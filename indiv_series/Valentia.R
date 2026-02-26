rm(list=ls())
library(dataresqc)
library(readxl)
library(tidyr)
library(dplyr)

source('/scratch2/ccorbella/code/dataprep/helpfun.R')

indir <- "/scratch3/PALAEO-RA/daily_data/original/"
name <- "Valentia"

outdir <- paste0("/scratch3/PALAEO-RA/daily_data/final/",name)

source <- "Mateus, C., Potito, A., & Curley, M. (2020). Reconstruction of a long‐term historical daily maximum and minimum air temperature network dataset for Ireland (1831‐1968). Geoscience Data Journal, 7(2), 102-115."
link   <- "https://www.edepositireland.ie/entities/publication/7271837d-dbe6-463f-b103-422f42a15446"

raw <- read.csv(paste0(indir, name, "/Valentia Observatory telegraphic reporting station_1850-1920.csv"),
                col.names=c("year", "month", "day", "maxF", "minF", "maxC", "minC"))

# Create a Date column for easy filtering
raw$Date <- as.Date(paste(raw$year, raw$month, raw$day, sep="-"))

# Define the split date
relocation_date <- as.Date("1892-03-12")

head(raw)
tail(raw)

# ---------------------------------------------------------
# PERIOD 1: Valentia Island (Knightstown) 1850 - 1892-03-11
# ---------------------------------------------------------

lat1 <- 51.925
lon1 <- -10.285
alt1 <- 9
name1 <- "Valentia_Island"
metaHead1 <- paste0(
  "Observer=Revenue House (1850-1892) | ",
  "Instrument=North Wall Screen",
  "Location=Knightstwon, Valentia Island | Note=Max temp to PRECEDING day until Feb 1912"
)

df1 <- raw %>%
  filter(Date < relocation_date) %>%
  mutate(
    hour=NA,
    minute=NA,
    meta.Tx=paste0("orig=",maxF, "F | obs.time=9a.m."),
    meta.Tn=paste0("orig=",minF, "F | obs.time=9a.m."),
  ) %>% select(year, month, day, hour, minute, Tx=maxC, Tn=minC, meta.Tx, meta.Tn) %>%
  filter(
    !is.na(Tx) | !is.na(Tn)
  )

head(df1)
tail(df1)

for (var in c("Tn", "Tx")) {
  stat_val <- ifelse(var == "Tn", "minimum", "maximum")
  meta_val <- if(var == "Tn") df1$meta.Tn else df1$meta.Tx
  
  write_sef_f(
    as.data.frame(df1[c("year","month", "day", "hour", "minute", var)]),
    outfile = outfile.name(paste0("ILMMT_", name1), var, df1, FALSE),
    outpath = outdir,
    cod     = paste0("ILMMT-", name1),
    lat     = lat1, lon = lon1, alt = alt1,
    sou     = source, link = link, nam = name1,
    var     = var, stat = stat_val, period = "day",
    units   = units(var), meta = meta_val, metaHead = metaHead1,
    keep_na = TRUE
  )
}

# ---------------------------------------------------------
# PERIOD 2: Caherciveen (Westwood House) 1892-03-12 - 1920
# ---------------------------------------------------------
lat2 <- 51.939
lon2 <- -10.244
alt2 <- 14
name2 <- "Valentia_Caherciveen"

metaHead2 <- paste0(
  "Observer=Westwood House Staff (1892-1920) | ",
  "Instrument=Stevenson Screen | ",
  "Location=Station moved to mainland (Caherciveen) | Note=Max temp to PRECEDING day until Feb 1912"
)

df2 <- raw %>%
  filter(Date >= relocation_date) %>%
  mutate(
    hour=NA, minute=NA,
    # Note: Observation time changed to 7 a.m. in July 1908
    meta.Tx=paste0("orig=", maxF, "F"),
    meta.Tn=paste0("orig=", minF, "F")
  ) %>% 
  select(year, month, day, hour, minute, Tx=maxC, Tn=minC, meta.Tx, meta.Tn) %>%
  filter(!is.na(Tx) | !is.na(Tn))

for (var in c("Tn", "Tx")) {
  stat_val <- ifelse(var == "Tn", "minimum", "maximum")
  meta_val <- if(var == "Tn") df2$meta.Tn else df2$meta.Tx
  
  write_sef_f(
    as.data.frame(df2[c("year","month", "day", "hour", "minute", var)]),
    outfile = outfile.name(paste0("ILMMT_", name2), var, df2, FALSE),
    outpath = outdir,
    cod     = paste0("ILMMT-", name2),
    lat     = lat2, lon = lon2, alt = alt2,
    sou     = source, link = link, nam = name2,
    var     = var, stat = stat_val, period = "day",
    units   = units(var), meta = meta_val, metaHead = metaHead2,
    keep_na = TRUE
  )
}