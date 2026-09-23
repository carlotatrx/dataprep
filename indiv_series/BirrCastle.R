rm(list=ls())
library(dataresqc)
library(readxl)
library(tidyr)
library(dplyr)

source('/scratch2/ccorbella/code/dataprep/helpfun.R')



indir <- "/scratch3/PALAEO-RA/daily_data/original/"
name <- "BirrCastle"

outdir <- paste0("/scratch3/PALAEO-RA/daily_data/SEF/",name)

source <- "Mateus, C., Potito, A., & Curley, M. (2020). Reconstruction of a long‐term historical daily maximum and minimum air temperature network dataset for Ireland (1831‐1968). Geoscience Data Journal, 7(2), 102-115."
link   <- "https://www.edepositireland.ie/entities/publication/7271837d-dbe6-463f-b103-422f42a15446"
lat <- 53.095 
lon <- -7.922
alt <- 50
metaHead <- "Observer=William Harding and G. Keating (1880), William Harding (1880-1881, 1887-1888), George Phillips (1881-1885), Benjamin Budds (1885-1887), Edmund Haines (1888-1892), James Perry (1892-1893), John Perry (1893-1894), Thomas Haines (1894-1898), A. Colvin (1899-1901), John L. Roe (1901-1903), George A. Roe (1903-1904), Wellesley J. Roe (1904-1911), Dr. Otto Boeddicker (1881-1916)  | Instrument=Stevenson screen (1880-1911), Double thermometer screen (1880-1911), Max thermometer B.T.492 (1887-1911), Min thermometer B.T.497 (1887-1911), Thermometer screen (1912-1954), Max thermometer M.O.2001 (1912-1925), Min thermometer M.O.2015 (1912-1925), Max thermometer M.O.11075 (1936-1950), Min thermometer M.O.1238 (1936-1941), Min thermometer M.O.10860 (1942-1950)"


# Birr Castle telegraphic reporting station -------------------------------
raw <- read.csv(paste0(indir, name, "/Birr Castle telegraphic reporting station_1880-1920.csv"),
                col.names=c("year", "month", "day", "maxF", "minF", "maxC", "minC"))
head(raw)

df <- raw %>%
  mutate(
    hour=NA,
    minute=NA,
    meta.Tx=paste0("orig=",maxF, "F | obs.time=8a.m."),
    meta.Tn=paste0("orig=",minF, "F | obs.time=8a.m.")
  ) %>% select(year, month, day, hour, minute, Tx=maxC, Tn=minC, meta.Tx, meta.Tn)

head(df)



# Birr Castle Second order station ----------------------------------------

raw <- read.csv(paste0(indir, name, "/Birr Castle second order station_1872-1911.csv"),
                col.names=c("year", "month", "day", "maxF", "minF", "maxC", "minC"))
head(raw)

df <- raw %>%
  mutate(
    hour=NA,
    minute=NA,
    meta.Tx=paste0("orig=",maxF, "F | obs.time=8a.m."),
    meta.Tn=paste0("orig=",minF, "F | obs.time=8a.m.")
  ) %>% select(year, month, day, hour, minute, Tx=maxC, Tn=minC, meta.Tx, meta.Tn)

head(df)

var <- "Tn"

write_sef_f(
  as.data.frame(df[c("year","month", "day", "hour", "minute","Tn")]),
  outfile = outfile.name(paste0("ILMMT_",name), var, df, FALSE),
  outpath = outdir,
  cod     = paste0("ILMMT-",name),
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "minimum",
  period    = "day",
  units   = units(var),
  metaHead = metaHead,
  meta = df$meta.Tn
)


var <- "Tx"

write_sef_f(
  as.data.frame(df[c("year","month", "day", "hour", "minute",var)]),
  outfile = outfile.name(paste0("ILMMT_",name), var, df, FALSE),
  outpath = outdir,
  cod     = paste0("ILMMT-",name),
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "maximum",
  period    = "day",
  units   = units(var),
  metaHead = metaHead, 
  meta=df$meta.Tx
)





