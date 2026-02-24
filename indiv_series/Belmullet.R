rm(list=ls())
library(dataresqc)
library(readxl)
library(tidyr)
library(dplyr)

source('/scratch2/ccorbella/code/dataprep/helpfun.R')

indir <- "/scratch3/PALAEO-RA/daily_data/original/"
name <- "Belmullet"

outdir <- paste0("/scratch3/PALAEO-RA/daily_data/final/",name)

source <- "Mateus, C., Potito, A., & Curley, M. (2020). Reconstruction of a long‐term historical daily maximum and minimum air temperature network dataset for Ireland (1831‐1968). Geoscience Data Journal, 7(2), 102-115."
link   <- "https://www.edepositireland.ie/entities/publication/7271837d-dbe6-463f-b103-422f42a15446"
lat <- 54.223
lon <- -9.988
alt <- 5
metaHead <- paste0(
  "Observer=Mary J. Tolan (1884-1896), Edward Tolan (1890-1891), Emily Tolan (1896-1899), ",
  "Joseph Hodge (1899), A. Marshall (1899-1904), James Heddon (1904-1906), Esceter (1906-1908), ",
  "William Barry (1909-1913), Charles Sammels (1913-1914), Arthur J. Read (1914-1918), ",
  "James Cole (1918-1919), A. Stone (1919), C. Kelly (1919-1920), S. J. Palmer (1920) | ",
  "Instrument=Stevenson thermometer screen (1884-1920), Max thermometer (1884-1920), ",
  "Min thermometer (1884-1920). | Location=relocation on 25 September 1899 to 54°6’N 10°4’W"
)

# To verify it is one single line:
print(metaHead)

raw <- read.csv(paste0(indir, name, "/Blacksod Point_Belmullet_1872-1920.csv"),
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





