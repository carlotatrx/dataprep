library(dplyr)
library(readxl)
library(purrr)
library(tidyr)
source('/scratch2/ccorbella/code/dataprep/helpfun.R')

outdir <- "/scratch3/PALAEO-RA/daily_data/final/Barcelona"


lat <- 41.382751
lon <- 2.172730
alt <- 20
name <- "Barriendos_Barcelona"
code <- "Barcelona"
source <- "Rodríguez, R., et al. `Long pressure series for Barcelona (Spain). Daily reconstruction and monthly homogenization.` International Journal of Climatology: A Journal of the Royal Meteorological Society 21.13 (2001): 1693-1704."
metaHead <- "Observer=Dr. F. Salva | location=Petritxol Street, nr. 11, top of buildiing"

raw <- read_excel('/scratch3/PALAEO-RA/daily_data/original/Barcelona/BAR1811def.XLS')

head(raw)

df <- raw %>%
  mutate(
    Hour=NA,
    Minute=NA,
    p=round(P2/10,2)
  )

var<-"p"  

write_sef_f(
  as.data.frame(df[,c("year","month", "day", "Hour", "Minute", "p")]),
  outfile = outfile.name(name, var, df, FALSE),
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
  units   = units(var),
  metaHead = paste0(metaHead, " | PTC=? | PGC=?"),
  keep_na = TRUE
)

