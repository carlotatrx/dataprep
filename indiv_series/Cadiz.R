# Cadiz

rm(list=ls())
library(dataresqc)
library(lubridate)
library(dplyr)

source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')
outdir <- "/scratch3/PALAEO-RA/daily_data/final"

code <-	"Cadiz"
name <-	"Cadiz"
lat  <-	36.5166666666667
lon	 <- -6.28333333333333
alt	 <- 13
source <- "Barriendos, Mariano et al. “DAILY METEOROLOGICAL OBSERVATIONS IN CADIZ­ SAN FERNANDO. ANALYSIS OF THE DOCUMENTARY SOURCES AND THE INSTRUMENTAL DATA CONTENT (1786-1996).” Climatic Change 53 (2002): 151-170."
link <- "https://doi.org/10.1007/978-94-010-0371-1_6"

# Read list of flagged dates
filepath <- "/scratch3/PALAEO-RA/daily_data/original/Cadiz/"
files <- list.files(filepath, full.names=TRUE)

df <- files %>%
  lapply(function(f) {
    read.table(f, header = FALSE, fill = TRUE, blank.lines.skip = TRUE) %>%
      rename(Date = V1, Tn = V2, Tx = V3, ta = V4, p = V5) %>%
      mutate(
        Date = as.Date(Date, format = "%d/%m/%Y"),
        Year = as.integer(format(Date, "%Y")),
        Month = as.integer(format(Date, "%m")),
        Day = as.integer(format(Date, "%d")),
        Hour = NA,
        Minute = NA
      )
  }) %>%
  bind_rows()

df[df == -999] <- NA
df <- subset(df, Year < 1950)

head(df)

vars <- c("Tn", "Tx", "ta", "p")

for (var in vars) {
  dat <- df[, c("Year", "Month", "Day", "Hour", "Minute", var)]
  
  if (var=="p") {
    # has diffferent dates and we have to remove NA rows before
    dat <- dat %>% drop_na(p)
  }
  
  write_sef_f(
    as.data.frame(dat),
    outfile  = outfile.name(code, var, dat, FALSE),
    outpath  = file.path(outdir, code),
    cod      = code,
    lat      = lat,
    lon      = lon,
    alt      = alt,
    sou      = source,
    link     = link,
    nam      = name,
    var      = var,
    stat     = case_when(
      var=="ta" ~ "mean",
      var=="Tn" ~ "min",
      var=="Tx" ~ "max",
      var=="p"  ~ "mean"
    ),
    period   = "day",
    units    = ifelse(var=="p", "hPa", "C"),
    keep_na  = FALSE
  )
}

