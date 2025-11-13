# Milan

rm(list=ls())
library(dataresqc)
library(lubridate)
library(dplyr)

source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')
outdir <- "/scratch3/PALAEO-RA/daily_data/final"

code <-	"Milan"
name <-	"Milan"
lat  <-	45.47
lon	 <- 9.19
alt	 <- 150
source <- "Maugeri, Maurizio et al. “Daily Milan temperature and pressure series (1763-1998): Completing and homogenising the data.” Climatic Change 53 (2002): 119-149."
link <- "https://doi.org/10.1023/A:1014923027396"


# Read list of flagged dates
filepath <- "/scratch3/PALAEO-RA/daily_data/original/Milan"
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
  
  if (var=="ta") {
    metaHead <- "ta=(Tx+Tn)/2"
  } else {
    metaHead <- NA_character_
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

