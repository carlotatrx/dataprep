# CBT

rm(list=ls())
library(dataresqc)
library(lubridate)
library(dplyr)

source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')
outdir <- "/scratch3/PALAEO-RA/daily_data/final"
  
code <-	"CBT"
name <-	"Central Belgium Temperature"
lat  <-	50.85045
lon	 <- 4.34878
alt	 <- 13
source <- "Demarée, Gaston R. et al. “The Long-Term Daily Central Belgium Temperature (CBT) Series (1767–1998) and Early Instrumental Meteorological Observations in Belgium.” Climatic Change 53 (2002): 269-293."
link <- "https://doi.org/10.1023/A:1014931211466"

# Read list of flagged dates
flagged <- read.table("/scratch3/PALAEO-RA/daily_data/original/CBT/Cbt-error.txt", header=FALSE)
flagged_dates <- as.Date(flagged$V1, format="%d/%m/%Y")


df <- read.table("/scratch3/PALAEO-RA/daily_data/original/CBT/CBT_T767_899.txt", header=F) %>%
  rename(
    Date = V1,
    Tn = V2,
    Tx = V3,
    ta = V4
  )

df$Date <- as.Date(df[,1], format="%d/%m/%Y")
df$Year <- format(df$Date, "%Y")
df$Month <- as.integer(format(df$Date, "%m"))
df$Day <- as.integer(format(df$Date, "%d"))
df$Hour <- NA
df$Minute <- NA

df[df == -999] <- NA

df$meta <- ifelse(df$Date %in% flagged_dates, "qc=Tn>Tx", "")

head(df)

vars <- c("Tn", "Tx", "ta")

for (var in vars) {
  dat <- df[, c("Year", "Month", "Day", "Hour", "Minute", var)]
  

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
      var=="Tx" ~ "max"
    ),
    period   = "day",
    units    = "C",
    metaHead = ifelse(var == "ta", "ta=(Tx+Tn)/2", ""),
    meta     = df$meta,
    keep_na  = FALSE
  )
}

