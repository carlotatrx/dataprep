# CET

rm(list=ls())
library(dataresqc)
library(lubridate)

source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')

code <-	"HadCET"
name <-	"Central England Temperature"
lat  <-	52.5
lon	 <- -1.9
alt	 <- 44
source <- "Parker, D.E., T.P. Legg, and C.K. Folland. 1992. A new daily Central England Temperature Series, 1772-1991. Int. J. Clim., Vol 12, pp 317-342"
link <- "https://www.metoffice.gov.uk/hadobs/hadcet/data/download.html"

df <- read.table("/scratch3/PALAEO-RA/daily_data/original/CET/meantemp_daily_totals.txt", header=T)

df$Date <- as.Date(df$Date)
df$Year <- format(df$Date, "%Y")
df$Month <- format(df$Date, "%m")
df$Day <- format(df$Date, "%d")
df$Hour <- NA
df$Minute <- NA

df <- df[, c("Year", "Month", "Day", "Hour", "Minute", "Value")]

head(df)

write_sef_f(df,
            outfile="CET_ta_subdaily.tsv",
            outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
            variable = "ta",
            cod = code,
            nam = name,
            lat = lat,
            lon = lon,
            alt = alt,
            sou = source,
            link = link,
            units = "C",
            stat = "point",
            period = 0,
            keep_na = F)
