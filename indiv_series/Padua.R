# Padua.R 

rm(list=ls())
library(dataresqc)
library(lubridate)

source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')

# Padua_Poleni ------------------------------------------------------------

# PD_Poleni contains daily temperatures at noon derived from the observations 
# of Giovanni Poleni and his son Francesco between 12/01/1725 and 31/12/1764

code <-	"Padua_Poleni"
name <-	"Padua"
lat  <-	45.4
lon	 <- 11.87
alt	 <- 18
source <- "Stefanini, Claudio; Becherini, Francesca; della Valle, Antonio; Camuffo, Dario (2024). Padua daily temperature 1725-2023. figshare. Dataset."
link <- "https://doi.org/10.6084/m9.figshare.25471507.v1"
metahead <- "Observer=Giovanni Poleni, Francesco Poleni"

df <- read.table("/scratch3/PALAEO-RA/daily_data/original/Padua/PD_Poleni.csv", header=T, sep=";")
head(df)

df$Date <- as.Date(df$Date, format="%d/%m/%Y")
df$Year <- as.numeric(format(df$Date, "%Y"))
df$Month <- as.numeric(format(df$Date, "%m"))
df$Day <- as.numeric(format(df$Date, "%d"))
df$Hour <- NA
df$Minute <- NA

df <- df[, c("Year", "Month", "Day", "Hour", "Minute", "t_h12")]

head(df)

write_sef_f(df,
            outfile=paste0(name,"_",get_date_range(df), "_ta_daily.tsv"),
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
            metaHead = metahead,
            keep_na = F)



# PD_Morgagni -------------------------------------------------------------

# PD_Morgagni contains daily temperatures at noon derived from the observations 
# of Giambattista Morgagni between 01/01/1740 and 31/12/1768;

metahead <- "Observer=Giovanni Poleni, Francesco Poleni"

df <- read.table("/scratch3/PALAEO-RA/daily_data/original/Padua/PD_Poleni.csv", header=T, sep=";")
head(df)

df$Date <- as.Date(df$Date, format="%d/%m/%Y")
df$Year <- as.numeric(format(df$Date, "%Y"))
df$Month <- as.numeric(format(df$Date, "%m"))
df$Day <- as.numeric(format(df$Date, "%d"))
df$Hour <- NA
df$Minute <- NA

df <- df[, c("Year", "Month", "Day", "Hour", "Minute", "t_h12")]

head(df)

write_sef_f(df,
            outfile=paste0(name,"_",get_date_range(df), "_ta_daily.tsv"),
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
            metaHead = metahead,
            keep_na = F)

