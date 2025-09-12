rm(list=ls())
library(dataresqc)
library(XLConnect)

setwd('~/scratch2_symboliclink/code/KF_assimilation/dataprep/')
source('helpfun.R')

name <- "Oxford"
lat  <-	51.77
lon	 <- -1.27
alt  <- 63
metahead <- "Observer=Locke"
source   <- "Boyle, General history of the air, London, 1692, pag. 104"

time.offset <- as.numeric(lon)*12/180

get_date_range <- function(df) {
  start.date <- paste0(df$Year[1],
                       sprintf("%02d", df$Month[1]),
                       sprintf("%02d", df$Day[1]))
  n <- nrow(df)
  end.date   <- paste0(df$Year[n],
                       sprintf("%02d", df$Month[n]),
                       sprintf("%02d", df$Day[n]))
  paste0(start.date, "-", end.date)
}

clean_num  <- function(x) suppressWarnings(as.numeric(trimws(x)))
clean_hour <- function(x) suppressWarnings(as.integer(trimws(x)))


## Define function to fill missing values with the value on the previous row
fill_variable <- function(timeseries, strings) {
  for (i in 2:length(timeseries)) {
    if (timeseries[i] %in% strings) {
      if(!timeseries[i-1] %in% strings) {
        timeseries[i] <- timeseries[i-1]
      }
    }
  }
  return(timeseries)
}

raw <- readWorksheetFromFile("/scratch3/PALAEO-RA/daily_data/original/Oxford/Europe_T3_GB_Oxford_Locke_1666-1683_subdaily.xls",
                             sheet=1, startRow=9, header=FALSE, colTypes="character")

raw[,1] <- fill_variable(raw[,1], NA)
raw[,2] <- fill_variable(raw[,2], NA)
raw[,3] <- fill_variable(raw[,3], NA)

## Time
dates <- as.Date(paste(raw[,1], raw[,2], raw[,3], sep="-"))
dates <- dates + 10 # convert from Julian to Gregorian calendar
years <- as.integer(format(dates,"%Y"))
months <- as.integer(format(dates,"%m"))
days <- as.integer(format(dates,"%d"))
hours <- as.integer(substr(raw[,4], 12, 13))
minutes <- as.integer(substr(raw[,4], 15, 16))
meta.ta <- paste0("orig.date=", dates-10, " | orig.time=", substr(raw[,4],12,16), " | orig.ta=", raw[,5])

df <- data.frame(Year = years,
                 Month = months,
                 Day = days,
                 Hour = hours,
                 Minute = minutes,
                 Value = as.numeric(raw[,5]),
                 Meta = meta.ta,
                 stringsAsFactors = FALSE)

head(df)

# save
write_sef_f(
  Data=as.data.frame(df[, c("Year", "Month", "Day", "Hour", "Minute", "Value")]),
  outfile=paste0(name,"_",get_date_range(df), "_ta_subdaily.tsv"),
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod="Kanold_Nuernberg",
  lat=lat, lon=lon, alt=alt,
  variable="ta",
  nam=name,
  meta=df$Meta,
  metaHead=metahead,
  units="unknown", sou=source, stat="point",keep_na = F,
  time_offset = time.offset
)
