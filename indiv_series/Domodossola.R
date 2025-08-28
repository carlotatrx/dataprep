rm(list=ls())

library(dataresqc)
library(dplyr)
library(readxl)

source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')
outdir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/'

lat <- 46.1134
lon <- 8.2874
ele <- 280

file.name <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_orig/Domodossola-Rosmini_1871-2014_giorn.xls'
df <- read_excel(file.name, sheet=1, skip=1,
                 col_types=c('date','numeric','numeric','numeric','numeric','numeric','numeric','numeric'),
                 col_names=c('Date_orig','Tmin','Tmax','Prcp','Hn','Nuv1','Nuv2','Nuv3'))

df <- df %>%
  mutate(ta = ifelse(Tmin==999.9 | Tmax==999.9, NA,
                     round((as.numeric(Tmin) + as.numeric(Tmax))/2,2)),
         hour=24,
         minute=0,
         date = seq(as.Date("1871/12/1"), as.Date("2014/06/30"), "day"),
         meta=paste0('orig.ta.min=',Tmin, 'C | orig.ta.max=',Tmax, "C")
         )

# columns for Year Month Day
df <- df %>%
  mutate(year = format(date,"%Y"),
         month = format(date, "%m"),
         day = format(date, "%d")
  )
head(df)

write_sef_f(Data = as.data.frame(df[, c('year','month','day', 'hour','minute', 'ta')]),
            outfile="Domodossola_ta_daily.tsv",
            outpath='/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/',
            cod="Domodossola",
            variable="ta",
            nam="Domodossola",
            lat=lat, lon=lon, alt=ele, stat='mean', period='day', meta=df$meta,
            sou="Collegio Alberoni --  Società Meteorologica Italiana", units="C", keep_na = F
            )
