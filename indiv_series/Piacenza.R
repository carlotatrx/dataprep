rm(list=ls())

library(dataresqc)
library(dplyr)
library(readxl)

source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')
outdir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/'

lat <- 45.035278
lon <- 9.725000
ele <- 78

file.name <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_orig/Piacenza-Alberoni_serie_1871-2022.xls'
df <- read_excel(file.name, sheet=2, skip=2,
                 col_types=c('numeric','numeric','numeric','numeric','numeric','numeric','numeric','numeric','numeric'),
                 col_names=c('day','month','year','Tmin','Tmax','Tmin_floor','Tmax_floor','precip','snow'))

df <- df %>% mutate(across(c(year, month, day), as.integer))

df <- df %>%
  mutate(ta=(as.numeric(Tmin) + as.numeric(Tmax))/2,
         hour=NA_integer_,
         minute=NA_integer_)

head(df)

write_sef_f(Data = df[, c('year','month','day', 'hour','minute', 'ta')],
            outfile="Piacenza_ta_daily.tsv",
            outpath='/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/',
            cod="Piacenza",
            variable="ta",
            nam="Piacenza",
            lat=lat, lon=lon, alt=ele, stat='mean',
            sou="Collegio Alberoni --  Società Meteorologica Italiana", units="C", keep_na = F
            )
