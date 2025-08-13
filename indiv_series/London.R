library(dataresqc)
source("/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/write_sef_f.R")

outdir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/'
library(readxl)
library(lubridate)

# earliest: Locke 1669-1675 -----------------------------------------------

## NOTE London_1 (Locke) was not formatted (unknown temperature scale)

rm(list=ls())

file.name <- '/scratch3/PALAEO-RA/DataRescue/BORIS/XLS/Europe/London/Europe_T3_GB_London_Locke_1669-1675_subdaily.xls'

df <- read_excel(file.name, skip=6)
colnames(df) <- c('Year','Month','Day','Date','ta','p','dd')
df <- df[,1:7]

df <- df %>%
  mutate(
    Year = suppressWarnings(as.integer(Year)),
    Month = suppressWarnings(as.integer(Month)),
    Day = suppressWarnings(as.integer(Day)),
    Hour = hour(Date),
    Miute = 0
  ) %>%
  tidyr::fill(Year, Month, Day, .direction = "down")

head(df)
