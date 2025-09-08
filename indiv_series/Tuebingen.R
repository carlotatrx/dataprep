rm(list=ls())
library(dataresqc)
library(readxl)
library(dplyr)
library(tidyr)
source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')

name <- "Tuebingen"
cod <- "Europe_Tuebingen_1"

file <- paste0('/scratch3/PALAEO-RA/DataRescue/Projects/Europe/3_Revised/',name,'/Europe_T3_DE_',name,'_1691-1694_subdaily.xlsx')
raw <- read_excel(file, sheet = 1, skip  = 6)

lat	<- 48.5216364
lon	<- 9.0576448
alt	<- NA
metaHead <- "Observer=Camerarius"

# pressure conversion
conversions <- list(p  = function(x, ta) round(convert_pressure(as.numeric(x), f = 27.07, lat = lat, alt = alt, atb=ta), 1))
time.offset <- as.numeric(lon)*12/180

# Keep only the first 10 columns and rename
df <- raw %>%
  select(c(1:4,7) %>%
  rename(
    Year = 1, Month = 2, Day = 3, Hour = 4, p.zoll = 4, p.in = 5, ta.R = 6, dd.orig = 7
  ) %>%
 ...
 # pressure units unknown