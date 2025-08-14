rm(list=ls())

library(dataresqc)
library(readxl)
library(lubridate)
library(dplyr)

source("/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/write_sef_f.R")

outdir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/'

# earliest: Locke 1669-1675 -----------------------------------------------

## NOTE London_1 (Locke) was not formatted (unknown temperature scale)


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

meta <- list(
  name='London',
  metaHead= 'obs=Locke',
  sou = 'Boyle, General history of the air, London, 1692, pag. 116'
)

#######################################################
# continue here! We don't know the lat lon
#######################################################

file.name.p <- '/scratch3/PALAEO-RA/Instrumental/SEF/Daily/QCed/London/Cornes_Boyle_London_16841211-16860210_p.tsv'
file.name.ta <- '/scratch3/PALAEO-RA/Instrumental/SEF/Daily/QCed/London/Cornes_Boyle_London_16841211-16860210_ta.tsv'

df.p <- read_sef(file.name.p)
df.ta <- read_sef(file.name.ta)
meta.p <- read_meta(file.name.p)

df.all <- merge(df.ta, df.p, by=c('Year', 'Month','Day', 'Hour', 'Minute'), suffixes=c('.ta','.p'))

df.all <- df.all %>%
  mutate(
    p.corr = round(convert_pressure(Value.p, f=0.75006156130264, lat=as.numeric(meta.p[['lat']]), atb=Value.ta),1),
    meta.p = paste0("atb=", Value.ta,"C | p_corr=", round(p.corr-Value.p,2))
  )


head(df.all)


write_sef_f(Data=df.all[, c("Year","Month","Day","Hour","Minute","p.corr")], outfile="London_16841211-16860210_p_subdaily.tsv",
            outpath='/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/',
            cod=meta[["id"]],
            variable=meta[['var']],
            nam=meta[["name"]],
            lat=meta[["lat"]],
            lon=meta[["lon"]], alt=meta[["alt"]], link=meta[['link']], stat='point',
            sou=meta[["source"]], units=meta[["units"]],keep_na = F, meta=df.all$meta.p,
            metaHead="Observer=Boyle | PTC=Y | PGC=Y"
)






