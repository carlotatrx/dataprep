rm(list=ls())

library(dataresqc)
library(readxl)
library(lubridate)
library(dplyr)

source("/home/ccorbella/scratch2_symboliclink/code/dataprep/helpfun.R")

outdir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/'


# Cornes pressure in one series -------------------------------------------

file <- '/scratch3/PALAEO-RA/daily_data/original/London/london_MSLP_1692_2007.txt'

raw <- suppressWarnings(read.table(
  file = file,
  skip = 31,
  #col_types = cols(.default = col_guess()),
  col.names = c("Year", "Month", "Day", "p.orig", "p.hom", "diff", "code", "holborn_code", "manley_code"),
  header = FALSE,
  na.strings = NA
))

head(raw)

df <- raw %>%
  mutate(
    Hour = NA_integer_,
    Minute = NA_integer_,
  ) %>% select(Year, Month, Day, Hour, Minute, p.orig)

write_sef_f(as.data.frame(df),
            outfile=paste0('London_Cornes_', get_date_range(df), "_p_daily.tsv"),
            outpath="/scratch3/PALAEO-RA/daily_data/final/London/",
            variable = "p",
            cod = "London",
            nam = "London",
            lat = 51.5,
            lon = -0.15,
            alt = 20,
            sou = " Cornes, R. C., Jones, P. D., Briffa, K. R. and Osborn, T. J. (2012) A daily series of mean sea-level pressure for London, 1692-2007. International Journal of Climatology 32, 641-656. doi: 10.1002/joc.2301",
            units = "hPa",
            period = "day",
            stat = "mean",
            metaHead = "PTC=Y | PGC=Y | a homogenized version exists, see `raw` file",
            keep_na = TRUE)


## qc
qc("/scratch3/PALAEO-RA/daily_data/final/London/London_Cornes_16920101-20071231_p_daily.tsv",
   outpath="/scratch3/PALAEO-RA/daily_data/final/London/")
write_flags_f(infile="/scratch3/PALAEO-RA/daily_data/final/London/London_Cornes_16920101-20071231_p_daily.tsv",
              qcfile="/scratch3/PALAEO-RA/daily_data/final/London/qc_London_p_daily.txt", 
              outpath="/scratch3/PALAEO-RA/daily_data/final/London/",
              match=FALSE)

# earliest: Locke 1669-1675 -----------------------------------------------

## NOTE London_1 (Locke) was not formatted (unknown temperature scale)

file.name <- '/scratch3/PALAEO-RA/DataRescue/BORIS/XLS/Europe/London/Europe_T3_GB_London_Locke_1669-1675_subdaily.xls'

raw <- read_excel(file.name, skip=6)
colnames(raw) <- c('Year','Month','Day','Date','ta','p','dd.orig', 'ddextra', 'comments')
raw <- raw[,1:7]

df <- raw %>%
  mutate(
    Year = suppressWarnings(as.integer(Year)),
    Month = suppressWarnings(as.integer(Month)),
    Day = suppressWarnings(as.integer(Day)),
    Hour = hour(Date),
    Minute = 0,
    
    dd_norm = dd_normalize(dd.orig),
    dd = dd2deg(dd_norm),
    
    meta = paste0(meta_time(Hour, Minute), " | orig.dd=", dd.orig)
  ) %>%
  tidyr::fill(Year, Month, Day, .direction = "down") %>%
  select(Year, Month, Day, Hour, Minute, dd, meta)

write_sef_f(as.data.frame(df),
            outfile="London_Locke_16720131-16750613_dd_subdaily.tsv",
            outpath="/scratch3/PALAEO-RA/daily_data/final/London/",
            variable = "dd",
            cod = "London_Locke",
            nam = "London",
            lat = 51.76,
            lon = 0.21,
            alt = 55,
            sou = "Cornes",
            link = "doi: 10.6084/m9.figshare.24242302.v1",
            units = "dd",
            stat = "point",
            time_offset = time.offset(0.21),
            keep_na = FALSE
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






