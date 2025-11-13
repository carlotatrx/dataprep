# Padua.R 

rm(list=ls())
library(dataresqc)
library(lubridate)
library(stringr)
library(suncalc) # to calculate sunrise hours
library(dplyr)   # to concatenate dfs
library(readr)   # to read text files

source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')

# Padua_Poleni ------------------------------------------------------------

# PD_Poleni contains daily temperatures at noon derived from the observations 
# of Giovanni Poleni and his son Francesco between 12/01/1725 and 31/12/1764

code <-	"Padua"
name <-	"Padua"
lat  <-	45.4
lon	 <- 11.87
alt	 <- 18
source <- "Stefanini, Claudio; Becherini, Francesca; della Valle, Antonio; Camuffo, Dario (2024). Padua daily temperature 1725-2023. figshare. Dataset."
link <- "https://doi.org/10.6084/m9.figshare.25471507.v1"
metahead <- "Observer=Giovanni Poleni, Francesco Poleni (1725-1764); Giambattista Morgagni (1740-1768); Giuseppe Toaldo (1766-1773)"

df <- read.table("/scratch3/PALAEO-RA/daily_data/original/Padua/PD_Poleni.csv", header=T, sep=";")
head(df)

df$Date <- as.Date(df$Date, format="%d/%m/%Y")
df$Year <- as.numeric(format(df$Date, "%Y"))
df$Month <- as.numeric(format(df$Date, "%m"))
df$Day <- as.numeric(format(df$Date, "%d"))
df$Hour <- 12
df$Minute <- 0
df$meta <- "obs.time=noon"

df <- df[, c("Year", "Month", "Day", "Hour", "Minute", "t_h12", "meta")]

head(df)



# PD_Morgagni -------------------------------------------------------------

# PD_Morgagni contains daily temperatures at noon derived from the observations 
# of Giambattista Morgagni between 01/01/1740 and 31/12/1768;


df2 <- read.table("/scratch3/PALAEO-RA/daily_data/original/Padua/PD_Morgagni.csv", header=T, sep=";")
head(df2)

df2$Date <- as.Date(df2$Date, format="%d/%m/%Y")
df2$Year <- as.numeric(format(df2$Date, "%Y"))
df2$Month <- as.numeric(format(df2$Date, "%m"))
df2$Day <- as.numeric(format(df2$Date, "%d"))
df2$Hour <- 12
df2$Minute <- 0
df2$meta <- "obs.time=noon"

df2 <- df2[, c("Year", "Month", "Day", "Hour", "Minute", "t_h12", "meta")]

head(df2)


# PD_Toaldo -------------------------------------------------------------

# PD_Toaldo contains daily temperatures near the sunrise collected
# by Giuseppe Toaldo between 01/05/1766 and 31/12/1773


df3 <- read.table("/scratch3/PALAEO-RA/daily_data/original/Padua/PD_Toaldo.csv", header=T, sep=";")
head(df3)

df3$Date <- as.Date(df3$Date, format="%d/%m/%Y")
df3$Year <- as.numeric(format(df3$Date, "%Y"))
df3$Month <- as.numeric(format(df3$Date, "%m"))
df3$Day <- as.numeric(format(df3$Date, "%d"))
sun <- suncalc::getSunlightTimes(date = df3$Date, lat=lat, lon=lon, tz="UTC")
sunrise_lmt <- sun$sunrise + lon/15*3600 # calculate UTC to padua local mean time
df3$Hour <- as.integer(format(sunrise_lmt, "%H"))
df3$Minute <- as.integer(format(sunrise_lmt, "%M"))
df3$meta <- "obs.time=sunrise"

df3 <- df3[, c("Year", "Month", "Day", "Hour", "Minute", "t_sunrise", "meta")]

head(df3)


# combine dfs -------------------------------------------------------------

df.all <- bind_rows(df, df2, df3) %>%
  mutate(Value = coalesce(t_h12, t_sunrise)) %>%
  select(-any_of(c("t_h12", "t_sunrise"))) %>%
  relocate(Value, .before=meta)

head(df.all)

# save --------------------------------------------------------------------

write_sef_f(df.all,
            outfile=paste0(name,"_",get_date_range(df.all), "_ta_daily.tsv"),
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

# ----------------------------------------------------------------------------------------
# -- Padua Pressure ----------------------------------------------------------------------
# ----------------------------------------------------------------------------------------

code <- "Padua"
name <- "Padua"
source <- "IMPROVE"
lat  <-	45.4
lon	 <- 11.87
alt	 <- 18

read_padua <- function(file) {
  # Read with minimal guessing; everything as character first
  df_raw <- read_table(file, col_names = FALSE, col_types = cols(.default = "c"), na = "-999.0")
  ncol_df <- ncol(df_raw)
  
  if (ncol_df == 3) {
    message("3-column file: ", file)
    df <- df_raw %>%
      transmute(
        Date = as.Date(X1, format = "%d/%m/%Y"),
        Year = year(Date),
        Month = month(Date),
        Day = day(Date),
        Hour = NA_integer_,
        Minute = NA_integer_,
        Tn = NA_real_,
        Tx = NA_real_,
        ta = as.numeric(X2),
        p  = as.numeric(X3)
      )
  } else if (ncol_df >= 5) {
    message("5-column file: ", file)
    df <- df_raw %>%
      transmute(
        Date   = as.Date(X1, format = "%d/%m/%Y"),
        Year   = year(Date),
        Month  = month(Date),
        Day    = day(Date),
        Hour   = NA_integer_,
        Minute = NA_integer_,
        Tn = as.numeric(X2),
        Tx = as.numeric(X3),
        ta = as.numeric(X4),
        p  = as.numeric(X5)
      )
  } else {
    warning("Unexpected column count in file: ", file)
    return(NULL)
  }
  
  return(df)
}



indir <- '/scratch3/PALAEO-RA/daily_data/original/Padua'
outdir <- '/scratch3/PALAEO-RA/daily_data/final'
files <- list.files(indir, pattern="^PD_PT.*\\.txt$", full.names=T)
files  

df <- bind_rows(lapply(files, read_padua))


head(df)

vars <- c("Tn", "Tx", "ta", "p")

for (var in vars) {
  dat <- df[, c("Year", "Month", "Day", "Hour", "Minute", var)]
  
  # has diffferent dates and we have to remove NA rows before
  dat <- dat %>% drop_na(var)
  
  
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
      var=="Tx" ~ "max",
      var=="p"  ~ "mean"
    ),
    period   = "day",
    metaHead = metahead,
    units    = ifelse(var=="p", "hPa", "C"),
    keep_na  = FALSE
  )
}
