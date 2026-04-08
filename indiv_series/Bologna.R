# Bologna

rm(list=ls())
library(dplyr)
library(readxl)
library(purrr)
library(tidyr)
source('/scratch2/ccorbella/code/dataprep/helpfun.R')

library(dataresqc)
library(stringr)
library(ls)



# Bologna PALAEO-RA 1787-1788 daily ---------------------------------------

outdir <- "/scratch3/PALAEO-RA/daily_data/final/Bologna"


lat <- 44.4967
lon <- 11.3526
alt <- 74
name <- "Bologna"
code <- "Palatine-Society_Bologna"
source <- "PALAEO-RA"
link   <- ""
raw1 <- read_excel("/scratch3/PALAEO-RA/daily_data/original/Bologna/Palatine-Society_T2_IT_Bologna_1787_daily.xls",
                  skip=6)
raw2 <- read_excel("/scratch3/PALAEO-RA/daily_data/original/Bologna/Palatine-Society_T2_IT_Bologna_1788_daily.xls",
                   skip=6)
df1 <- raw %>%
  fill(
    c("Inches", "Degrees"), 
    .direction = "down"
  ) %>% mutate(
    Hour=NA_integer_,
    Minute = NA_integer_
  )

df2 <- raw2 %>%
  fill(
    c("Inches", "Degrees"), 
    .direction = "down"
  ) %>% mutate(
    Hour=NA_integer_,
    Minute = NA_integer_
  )

df <- df1 %>% rbind(df2)
head(df)

df <- df %>%
  mutate(
    rh= `...9`,
    R = as.double(R),
    ta = round(R*1.25,2),
    p = round(convert_pressure(p=(Inches+Lines/12), f=27.07, lat=lat, alt=alt,
                                  atb=ta),2),
    meta.p = paste0("orig=", Inches, "in", Lines, "l | atb=", round(R,2), "R"),
    meta.ta = paste0("orig=", round(R,2), "R")
  ) %>% select(Year, Month, Day, Hour, Minute, p, ta, meta.ta, meta.p,rh)

head(df)

# save
var <-"p"
write_sef_f(
  as.data.frame(df[,c("Year","Month", "Day", "Hour", "Minute", "p")]),
  outfile = outfile.name("PALAEO-RA_Bologna", var, df, FALSE),
  outpath = outdir,
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "day",
  period  = "day",
  units   = units(var),
  metaHead = "PTC=Y | PGC=Y",
  meta    = df$meta.p,
  keep_na = TRUE
)

var <-"ta"
write_sef_f(
  as.data.frame(df[,c("Year","Month", "Day", "Hour", "Minute", "ta")]),
  outfile = outfile.name("PALAEO-RA_Bologna", var, df, FALSE),
  outpath = outdir,
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "day",
  period  = "day",
  units   = units(var),
  meta    = df$meta.ta,
  keep_na = TRUE
)

var <-"rh"
write_sef_f(
  as.data.frame(df[,c("Year","Month", "Day", "Hour", "Minute", "rh")]),
  outfile = outfile.name("PALAEO-RA_Bologna", var, df, FALSE),
  outpath = outdir,
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "day",
  period  = "day",
  metaHead="Hygrometer (Siccitatis)",
  units   = "unknown",
  keep_na = TRUE
)




# Bologna -----------------------------------------------------------------

source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')
sef_test_path <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/sef_tests/'
outpath_preprocessed <- "/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed"


vars <- c('ta','p','rr','dd')
units <- c('C','hPa','mm','degree')
stats <- ifelse(vars == "rr", "sum", "point")
period <- ifelse(vars=='rr', 'day', 0)

indir <- "/scratch3/PALAEO-RA/DataRescue/Projects/Palatine-Society/4_Formatted/Bologna"

for (i in seq_along(vars)) {
  var <- vars[i]
  pattern <- paste0("PALAEO-RA_Palatine-Society_Bologna_.*_", var, ".tsv")
  
  file.name <- dir_ls(indir, regexp=pattern)
  cat("working on file ",file.name)
  
  df <- read.delim(file.name,
                   header=T, sep="\t", skip=12)
  meta <- read_meta(file.name)
  

  df <- df %>%
    mutate(
      Minute = 0,
      Period = NULL,
      # Construct YYYY-MM-DD from Year/Month/Day
      date_string = sprintf("%04d-%02d-%02d", Year, Month, Day),
  
      Meta = str_replace_all(Meta, "\\|", " | "),
      Meta = str_replace_all(Meta, "orig=", paste0("orig_", var, "=")),
      Meta = str_replace_all(Meta, "orig.time", "orig_time"),
      # Remove orig.date=YYYY-MM-DD  |  if it matches date_string
      Meta = if_else(
        str_detect(Meta, paste0("orig.date=", date_string, " \\| ")),
        str_remove(Meta, paste0("orig.date=", date_string, " \\| ")),
        Meta
      )
    ) %>%
    select(-date_string)  # Drop helper column

  write_sef_f(Data=df, outfile=paste0("Bologna_",var,"_subdaily.tsv"),
              outpath=outpath_preprocessed,
              cod=meta[["id"]],
              metaHead=meta[["meta"]],
              variable=var,
              nam=meta[["name"]], link=meta[['link']],
              lat=meta[["lat"]],
              lon=meta[["lon"]], alt=meta[["alt"]], 
              period=period[i],
              sou=meta[["source"]],
              units=units[i], stat=stats[i], meta=df$Meta, keep_na = F
  )
}
