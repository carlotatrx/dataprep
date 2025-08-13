# Bologna
library(dataresqc)
library(dplyr)
library(stringr)
library(ls)

source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')
sef_test_path <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/sef_tests/'
outpath_preprocessed <- "/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed"

# Bologna -----------------------------------------------------------------
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
