library(dataresqc)
library(dplyr)
library(lubridate)
library(readr)

source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')

indir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/'


sef_calendar_correction <- function(infile,
                                    outfile,
                                    cutoff_date, ## format YYYY-MM-DD
                                    skip_days,
                                    indir='/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/') {
  # Read SEF file
  df <- read.delim(paste0(indir, infile), header=T, sep='\t', stringsAsFactors = F,
                   skip=12)

  meta <- read_meta_nonofficial(paste0(indir,infile))

  # Ensure Meta column exists and is character
  if (!"Meta" %in% names(df)) {
    df$Meta <- NA_character_
  }
  
  cutoff_date <- as.Date(cutoff_date)


  # Create full date column
  df$orig_date <- as.Date(with(df, sprintf("%04d-%02d-%02d", Year, Month, Day)))
  df$apply_shift <- df$orig_date > cutoff_date
  df$full_date <- df$orig_date
  df$full_date[df$apply_shift] <- df$orig_date[df$apply_shift] + skip_days

  # Sanity check: full_date should be orig_date + 12 where shifted
  stopifnot(all(df$full_date[df$apply_shift] - df$orig_date[df$apply_shift] == skip_days))
  
  # Apply date shift where needed
  df <- df %>%
    mutate(Meta = if_else(apply_shift,
                          if_else(is.na(Meta) | Meta == "",
                                  paste0("orig_date=", format(orig_date, "%Y-%m-%d")),
                                  paste0(Meta, " | orig_date=", format(orig_date, "%Y-%m-%d"))),
                          Meta))
  

  # Update Year, Month, Day
  df$Year <- year(df$full_date)
  df$Month <- month(df$full_date)
  df$Day <- day(df$full_date)
  
  # Drop helper columns. Period col needs to be dropped too!
  df <- df %>% select(-orig_date, -apply_shift, -full_date, -Period)
  meta.head <- if_else(meta[['meta']] == "", "CC=Y", paste0(meta[['meta']]," | CC=Y"))

  write_sef_f(Data=df,
              outpath=indir, outfile=outfile,
              cod=meta[["id"]],
              variable=meta[['var']],
              nam=meta[["name"]],
              lat=meta[["lat"]],
              lon=meta[["lon"]], alt=meta[["alt"]], sou=meta[["source"]], 
              metaHead = meta.head,
              link=meta[["link"]], units=meta[['units']], stat=meta[['stat']],
              meta=df$Meta, keep_na = F)
}

# Dnipro
sef_calendar_correction(
  infile = "Dnipro_ta_subdaily.tsv",
  outfile = "Dnipro_ta_subdaily_corrected.tsv",
  cutoff_date = as.Date("1838-12-15"),
  skip_days = 12
)

sef_calendar_correction(
  infile = "Dnipro_p_subdaily.tsv",
  outfile = "Dnipro_p_subdaily_corrected.tsv",
  cutoff_date = as.Date("1838-12-15"),
  skip_days = 12
)

# Kharkiv
sef_calendar_correction(
  infile = "Kharkiv_ta_subdaily.tsv",
  outfile = "Kharkiv_ta_subdaily_corrected.tsv",
  cutoff_date = as.Date("1845-05-15"),
  skip_days = 12
)

sef_calendar_correction(
  infile = "Kharkiv_p_subdaily.tsv",
  outfile = "Kharkiv_p_subdaily_corrected.tsv",
  cutoff_date = as.Date("1845-05-15"),
  skip_days = 12
)

# Kherson
sef_calendar_correction(
  infile = "Kherson_ta_subdaily.tsv",
  outfile = "Kherson_ta_subdaily_corrected.tsv",
  cutoff_date = as.Date("1700-01-01"), # since the start
  skip_days = 12
)

sef_calendar_correction(
  infile = "Kherson_p_subdaily.tsv",
  outfile = "Kherson_p_subdaily_corrected.tsv",
  cutoff_date = as.Date("1700-01-01"),
  skip_days = 12
)

# Lugansk
sef_calendar_correction(
  infile = "Lugansk_ta_subdaily.tsv",
  outfile = "Lugansk_ta_subdaily_corrected.tsv",
  cutoff_date = as.Date("1840-01-15"),
  skip_days = 12
)

sef_calendar_correction(
  infile = "Lugansk_p_subdaily.tsv",
  outfile = "Lugansk_p_subdaily_corrected.tsv",
  cutoff_date = as.Date("1840-01-15"),
  skip_days = 12
)

# Poltava
sef_calendar_correction(
  infile = "Poltava_ta_subdaily.tsv",
  outfile = "Poltava_ta_subdaily_corrected.tsv",
  cutoff_date = as.Date("1839-01-15"),
  skip_days = 12
)

sef_calendar_correction(
  infile = "Poltava_p_subdaily.tsv",
  outfile = "Poltava_p_subdaily_corrected.tsv",
  cutoff_date = as.Date("1839-01-15"),
  skip_days = 12
)

