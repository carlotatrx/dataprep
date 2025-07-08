library(dataresqc)
library(dplyr)
library(lubridate)
library(readr)

source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')

indir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/'

infile <- "Dnipro_ta_subdaily.tsv"
output_file <- "Dnipro_ta_subdaily_corrected.tsv"

# Read SEF file
df <- read.delim(paste0(indir, infile), header=T, sep='\t', stringsAsFactors = F,
                 skip=12)

meta <- read_meta_nonofficial(paste0(indir,infile))

# Ensure Meta column exists and is character
if (!"Meta" %in% names(df)) {
  df$meta <- NA_character_
}

# Define cutoff date
cutoff_date <- as.Date("1838-12-15")

# Create full date column
df$orig_date <- as.Date(with(df, sprintf("%04d-%02d-%02d", Year, Month, Day)))
df$apply_shift <- df$orig_date > cutoff_date

# Apply date shift where needed
df <- df %>%
  mutate(
    full_date = if_else(apply_shift, orig_date + 12, orig_date),
    Meta = if_else(apply_shift,
                   if_else(
                     is.na(Meta) | Meta == "", paste0("orig_date=", orig_date), paste0(Meta, " | orig_date=", orig_date)
                     ),
                   Meta)
  )

# Update Year, Month, Day
df$Year <- year(df$full_date)
df$Month <- month(df$full_date)
df$Day <- day(df$full_date)

# Drop helper columns
df <- df %>% select(-orig_date, -apply_shift, -full_date)

write_sef_f(Data=df,
            outpath=indir, outfile="Dnipro_p_subdaily_dates.tsv",
            cod=meta[["id"]],
            variable=meta[['var']],
            nam=meta[["name"]],
            lat=meta[["lat"]],
            lon=meta[["lon"]], alt=meta[["alt"]], sou=meta[["source"]], 
            metaHead = if_else(meta[['meta']] == "", "CC=Y", paste0(meta[['meta']]," | CC=Y")),
            link=meta[["link"]], units=meta[['units']], stat=meta[['stat']],
            meta=df$Meta, keep_na = F)

