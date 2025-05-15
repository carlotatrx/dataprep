library(dataresqc)
library(dplyr)
library(tidyr)

source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/write_sef_f.R')

indir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_orig/Ukraine/'
outdir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/'

# Dnipro -----------------------------------------------------------------
infile <- 'Dnipro.tsv'
df <- read.delim(paste0(indir, infile), header=T, sep='\t', stringsAsFactors = F,
                 skip=12)

read_meta_nonofficial <- function (file = file.choose(), parameter = NULL) {
  header <- read.table(file = file, quote = "", comment.char = "", 
                       sep = "\t", nrows = 12, stringsAsFactors = FALSE, fill = TRUE)
  pars <- c("version", "id", "name", "lat", "lon", "alt", "source", 
            "link", "var", "stat", "units", "meta")
  if (is.null(parameter)) {
    out <- header[, 2]
    names(out) <- pars
  }
  else {
    out <- header[match(parameter, pars), 2]
    names(out) <- parameter
  }
  return(out)
}

meta <- read_meta_nonofficial(paste0(indir,infile))

# Helper: split time column into Hour and Minute
split_hour_minute <- function(time_str) {
  time_parts <- strsplit(time_str, ":")[[1]]
  list(
    hour = as.integer(time_parts[1]),
    minute = as.integer(time_parts[2])
  )
}

# Create temperature dataframe
df.ta.Dnipro <- df %>%
  separate(Hour, into = c("Hour", "Minute"), sep = ":", convert = TRUE) %>%
  mutate(value = T * 1.25, Minute='NA') %>%
  select(Year, Month, Day, Hour, Minute, value)

df.p.Dnipro <- df %>%
  separate(Hour, into = c("Hour", "Minute"), sep = ":", convert = TRUE) %>%
  mutate(
  value = case_when(
      Year >= 1833 & Year <= 1838 ~ P * 33.8639,        # inches → hPa
      Year %in% c(1839:1842, 1850) ~ P * 1.27/(750.06)*1000,       
      # R.s.l. → hPa ccording to (Shostin, 1975;Lamb, 1986) this unit can be converted into millimetres
      # based on the relation 1 R.s.l. = 1.27 mm. That is, 1,000 hPa  = 1,000 mbar = 750.06 mmHg = 590.60 R.s.l
      TRUE ~ NA_real_                                   # fallback if unexpected year
    )
  ) %>%
  select(Year, Month, Day, Hour, Minute, value)

write_sef_f(Data=df.ta.Dnipro,
            outpath=outdir, outfile="Dnipro_ta_subdaily.tsv",
            cod=meta[["id"]],
            variable='ta',
            nam=meta[["name"]],
            lat=meta[["lat"]],
            lon=meta[["lon"]], alt=meta[["alt"]], sou=meta[["source"]], metaHead = 'orig_ta=Reaumur',
            link=meta[["link"]], units='C', stat="point",
            meta="", keep_na = F)

write_sef_f(Data=df.p.Dnipro,
            outpath=outdir, outfile="Dnipro_p_subdaily.tsv",
            cod=meta[["id"]],
            variable='ta',
            nam=meta[["name"]],
            lat=meta[["lat"]],
            lon=meta[["lon"]], alt=meta[["alt"]], sou=meta[["source"]], metaHead = 'orig_p=in(1833-1838);R.s.l.(1839-1842,1850)',
            link=meta[["link"]], units='hPa', stat="point",
            meta="", keep_na = F)

# Kherson -----------------------------------------------------------------

infile <- 'Kherson_pre1825.tsv'
df <- read.delim(paste0(indir, infile), header=T, sep='\t', stringsAsFactors = F,
                 skip=12)

# convert Reaumur to C, in to hPa
df$T_on_N <- ifelse(df$T_on_N == -999.9, NA, df$T_on_N * 1.25)
df$T_on_S <- ifelse(df$T_on_S == -999.9, NA, df$T_on_S * 1.25)
df$P <-      ifelse(df$P == -999.9, NA, df$P * 33.8638)

df <- df %>% filter(P < 3000, P > 950)

# extract hour and minute
df$Minute <- as.integer(sub(".*:","", df$Hour))
df$Hour <- as.integer(sub(":.*", "", df$Hour))

### add 12 days to each date
df$date <- as.Date(sprintf("%04d-%02d-%02d", df$Year, df$Month, df$Day))
df$date <- df$date + 12 # add 12 days

# Update the Year, Month, and Day columns based on the new date
df$Year <- as.integer(format(df$date, "%Y"))
df$Month <- as.integer(format(df$date, "%m"))
df$Day <- as.integer(format(df$date, "%d"))

# read meta from lines
meta_lines <- strsplit(readLines(paste0(indir, infile), n=12), '\t')

meta <- list()
for (line in meta_lines) {
  if (length(line) >= 2 ** nzchar(line[2])) {
    key <- line[1]
    value <- line[2]
    meta[[key]] <- value
  } else { # for lines wo a value (e.g. SEF, ID)
    meta[[line[1]]] <- NA
  }
}

# create two different files for T on N and T on S
df_T_on_N <- df[, c("Year","Month","Day","Hour","Minute","T_on_N")]
colnames(df_T_on_N)[colnames(df_T_on_N)=="T_on_N"] <- "Value"
df_T_on_S <- df[, c("Year","Month","Day","Hour","Minute","T_on_S")]
colnames(df_T_on_S)[colnames(df_T_on_S)=="T_on_S"] <- "Value"

# generate TSV for T_on_N
meta[['Vbl']] <- 'ta'
meta[['Units']] <- 'C'
meta[['metaHead']] <- 'date shifted +12 days | Hours correspond to original day periods morning, midday, evening'

write_sef_f(Data=df_T_on_S,
            outpath=outdir, outfile="Kherson_ta_subdaily_TonS.tsv",
            cod=meta[["ID"]],
            variable=meta[["Vbl"]],
            nam=meta[["Name"]],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]], metaHead = meta[['metaHead']],
            link=meta[["Link"]], units=meta[["Units"]], stat="point",
            meta="orig_ta=Reaumur | T_on_S", keep_na = F)
write_sef_f(Data=df_T_on_N,
            outpath=outdir, outfile="Kherson_ta_subdaily_TonN.tsv",
            cod=meta[["ID"]],
            variable=meta[["Vbl"]],
            nam=meta[["Name"]],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]], metaHead = meta[['metaHead']],
            link=meta[["Link"]], units=meta[["Units"]], stat="point",
            meta="orig_ta=Reaumur | T_on_N", keep_na = F)

# create pressure file
df_p <- df[, c("Year","Month","Day","Hour","Minute","P")]
colnames(df_p)[colnames(df_p)=="P"] <- "Value"

# generate TSV for T_on_N
meta[['Vbl']] <- 'p'
meta[['Units']] <- 'hPa'

write_sef_f(Data=df_p,
            outpath=outdir, outfile="Kherson_p_subdaily",
            cod=meta[["ID"]],
            variable=meta[["Vbl"]],
            nam=meta[["Name"]],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]], metaHead = meta[['metaHead']],
            link=meta[["Link"]], units=meta[["Units"]], stat="point",
            meta="orig_p=inches", keep_na = F)

# Kyiv --------------------------------------------------------------------

infile <- 'Kyiv_ta_pre1837.tsv'
df <- read.delim(paste0(indir, infile), header=T, sep='\t', stringsAsFactors = F,
                 skip=12)

# convert Reamur to Kelvin
df$Value <- ifelse(df$Value == -999.9, NA, df$Value * 1.25)

# extract hour and minute
df$Minute <- as.integer(sub(".*:","", df$Hour))
df$Hour <- as.integer(sub(":.*", "", df$Hour))

# rearrange columns
df <- df[, c("Year","Month","Day","Hour","Minute","Value")]

# read meta for SEF
meta_lines <- strsplit(readLines(paste0(indir, infile), n=12), '\t')

meta <- list()
for (line in meta_lines) {
  if (length(line) >= 2 ** nzchar(line[2])) {
    key <- line[1]
    value <- line[2]
    meta[[key]] <- value
  } else { # for lines wo a value (e.g. SEF, ID)
    meta[[line[1]]] <- NA
  }
}
meta[['Vbl']] <- 'ta'
meta[['Units']] <- 'C'
meta[['metaHead']] <- 'Hours correspond to original day periods morning, midday, evening'

write_sef_f(Data=df,
            outpath=outdir, outfile='Kyiv_ta_subdaily.tsv',
            cod=meta[["ID"]],
            variable=meta[["Vbl"]],
            nam=meta[["Name"]],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]],
            link=meta[["Link"]], units=meta[["Units"]], stat="point",
            meta="orig_ta=Reaumur", keep_na = F)
 