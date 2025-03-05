##
indir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_orig/Ukraine/'
outdir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/Ukraine/'


# Kherson -----------------------------------------------------------------

infile <- 'Kherson_pre1825.tsv'
df <- read.delim(paste0(indir, infile), header=T, sep='\t', stringsAsFactors = F,
                 skip=12)

# convert Reaumur to K, in to hPa
df$T_on_N <- ifelse(df$T_on_N == -999.9, NA, df$T_on_N * 1.25 + 273.13)
df$T_on_S <- ifelse(df$T_on_S == -999.9, NA, df$T_on_S * 1.25 + 273.13)
df$P <-      ifelse(df$P == -999.9, NA, df$P * 33.8638)

# extract hour and minute
df$Minute <- as.integer(sub(".*:","", df$Hour))
df$Hour <- as.integer(sub(":.*", "", df$Hour))

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
meta[['ID']] <- 'unknown_code'

# create two different files for T on N and T on S
df_T_on_N <- df[, c("Year","Month","Day","Hour","Minute","T_on_N")]
colnames(df_T_on_N)[colnames(df_T_on_N)=="T_on_N"] <- "Value"
df_T_on_S <- df[, c("Year","Month","Day","Hour","Minute","T_on_S")]
colnames(df_T_on_S)[colnames(df_T_on_S)=="T_on_S"] <- "Value"

# generate TSV for T_on_N
meta[['Vbl']] <- 'ta'
meta[['Units']] <- 'K'
meta[['metaHead']] <- 'Hours correspond to original day periods morning, midday, evening'

write_sef(Data=df_T_on_S,
          outpath=outdir,
          cod=meta[["ID"]],
          variable=meta[["Vbl"]],
          nam=meta[["Name"]],
          lat=meta[["Lat"]],
          lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]], metaHead = meta[['metaHead']],
          link=meta[["Link"]], units=meta[["Units"]], stat="point",
          meta="orig_ta=Reaumur | T_on_S", keep_na = F)
write_sef(Data=df_T_on_N,
          outpath=outdir,
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

write_sef(Data=df_p,
          outpath=outdir,
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
df$value <- ifelse(df$value == -999.9, NA, df$value * 1.25 + 273.13)

# extract hour and minute
df$Minute <- as.integer(sub(".*:","", df$Hour))
df$Hour <- as.integer(sub(":.*", "", df$Hour))

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
meta[['ID']] <- 'unknown_code'
meta[['Vbl']] <- 'ta'
meta[['Units']] <- 'K'
meta[['metaHead']] <- 'Hours correspond to original day periods morning, midday, evening'

final_df <- df[, c("Year","Month","Day","Hour","Minute","Value")]
write_sef(Data=final_df,
          outpath=outdir,
          cod=meta[["ID"]],
          variable=meta[["Vbl"]],
          nam=meta[["Name"]],
          lat=meta[["Lat"]],
          lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]],
          link=meta[["Link"]], units=meta[["Units"]], stat="point",
          meta="orig_ta=Reaumur", keep_na = F)
 