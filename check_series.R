library(dataresqc)
library(hclim)
source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')

indir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/'
outdir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/sef_tests/'

files <- list.files(indir, pattern = "_corrected\\.tsv$", full.names = TRUE)

# Loop through files
for (file in files[5:6]) {
  # Get just the file name
  file.name <- basename(file)
  
  # Extract base name without _corrected.tsv
  base <- sub("_subdaily_corrected\\.tsv$", "", file.name)
  base <- sub("_corrected\\.tsv$", "", base)  # In case it's not subdaily
  
  cat("\n=== Checking:", file.name, "===\n")
  
  # Run check_sef and print results
  check_result <- check_sef(file)
  print(check_result)
  
  # Run qc and save results
  qc(file, outpath=outdir)

  # Plot decimals and save as PNG
  plot_decimals(file, outfile = paste0(outdir, paste0(base, "_decimals")))

}

# Poltava
infile <- paste0(indir, 'Poltava_ta_subdaily_corrected.tsv')

check_sef(infile)
climatic_outliers(infile, outpath=outdir)
duplicate_columns(infile, outpath=outdir)
duplicate_times(infile, outpath=outdir)
plot_daily(infile, outfile=paste0(outdir,'Poltava_daily'))
plot_decimals(infile, outfile=paste0(outdir,'Poltava_decimals'))
plot_subdaily(infile, outfile=paste0(outdir,'Poltava_subdaily'))
qc(infile, outpath=outdir)

# flags for Kherson


write_flags_f(infile='/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/Kherson_p_subdaily_corrected.tsv',
              qcfile='/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/sef_tests/qc_00033902_p_subdaily.txt',
              outpath='/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/',
              match=F)

# Find problematic rows
df <- read.delim(infile, sep = "\t", skip=12, header = TRUE, stringsAsFactors = FALSE)
df.dailydata <- df[,1:3]
df.dailydata$Value <- df$Value
df.dailydata$var <- 'ta'
daily_repetition(df.dailydata, meta=c("STA001", "40.7128", "-74.0060", "10", "TMAX", "degC"),outpath=outdir) ## dummy meta

str(df)

bad_hour <- which(is.na(df$Hour) | !grepl("^\\d+$", df$Hour))
bad_minute <- which(is.na(df$Minute) | !grepl("^\\d+$", df$Minute))

print(df[bad_hour, ])
print(df[bad_minute, ])
