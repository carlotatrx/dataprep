library(dataresqc)
library(hclim)
source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')
library(glue)

# bash script to change filenames
# for file in *_corrected.tsv; do [ -e "$file" ] || continue; new_name="${file/_corrected.tsv/.tsv}"; mv "$file" "$new_name"; done

indir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed'
outdir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/sef_tests/'

files <- list.files(indir, pattern = "_subdaily.tsv$", full.names = TRUE)
files <- files[!grepl("Barcelona_ta_subdaily.tsv|Bologna_rr_subdaily.tsv", files)] # Bcn we don't want now

files <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/Piacenza_ta_daily.tsv'
print(files)

# Loop through files
for (file in files) {
  # Get just the file name
  file.name <- sub(".tsv$", "", basename(file))

  cat("\n=== Checking:", file.name, "===\n")
  
  # Run check_sef and print results
  check_result <- check_sef(file)
  print(check_result)
  
  # Run qc and save results
  qc(file, outpath=outdir)

  # Plot decimals and save as PNG
  base <- gsub(".tsv", "", basename(file))
  plot_decimals(file, outfile = paste0(outdir, paste0(base, "_decimals")))

}

###############################################################################################
############################# write flags #####################################################
###############################################################################################

dir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed'

# Bologna
for (var in c('ta', 'p', 'dd')) {
  write_flags_f(infile=glue('{dir}/Bologna_{var}_subdaily.tsv'),
                qcfile=glue('{dir}/sef_tests/qc_Palatine-Society_Bologna_{var}_subdaily.txt'), # station ID 00034504
                outpath=dir,
                match=F)
}


# Delemont-Delsberg
for (var in c('ta', 'p')) {
  write_flags_f(infile=glue('{dir}/Delemont_{var}_subdaily.tsv'),
                qcfile=glue('{dir}/sef_tests/qc_JU01_Delemont_{var}_subdaily.txt'), # station ID 00034504
                outpath=dir,
                match=F)
}

# Dnipro
for (var in c('ta', 'p')) {
  write_flags_f(infile=glue('{dir}/Dnipro_{var}_subdaily.tsv'),
                qcfile=glue('{dir}/sef_tests/qc_Dnipro_{var}_subdaily.txt'), # station ID 00034504
                outpath=dir,
                match=F)
}

# Kamyantes-podilskyi
for (var in c('ta', 'p')) {
  write_flags_f(infile=glue('{dir}/Kamyanets_{var}_subdaily.tsv'),
                qcfile=glue('{dir}/sef_tests/qc_00033548_{var}_subdaily.txt'),
                outpath=dir,
                match=F)
}

# Kharkiv
for (var in c('ta', 'p')) {
  write_flags_f(infile=glue('{dir}/Kharkiv_{var}_subdaily.tsv'),
                qcfile=glue('{dir}/sef_tests/qc_Kharkiv_{var}_subdaily.txt'),
                outpath=dir,
                match=F)
}


# Kherson
for (var in c('ta', 'p')) {
  write_flags_f(infile=glue('{dir}/Kherson_{var}_subdaily.tsv'),
                qcfile=glue('{dir}/sef_tests/qc_00033902_{var}_subdaily.txt'),
                outpath=dir,
                match=F)
}

# Kyiv
for (var in c('ta', 'p')) {
  write_flags_f(infile=glue('{dir}/Kyiv_{var}_subdaily.tsv'),
                qcfile=glue('{dir}/sef_tests/qc_00033345_{var}_subdaily.txt'),
                outpath=dir,
                match=F)
}

# London
write_flags_f(infile=glue('{dir}/London_17870101-18221231_p_subdaily.tsv'),
              qcfile=glue('{dir}/sef_tests/qc_london-11_p_subdaily.txt'),
              outpath=dir,
              match=F)

write_flags_f(infile=glue('{dir}/London_18230101-18411231_p_subdaily.tsv'),
              qcfile=glue('{dir}/sef_tests/qc_london-12_p_subdaily.txt'),
              outpath=dir,
              match=F)

# Lugansk
for (var in c('ta', 'p')) {
  write_flags_f(infile=glue('{dir}/Lugansk_{var}_subdaily.tsv'),
                qcfile=glue('{dir}/sef_tests/qc_Lugansk_{var}_subdaily.txt'),
                outpath=dir,
                match=F)
}

# Mulhouse
for (var in c('ta', 'p')) {
  write_flags_f(infile=glue('{dir}/Mulhouse_{var}_subdaily.tsv'),
                qcfile=glue('{dir}/sef_tests/qc_Mulhouse_{var}_subdaily.txt'),
                outpath=dir,
                match=F)
}

# Odesa
for (var in c('ta', 'p')) {
  write_flags_f(infile=glue('{dir}/Odesa_{var}_subdaily.tsv'),
                qcfile=glue('{dir}/sef_tests/qc_Odesa_{var}_subdaily.txt'),
                outpath=dir,
                match=F)
}

# Piacenza
write_flags_f(infile=glue('{dir}/Piacenza_ta_daily.tsv'),
              qcfile=glue('{dir}/sef_tests/qc_Piacenza_ta_daily.txt'),
              outpath=dir,
              match=F)

# Poltava
for (var in c('ta', 'p')) {
  write_flags_f(infile=glue('{dir}/Poltava_{var}_subdaily.tsv'),
                qcfile=glue('{dir}/sef_tests/qc_Poltava_{var}_subdaily.txt'),
                outpath=dir,
                match=F)
}


###############################################################################################

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
