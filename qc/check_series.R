rm(list=ls())
library(dataresqc)
library(hclim)
source('/home/ccorbella/scratch2_symboliclink/code/dataprep/helpfun.R')
library(glue)

# bash script to change filenames
# for file in *_corrected.tsv; do [ -e "$file" ] || continue; new_name="${file/_corrected.tsv/.tsv}"; mv "$file" "$new_name"; done

indir <- '/scratch3/PALAEO-RA/daily_data/final/'

### quick one
dirname <- "Marseille/"
filename <- "Marseille_17090105-17090123_ta_daily.tsv"

qc(glue(indir,dirname,filename), outpath=glue(indir, dirname))

qcfilename <- "qc_ACRE_Bessested_p_subdaily.txt"
write_flags_f(infile=glue(indir,dirname,filename),
              qcfile=glue(indir,dirname, qcfilename),
              outpath=glue(indir,dirname),
              match=FALSE)

###

outdir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/sef_tests/'
dir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed'

files <- list.files(indir, pattern = "_daily.tsv$", full.names = TRUE)
files <- files[!grepl("Barcelona_ta_subdaily.tsv|Bologna_rr_subdaily.tsv", files)] # Bcn we don't want now

name <- "Padua"
files <- c(paste0('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/',name,'_17250112-19971130_p_daily.tsv'),
           paste0('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/',name,'_17250112-19970531_ta_daily.tsv'))
files <- paste0('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/Delft_Rijnsburg_17060101-17341231_ta_subdaily.tsv')
print(files)

# Loop through files
for (file in files) {
  # Get just the file name
  file.name <- sub(".tsv$", "", basename(file))

  cat("\n=== Checking:", file.name, "===\n")
  
  # Run check_sef and print results
  check_result <- check_sef(file)
  print(check_result)

  # Plot decimals and save as PNG
  base <- sub("\\.tsv$", "", basename(file))
  
  # Run qc and save results
  qc(file, outpath=outdir)
  
  plot_decimals(file, outfile = paste0(outdir, paste0(base, "_decimals")))

}

###############################################################################################
############################# write flags #####################################################
###############################################################################################

# Bologna
for (var in c('ta', 'p', 'dd')) {
  write_flags_f(infile=glue('{dir}/Bologna_{var}_subdaily.tsv'),
                qcfile=glue('{dir}/sef_tests/qc_Palatine-Society_Bologna_{var}_subdaily.txt'), 
                outpath=dir,
                match=F)
}

# CET
var<- 'ta'
name <- 'CET'
write_flags_f(infile=glue('{dir}/{name}_{var}_daily.tsv'),
              qcfile=glue('{dir}/sef_tests/qc_Had{name}_{var}_daily.txt'), # station ID HadCET
              outpath=dir,
              match=F)

# Delemont-Delsberg
for (var in c('ta', 'p')) {
  write_flags_f(infile=glue('{dir}/Delemont_{var}_subdaily.tsv'),
                qcfile=glue('{dir}/sef_tests/qc_JU01_Delemont_{var}_subdaily.txt'), 
                outpath=dir,
                match=F)
}

# Delft
write_flags_f(infile=glue('{dir}/Delft_Rijnsburg_17060101-17341231_ta_subdaily.tsv'),
              qcfile=glue('{dir}/sef_tests/qc_KNMI-49_Delfft_Rijnsburg_ta_subdaily.txt'), 
              outpath=dir,
              match=F)

# Domodossola
write_flags_f(infile=glue('{dir}/Domodossola_ta_daily.tsv'),
              qcfile=glue('{dir}/sef_tests/qc_Domodossola_ta_daily.txt'), 
              outpath=dir,
              match=F)

# Dnipro
for (var in c('ta', 'p')) {
  write_flags_f(infile=glue('{dir}/Dnipro_{var}_subdaily.tsv'),
                qcfile=glue('{dir}/sef_tests/qc_Dnipro_{var}_subdaily.txt'), # station ID 00034504
                outpath=dir,
                match=F)
}

# Florence
for (var in c('ta', 'p')) {
  write_flags_f(infile=glue('{dir}/Dnipro_{var}_subdaily.tsv'),
                qcfile=glue('{dir}/sef_tests/qc_Dnipro_{var}_subdaily.txt'),
                outpath=dir,
                match=F)
}

# Frankfurt
# molt trist, només sabem les unitats de la direcció del vent
for (var in c('dd','p', 'ta')) {
  write_flags_f(
    infile=glue('{dir}/Frankfurt_{var}_subdaily.tsv'),
    qcfile=glue('{dir}/sef_tests/qc_Frankfurt_{var}_subdaily.txt'),
    outpath=dir,
    match=F
  )
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

# Montpellier
for (var in c('ta', 'p', 'dd')) {
  write_flags_f(infile=glue('{dir}/Montpellier_{var}_subdaily.tsv'),
                qcfile=glue('{dir}/sef_tests/qc_Montpellier_{var}_subdaily.txt'),
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

# Padua
var<- 'ta'
name <- 'Padua'
write_flags_f(infile=glue('{dir}/{name}_17250112-17731231_ta_daily.tsv'),
              qcfile=glue('{dir}/sef_tests/qc_{name}_{var}_daily.txt'),
              outpath=dir,
              match=F)

write_flags_f(infile=glue('{dir}/{name}_17250112-19970531_ta_daily.tsv'),
              qcfile=glue('{dir}/sef_tests/qc_IMPROVE_{name}_ta_daily.txt'),
              outpath=dir,
              match=F)

write_flags_f(infile=glue('{dir}/{name}_17250112-19971130_p_daily.tsv'),
              qcfile=glue('{dir}/sef_tests/qc_IMPROVE_{name}_p_daily.txt'),
              outpath=dir,
              match=F)


# Piacenza
for (var in c('ta','rr')) {
  write_flags_f(infile=glue('{dir}/Piacenza_{var}_daily.tsv'),
                qcfile=glue('{dir}/sef_tests/qc_Piacenza_{var}_daily.txt'),
                outpath=dir,
                match=F)
}

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
