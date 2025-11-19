library(glue)
library(dataresqc)
source("/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R")

# Define paths
indir <- "/scratch3/PALAEO-RA/daily_data/final"
outdir <- "/scratch3/PALAEO-RA/daily_data/tmp/sef_tests"

# Read file list
files <- readLines("/scratch3/PALAEO-RA/daily_data/tmp/files_no_qc.txt")

# Open connection
log_conn <- file("/scratch3/PALAEO-RA/daily_data/tmp/qc_log.txt", open = "wt")

# Redirect both output and messages
sink(log_conn, type = "output")
sink(log_conn, type = "message")

# Loop over all files
for (f in files[1:length(files)]) {
  infile <- glue("{indir}/{f}")
  basename_noext <- tools::file_path_sans_ext(basename(infile))
  qcfile <- glue("{outdir}/qc_{basename_noext}.txt")
  
  message("\n[INFO] Running QC for: ", basename_noext)
  
  # Step 1: run QC to generate qcfile
  qc(infile, outpath = outdir)

  # Step 2: check if expected QC file was created
  # otherwise make an empty one
  if (!file.exists(qcfile)) {
    message("[INFO] writing empty QC file for: ", basename_noext)
    placeholder <- data.frame(
      Var = character(),
      Year = numeric(),
      Month = numeric(),
      Day = numeric(),
      Hour = numeric(),
      Minute = numeric(),
      Value = numeric(),
      Test = character()
    )
    write.table(placeholder, qcfile, sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  # Step 3: apply flags
  write_flags_f(
    infile = infile,
    qcfile = qcfile,
    outpath = outdir,
    match = FALSE
  )
  
  # Flush after each file to force output to appear in log
  flush(log_conn)
}

sink(type = "message")
sink(type = "output")
close(log_conn)

