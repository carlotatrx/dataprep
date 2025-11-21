rm(list = ls())
library(glue)
library(readr)   # for read_lines()

# -----------------------------
# PATHS
# -----------------------------
orig_root <- "/scratch3/PALAEO-RA/daily_data/final"
qc_root   <- "/scratch3/PALAEO-RA/daily_data/tmp/sef_tests"

# -----------------------------
# VARIABLES THAT DO NOT HAVE QC FILES
# -----------------------------
no_qc_vars <- c("w", "eee")

# -----------------------------
# OUTPUT LIST OF SUCCESSFUL MATCHES
# -----------------------------
successful_matches <- list()

# -----------------------------
# LIST original FILES
# -----------------------------
orig_files <- list.files(orig_root, pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE)

# -----------------------------
# LOOP THROUGH ORIGINAL FILES
# -----------------------------
for (orig in orig_files) {
  
  fname <- basename(orig)
  
  # Skip QC files
  if (grepl("_qc\\.tsv$", fname)) {
    next
  }
  
  # check variable and skip those that don't have qc's done on them.
  header_lines <- read_lines(orig, n_max=20)
  var_line     <- header_lines[grepl("^Vbl", header_lines)]
  var          <- sub("^Vbl\\s+", "", var_line)
    
  # Construct expected QC filename
  qc_fname <- sub("\\.tsv$", "_qc.tsv", fname)
  qc_path <- file.path(qc_root, qc_fname)
  
  # ----------------------------------
  # SKIP CHECKS OF NON-QC VARS
  # ----------------------------------
  if (var %in% no_qc_vars && !file.exists(qc_path)) {
    # silent skip
    next
  }
  
  
  # ----------------------------------
  # CHECK IF QC FILE IS MISSING
  # ----------------------------------
  if (!file.exists(qc_path)) {
    warning(glue("QC file missing for: {fname}"))
    next
  }
  
  # ----------------------------------
  # READ HEADER LINES FROM BOTH FILES
  # ----------------------------------
  orig_lines <- read_lines(orig, n_max = 50)
  qc_lines   <- read_lines(qc_path, n_max = 50)
  
  # extract SEF header
  orig_header <- orig_lines[1:(grep("^Year", orig_lines)[1] - 1)]
  qc_header   <- qc_lines[1:(grep("^Year", qc_lines)[1] - 1)]
  
  # ----------------------------------
  # CHECK HEADERS MATCH (except META)
  # ----------------------------------
  # Copy headers so we can modify them
  oh <- orig_header[1:length(orig_header)-1]
  qh <- qc_header[1:length(qc_header)-1]
  
  
  # Compare headers
  if (!identical(oh, qh)) {
    warning(glue("Header mismatch for {fname}"))
    next
  }
  
  # ----------------------------------
  # CHECK BODY DIMENSIONS MATCH
  # ----------------------------------
  orig_df <- try(read_tsv(orig, skip = length(orig_header), show_col_types = FALSE), silent = TRUE)
  qc_df   <- try(read_tsv(qc_path, skip = length(qc_header), show_col_types = FALSE), silent = TRUE)
  
  if (inherits(orig_df, "try-error") || inherits(qc_df, "try-error")) {
    warning(glue("Error reading dataframes for {fname}"))
    next
  }
  
  if (!all(dim(orig_df) == dim(qc_df))) {
    warning(glue("Dimension mismatch for {fname}"))
    next
  }
  
  # ----------------------------------
  # IF WE GET HERE → SUCCESS
  # ----------------------------------
  successful_matches[[length(successful_matches) + 1]] <- list(
    file = fname,
    original_path = orig,
    qc_path = qc_path
  )
}

# -----------------------------
# PRINT SUMMARY
# -----------------------------
cat("\n\n### Successfully matched files:\n")
print(succesful_matches)
