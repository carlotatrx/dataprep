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
no_qc_vars <- c("w", "eee", "rrt")

# -----------------------------
# OUTPUT LIST OF SUCCESSFUL MATCHES
# -----------------------------
successful_matches <- list()

# -----------------------------
# LIST original FILES
# -----------------------------
orig_files <- list.files(orig_root, pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE)

# OUTPUT
warnings_file <- "/scratch3/PALAEO-RA/daily_data/tmp/qc_warnings.txt"
warnings_list <- character()


# ----------------------------------
# function to NORMALIZE HEADERS BEFORE COMPARISON
# ----------------------------------
# e.g. Link NA and Link "" should not count as a mismatch

normalize_header <- function(lines) {
  out <- lines
  
  # Trim trailing whitespace on the entire line first
  out <- trimws(out, which = "right")
  
  # Normalize any header of format: Field<TAB>value
  for (i in seq_along(out)) {
    parts <- strsplit(out[i], "\t")[[1]]
    
    # skip if it's not key-value style (e.g., SEF 1.0.0)
    if (length(parts) == 2) {
      key <- trimws(parts[1], which="both")
      val <- trimws(parts[2], which="both")
      
      # Treat NA / N/A / empty as empty
      if (val %in% c("", " ", "NA", "N/A", "na", "Na")) {
        val <- ""
      }
      
      # rebuild line
      out[i] <- paste0(key, "\t", val)
    }
  }
  
  # trim once more for safety
  out <- trimws(out, which="both")
  
  return(out)
}

# function to capture the warnings
withCallingHandlers({
  
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
    # CHECK HEADERS MATCH (except META) and apply normalization
    # ----------------------------------
    # Copy headers so we can modify them
    oh <- normalize_header(orig_header[1:length(orig_header)-1])
    qh <- normalize_header(qc_header[1:length(qc_header)-1])
    
    
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
  
}, warning = function(w) {
  
  # store warnings
  warnings_list <<- c(warnings_list, w$message)
  
  # so it doesn't print to console
  invokeRestart("muffleWarning")
})

if (length(warnings_list) > 0){
  writeLines(warnings_list, warnings_file)
  cat(glue("Saved {length(warnings_list)} warnings to {warnings_file}\n"))
} else {
  cat("No warnings generated.\n")
}

# -----------------------------
# PRINT and SAVE SUMMARY
# -----------------------------

outfile <- "/scratch3/PALAEO-RA/daily_data/tmp/qcd_files.txt"

if (length(successful_matches) > 0) {
  
  # convert list to df
  df <- do.call(rbind, lapply(successful_matches, as.data.frame))
  
  # write to file
  write.table(
    df,
    file = outfile,
    sep  = '\t',
    quote = FALSE,
    row.names = FALSE
  )
  cat(glue("\nSaved {nrow(df)} successful qc'd files to: {outfile}\n"))
  
} else {
  cat("\nNo files were qc'd.\n")
}


# --------------------------------------------------------------------------
# safely replace final files w final qc'd files ----------------------------
# --------------------------------------------------------------------------

non_qcd_dir <- "/scratch3/PALAEO-RA/daily_data/tmp/non_qcd_files"

# write everything down
log_file <- "/scratch3/PALAEO-RA/daily_data/tmp/replacement_log.txt"
log_conn <- file(log_file, open="wt")

log_message <- function(...) {
  msg <- paste(...)
  writeLines(msg, log_conn)
  message(msg)
}

log_message("=== Starting replacement process ===")

# loop through successful matches
for (i in seq_len(nrow(df))) {
  
  orig_file <- df$original_path[i]
  qc_file   <- df$qc_path[i]
  
  orig_dir   <- dirname(orig_file)
  orig_fname <- basename(orig_file)
  qc_fname   <- basename(qc_file)
  
  # new locations
  backup_path <- file.path(non_qcd_dir, orig_fname)
  new_path    <- file.path(orig_dir, qc_fname)
  
  log_message(glue::glue("Processing {fname}"))
  
  # Safety check: original exists
  if (!file.exists(orig_file)) {
    log_message(glue::glue("  ERROR: Original file missing ({orig_file}). Skipping."))
    next
  }
  
  # Move original to backup folder
  ok1 <- file.rename(orig_file, backup_path)
  if (!ok1) {
    log_message(glue::glue("  ERROR: Could not move original to backup: {orig}"))
    next
  }
  
  # Move QC file to original location
  ok2 <- file.rename(qc_file, new_path)
  if (!ok2) {
    log_message(glue::glue("  ERROR: Could not move QC file into place: {qc_file}"))
    # try to restore original back (safety fallback)
    file.rename(backup_path, orig)
    next
  }
  
  log_message(glue::glue("  SUCCESS: Replaced {orig_fname}"))
}



log_message("=== Replacement process completed ===")
close(log_con)

message(glue::glue("Replacement log written to: {log_file}"))








