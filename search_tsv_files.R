library(dataresqc)

# Define directories
dirs <- c("/home/ccorbella/scratch2_symboliclink/files/1807_USBstick",
          "/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed")

# Initialize lists
ta_files <- c()
p_files <- c()
other_files <- c()

# Loop over directories
for (d in dirs) {
  files <- list.files(d, pattern="\\.tsv$", full.names=TRUE)
  for (f in files) {
    meta <- tryCatch(read_meta(f), error = function(e) NULL)
    if (!is.null(meta)) {
      if (meta["var"] == "ta") {
        ta_files <- c(ta_files, f)
      } else if (meta["var"] == "p") {
        p_files <- c(p_files, f)
      } else {
        other_files <- c(other_files, f)
      }
    }
  }
}

# Result: ta_files and p_files contain full paths
print(ta_files)
print(p_files)
print(other_files)



# extract lat lon from files ----------------------------------------------

# Helper function to extract metadata
get_meta_table <- function(file_list) {
  table <- data.frame()
  for (f in file_list) {
    meta <- tryCatch(read_meta(f), error = function(e) NULL)
    if (!is.null(meta)) {
      row <- data.frame(
        ID = meta["id"],
        Name = meta["name"],
        lat = as.numeric(meta["lat"]),
        lon = as.numeric(meta["lon"]),
        alt = as.numeric(meta["alt"]),
        filepath = f,
        stringsAsFactors = FALSE
      )
      table <- rbind(table, row)
    }
  }
  return(table)
}

ta_info <- get_meta_table(ta_files)
p_info  <- get_meta_table(p_files)

# save
write.csv(ta_info, "/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/ta_latlon.csv", row.names=FALSE, quote=FALSE) 
write.csv(p_info, "/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/p_latlon.csv", row.names=FALSE, quote=FALSE) 


