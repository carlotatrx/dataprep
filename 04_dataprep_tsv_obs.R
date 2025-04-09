library(dataresqc)
library(dplyr)
library(tidyr)

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
      id <- meta["id"]
      name <- meta["name"]
      print(name)
      print("start")
      
      # Custom rename logic for stupidly named stations
      if (name == "Royal Society - Somerset House") {
        id <- paste0(id, " | Royal Society - Somerset House")
        name <- "London"
      } else if (name == "Stockholm Old Astronomical Observatory") {
        id <- paste0(id, " | Stockholm Old Astronomical Observatory")
        name <- "Stockholm"
      } else if (name == "Observatoire") {
        id <- paste0(id, "| Observatoire")
        name <- "Paris"
      } else if (name=='Geneve') {
        name <- 'Geneva'
      }
      
      print(name)
      print("end")
      
      row <- data.frame(
        ID = id,
        Name = name,
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

# get_meta_table <- function(file_list) {
#   table <- data.frame()
#   for (f in file_list) {
#     meta <- tryCatch(read_meta(f), error = function(e) NULL)
#     if (!is.null(meta)) {
#       row <- data.frame(
#         ID = meta["id"],
#         Name = meta["name"],
#         lat = as.numeric(meta["lat"]),
#         lon = as.numeric(meta["lon"]),
#         alt = as.numeric(meta["alt"]),
#         filepath = f,
#         stringsAsFactors = FALSE
#       )
#       table <- rbind(table, row)
#     }
#   }
#   return(table)
# }

ta_info <- get_meta_table(ta_files)
p_info  <- get_meta_table(p_files)

# save
write.csv(ta_info, "/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/ta_latlon.csv", row.names=FALSE, quote=FALSE) 
write.csv(p_info, "/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/p_latlon.csv", row.names=FALSE, quote=FALSE) 


# create df w values ------------------------------------------------------

build_obs_df <- function(file_list) {
  dfs <- list()

  for (f in file_list) {
    meta <- tryCatch(read_meta(f), error = function(e) NULL)
    if (is.null(meta)) next
    
    x <- tryCatch(read_sef(f, all=TRUE), error = function(e) NULL)
    if (is.null(x) || nrow(x) == 0) next
    
    x$Year  <- as.integer(x$Year)
    x$Month <- as.integer(x$Month)
    x$Day   <- as.integer(x$Day)
    
    # Remove flagged values
    ko <- grep("qc", x$Meta)
    # if (length(ko) > 0) x <- x[-ko, ]
    
    # Add Time column and order just in case
    if (all(c("Hour", "Minute") %in% names(x))) {
      x$Time <- x$Hour + x$Minute / 60
    } else {
      x$Time <- NA
    }
    
    x <- x[order(x$Year, x$Month, x$Day, x$Time), ]
    
    # Decide if subdaily
    is_subdaily <- grepl("subdaily", f, ignore.case = TRUE)
    
    if (is_subdaily) {
      # Take first subdaily observation per day
      df_val <- x %>%
        group_by(Year, Month, Day) %>%
        summarise(Value = first(Value), .groups = "drop")
    } else {
      # Assume already daily — just keep columns
      # df_val <- x[, c("Year", "Month", "Day", "Value")]
      df_val <- x %>%
        group_by(Year, Month, Day) %>%
        mutate(n = n()) %>%
        filter(row_number() == 1) %>%
        ungroup() %>%
        select(Year, Month, Day, Value)
      
      dropped_vals <- x %>%
        group_by(Year, Month, Day) %>%
        filter(n() > 1 & row_number() > 1) %>%
        select(Year, Month, Day, Value)
      
      if (nrow(dropped_vals) > 0) {
        cat("\n🗑️  Dropping repeated values in:", colname, "\n")
        print(dropped_vals)
      }
      
    }
    
    # Get column name
    colname <- meta["name"]

    # Apply the same renaming logic as above for stupidly named stations
    if (colname == "Royal Society - Somerset House") {
      colname <- "London"
    } else if (colname == "Stockholm Old Astronomical Observatory") {
      colname <- "Stockholm"
    } else if (colname == "Observatoire") {
      colname <- "Paris"
    } else if (colname=='Geneve') {
      colname <- 'Geneva'
    }
    
    if (is_subdaily) colname <- paste0(colname, "_SUBs")
    
    dupes <- x %>%
      mutate(row_num = row_number()) %>%
      group_by(Year, Month, Day) %>%
      filter(n() > 1) %>%
      arrange(Year, Month, Day, row_num)
    
    if (nrow(dupes) > 0) {
      cat("\n⚠️  Duplicates found in:", colname, "\n")
      print(dupes %>% select(row_num, Year, Month, Day, Value, Meta))
    }

    # Rename and store
    df <- df_val %>%
      rename(!!colname := Value)
    
    dfs[[length(dfs) + 1]] <- df
  }
  
  # Merge all data frames by Year, Month, Day
  if (length(dfs) == 0) return(NULL)
  
  combined <- purrr::reduce(dfs, full_join, by = c("Year", "Month", "Day"))
  combined <- arrange(combined, Year, Month, Day)
  return(combined)
}

# Now apply it
ta_df <- build_obs_df(ta_files)
p_df  <- build_obs_df(p_files)

ta_df <- ta_df %>% filter(Year>=1806, Year<=1850)
p_df <- p_df %>% filter(Year>=1806, Year<=1850)

head(ta_df)
head(p_df)

write.csv(ta_df, "/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/ta_obs.csv", row.names=FALSE, quote=FALSE) 
write.csv(p_df, "/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/data/p_obs.csv", row.names=FALSE, quote=FALSE) 
