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

#' Write Data in SEF Format
#'
#' This function writes meteorological or climate data in the `SEF` (Structured Environmental Format) version 1.0.0.
#'
#' @param Data A dataframe containing the data to be written. It must have columns: Year, Month, Day, Hour, Minute, and Value.
#' @param outpath A character string specifying the directory where the file will be saved.
#' @param outfile A character string specifying the filename (without extension) under which the file will be saved.
#' @param variable A character string indicating the variable name (e.g., "temperature").
#' @param cod A character string representing the station ID.
#' @param nam A character string for the station name (default: "").
#' @param lat A numeric value specifying latitude (default: "").
#' @param lon A numeric value specifying longitude (default: "").
#' @param alt A numeric value specifying altitude (default: "").
#' @param sou A character string for the data source (default: "").
#' @param link A character string for a related URL or link (default: "").
#' @param units A character string indicating the measurement units (e.g., "Celsius").
#' @param stat A character string indicating the statistic type (e.g., "mean").
#' @param metaHead A character string for metadata information (default: "").
#' @param meta A character string for additional metadata (default: "").
#' @param period A character string indicating the data period (default: "").
#' @param time_offset A numeric value representing time offset in hours (default: 0).
#' @param note A character string for any additional notes (default: "").
#' @param keep_na A logical value indicating whether to keep `NA` values in the output file (default: FALSE).
#'
#' @return Writes a `.tsv` file in the specified `outpath` and prints a message confirming the file was written.
#' @export
#'
#' @examples
#' Data <- data.frame(
#'     Year = c(2023, 2023),
#'     Month = c(3, 3),
#'     Day = c(15, 16),
#'     Hour = c(12, 14),
#'     Minute = c(30, 45),
#'     Value = c(25.5, 26.0)
#' )
#' 
#' write_sef(Data, "output_directory/", "output_file", "temperature", "001", 
#'           "Sample Station", "40.7128", "-74.0060", "10", "ExampleSource", 
#'           "http://example.com", "Celsius", "mean", period="day", meta="Metadata", 
#'           note="sample", keep_na=TRUE)

write_sef_f <- function(Data, outpath, outfile, variable, cod, nam = "", lat = "", lon = "", 
                        alt = "", sou = "", link = "", units, stat, metaHead = "", 
                        meta = "", period = "", time_offset = 0, note = "", keep_na = FALSE) 
{
  for (i in 1:ncol(Data)) Data[, i] <- as.character(Data[, i])
  
  header <- array(dim = c(12, 2), data = "")
  header[1, ] <- c("SEF", "1.0.0")
  header[2, ] <- c("ID", trimws(as.character(cod)))
  header[3, ] <- c("Name", trimws(as.character(nam)))
  header[4, ] <- c("Lat", trimws(as.character(lat)))
  header[5, ] <- c("Lon", trimws(as.character(lon)))
  header[6, ] <- c("Alt", trimws(as.character(alt)))
  header[7, ] <- c("Source", trimws(as.character(sou)))
  header[8, ] <- c("Link", trimws(as.character(link)))
  header[9, ] <- c("Vbl", trimws(as.character(variable)))
  header[10, ] <- c("Stat", trimws(as.character(stat)))
  header[11, ] <- c("Units", trimws(as.character(units)))
  header[12, ] <- c("Meta", trimws(as.character(metaHead)))
  
  if (stat == "point" & !all(as.character(period) == "0")) {
    period <- "0"
    warning("Period forced to 0 because of 'stat'")
  }
  
  if (!all(time_offset == 0) & !all(is.na(as.integer(Data[, 4]) + as.integer(Data[, 5])))) {
    times <- ISOdate(Data[, 1], Data[, 2], Data[, 3], Data[, 4], Data[, 5])
    times <- times - time_offset * 3600
    Data[which(!is.na(times)), 1] <- as.integer(substr(times[which(!is.na(times))], 1, 4))
    Data[which(!is.na(times)), 2] <- as.integer(substr(times[which(!is.na(times))], 6, 7))
    Data[which(!is.na(times)), 3] <- as.integer(substr(times[which(!is.na(times))], 9, 10))
    Data[which(!is.na(times)), 4] <- as.integer(substr(times[which(!is.na(times))], 12, 13))
    Data[which(!is.na(times)), 5] <- as.integer(substr(times[which(!is.na(times))], 15, 16))
  }
  
  DataNew <- data.frame(
    Year = Data[, 1], Month = Data[, 2], Day = Data[, 3], 
    Hour = Data[, 4], Minute = Data[, 5], Period = as.character(period), 
    Value = Data[, 6], Meta = as.character(meta), stringsAsFactors = FALSE
  )
  
  if (!keep_na) {
    DataNew <- DataNew[which(!is.na(DataNew$Value)), ]
  }
  
  if (substr(outpath, nchar(outpath), nchar(outpath)) != "/") {
    outpath <- paste0(outpath, "/")
  }
  
  # Ensure outfile is specified and formatted correctly
  if (missing(outfile) || is.na(outfile) || outfile == "") {
    stop("You must specify a valid filename in 'outfile'.")
  }
  
  if (!grepl("\\.tsv$", outfile)) {
    outfile <- paste0(outfile, ".tsv")
  }
  
  filename <- paste0(outpath, outfile)
  
  write.table(header, file = filename, quote = FALSE, row.names = FALSE, 
              col.names = FALSE, sep = "\t", dec = ".", fileEncoding = "UTF-8")
  
  write.table(t(names(DataNew)), file = filename, quote = FALSE, 
              row.names = FALSE, col.names = FALSE, sep = "\t", fileEncoding = "UTF-8", append = TRUE)
  
  write.table(DataNew, file = filename, quote = FALSE, row.names = FALSE, 
              col.names = FALSE, sep = "\t", dec = ".", fileEncoding = "UTF-8", append = TRUE)
  
  message(paste("Data written to file", filename))
}



write_flags_f <- function (infile, qcfile, outpath, note = "", match = TRUE)
{
  library(dataresqc)
  Data <- read_sef(infile, all = TRUE)
  header <- read.table(file = infile, quote = "", comment.char = "", 
                       sep = "\t", nrows = 12, stringsAsFactors = FALSE, fill = TRUE)
  header[which(is.na(header[, 2]) & !header[, 1] %in% c("Lat", 
                                                        "Lon", "Alt")), 2] <- ""
  vbl <- read_meta(infile, "var")
  uts <- read_meta(infile, "units")
  Data$Value <- dataresqc:::check_units(Data$Value, vbl, uts)
  if (vbl %in% c("ta", "tb", "td", "t_air", "t_wet", "t_dew", 
                 "Tx", "Tn", "dep_dew", "ibt", "atb", "Txs", "TGs", "Tns", 
                 "TGn", "t_snow", "Ts", "t_water")) {
    uts <- "C"
  }
  else if (vbl %in% c("p", "mslp", "pppp")) {
    uts <- "hPa"
  }
  else if (vbl %in% c("rr", "sw", "rrls")) {
    uts <- "mm"
  }
  else if (vbl %in% c("sd", "fs")) {
    uts <- "cm"
  }
  else if (vbl == "w") {
    uts <- "m/s"
  }
  else if (vbl == "rh") {
    uts <- "%"
  }

  flags <- read.table(qcfile, stringsAsFactors = FALSE, header = TRUE, sep = "\t")
  colnames(flags) <- c("Var", "Year", "Month", "Day", "Hour", "Minute", "Value", "Test")
  if (flags$Var[1] != Data$Var[1]) stop("Variable mismatch")
  
  if (match) {
    dates <- paste(flags$Year, flags$Month, flags$Day, flags$Hour, flags$Minute, flags$Value)
    i <- which(paste(Data$Year, Data$Month, Data$Day, Data$Hour, Data$Minute, Data$Value) %in% dates)
  } else {
    dates <- paste(flags$Year, flags$Month, flags$Day, flags$Hour, flags$Minute)
    i <- which(paste(Data$Year, Data$Month, Data$Day, Data$Hour, Data$Minute) %in% dates)
  }

  if (length(i) > 0) {
    if (length(i) < nrow(flags)) {
      if (dim(flags)[2] == 6) {
        if (match) {
          k <- which(!dates %in% paste(Data$Year, Data$Month, 
                                       Data$Day, Data$Value))
        }
        else {
          k <- which(!dates %in% paste(Data$Year, Data$Month, 
                                       Data$Day))
        }
      }
      else if (dim(flags)[2] == 8) {
        if (match) {
          k <- which(!dates %in% paste(Data$Year, Data$Month, 
                                       Data$Day, Data$Hour, Data$Minute, Data$Value))
        }
        else {
          k <- which(!dates %in% paste(Data$Year, Data$Month, 
                                       Data$Day, Data$Hour, Data$Minute))
        }
        warning("Flags for some observations could not be written.")
        flags <- flags[-k, ]
        flags <- flags[seq_len(nrow(flags)), , drop = FALSE]  # RESET ROW INDEXES
      }
    }
    Data$Meta[i] <- paste0(Data$Meta[i], " | qc=", flags$Test)
    Data$Meta <- gsub("^\\|", "", Data$Meta)
    meta_string <- paste0("QC software=dataresqc v", packageVersion("dataresqc"))
    if (header[12, 2] == "") {
      header[12, 2] <- meta_string
    }
    else {
      header[12, 2] <- paste(header[12, 2], meta_string, 
                             sep = " | ")
    }
    filename <- paste0(sub("\\.tsv","",basename(infile)), "_qc")
    if (note != "") {
      filename <- paste0(tools::file_path_sans_ext(filename), "_", gsub(" ", "_", note), ".tsv")
    }
    
    write_sef(Data = Data[, c("Year", "Month", "Day", "Hour", "Minute", "Period", "Value")], outpath = outpath, 
              variable = vbl, cod = header[2, 2], nam = header[3, 2], lat = header[4, 2], lon = header[5, 2],
              alt = header[6, 2], sou = header[7, 2], link = header[8, 2], 
              stat = header[10, 2], units = uts, metaHead = header[12, 2],
              meta = Data[, 9], period = Data[, 7], outfile = filename, 
              keep_na = TRUE)
  }
  else warning("No matches found: possibly incorrect input files")
}
