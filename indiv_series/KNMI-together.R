###############################################################################
#  KNMI Historical Subdaily Observations → SEF Conversion
#  Author:   Carlota Corbella
#  Date:     2025-10-13
#
#    This script converts early KNMI “antieke waarnemingen” station records
#    as obtained in https://www.knmi.nl/nederland-nu/klimatologie/daggegevens/antieke-waarnemingen
#    into the WMO Station Exchange Format (SEF 1.0.0).
#
#  Overview of workflow:
#    1.  Load libraries and helper functions.
#    2.  Define generic utilities for:
#           – Wind-direction normalization (Dutch → English 16-point).
#           – Unit conversions (°F →°C, pressure (Rijnlandse inches, mmHg, etc.) → hPa,
#             precipitation → mm).
#           – Time stamping from either observation number or clock time.
#    3.  Read and standardize raw KNMI files using station-specific “skip”
#        line counts and file orders.
#    4.  Compute physical variables (ta, p, dd, rr, rh) and associated
#        metadata for each record.
#    5.  Write results to SEF files, one per variable, per station, with
#        automatic header metadata and station-level notes.
#
#  Features:
#    • Handles both Fahrenheit and Celsius temperature scales.
#    • Converts historical pressure units to modern hPa using station altitude 
#      and temperature corrections (convert_pressure function from dataresqc).
#    • Supports subdaily data either by “obsnum” (1=08h, 2=13h, 3=22h)
#      or explicit hour/minute columns (hh/hhmm).
#    • Normalizes Dutch wind-direction codes (e.g. “OZO” → “ESE”).
#    • Allows station-specific metaHead fields:
#          meta_head_p   → added only for pressure (“p”)
#          meta_head_all → added for all variables (e.g. altitude, observer)
#    • Includes time_offset only when the source file provides explicit obs
#      times (hh or hhmm).
#
#  Inputs:
#    – Station registry (tibble “stations”) defining per-station metadata:
#        name, code, lat, lon, alt, path, skip lines, file pattern/order,
#        variable scaling and pressure mode, meta_head_p, meta_head_all.
#    – Raw KNMI text files in /scratch3/PALAEO-RA/daily_data/original/<station>
#
#  Outputs:
#    – SEF 1.0.0 files in OUTDIR (/scratch3/PALAEO-RA/daily_data/final/)
#      organized by station and variable, e.g.:
#         KNMI2/Leiden/Leiden_1740-1742_ta_subdaily.sef
#         KNMI2/Haarlem/Haarlem_1735-1742_p_subdaily.sef
#
#  Dependencies:
#    Packages: lubridate, dplyr, readr, tidyr, purrr, dataresqc
#    Dependencies: helpfun.R
#
#  Notes:
#    Some stations have to be processed manually, e.g. Utrecht1836_46.txt
#    and other Utrechts (I already have the code and don't want to automatize
#    it here because it's time-consuming as it's the only station with this
#    much information and daily data such as tmin, tmax).
#
###############################################################################

rm(list = ls())
library(lubridate)
library(dplyr)
library(readr)
library(tidyr)
library(dataresqc)
library(purrr)
source("/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R")

source <- "KNMI"
link   <- "https://www.knmi.nl/nederland-nu/klimatologie/daggegevens/antieke-waarnemingen"

STATION_NAME <- "Maastricht" # in case you only want to process 1 station
OUTDIR <- "/scratch3/PALAEO-RA/daily_data/final/"

time.offset <- function(lon) {as.numeric(lon)*12/180}

outfile.name <- function(name, var, df, subdaily=TRUE) {
  subdaily.str <- ifelse(subdaily, "subdaily", "daily")
  df.short <- df[!is.na(df[[var]]), , drop=FALSE]
  paste0(name,"_",
         get_date_range(if (nrow(df.short)) df.short else df),
         "_", var,"_", subdaily.str)
}

units <- function(var) {
  case_when(
    var == "ta" ~ "C",
    var == "p"  ~ "hPa",
    var == "dd" ~ "deg",
    var == "rr" ~ "mm",
    var == "rh" ~ "perc",
    TRUE ~ "unknown"
  )
}

# wind conversion
directions <- c("N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
                "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW")

dd_normalize <- function(x) {
  y <- toupper(trimws(x))
  y <- recode(y,
              "NWE"  = "WNW",
              "NWW"  = "WNW",
              "SWW"  = "WSW",
              "WENW" = "WNW",
              "ENW"  = "ENE",
              "NEE"  = "ENE",
              "NN0"  = "NNW",   # OCR '0' -> 'W'
              "NNO"  = "NNW",   # 'O' (Oeste) -> W
              "SWS"  = "SSW",
              "ES"   = "SE",
              "SES"  = "SE",
              "NSW"  = "SSW",
              "C"    = "calm",
              .default = y,
              .missing = NA_character_
  )
  ifelse(y %in% c(directions, "calm"), y, NA_character_)
}
dd2deg <- function(x) {
  deg <- rep(NA_character_, length(x))
  i_dir <- !is.na(x) & x %in% directions
  i_calm <- !is.na(x) & x == "calm"
  deg[i_dir] <- as.character(22.5 * (match(x[i_dir], directions) - 1))
  deg[i_calm] <- "calm"
  deg
}

dir_deg <- setNames(seq(0, 337.5, by=22.5), directions)
# base Dutch -> English mapping (no 'T')
dutch_map <- c(
  "N"="N", "NO"="NE", "NNO"="NNE",
  "ONO"="ENE", "O"="E", "OZO"="ESE",
  "ZO"="SE", "ZZO"="SSE", "Z"="S",
  "ZZW"="SSW", "ZW"="SW", "WZW"="WSW",
  "W"="W", "WNW"="WNW", "NW"="NW", "NNW"="NNW",
  # Iberian variants
  "O"="E", "NO"="NE", "SO"="SE", "Z0"="S", "N0"="N",
  # common typos
  "WZW"="WSW", "WZW"="WSW"
)

# helper: English 16-pt to degrees and back
to_deg <- function(code) dir_deg[[code]]
nearest_dir <- function(deg) names(dir_deg)[which.min(abs(((deg - unname(dir_deg) + 180) %% 360) - 180))]

# Normalize one token (Dutch → English 16-pt or "calm"/NA), handling "T"
norm_dutch_token <- function(tok) {
  if (is.na(tok) || !nzchar(tok)) return(NA_character_)
  y <- toupper(trimws(tok))
  if (y %in% c("C","CALM","V")) return("calm")   # C=Calm (Calmo/Calma), V=Variable → keep as calm label
  # direct map first
  if (y %in% names(dutch_map)) return(dutch_map[[y]])
  # if contains 'T' (range like "ZWTZ" = SW to S, "NOTN" = NE to N)
  # or also lower case (e.g. den Helder, "ZWtZ")
  if (grepl("T", y, ignore.case = TRUE)) {
    parts <- strsplit(toupper(y), "T", fixed = TRUE)[[1]]
    if (length(parts) != 2) return(NA_character_)
    left  <- dutch_map[toupper(parts[1])]
    right <- dutch_map[toupper(parts[2])]
    if (any(is.na(c(left, right)))) return(NA_character_)
    la <- to_deg(left); ra <- to_deg(right)
    # angular midpoint (shortest arc)
    diff <- ((ra - la + 540) %% 360) - 180
    mid  <- (la + diff/2) %% 360
    return(nearest_dir(mid))
  }
  # other single Dutch combos like "ONO","OZO","ZZO","ZZW","NNO" covered above;
  # handle simple two-letter Dutch not in map by translating letters and retry
  y2 <- chartr("OZ", "ES", y) # O->E, Z->S
  if (y2 %in% directions) return(y2)
  NA_character_
}

# Vectorized normalizer
dd_normalize_nl <- function(x) vapply(x, norm_dutch_token, character(1))

# check that array is sorted by date
date_sorted_check <- function(x, nm = "Date") stopifnot(all(x[[nm]] == sort(x[[nm]])))


# ---- 1) A single reader that normalizes column names -------------------------
read_knmi_files <- function(files, skip, name_repair = TRUE) {
  cat("\n  • Reading", length(files), "file(s) with skip =", skip, "\n")
  read_one <- function(f) {
    cat("    →", basename(f), "\n")
    df <- tryCatch(
      read.csv(f, skip = skip, header = TRUE, fill = TRUE,
               strip.white = TRUE, stringsAsFactors = FALSE),
      error = function(e) {
        cat("      !! Failed to read", f, ":", e$message, "\n")
        return(data.frame())
      }
    )

    
    if (nrow(df) == 0) return(df)
    
    orig_names <- names(df) # keep just in case
    
    up <- toupper(names(df))
    rename_map <- c(
      "STN"     = "station",
      "YYYYMMDD"= "Date", "DATE"="Date", "DATUM"="Date",
      
      # time columns
      "M"       = "obsnum",
      "HH"      = "hh",
      "HHMM"    = "hhmm",
      
      # core variables
      "P"       = "porig",
      "T"       = "taorig",
      "DD"      = "ddorig",
      "R"       = "rrorig",
      "U"       = "rh",
      "RH"      = "rh",
      "W2"      = "w2",
      "WW"      = "ww",
      "FH"      = "fh",
      
      # present in long format only
      "TN"      = "tmin",
      "TX"      = "tmax",
      "N"       = "n",
      "FD"      = "fd",
      "FK"      = "fk",
      "FX"      = "fx"
    )
    
    # map original column names to new standard
    new_names <- names(df)
    for (i in seq_along(up)) {
      key <- up[i]
      if (key %in% names(rename_map)) new_names[i] <- rename_map[[key]]
    }
    names(df) <- new_names
    
    # make sure's there's obs time col
    for (nm in c("obsnum","hh","hhmm")) if (!nm %in% names(df)) df[[nm]] <- NA

    df
  }
  
  out <- purrr::map(files, read_one) %>% list_rbind()
  # cat("  → Combined rows:", nrow(out), " cols:", ncol(out), "\n")
  out
}

# ---- 2) Time stamping that works for obsnum OR hh/hhmm -----------------------
stamp_time <- function(df) {
  
  # ensure there're always to columns to pick from
  for (nm in c("obsnum","hh","hhmm")) if (!nm %in% names(df)) df[[nm]] <- NA_integer_
  
  df %>%
    mutate(
      Date  = ymd(as.numeric(Date)),
      Year  = year(Date),
      Month = month(Date),
      Day   = day(Date),
      Hour  = dplyr::coalesce(
        # classic 3 obs/day via obsnum
        case_when(
          obsnum == 1 ~ 8L,
          obsnum == 2 ~ 13L,
          obsnum == 3 ~ 22L,
          TRUE ~ NA_integer_
        ),
        # hh like 8, 13, 22
        suppressWarnings(as.integer(hh) %% 24L),
        # hhmm like 800, 1300, 2200
        suppressWarnings(as.integer(hhmm) %/% 100L)
      ),
      Minute = 0L,
      meta.time = ifelse(!is.na(obsnum),
                         paste0("obs.num=", obsnum),
                         paste0("orig.time=", sprintf("%02d", Hour), ":00"))
    ) %>%
    arrange(Date, Hour, Minute) %>%
    { date_sorted_check(., "Date"); . }
}

# ---- 3) Variable computations with station-level knobs -----------------------
# knobs:
#   ta_scale: "F10" (tenths of Fahrenheit), "C10" (tenths of Celsius), "C1" (already °C), "MIX_1836C10_else_C10" (your Utrecht special)
#   p_mode: "rijnlandse_5digits" (ii.ll.q packed), "mmHg10", "mmHg100", NA (no p)
#   rr_factor: multiply rrorig by this (e.g., 0.22, 0.1) or NA (no rr)
compute_vars <- function(df, lat, lon, alt,
                         ta_scale = "F10",
                         p_mode   = NA_character_,
                         rr_factor = NA_real_) {
  
  # temperature
  ta <- case_when(
    ta_scale == "F10" ~ round((as.numeric(df$taorig)/10 - 32) * 5/9, 1),
    ta_scale == "F1"  ~ round((as.numeric(df$taorig) - 32) * 5/9, 1),
    ta_scale == "C10" ~ round(as.numeric(df$taorig)/10, 1),
    ta_scale == "C1"  ~ round(as.numeric(df$taorig), 1),
    ta_scale == "MIX_1836F10_else_C10" ~ round(if_else(df$Year == 1836,
                                                       (as.numeric(df$taorig)/10 - 32) * 5/9,
                                                       as.numeric(df$taorig)/10), 1),
    TRUE ~ NA_real_
  )
  
  # wind
  dd_norm <- dd_normalize_nl(df$ddorig)
  stopifnot(all(dd_norm[!is.na(dd_norm)] %in% c(directions, "calm")))
  dd <- dd2deg(dd_norm)
  
  # precipitation
  rr <- if (!is.na(rr_factor) && "rrorig" %in% names(df)) round(as.numeric(df$rrorig) * rr_factor, 1) else NULL
  
  # relative humidity
  rh <- if ("rh" %in% names(df)) as.numeric(df$rh) else NULL
  
  # pressure
  p <- NULL; meta.p <- NULL
  if (!is.na(p_mode) && "porig" %in% names(df)) {
    if (p_mode == "rijnlandse_5digits") {
      s  <- sprintf("%05d", as.integer(df$porig))
      ii <- as.numeric(substr(s, 1, 2))
      ll <- as.numeric(substr(s, 3, 4))
      q  <- as.numeric(substr(s, 5, 5))
      pinch <- ii + ll/12 + q/48
      
      # default to Rijnlandse inch (2.62 cm -> f=26.2)
      p  <- round(convert_pressure(pinch, f = 26.2, lat = lat, alt = alt, atb = ta), 1)
      meta.p <- paste0(df$meta.time, " | orig.p=", ii, ".", ll, ".", q,
                       " Rijnlandse_inch.line.quarter | atb=", ta, "C")
    } else if (p_mode == "rijnlandse_3digits") {  # e.g. Breda
      s  <- sprintf("%03d", as.integer(df$porig))
      ii <- as.numeric(substr(s, 1, 2))
      ll <- as.numeric(substr(s, 3, 3))
      pinch <- ii + ll/12
      pinch[pinch == 0] <- NA # only for Breda, there's 0 here and there
      p  <- round(convert_pressure(pinch, f = 25.73, lat = lat, alt = alt, atb = ta), 1)
      meta.p <- paste0(df$meta.time, " | orig.p=", ii, ".", ll, 
                       " Rijnlandse_inch.line | atb=", ta, "C")
    } else if (p_mode == "mmHg10") {
      p <- round(convert_pressure(as.numeric(df$porig)/10, lat = lat, alt = alt, atb = ta), 1)
      meta.p <- paste0(df$meta.time, " | orig.p=", as.numeric(df$porig)/10, "mmHg | atb=", ta, "C")
    } else if (p_mode == "mmHg100") {
      p <- round(convert_pressure(as.numeric(df$porig)/100, lat = lat, alt = alt, atb = ta), 1)
      meta.p <- paste0(df$meta.time, " | orig.p=", as.numeric(df$porig)/100, "mmHg | atb=", ta, "C")
    }
  }
  
  # debug for wrong data
  if (!is.null(p)) {
    print(df[!is.na(p) & p < 900, c("Year", "Month", "Day")])
  } else {
    print("No pressure recorded in this station.")
  }
  
  rr_origunit <- ifelse(!is.na(rr_factor) & rr_factor == 0.1, "mm", "l")
  
  tibble(
    Year = df$Year, Month = df$Month, Day = df$Day, Hour = df$Hour, Minute = df$Minute,
    dd = dd, ta = ta, p = p, rr = rr, rh = rh,
    meta.dd = paste0(df$meta.time, " | orig.dd=", df$ddorig),
    meta.rh = df$meta.time,
    meta.ta = paste0(df$meta.time, " | orig.ta=",
                     dplyr::case_when(
                       ta_scale == "F10" ~ as.numeric(df$taorig)/10,
                       ta_scale == "F1"  ~ as.numeric(df$taorig),
                       ta_scale == "C10" ~ as.numeric(df$taorig)/10,
                       ta_scale == "C1"  ~ df$taorig,
                       ta_scale == "MIX_1836F10_else_C10" & df$Year == 1836 ~ as.numeric(df$taorig)/10,
                       TRUE ~ as.numeric(df$taorig)/10
                     ),
                     ifelse(ta_scale == "F10" | ta_scale =="F1" | (ta_scale == "MIX_1836F10_else_C10" & df$Year == 1836),
                            "F", "C")),
    meta.rr = if (!is.null(rr)) paste0(df$meta.time, " | orig.rr=", df$rrorig/10, rr_origunit) else NA_character_,
    meta.p  = meta.p
  )
}



# ---- 4) A single writer for multiple variables ------------------------------
write_sef_all_vars <- function(df_vars, name, code, lat, lon, alt, outdir,
                               meta_head_p = "", meta_head_all = "",
                               source = "KNMI", time_offset_bool = FALSE,
                               link = "https://www.knmi.nl/nederland-nu/klimatologie/daggegevens/antieke-waarnemingen") {
  base_cols <- df_vars %>% select(Year, Month, Day, Hour, Minute)
  
  vars_possibilities <-  c("dd","p","ta","rr", "rh")
  
  # select the actual variables found in the df
  vars <- vars_possibilities[vars_possibilities %in% names(df_vars)]
  
  if (!length(vars)) {
    cat("  (no writeable variables found)\n")
    return(invisible(NULL))
  }
  
  # make directory if it doesn't exist
  station_dir <- file.path(outdir, name)
  if (!dir.exists(station_dir)) {
    dir.create(station_dir, recursive = TRUE)
    cat("  Created directory:", station_dir, "\n")
  }
  
  # calculate time_offset if there is
  time_offset <- if (time_offset_bool) time.offset(lon) else 0L

  for (var in vars) {
    
    # join metaHeads if more than one
    metaheads <- c(
      if (!is.na(meta_head_all) && nzchar(meta_head_all)) meta_head_all else NULL,
      if (var=="p" && !is.na(meta_head_p) && nzchar(meta_head_p)) meta_head_p else NULL
    )
    metaHead <- if (length(metaheads)) paste(metaheads, collapse = " | ") else ""
    
    # select proper meta column
    metaname <- paste0("meta.", var)
    metavec  <- if(metaname %in% names(df_vars)) {
      df_vars[[metaname]]
    } else {
      rep("", nrow(df_vars))
    }
    
    dat <- bind_cols(base_cols, setNames(df_vars[var], var))
    
    cat("  → Writing", var, "rows:", sum(!is.na(dat[[var]])), "/", nrow(dat), "\n")
    write_sef_f(
      as.data.frame(dat),
      outfile  = outfile.name(name, var, dat, TRUE),
      outpath  = file.path(outdir, name),
      cod      = code,
      lat      = lat,
      lon      = lon,
      alt      = alt,
      sou      = source,
      link     = link,
      nam      = name,
      var      = var,
      stat     = "point",
      period   = 0L,
      units    = units(var),
      meta     = metavec,
      metaHead = metaHead,
      time_offset = time_offset,
      keep_na  = FALSE
    )
  }
}


# ---- 5) Station registry -------------------------------------------
stations <- read_csv("/scratch3/PALAEO-RA/daily_data/original/knmi-inventory.csv")

# convert column of file order:
stations$file_order <- lapply(stations$file_order, function(x) {
  if (is.null(x) || is.na(x)) return(NULL)
  as.numeric(unlist(strsplit(gsub(" ", "", x), ",")))
})

# filter stations to use now
if (!is.na(STATION_NAME)) {
  stations <- stations[grepl("Maastricht", stations$name), ]
}
print(stations)

# ---- 6) One driver to process a station row ---------------------------------
process_station <- function(st_row) {
  cat("\n====================\nProcessing station:", st_row$name, "\n")
  
  all_files <- list.files(st_row$dir, full.names = TRUE)

  # Filter by pattern if specified
  if (!is.na(st_row$file_pattern[[1]])) {
    all_files <- all_files[grepl(st_row$file_pattern[[1]], basename(all_files))]
  }
  
  # Reorder if specified
  if (!is.null(st_row$file_order[[1]])) {
    all_files <- all_files[st_row$file_order[[1]]]
  }
  
  if (length(all_files) == 0) {
    cat("  !! No files found for", st_row$name, "→ skipping\n")
    return(invisible(NULL))
  }
  
  raw  <- read_knmi_files(all_files, skip = st_row$skip)
  
  df   <- stamp_time(raw)
  has_hour <- (("hh"   %in% names(raw)) && any(!is.na(raw$hh))) ||
              (("hhmm" %in% names(raw)) && any(!is.na(raw$hhmm)))
  
  df_vars <- compute_vars(
    df, 
    lat = st_row$lat, 
    lon = st_row$lon, 
    alt = st_row$alt,
    ta_scale  = st_row$ta_scale %||% NA_character_,
    p_mode    = st_row$p_mode %||% NA_character_,
    rr_factor = st_row$rr_factor %||% NA_real_
  )
  
  summarise_station(st_row$name, OUTDIR, df_vars)
  
  write_sef_all_vars(
    df_vars, 
    name = st_row$name, 
    code = st_row$code,
    lat = st_row$lat, 
    lon = st_row$lon, 
    alt = st_row$alt,
    outdir = OUTDIR, 
    meta_head_p   = st_row$meta_head_p,
    meta_head_all = st_row$meta_head_all,
    time_offset_bool = has_hour
  )
  cat("  ✓ Completed", st_row$name, "\n")
}


summarise_station <- function(name, outdir, df_vars) {
  cat("\n--- Summary for", name, "---\n")
  vars <- intersect(c("dd","p","ta","rr","rh"), names(df_vars))
  for (v in vars) {
    n_non_na <- sum(!is.na(df_vars[[v]]))
    n_total  <- nrow(df_vars)
    rng <- if (is.numeric(df_vars[[v]])) {
      paste0("[", suppressWarnings(min(df_vars[[v]], na.rm=TRUE)), ", ",
             suppressWarnings(max(df_vars[[v]], na.rm=TRUE)), "]")
    } else {
      # e.g. 'dd' can be "calm" or "135"
      paste0("examples: ", paste(head(unique(df_vars[[v]]), 3), collapse=", "))
    }
    cat(sprintf("  %s: %d/%d non-NA; %s\n", v, n_non_na, n_total, rng))
  }
  dest <- file.path(outdir, name)
  if (dir.exists(dest)) {
    wrote <- list.files(dest, full.names = FALSE)
    cat("  Files in output dir:", length(wrote), "\n")
    if (length(wrote)) cat("   ", paste(head(wrote, 8), collapse="  "), if (length(wrote)>8) " ..." else "", "\n")
  } else {
    cat("  (Output directory does not exist yet)\n")
  }
}


# ---- 7) Run everything -------------------------------------------------------
for (i in seq_len(nrow(stations))) {
  st_row <- stations[i, ]
  process_station(st_row)
}


