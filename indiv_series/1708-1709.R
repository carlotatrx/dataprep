# Parse Lenke (1709) pasted table into per-station SEF-like tables:
# Year, Month, Day, Value
#
# Usage:
#   raw_txt <- readLines("lenke_table.txt", encoding = "UTF-8")
#   out <- lenke_to_sef(raw_txt)
#   names(out)            # station names
#   out[["Berlin"]]       # one station table
#   # optional: write each station to a TSV
#   purrr::iwalk(out, ~ readr::write_tsv(.x, paste0("SEF_", .y, "_ta_1709.tsv")))

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(readxl)
})

x <- read_xlsx("/scratch3/PALAEO-RA/daily_data/original/Berlin/Lenke_Tab3_1708_1709.xlsx",
               skip=3, col_names=FALSE)

# Month mapping (German)
month_map <- c(
  januar = 1, februar = 2, maerz = 3, maerz_uml = 3, april = 4,
  mai = 5, juni = 6, juli = 7, august = 8, september = 9,
  oktober = 10, november = 11, dezember = 12
)

norm_month <- function(s) {
  s %>%
    str_to_lower() %>%
    str_replace_all("ä", "ae") %>%
    str_replace_all("ö", "oe") %>%
    str_replace_all("ü", "ue") %>%
    str_replace_all("ß", "ss")
}

is_month_row <- function(v) {
  if (length(v) == 0 || is.na(v)) return(FALSE)
  v <- as.character(v)
  str_detect(str_trim(v), "^\\p{L}+\\s+\\d{4}$")
}

parse_month_year <- function(v) {
  m <- str_match(str_trim(v), "^(\\p{L}+)\\s+(\\d{4})$")
  mon <- norm_month(m[,2])
  yr  <- as.integer(m[,3])
  mo  <- unname(month_map[mon])
  if (is.na(mo)) stop("Unknown month name in sheet: ", m[,2])
  list(year = yr, month = mo)
}

# Helper to coerce values, stripping "*" (open-window mark)
to_num <- function(v) {
  v <- as.character(v)
  v <- str_replace_all(v, ",", ".")
  v <- str_replace(v, "\\*$", "")
  suppressWarnings(as.numeric(v))
}

# Iterate blocks
blocks <- list()
i <- 1
n <- nrow(x)

while (i <= n) {
  month_cell <- x[[2]][[i]]  # column B
  if (is_month_row(month_cell)) {
    ym <- parse_month_year(month_cell)
    
    # day header is next row, columns B..end
    if (i + 1 > n) break
    days <- suppressWarnings(as.integer(unlist(x[i + 1, 2:ncol(x)])))
    
    # station rows start at i+2 until next month row or empty row
    j <- i + 2
    while (j <= n) {
      # stop if next month starts
      if (is_month_row(x[[2]][[j]])) break
      
      station <- x[[1]][j]  # column A
      # stop if fully empty line
      if (all(is.na(unlist(x[j, ])))) { j <- j + 1; next }
      
      if (!is.na(station) && str_trim(as.character(station)) != "") {
        vals_raw <- unlist(x[j, 2:ncol(x)])
        vals <- to_num(vals_raw)
        
        df <- tibble(
          station = as.character(station),
          Year = ym$year,
          Month = ym$month,
          Day = days,
          Value = vals
        )
        
        blocks[[length(blocks) + 1]] <- df
      }
      j <- j + 1
    }
    i <- j
  } else {
    i <- i + 1
  }
}

all_df <- bind_rows(blocks) %>%
  arrange(station, Year, Month, Day)

# One table per station (SEF-like)
out <- all_df %>%
  group_by(station) %>%
  group_split(.keep = TRUE) %>%
  set_names(map_chr(., ~ .x$station[1])) %>%
  map(~ select(.x, Year, Month, Day, Value))

# Example:
out[["Berlin"]]
station_names <- names(out)  # stations found

station_meta <- tibble::tribble(
  ~station,       ~lat,     ~lon,      ~alt,
  "Berlin",        52.5656, 13.3106,  36,
  "Danzig",        54.35, 18.53,  18,
  "Delft",         52.1014, 5.1867,  0,
  "Halle",         51.483056, 11.966667,  102,
  "Hamburg",       53.550, 10.000,  116,
  "Jena",          50.927054, 11.5892372,  150,
  "Kiel",          54.31667, 10.15, 5,
  "Kopenhagen",    55.6760968, 12.5683372,  40,
  "Königsberg",    54.7167, 20.5,  23,
  "Montpellier",   43.61077, 3.876716,  57,
  "Paris",         48.85, 2.34, 57,
  "Upminster",     51.55591, 0.248894,  19,
  "Zeitz",         51.04778, 12.13833, 210
)

missing <- station_meta %>%
  dplyr::filter(station %in% station_names) %>%
  dplyr::filter(is.na(lat) | is.na(lon) | is.na(alt))

if (nrow(missing) > 0) {
  stop("Missing lat/lon/alt for: ", paste(missing$station, collapse = ", "))
}

var<-"ta"
for (station in station_names) {
  m <- dplyr::filter(station_meta, station == !!station)
  if (nrow(m) != 1) stop("Metadata not unique for station: ", station)
  dirname <- case_when(
    station == "Danzig" ~ "Gdansk",
    station == "Kopenhagen" ~ "Copenhagen",
    station == "Königsberg" ~ "Kaliningrad",
    TRUE ~ station
    )
  
  df <- out[[station]] %>%
    dplyr::mutate(
      Hour   = NA_integer_,
      Minute = NA_integer_,
      .before = Value
    )
  
  
  write_sef_f(
    as.data.frame(df),
    outfile = outfile.name(dirname, var, df, subdaily=FALSE),
    outpath = file.path("/scratch3/PALAEO-RA/daily_data/final", dirname),
    cod     = station,
    lat     = lat,
    lon     = lon,
    alt     = alt,
    sou     = "Untersuchung der ältesten Temperaturmessungen mit Hilfe des strengen Winters 1708-1709, Walter Lenke",
    link    = "",
    nam     = paste0("Lenke_", dirname),
    var     = var,
    stat    = "point",
    units   = units(var),
    metaHead = "coords=approx",
    keep_na = TRUE
  )
}
