rm(list=ls())
library(dataresqc)
library(tidyr)
library(dplyr)
library(lubridate)
source('/scratch2/ccorbella/code/dataprep/helpfun.R')

outdir <- "/scratch3/PALAEO-RA/daily_data/final/Zurich"



# read SEF/TSV, keep strings as-is
f <- "/scratch3/PALAEO-RA/daily_data/original/Zurich/ZH01_Zurich_Scheuchzer.tsv"
raw <- read.table(f, sep = "\t", header = TRUE,stringsAsFactors = FALSE,
                 na.strings = c("NA", "", "NaN"))

head(raw)
# columns with at least one non-NA value
non_empty_cols <- names(raw)[colSums(!is.na(raw)) > 0]

name <- unique(na.omit(raw$station))
code <- unique(na.omit(raw$place))
lat	<- unique(na.omit(raw$latitude))
lon	<- unique(na.omit(raw$longitude))
alt	<- unique(na.omit(raw$altitude))
metaHead <- paste0("Observer=", unique(na.omit(raw$observer)))
source <- "CHIMES"
link   <- "https://doi.org/10.5194/cp-15-1345-2019"


# print them
print(non_empty_cols)
colnames(raw)

keep_cols <- c("year","month","day","hour","minutes","time", "time_orig","air_pressure_orig","air_pressure_unit",
               "air_temperature_orig", "wind_from_direction",
               "wind_from_direction_V2",           
               "air_pressure_V2",                "air_pressure_orig_V2"  ,           
               "air_pressure_unit_V2",              "air_pressure2"       ,             
               "air_pressure2_orig"   ,             "air_pressure2_unit"   ,            
               "accumulated_precipitation_V2" ,    "accumulated_precipitation_orig_V2",
               "air_pressure_hPa"              ,    "air_pressure2_hPa"   )
raw2 <- raw %>%
  select(all_of(keep_cols)) 

head(raw2)

# wind direction function for Zurich Schuechzer
dd_normalize_zh01 <- function(x) {
  y <- toupper(trimws(x))
  y[y %in% c("", "-", "NA")] <- NA_character_
  
  # clean: keep only letters and spaces (turn everything else into space)
  y <- gsub("[^A-Z]+", " ", y)
  y <- gsub("\\b(INT|INTER|FORT|FORTIS|NUB|NUBIB|NUBIBUS|IMO|INFRA|AERE|A\\b|AD\\b|VAR|VARIANT|VARIANTES|PLAGA|OCCIDENTALI|VENTI|IN|I)\\b",
            " ", y)
  y <- gsub("\\s+", " ", trimws(y))
  
  # split into tokens
  toks <- strsplit(y, " ", fixed = TRUE)
  
  norm_one <- function(tt) {
    if (length(tt) == 0) return(NA_character_)
    tt <- tt[tt != ""]
    if (length(tt) == 0) return(NA_character_)
    
    # take the first token that contains any direction letters
    t1 <- tt[1]
    if (is.na(t1) || t1 == "") return(NA_character_)
    
    # some tokens are glued like "NNOO" "NOO" "SOO" "SSSO" "SSWW" "NWW" "NWN" "NON"
    t1 <- dplyr::recode(t1,
                        # German notation: O=E system
                        "O"    = "E",
                        "NO"   = "NE",
                        "SO"   = "SE",
                        "ONO"  = "ENE",
                        "OSO"  = "ESE",
                        "NOO"  = "ENE",
                        "SOO"  = "ESE",
                        "NNO"  = "NNE",
                        "SSO"  = "SSE",
                        
                        # weird combined forms seen in your list
                        "NNOO" = "NNE",  # treat as NNO
                        "SSSO" = "SSE",
                        "SSWW" = "WSW",
                        "NWW"  = "WNW",
                        "NWN"  = "NW",
                        "NON"  = "NE",
                        "WO"   = "WSW",
                        "WS"   = "SW",
                        "SWS"  = "SSW",
                        "SWW"  = "WSW",
                        "WWS"  = "WSW",
                        "OS"   = "ESE",
                        "ON"   = "ENE",
                        "SOS"  = "SE",   # appears in your list
                        "OSO"  = "ESE",
                        .default = t1
    )
    
    # if still not a standard 16-wind, try to salvage by picking first recognizable piece
    # e.g. "NNO" would already be mapped; but "NN O" etc might appear after cleaning
    if (t1 %in% directions) return(t1)
    
    # try to map single-letter O/W/N/S if present
    if (t1 %in% c("O","W","N","S")) return(dplyr::recode(t1, "O"="E", .default=t1))
    NA_character_
  }
  
  out <- vapply(toks, norm_one, FUN.VALUE = NA_character_)
  ifelse(out %in% directions, out, NA_character_)
}

# check for unique values
unique(na.omit(raw2$wind_from_direction_V2))

# treu els factors de després, e.g. 26.1.33333333333 -> 26.1.33
trim_frac2 <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    grepl("^\\d+\\.\\d+\\.\\d+\\.\\d+", x) ~ sub("(^\\d+\\.\\d+\\.\\d+\\.)(\\d{2}).*", "\\1\\2", x),
    TRUE ~ x
  )
}


df <- raw2 %>%
  mutate(
    year    = as.integer(year),
    month   = as.integer(month),
    day     = as.integer(day),
    hour = case_when(
      !is.na(hour) ~ as.integer(hour),
      time_orig == "vormittag" ~ 11L,
      time_orig == "nachmittag" ~ 15L,
      TRUE ~ NA_integer_
    ),
    minute = case_when(
      !is.na(minutes) ~ as.integer(minutes),
      time_orig %in% c("vormittag", "nachmittag") ~ 0L,
      TRUE ~ NA_integer_
    ),
    
    p1 = round(air_pressure_hPa,2),
    p2 = round(air_pressure2_hPa,2),
    
    dd.norm1 = dd_normalize_zh01(wind_from_direction),
    dd1 = as.numeric(dd2deg(dd.norm1)),
    dd.norm2 = dd_normalize_zh01(wind_from_direction_V2),
    dd2 = as.numeric(dd2deg(dd.norm2)),
    
    rr = round(accumulated_precipitation_V2,2),
    air_pressure_orig  = trim_frac2(air_pressure_orig),
    air_pressure2_orig = trim_frac2(air_pressure2_orig),
    
    ta = trim_frac2(air_temperature_orig),

    meta.time = case_when(
      is.na(time_orig) ~ NA_character_,
      grepl("^[0-9]", time_orig) ~ paste0(
        "orig.time=",
        sub(":[0-9]{2}$", "", time_orig)
      ),
      TRUE ~ paste0("orig.time=", time_orig)
    ),
    meta.ta  = paste0(meta.time, " | orig=", ta),
    
    meta.p1 = paste0(meta.time, " | orig=", air_pressure_orig, air_pressure_unit),
    meta.p2 = paste0(meta.time, " | orig=", air_pressure2_orig, air_pressure2_unit),
    
    meta.dd1 = paste0(meta.time, " | orig=", wind_from_direction),
    meta.dd2 = paste0(meta.time, " | orig=", wind_from_direction_V2),
    meta.rr  = paste0(meta.time, " | orig=", accumulated_precipitation_orig_V2),
  ) %>% select(year, month, day, hour, minute, p1, p2, dd1, dd2, rr, ta,
               meta.ta, meta.p1, meta.p2, meta.dd1, meta.dd2, meta.rr)
head(df)


write_Zurich_Scheuchzer_sef <- function(df, var_name, value_col, meta_col, units, stat = "point") {
  keep <- !grepl("orig\\s*=\\s*(<NA>|NA)\\w*", df[[meta_col]])
  d <- df[keep, c("year","month","day","hour","minute", value_col, meta_col)]
  
  dat  <- d[, c("year","month","day","hour","minute", value_col)]
  names(dat)[names(dat) == value_col] <- var_name
  
  write_sef_f(
    dat,
    outfile  = outfile.name(name, var_name, dat, subdaily = TRUE),
    outpath  = outdir,
    cod      = code,
    lat      = lat,
    lon      = lon,
    alt      = alt,
    sou      = source,
    link     = link,
    nam      = name,
    var      = var_name,
    stat     = stat,
    units    = units,
    metaHead = metaHead,
    meta     = d[[meta_col]],
    keep_na  = TRUE
  )
}

write_Zurich_Scheuchzer_sef(df, "ta", "ta",  "meta.ta",  units = "unknown")
write_Zurich_Scheuchzer_sef(df, "p",  "p1",  "meta.p1",  units = "hPa")
write_Zurich_Scheuchzer_sef(df, "p_bis",  "p2",  "meta.p2",  units = "hPa")
write_Zurich_Scheuchzer_sef(df, "rr", "rr",  "meta.rr",  units = "unknown")
write_Zurich_Scheuchzer_sef(df, "dd", "dd1", "meta.dd1", units = "degree")
write_Zurich_Scheuchzer_sef(df, "dd_bis", "dd2", "meta.dd2", units = "degree")

