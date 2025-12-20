rm(list=ls())
library(dataresqc)
library(readxl)
library(tidyr)
library(dplyr)
library(lubridate)
library(tibble)
library(stringr)
source('/scratch2/ccorbella/code/dataprep/helpfun.R')

outdir <- "/scratch3/PALAEO-RA/daily_data/final/Loebau"
name <- "Lausitz"
code <- "Loebau"
lat	<- 51.08333333
lon	<- 14.66666667
alt	<- 260
source <- "PALAEO-RA"
link   <- ""

# account for winds -------------------------------------------------------

dd_normalize_Loebau <- function(x) {
  
  y <- toupper(trimws(x))
  
  # 1. hard NA / calm
  y[y %in% c("NA", "M", "0", "-", "--", "- - -")] <- NA_character_
  y[grepl("CIRCUM|CALM", y)] <- "CALM"
  
  # 2. remove junk (quotes, words, line breaks)
  y <- gsub("[^A-Z&\\.]", "", y)
  
  # 3. normalize separators
  y <- gsub("\\.|UND|U", "&", y)
  y <- gsub("&+", "&", y)
  
  # 4. historical German: O = East
  y <- gsub("O", "E", y)
  
  # 5. split composites and keep first valid direction
  pick_dir <- function(s) {
    
    parts <- unique(unlist(strsplit(s, "&")))
    parts <- parts[parts %in% c("N","E","S","W",
                                "NE","SE","SW","NW")]
    
    if (length(parts) == 0) return(NA_character_)
    
    # exact cardinal pairs → intercardinal
    if (length(parts) == 2 && all(parts %in% c("N","E","S","W"))) {
      pair <- paste(sort(parts), collapse = "")
      return(switch(pair,
                    "ES" = "SE",
                    "EN" = "NE",
                    "NW" = "NW",
                    "SW" = "SW",
                    "NS" = NA_character_,  # opposite directions → undefined
                    "EW" = NA_character_,
                    NA_character_))
    }
    
    # diagonal + cardinal → secondary intercardinal
    if (length(parts) == 2) {
      pair <- paste(sort(parts), collapse = "")
      return(switch(pair,
                    "SW" = "WSW",
                    "NW" = "WNW",
                    "EN" = "ENE",
                    "ES" = "ESE",
                    NA_character_))
    }
    
    # everything else → first valid
    parts[1]
  }
  
  out <- vapply(y, pick_dir, character(1))
  
  # 6. allow calm explicitly
  out[out == "CALM"] <- "calm"
  
  out
}


# 1720-1726  --------------------------------------------------------------
raw <- read_excel("/scratch3/PALAEO-RA/daily_data/original/Loebau/Kanold_T3_DE_Loebau_1720-1726_subdaily.xls",
                  na="-", skip=8,
                  col_names=c("Year", "Month", "Day", "pgr", "plin", "ta", "ta_string", "dd.orig", "notes", "notes2", "notes3")
                 )

raw2 <- read_excel("/scratch3/PALAEO-RA/daily_data/original/Loebau/Kanold_T3_DE_Loebau_1727-1730_subdaily.xls",
                   na="-", skip=8,
                   col_names=c("Year", "Month", "Day", "pgr", "plin", "ta", "ta_string", "dd.orig", "notes")
                  )
head(raw)
head(raw2)

df <- raw %>%
  mutate(
    Hour =NA_integer_,
    Minute=NA_integer_,
    dd = dd2deg(dd_normalize_Loebau(dd.orig)),
    ta = ifelse(!is.na(ta_string), paste0(ta, ta_string), ta),
    meta.dd = paste0("orig.dd=",dd.orig),
  ) %>% select(Year, Month, Day, Hour, Minute,ta,dd, meta.dd)

df2 <- raw2 %>%
  mutate(
    Hour =NA_integer_,
    Minute=NA_integer_,
    dd = dd2deg(dd_normalize_Loebau(dd.orig)),
    ta = ifelse(!is.na(ta_string), paste0(ta, ta_string), ta),
    meta.dd = paste0("orig.dd=",dd.orig),
  ) %>% select(Year, Month, Day, Hour, Minute,ta,dd, meta.dd)
head(df)
head(df2)

df_all <- bind_rows(df, df2)

base_cols <- df_all[c("Year","Month","Day","Hour","Minute")]

# dd
var <-"dd"
dat <- cbind(base_cols, setNames(df_all[var], var))
write_sef_f(
  dat,
  outfile = outfile.name(name, var, dat, subdaily=TRUE),
  outpath = outdir,
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "point",
  units   = units(var),
  meta    = df_all$meta.dd,
  metaHead = "Observer=C.Trautmann | obs.time=unknown",
  keep_na = TRUE
)

# ta
var <-"ta"
dat <- cbind(base_cols, setNames(df_all[var], var))
write_sef_f(
  dat,
  outfile = outfile.name(name, var, dat, subdaily=TRUE),
  outpath = outdir,
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "point",
  units   = "unknown",
  metaHead = "Observer=C.Trautmann | obs.time=unknown",
  keep_na = TRUE
)
