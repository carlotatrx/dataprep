rm(list=ls())
library(dataresqc)
library(readxl)
library(dplyr)
library(stringr)
source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')

file <- '/scratch3/PALAEO-RA/daily_data/original/Europe_T3_PL_Warsaw_1725-1728_subdaily.xlsx'
raw <- read_excel(file, sheet = 1, skip  = 6)

lat	<- 52.2296756
lon	<- 21.0122287
alt	<- NA
name <- "Warsaw"
metaHead <- "Observer=Erndtel"


time.offset <- as.numeric(lon)*12/180

###### Wind direction ########
# Convert Warsaw-style wind codes (O=East) to daily degrees
# Strategy: split on '.'/' ' separators, map each token to degrees, circular mean per row.

# Base 16-point (with O=East) + 32-point extensions + a few forgiving shortcuts
dir_deg <- c(
  "N"= 0, "NNO"= 22.5, "NO"= 45, "ONO"= 67.5, "O"= 90, "OSO"=112.5,
  "SO"=135, "SSO"=157.5, "S"=180, "SSW"=202.5, "SW"=225, "WSW"=247.5,
  "W"=270, "WNW"=292.5, "NW"=315, "NNW"=337.5,
  # Forgiving shorthands occasionally seen in manuscripts:
  # 'WN' ~ WNW (292.5), 'WS' ~ WSW (247.5), 'NO' and 'SO' already in table.
  "WN"=292.5, "WS"=247.5
)

# helper: map a single token to degrees
map_token <- function(tok){
  if (is.na(tok)) return(NA_real_)
  t0 <- toupper(tok)
  t0 <- gsub("[^A-Z]", "", t0, perl=T)  # keep letters only
  if (!nzchar(t0)) return(NA_real_)     # if t0 is empty
  
  # Direct hit
  hit <- unname(dir_deg[t0])
  if (!is.na(hit)) return(hit)
  
  # Try collapsing repeated letters like 'NNOO' -> 'NNO' (rare OCR oddity)
  t2 <- gsub("(.)\\1+", "\\1", t0, perl=T)
  hit <- unname(dir_deg[t2])
  if (!is.na(hit)) return(hit)
  
  # else return NA  
  NA_real_
}

# circular mean in degrees (0-360)
circ_mean_deg <- function(deg){
  deg <- deg[!is.na(deg)]
  if (!length(deg)) return(NA_real_)
  rad <- deg * pi/180
  mX <- mean(cos(rad)); mY <- mean(sin(rad))
  out <- (atan2(mY, mX) * 180/pi) %% 360
  out
}

# Main function: vectorized over a character vector of codes
wind_to_dd <- function(x){
  unknown <- character(0)
  out <- vapply(x, function(s){
    if (is.na(s)) return(NA_real_)
    s2 <- toupper(trimws(as.character(s)))
    if (!nzchar(s2)) return(NA_real_)
    # standardize separators '.' and spaces
    parts <- unlist(strsplit(s2, "[\\.[:space:]]+", perl=T))
    parts <- parts[nzchar(parts)]
    if (!length(parts)) return(NA_real_)
    vals <- vapply(parts, function(p){
      d <- map_token(p)
      if (is.na(d) && nzchar(p)) unknown <<- unique(c(unknown, p))
      d
    }, numeric(1))
    circ_mean_deg(vals)
  }, numeric(1))
  attr(out, "unknown_tokens") <- sort(unique(unknown))
  out
}


# Keep only the first 10 columns and rename
df <- raw %>%
  select(1:8) %>%
  rename(
    Year = 1, Month = 2, Day = 3, p.morning.zoll = 4, p.morning.in = 5, p.afternoon.zoll = 6, p.afternoon.in = 7, dd.orig = 8
  ) %>%
  # Fill down Year/Month/Day
  mutate(across(where(is.character), ~ na_if(trimws(.x), ""))) %>%
  mutate(
    Year = as.integer(Year),
    Month = as.integer(Month),
    Day = as.integer(Day),
    Hour = NA,
    Minute = NA,
    
    dd = round(wind_to_dd(dd.orig),2),
    meta.dd = paste0("orig.dd=",dd.orig)
  )
head(df)

# save
# dd
write_sef_f(
  Data=as.data.frame(df[, c("Year", "Month", "Day", "Hour", "Minute", "dd")]),
  outfile=paste0(name,"_dd_subdaily.tsv"),
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod=name,
  lat=lat, lon=lon, alt=alt,
  variable="dd",
  nam=name,
  meta=df$meta.dd,
  units="deg", sou="PALAEO-RA", stat="point",keep_na = F
)
