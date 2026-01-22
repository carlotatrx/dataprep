rm(list=ls())
library(dataresqc)
library(readxl)
library(dplyr)
library(tidyr)
library(hms)
source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')

# this won't work bc I slightly modified the cell structure
# file1 <- '/scratch3/PALAEO-RA/DataRescue/BORIS/XLS/Europe/Frankfurt/Europe_T3_DE_Frankfurt_1749_subdaily.xlsx'
# file2 <- '/scratch3/PALAEO-RA/DataRescue/BORIS/XLS/Europe/Frankfurt/Europe_T3_DE_Frankfurt_1750-1752_subdaily.xlsx'
# file3 <- '/scratch3/PALAEO-RA/DataRescue/BORIS/XLS/Europe/Frankfurt/Europe_T3_DE_Frankfurt_1753-1755_subdaily.xlsx'

# do this instead:
file1 <- '/scratch3/PALAEO-RA/daily_data/original/Frankfurt/Europe_T3_DE_Frankfurt_1749_subdaily.xlsx'
file2 <- '/scratch3/PALAEO-RA/daily_data/original/Frankfurt/Europe_T3_DE_Frankfurt_1750-1752_subdaily.xlsx'
file3 <- '/scratch3/PALAEO-RA/daily_data/original/Frankfurt/Europe_T3_DE_Frankfurt_1753-1755_subdaily.xlsx'

raw1 <- read_excel(file1, sheet = 1, skip  = 6)
raw2 <- read_excel(file2, sheet = 1, skip  = 6)
raw3 <- read_excel(file3, sheet = 1, skip = 6)

lat <- 50.11180
lon <- 8.68131
ele <- 107

# Keep only the first 10 columns and rename
df1 <- raw1 %>%
  select(1:9) %>%
  rename(
    Year = 1, Month = 2, Day = 3, Hour = 4, p = 5,
    ta = 6, ta.corrected = 7, direction = 8, force=9
  ) %>%
  # Fill down Year/Month/Day
  mutate(across(where(is.character), ~ na_if(trimws(.x), ""))) %>%
  fill(Year, Month, Day, .direction = "down") %>%
  mutate(
    Year = as.integer(Year),
    Month = as.integer(Month),
    Day = as.integer(Day),
    Hour = round(suppressWarnings(as.numeric(Hour))* 86400/ 3600, 0), # convert stupidity
    p = suppressWarnings(as.numeric(p)),
    p.zoll = as.integer(p%/%10), # first 2 digits
    p.in = as.integer(round((p%%10) * 10)), # everything after decimal * 10
    ta = suppressWarnings(as.numeric(ta)),
    ta.corrected = suppressWarnings(as.numeric(ta.corrected)),
    force = suppressWarnings(as.numeric(force))
  ) %>% 
  select(-any_of("p"))
head(df1)

# Keep only the first 10 columns and rename
df2 <- raw2 %>%
  select(1:10) %>%
  rename(
    Year = 1, Month = 2, Day = 3, Hour = 4, direction = 5, force=6,
    ta = 7, p.zoll = 8, p.in = 9, hygrom = 10 
  ) %>%
  # Fill down Year/Month/Day
  mutate(across(where(is.character), ~ na_if(trimws(.x), ""))) %>%
  fill(Year, Month, Day, p.zoll, .direction = "down") %>%
  mutate(
    Year = as.integer(Year),
    Month = as.integer(Month),
    Day = as.integer(Day),
    Hour = suppressWarnings(as.integer(Hour)),
    p.zoll = suppressWarnings(as.numeric(p.zoll)),
    p.in = suppressWarnings(as.numeric(p.in)),
    ta = suppressWarnings(as.numeric(ta)),
    force = suppressWarnings(as.numeric(force)),
    direction = ifelse(direction=="-", NA, direction)
  )
head(df2)

df3 <- raw3 %>%
  select(1:11) %>%
  rename(
    Year = 1, Month = 2, Day = 3, Hour = 4, direction = 5, force=6,
    ta = 7, p.zoll = 8, p.in = 9, p.zoll.min = 10, p.in.min = 11
  ) %>%
  # Fill down Year/Month/Day
  mutate(across(where(is.character), ~ na_if(trimws(.x), ""))) %>%
  fill(Year, Month, Day, p.zoll, p.zoll.min, .direction = "down") %>%
  mutate(
    Year = as.integer(Year),
    Month = as.integer(Month),
    Day = as.integer(Day),
    Hour = suppressWarnings(as.integer(Hour)),
    p.zoll = suppressWarnings(as.numeric(p.zoll)),
    p.in = suppressWarnings(as.numeric(p.in)),
    p.zoll.min = suppressWarnings(as.numeric(p.zoll.min)),
    p.in.min = suppressWarnings(as.numeric(p.in.min)),
    ta = suppressWarnings(as.numeric(ta)),
    force = suppressWarnings(as.numeric(force)),
    direction = ifelse(direction=="-", NA, direction)
  )
head(df3)

# concatenate
all.cols <- Reduce(union, list(names(df1),names(df2),names(df3)))

# add missing cols as NA
add.cols.NA <- function(x, cols) {
  missing.cols <- setdiff(cols, names(x))
  if (length(missing.cols)) {
    for (m in missing.cols) x[[m]] <- NA
  }
  x < x[, cols]   # put cols in correct order (in case there was a missing one, it's already in cols)
  x
}

df1.all <- add.cols.NA(df1, all.cols)
df2.all <- add.cols.NA(df2, all.cols)
df3.all <- add.cols.NA(df3, all.cols)

df.all <- bind_rows(df1.all, df2.all, df3.all) %>%
  arrange(Year, Month, Day, Hour)

## Define wind directions
dirs <- c("N", "NNO", "NO", "ONO", "O", "OSO", "SO", "SSO", "S", 
          "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW")

library(stringr)
# Build regex: longest-first, match as a standalone token (no letters around)
dir.ordered <- dirs[order(nchar(dirs), decreasing = TRUE)]
dir.pattern <- paste0("(?<![A-Z])(", paste(dir.ordered, collapse="|"), ")(?![A-Z])")

df.all <- df.all %>%
  mutate(
    dir.clean = str_extract(toupper(trimws(direction)),  # column direction
                            regex(dir.pattern, ignore_case=F)),
    dd = ifelse(
      is.na(dir.clean),
      NA_real_,
      round(22.5 * (match(dir.clean, dirs) - 1), 0)
    ),
    
    ta.C = round((ta-32) * 5/9, 1),
    
    p.together = if_else(
      is.na(p.in),
      NA_real_,
      round(convert_pressure(p.zoll + p.in/100, f=27.07, lat=lat, atb=ta.C, alt=ele),1)  # consider the offset for realistic values
    ),
    
    meta.dd = if_else(
      !is.na(force),
      paste0("orig.dd=", direction, " | strength=",force),
      paste0("orig.dd=", direction)
    ),
    meta.p = if_else(
      !is.na(p.zoll.min) & !is.na(p.in.min),
      paste0("orig.p=", p.zoll, "|", p.in, " | alt.p.min=", p.zoll.min, "|", p.in.min),
      paste0("orig.p=", p.zoll, "|", p.in)
    ),
    meta.ta = if_else(
      !is.na(ta.corrected),
      paste0("orig.ta=", ta, "F | orig.ta.corrected=", ta.corrected, "F"),
      paste0("orig.ta=", ta)
    ),
    Minute = NA
  )

head(df.all)

# save

# ta
write_sef_f(
  Data=as.data.frame(df.all[, c("Year", "Month", "Day", "Hour", "Minute", "ta.C")]),
  outfile="Frankfurt_ta_subdaily.tsv",
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod="Frankfurt",
  lat=lat, lon=lon, alt=ele,
  variable="ta",
  nam="Frankfurt",
  meta=df.all$meta.ta,
  metaHead = "assumed Fahrenheit as original units | approx.coords",
  units="C", sou="PALAEO-RA", stat="point",keep_na = F
)

# p
write_sef_f(
  Data=as.data.frame(df.all[, c("Year", "Month", "Day", "Hour", "Minute", "p.together")]),
  outfile="Frankfurt_p_subdaily.tsv",
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod="Frankfurt",
  variable="p",
  nam="Frankfurt",
  lat=lat, lon=lon, alt=ele,
  meta=df.all$meta.p,
  metaHead="assumed Paris inch as original units with a +4 inch offset | approx. coords | PTC=Y | PGC=Y",
  units="hPa", sou="PALAEO-RA", stat="point", keep_na = F
)

# dd
write_sef_f(
  Data=as.data.frame(df.all[, c("Year", "Month", "Day", "Hour", "Minute", "dd")]),
  outfile="Frankfurt_dd_subdaily.tsv",
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod="Frankfurt",
  variable="dd",
  lat=lat, lon=lon, alt=ele,
  nam="Frankfurt",
  meta=df.all$meta.dd,
  metaHead = "approx. coords",
  units="deg", sou="PALAEO-RA", stat="point", keep_na = F
)
