## NOTES:
## Temperature (unknown Fahrenheit's thermometer) was not formatted

library(dataresqc)
library(XLConnect)
library(readxl)
library(tidyr)
library(dplyr)
source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')

name <- "Nürnberg"
lat <-	49.45
lon	<- 11.08333333
alt	<- 316
metahead <- "Observer=Rost"

time.offset <- as.numeric(lon)*12/180

## Define function to fill missing values with the value on the previous row
fill_variable <- function(timeseries, strings) {
  for (i in 2:length(timeseries)) {
    if (timeseries[i] %in% strings) {
      if(!timeseries[i-1] %in% strings) {
        timeseries[i] <- timeseries[i-1]
      }
    }
  }
  return(timeseries)
}

get_date_range <- function(df) {
  start.date <- paste0(df$Year[1],
                       sprintf("%02d", df$Month[1]),
                       sprintf("%02d", df$Day[1]))
  n <- nrow(df)
  end.date   <- paste0(df$Year[n],
                       sprintf("%02d", df$Month[n]),
                       sprintf("%02d", df$Day[n]))
  paste0(start.date, "-", end.date)
}

clean_num  <- function(x) suppressWarnings(as.numeric(trimws(x)))
clean_hour <- function(x) suppressWarnings(as.integer(trimws(x)))


meta_frac <- function(g, b1, b2) {
  if (is.na(b1) | is.na(b2)) {
    as.character(g)  # just return the integer part
  } else {
    paste0(g, ".", b1, "/", b2)
  }
}

ta_frac <- function(g, b1, b2) {
  if (is.na(b1) | is.na(b2)) {
    as.integer(g)  # just return the integer part
  } else {
    round(as.numeric(g) + as.numeric(b1)/as.numeric(b2)/10 ,2)
  }
}


## SUBDAILY

# 17181221-17241231 -----------------------------------------------------
files <- list.files("/scratch3/PALAEO-RA/daily_data/original/Nuernberg", pattern="subdaily", full.names=TRUE)

# select only files below 1729
yrs <- regmatches(basename(files), gregexpr("\\d{4}", basename(files)))
files <- files[sapply(yrs, function(x) max(as.integer(x), na.rm = TRUE)) < 1729]

out <- list()
for (i in 1:3) {
  print(i)
  ## Read data
  raw <- readWorksheetFromFile(files[i], sheet=1, startRow=8, header=T, colTypes="character")
  raw[,1] <- fill_variable(raw[,1], NA)
  raw[,2] <- fill_variable(raw[,2], NA)
  raw[,3] <- fill_variable(raw[,3], NA)
  ## Create output data frame
  out[[i]] <- data.frame(Year = as.integer(raw[,1]),
                         Month = as.integer(raw[,2]),
                         Day = as.integer(raw[,3]),
                         Hour = as.integer(raw[,4]),
                         Minute = NA_integer_,
                         ta = mapply(ta_frac, raw$Grad, raw$bruch1, raw$bruch2),
                         ta.meta = paste0("orig.ta.A=", mapply(meta_frac, raw$Grad, raw$bruch1, raw$bruch2),
                                          " | orig.ta.B=", mapply(meta_frac, raw$Grad, raw$bruch1, raw$bruch2)),
                         stringsAsFactors = FALSE)
}

## Merge data frames
out <- do.call(rbind, out)

start.date <- paste0(out$Year[1], out$Month[1], out$Day[1])
n <- nrow(out)
end.date <- paste0(out$Year[n], out$Month[n], out$Day[n])
  
## Write SEF files
write_sef_f(out[, 1:6],
            outfile=paste0("Nuernberg_",start.date, "-", end.date, "_ta_subdaily.tsv"),
            outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
            variable = "ta",
            cod = "Kanold_Nuernberg",
            nam = name,
            lat = lat,
            lon = lon,
            alt = alt,
            sou = "PALAEO-RA",
            units = "unknown",
            stat = "point",
            meta = out$ta.meta,
            metaHead = metahead,
            period = 0,
            time_offset = time.offset,
            keep_na = F)


# 1725-1726 ---------------------------------------------------------------

file <- "/scratch3/PALAEO-RA/daily_data/original/Nuernberg/Kanold_T3_DE_Nuernberg_1725-1726_subdaily.xls"
raw <- read_excel(file, sheet=1, skip=7)

# need to redefine for this data
make_frac <- function(b1, b2) {
  ifelse(is.na(b1) | is.na(b2), "", paste0(b1, "/", b2))
}

ta_frac <- function(z, l, b1, b2) {
  frac  <- make_frac(b1, b2)
  parts <- c(
    ifelse(is.na(z), NA, as.integer(z)),
    ifelse(is.na(l), NA, as.integer(l)),
    ifelse(frac == "", NA, frac)
  )
  paste(na.omit(parts), collapse = ".")
}

make_ta <- Vectorize(ta_frac)

# Keep only the first 10 columns and rename
df <- raw %>%
  select(1:4,6:10) %>%
  rename(
    Year = 1, Month = 2, Day = 3, Hour = 4, ta.zoll = 5,
    ta.l = 6, ta.bruch1 = 7, ta.bruch2 = 8, sign = 9
  ) %>%
  # Fill down Year/Month/Day
  mutate(across(where(is.character), ~ na_if(trimws(.x), ""))) %>%
  fill(Day, .direction = "down") %>%
  mutate(
    Year = as.integer(Year),
    Month = as.integer(Month),
    Day = as.integer(Day),
    Hour = as.integer(Hour),
    Minute = NA_integer_,
    ta.zoll = suppressWarnings(as.numeric(ta.zoll)),
    ta.l = suppressWarnings(as.numeric(ta.l)),
    ta.bruch1 = suppressWarnings(as.numeric(ta.bruch1)),
    ta.bruch2 = suppressWarnings(as.numeric(ta.bruch2)),
    ta = make_ta(ta.zoll, ta.l, ta.bruch1, ta.bruch2),
    meta.ta = paste0("orig.time=", Hour,
      " | orig.ta=", ta,
      ifelse(!is.na(sign), paste0(" | sign=", sign), "")
    )
  )
head(df)

# save
write_sef_f(
  Data=as.data.frame(df[, c("Year", "Month", "Day", "Hour", "Minute", "ta")]),
  outfile=paste0("Nuernberg_",get_date_range(df), "_ta_subdaily.tsv"),
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod="Kanold_Nuernberg",
  lat=lat, lon=lon, alt=alt,
  variable="ta",
  nam=name,
  meta=df$meta.ta,
  metaHead=metahead,
  units="unknown", sou="PALAEO-RA", stat="point",keep_na = F,
  time_offset = time.offset
)



# 1727_part1 ---------------------------------------------------------------

file <- "/scratch3/PALAEO-RA/daily_data/original/Nuernberg/Kanold_T3_DE_Nuernberg_1727_subdaily_part1.xls"
raw <- read_excel(file, sheet=1, skip=7)

make_ratio <- function(g1, g2) {
  # vectorized over g1, g2
  ifelse(is.na(g1) & is.na(g2), NA_character_,
         ifelse(is.na(g2), as.character(g1),
                ifelse(is.na(g1), as.character(g2),
                       paste0(g1, "/", g2))))
}


# Keep only the first 10 columns and rename
df <- raw %>%
  select(1:9) %>%
  rename(
    Year   = 1,
    Month  = 2,
    Day    = 3,
    Hour   = 4,
    Grad1  = 5,          # original Grad1
    Grad2  = 6,          # original Grad2
    BuchA  = 7,          # Buchstabe...7 (sign A)
    GradB  = 8,          # "Grad" alt column
    BuchB  = 9           # Buchstabe...9 (sign B)
  ) %>%
  mutate(across(where(is.character), ~ na_if(trimws(.x), ""))) %>%
  fill(Year, Month, Day, .direction = "down") %>%
  mutate(
    Year = as.integer(Year),
    Month = as.integer(Month),
    Day = as.integer(Day),
    Hour = as.integer(Hour),
    Minute = NA_integer_,
    
    # numeric versions of the grads
    Grad1_num = clean_num(Grad1),
    Grad2_num = clean_num(Grad2),
    GradB_num = clean_num(GradB),
    
    ta = make_ratio(Grad1_num, Grad2_num),
   
    meta.ta = paste0(
      "orig.time=", Hour,
      " | orig.ta.A=", ta,
      ifelse(!is.na(BuchA), paste0(" | sign.A=", BuchA), ""),
      ifelse(!is.na(GradB_num), paste0(" | orig.ta.B=", GradB_num), ""),
      ifelse(!is.na(BuchB), paste0(" | sign.B=", BuchB), "")
    )
  )
head(df)

# save
write_sef_f(
  Data=as.data.frame(df[, c("Year", "Month", "Day", "Hour", "Minute", "ta")]),
  outfile=paste0("Nuernberg_",get_date_range(df), "_ta_subdaily.tsv"),
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod="Kanold_Nuernberg",
  lat=lat, lon=lon, alt=alt,
  variable="ta",
  nam=name,
  meta=df$meta.ta,
  metaHead=metahead,
  units="unknown", sou="PALAEO-RA", stat="point",keep_na = F,
  time_offset = time.offset
)



# 1727_part2 ---------------------------------------------------------------

file <- "/scratch3/PALAEO-RA/daily_data/original/Nuernberg/Kanold_T3_DE_Nuernberg_1727_subdaily_part2.xls"
raw <- read_excel(file, sheet=1, skip=7)

# Keep only the first 10 columns and rename
df <- raw %>%
  select(1:9) %>%
  rename(
    Year   = 1,
    Month  = 2,
    Day    = 3,
    Hour   = 4,
    Grad1  = 5,          # original Grad1
    Grad2  = 6,          # original Grad2
    BuchA  = 7,          # Buchstabe...7 (sign A)
    GradB  = 8,          # "Grad" alt column
    BuchB  = 9           # Buchstabe...9 (sign B)
  ) %>%
  mutate(across(where(is.character), ~ na_if(trimws(.x), ""))) %>%
  fill(Year, Month, Day, .direction = "down") %>%
  mutate(
    Year = as.integer(Year),
    Month = as.integer(Month),
    Day = as.integer(Day),
    Hour = as.integer(Hour),
    Minute = NA_integer_,
    
    # numeric versions of the grads
    Grad1_num = clean_num(Grad1),
    Grad2_num = clean_num(Grad2),
    GradB_num = clean_num(GradB),
    
    ta = make_ratio(Grad1_num, Grad2_num),
    
    meta.ta = paste0(
      "orig.time=", Hour,
      " | orig.ta.A=", ta,
      ifelse(!is.na(BuchA), paste0(" | sign.A=", BuchA), ""),
      ifelse(!is.na(GradB_num), paste0(" | orig.ta.B=", GradB_num), ""),
      ifelse(!is.na(BuchB), paste0(" | sign.B=", BuchB), "")
    )
  )
head(df)

# save
write_sef_f(
  Data=as.data.frame(df[, c("Year", "Month", "Day", "Hour", "Minute", "ta")]),
  outfile=paste0("Nuernberg_",get_date_range(df), "_ta_subdaily.tsv"),
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod="Kanold_Nuernberg",
  lat=lat, lon=lon, alt=alt,
  variable="ta",
  nam=name,
  meta=df$meta.ta,
  metaHead=metahead,
  units="unknown", sou="PALAEO-RA", stat="point",keep_na = F,
  time_offset = time.offset
)


##############
# 1728 & 1729-1730 ---------------------------------------------------------------
##############

files <- c(
  "/scratch3/PALAEO-RA/daily_data/original/Nuernberg/Kanold_T3_DE_Nuernberg_1728_subdaily.xls",
  "/scratch3/PALAEO-RA/daily_data/original/Nuernberg/Kanold_T3_DE_Nuernberg_1729-1730_subdaily.xls"
)

# read + process all files into one df
dfs <- lapply(files, function(file) {
  raw <- read_excel(file, sheet=1, skip=7)
  
  raw %>%
    select(1:6) %>% # only the relevant first 6 cols
    rename(
      Year   = 1,
      Month  = 2,
      Day    = 3,
      Hour   = 4,
      Grad   = 5,
      sign   = 6
    ) %>%
    mutate(across(where(is.character), ~ na_if(trimws(.x), ""))) %>%
    fill(Year, Month, Day, .direction = "down") %>%
    mutate(
      Year   = as.integer(Year),
      Month  = as.integer(Month),
      Day    = as.integer(Day),
      Hour   = clean_hour(Hour),
      Minute = NA_integer_,
      ta     = round(clean_num(Grad), 1),
      meta.ta = paste0(
        "orig.time=", Hour,
        " | orig.ta=", ta,
        ifelse(!is.na(sign), paste0(" | sign=", sign), "")
      )
    )
})

# bind and sort
df <- bind_rows(dfs) %>%
  arrange(Year, Month, Day, Hour)

head(df)
tail(df)

# save
write_sef_f(
  Data=as.data.frame(df[, c("Year", "Month", "Day", "Hour", "Minute", "ta")]),
  outfile=paste0("Nuernberg_",get_date_range(df), "_ta_subdaily.tsv"),
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod="Kanold_Nuernberg",
  lat=lat, lon=lon, alt=alt,
  variable="ta",
  nam=name,
  meta=df$meta.ta,
  metaHead=metahead,
  units="unknown", sou="PALAEO-RA", stat="point",keep_na = F,
  time_offset = time.offset
)




##############
# 1732-1736 & 1737-1740 & 1741-1743 ---------------------------------------------------------------
##############

files <- c(
  "/scratch3/PALAEO-RA/daily_data/original/Nuernberg/Europe_T3_DE_Nuernberg_1732-1736_subdaily.xlsx",
  "/scratch3/PALAEO-RA/daily_data/original/Nuernberg/Europe_T3_DE_Nuernberg_1737-1740_subdaily.xls",
  "/scratch3/PALAEO-RA/daily_data/original/Nuernberg/Europe_T3_DE_Nuernberg_1741-1743_subdaily.xls"
)

# read + process all files into one df
dfs <- lapply(files, function(file) {
  raw <- read_excel(file, sheet=1, skip=6)
  
  raw %>%
    select(1:3, 13:15) %>% # only the relevant first 6 cols
    rename(
      Year   = 1,
      Month  = 2,
      Day    = 3,
      Grad   = 4,
      part   = 5,
      sign   = 6
    ) %>%
    mutate(across(where(is.character), ~ na_if(trimws(.x), ""))) %>%
    fill(Year, Month, Day, .direction = "down") %>%
    mutate(
      Year   = as.integer(Year),
      Month  = as.integer(Month),
      Day    = as.integer(Day),
      Hour   = NA_integer_,
      Minute = NA_integer_,
      ta     = paste0(clean_num(Grad), "|",part),
      meta.ta = paste0(
        "orig.ta=", ta,
        ifelse(!is.na(sign), paste0(" | sign=", sign), "")
      )
    ) %>%
    select("Year", "Month", "Day", "Hour", "Minute", "ta", "meta.ta")
})

# bind and sort
df <- bind_rows(dfs) %>%
  arrange(Year, Month, Day)

head(df)
tail(df)

# save
write_sef_f(
  Data=as.data.frame(df[, c("Year", "Month", "Day", "Hour", "Minute", "ta")]),
  outfile=paste0("Nuernberg_",get_date_range(df), "_ta_daily.tsv"),
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod="Europe_Nuernberg",
  lat=lat, lon=lon, alt=alt,
  variable="ta",
  nam=name,
  link="https://www.astronomie-nuernberg.de/index.php?category=doppelmayr&page=doppelmayr-observationes-meteorologicas",
  meta=df$meta.ta,
  metaHead="Observer=Doppelmayr",
  units="unknown", sou="PALAEO-RA", stat="point",keep_na = F,
  time_offset = time.offset
)





