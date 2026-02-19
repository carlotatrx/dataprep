rm(list=ls())
library(dataresqc)
library(readxl)
library(tidyr)
library(dplyr)
source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')


# Vienna from Lucas -------------------------------------------------------

name <- "Vienna"
code <-"WIE"

outdir <- paste0("/scratch3/PALAEO-RA/daily_data/final/", name)

raw <- read.csv(paste0("/scratch3/PALAEO-RA/daily_data/original/", name, "/", code,"_SFP.csv"), na=c("","NA"))

head(raw)

lat <- 48.249
lon <- 16.356
alt <- 198

df <- raw %>%
  mutate(
    year=year(dates),
    month=month(dates),
    day=day(dates),
    Hour=NA,
    minute=NA,
  )

df <- df %>%
  filter(!is.na(SFP))
head(df)


var<-"p"
write_sef_f(
  as.data.frame(df[c("year","month", "day", "Hour", "minute","SFP")]),
  outfile = outfile.name(paste0("Pfister_",name), var, df, FALSE),
  outpath = outdir,
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = "Pfister, L., Wilhelm, L., Brugnara, Y., Imfeld, N., & Brönnimann, S. (2024). Weather type reconstruction using machine learning approaches. EGUsphere, 2024, 1-33.",
  link    = "https://doi.org/10.5194/wcd-6-571-2025",
  nam     = name,
  var     = var,
  stat    = "mean",
  period    = "day",
  units   = units(var),
  metaHead = "PTC=Y | PGC=Y | homogenized=Y | QC=Y",
)


# Vienna from PALAEO-RA ---------------------------------------------------

name <- "Vienna"
lat	<- 48.2081743
lon	<- 16.3738189
alt	<- NA
metaHead <- "Observer=Beintema van Peima"

file <- "/scratch3/PALAEO-RA/daily_data/original/Vienna/Europe_T3_AT_Vienna_1709-1715_subdaily.xls"
raw <- read_excel(file, sheet=3, skip=6)

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


# sheet number 3 ----------------------------------------------------------

# Keep only the first columns and rename
df <- raw %>%
  select(1:3, 6:9) %>%
  rename(
    Year = 1, Month = 2, Day = 3, ta.sign.morn = 4, ta.sign.eve = 6, ta.morn = 5, ta.eve = 7
  ) %>%
  # Fill down Year/Month/Day
  mutate(across(where(is.character), ~ na_if(trimws(.x), ""))) %>%
  fill(Year, Month, .direction = "down") %>%
  mutate(
    Year = as.integer(Year),
    Month = as.integer(Month),
    Day = as.integer(Day),
    Minute = NA_integer_,
    ta.sign.morn = as.character(ta.sign.morn),
    ta.sign.eve  = as.character(ta.sign.eve)
    
  )
head(df)

df.long <- df %>%
  pivot_longer(
    cols = c(ta.morn, ta.eve),
    names_to = "tod",
    names_prefix = "ta\\.",
    values_to = "ta"
  ) %>%
  mutate(
    ta = round(ta,1),
    Hour = case_when(
      tod=="morn" ~ 8L,
      tod=="eve"  ~ 20L,
      TRUE        ~ NA_integer_
    ),
    orig.hour = case_when(
      tod=="morn" ~ "a.m.",
      tod=="eve"  ~ "p.m.",
      TRUE        ~ NA_character_
    ),
    ta.orig = case_when(
      !is.na(ta.sign.morn) ~ paste0(ta.sign.morn, ta),
      !is.na(ta.sign.eve)  ~ paste0(ta.sign.eve,  ta),
      TRUE                 ~ as.character(ta)
    ),
    ta = case_when(
      tod == "morn" & grepl("f", tolower(ta.sign.morn)) ~ -ta,
      tod == "eve"  & grepl("f", tolower(ta.sign.eve))  ~ -ta,
      TRUE ~ ta
    ),
    meta = paste0("orig.hour=", orig.hour, " | orig.ta=", ta.orig)
  ) %>%
  select(Year, Month, Day, Hour, Minute, Value = ta, meta) %>%
  filter(!is.na(Value))
head(df.long)


# sheet number 4 ----------------------------------------------------------

raw2 <- read_excel(file, sheet=4, skip=6)

df2 <- raw2 %>%
  select(1:4, 7:10) %>%
  rename(
    Year = 1, Month = 2, Day = 3, orig.hour = 4, ta.sign.morn = 5, ta.sign.eve = 7, ta.morn = 6, ta.eve = 8
  ) %>%
  # Fill down Year/Month/Day
  mutate(across(where(is.character), ~ na_if(trimws(.x), ""))) %>%
  fill(Year, Month, Day, .direction = "down") %>%
  mutate(
    Year = as.integer(Year),
    Month = as.integer(Month),
    Day = as.integer(Day),
    orig.hour = suppressWarnings(as.integer(orig.hour)),
    Minute = NA_integer_,
    ta.sign.morn = as.character(ta.sign.morn),
    ta.sign.eve  = as.character(ta.sign.eve),
    ta.morn = round(suppressWarnings(as.numeric(ta.morn)),1),
    ta.eve  = round(suppressWarnings(as.numeric(ta.eve)), 1)
  )
head(df2)

df2.long <- df2 %>%
  pivot_longer(cols = c(ta.morn, ta.eve), names_to = "tod", values_to = "ta") %>%
  mutate(
    rounding = if_else(tod == "ta.morn", ta.sign.morn, ta.sign.eve),
    Hour = case_when(
      tod == "ta.morn" & orig.hour <= 12 ~ orig.hour,        # keep as-is
      tod == "ta.eve"  & orig.hour < 12  ~ orig.hour + 12,    # convert to PM
      tod == "ta.eve"  & orig.hour >= 12 ~ orig.hour,         # already PM
      TRUE                          ~ NA_integer_
    ),

    meta.hour = if_else(tod=="ta.morn","g.a.","p.m."),
    orig.ta = if_else(!is.na(rounding), paste0(rounding, ta), as.character(abs(ta))),
    ta = case_when(
      tod == "ta.morn" & grepl("f", tolower(ta.sign.morn)) ~ -ta,
      tod == "ta.eve"  & grepl("f", tolower(ta.sign.eve))  ~ -ta,
      TRUE ~ ta
    ),
    meta = paste0("orig.hour=", orig.hour, meta.hour, " | orig.ta=", orig.ta)
  ) %>%
  select(Year, Month, Day, Hour, Minute, Value = ta, meta) %>%
  filter(!is.na(Value))

head(df2.long)


# sheet number 2 ----------------------------------------------------------

raw3 <- read_excel(file, sheet=2, skip=6)

df3 <- raw3 %>%
  select(1:7) %>%
  rename(
    Year = 1, Month = 2, Day = 3, F.morn = 4, F.eve = 5, C.morn = 6, C.eve = 7
  ) %>%
  # Fill down Year/Month/Day
  mutate(across(where(is.character), ~ na_if(trimws(.x), ""))) %>%
  fill(Year, Month, .direction = "down") %>%
  mutate(
    
    Year = as.integer(Year),
    Month = as.integer(Month),
    Day = as.integer(Day),
    C.eve = as.numeric(C.eve),
    ta.morn = case_when(
      F.morn > 0 ~ -F.morn,
      C.morn > 0 ~ C.morn,
      TRUE ~ NA_real_
    ),
    ta.eve = case_when(
      F.eve > 0 ~ -F.eve,
      C.eve > 0 ~ C.eve,
      TRUE ~ NA_real_
    ),
    ta.sign.morn = case_when(
      F.morn > 0 ~ "F",
      C.morn > 0 ~ "C",
      TRUE ~ NA_character_
    ),
    ta.sign.eve = case_when(
      F.eve > 0 ~ "F",
      C.eve > 0 ~ "C",
      TRUE ~ NA_character_
    ),
    
    Minute = NA_integer_
  )
head(df3)

df3.long <- df3 %>%
  pivot_longer(cols = c(ta.morn, ta.eve), names_to = "tod", values_to = "ta") %>%
  mutate(
    rounding = if_else(tod == "ta.morn", ta.sign.morn, ta.sign.eve),
    Hour = case_when(
      tod == "ta.morn" ~ 8L,       
      tod == "ta.eve"  ~ 20L,
      TRUE             ~ NA_integer_
    ),
    
    meta.hour = if_else(tod=="ta.morn","a.m.","p.m."),
    orig.ta = if_else(!is.na(rounding), paste0(rounding, ".", abs(ta)), as.character(abs(ta))),

    meta = paste0("orig.hour=", meta.hour, " | orig.ta=", orig.ta)
  ) %>%
  select(Year, Month, Day, Hour, Minute, Value = ta, meta) %>%
  filter(!is.na(Value))

head(df3.long)


# combine and save --------------------------------------------------------

df.all <- bind_rows(df.long, df2.long, df3.long) %>%
  mutate(Value=round(Value, 1)) %>%
  arrange(Year, Month, Day, Hour)

# save
write_sef_f(
  Data=as.data.frame(df.all[, c("Year", "Month", "Day", "Hour", "Minute", "Value")]),
  outfile=paste0("Vienna_",get_date_range(df.all), "_ta_subdaily.tsv"),
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod="Europe_Vienna_1",
  lat=lat, lon=lon, alt=alt,
  variable="ta",
  nam=name,
  meta=df.all$meta,
  metaHead=metaHead,
  units="unknown", sou="PALAEO-RA", stat="point",keep_na = F,
)
