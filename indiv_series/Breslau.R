rm(list=ls())
library(dataresqc)
library(readxl)
library(tidyr)
library(dplyr)
library(lubridate)
library(tibble)
library(stringr)
source('/scratch2/ccorbella/code/dataprep/helpfun.R')

outdir <- "/scratch3/PALAEO-RA/daily_data/final/Breslau"


# new data from Rajmund email ---------------------------------------------

name <- "Przybylak_Breslau"
code <- "Wroclaw"
lat <- 51.11
lon <- 17.04
alt	<- 117
source <- "Galle J.G., 1857, Grundzüge der Schlesischen Klimatologie, Josef Max & Komp., Breslau, XXIII+128 pp."
link <- ""
metaHead <- "coords=approx"

file <- "/scratch3/PALAEO-RA/daily_data/original/Breslau/Wroclaw_Daily_air_temperature_1791-1854.xlsx"

raw <- read_excel(file,
                  sheet=1,
                  skip =1,
                  na="*")

df <- raw %>%
  mutate(Year=as.integer(Year)) %>%
  fill(Year, .direction="down")  %>%
  pivot_longer(
    cols = c(Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec),
    names_to  = "Month",
    values_to = "Value"
  ) %>%
  mutate(
    Month = match(Month, month.abb),
    Day = as.integer(Day),
    Hour=NA_integer_,
    Minute=NA_integer_,
    Date = make_date(Year, Month, Day),
    Value = as.numeric(Value),
    
  ) %>% arrange(Year, Date, Month)

df

# sanity checks
# all years should be valid
stopifnot(nrow(df %>% filter(is.na(Year)))==0)

# check false dates
false_dates <- df %>%
  filter(is.na(Date)) %>%
  select(Year, Month, Day) %>%
  distinct() %>%
  arrange(Year, Month, Day)
false_dates

false_dates_with_values <- df %>%
  filter(is.na(Date) & !is.na(Value)) %>%
  select(Year, Month, Day, Value) %>%
  arrange(Year, Month, Day)

# sanity check
# there shouldn't be any false date with a non-NA value
stopifnot(nrow(false_dates_with_values)==0)


# once all is good, drop non-valid Dates
df2 <- df %>%
  filter(!is.na(Date)) %>%
  select(!Date)
df2


# save
var <-"ta"
write_sef_f(
  as.data.frame(df2[,c("Year","Month", "Day", "Hour", "Minute", "Value"),]),
  outfile = outfile.name(name, var, df2, FALSE),
  outpath = outdir,
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "mean",
  period  = "day",
  units   = units(var),
  metaHead = metaHead,
  keep_na = TRUE
)


# 1773-1781 (Rajmund) -----------------------------------------------------

# Air temperature data in Wroc³aw (SW Poland) in 1773-1781 based on a newly discovered series of meteorological measurements
name <- "Przybylak_Wroclaw"
code <- "Breslau"
lat <- 51.11684
lon <- 17.03042
alt	<- 117
source <- "Przybylak R., Pospieszyñska A., Oliñski P., Air temperature changes in Wroc³aw (SW Poland) in 1773-1781 based on a newly discovered series of meteorological measurements, Climate of the Past."
link <- "https://doi.org/10.18150/PYVVWU"
metaHead <- "Observer=Johann Ephraim Scheibel | Location=area of the Gymnasium, by Church of St. Elizabeth | orig.unit=degF"

raw <- read.csv("/scratch3/PALAEO-RA/daily_data/original/Breslau/Wroclaw_air_temp_1773-1781.csv",
                header=TRUE, sep=";")

df <- raw %>%
  pivot_longer(
    cols = c(Morning, Midday, Evening),
    names_to = "tod",
    values_to = "ta"
  ) %>%
  mutate(
    Hour=case_when(
      tod=="Morning" ~ 7L,
      tod=="Midday" ~ 12L,
      tod=="Evening" ~ 20L,
      TRUE ~ NA
    ),
    Minute=0L,
    meta=paste0("orig.time=", tod)
  ) %>%
  select(Year, Month, Day, Hour, Minute, ta, meta)

# check hours are properly converted
stopifnot(sum(is.na(df$ta))==0)

head(df)

# save
var <-"ta"
write_sef_f(
  as.data.frame(df),
  outfile = outfile.name(name, var, df, TRUE),
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
  meta    = df$meta,
  metaHead = metaHead,
  keep_na = TRUE
)

# 1717-1726  --------------------------------------------------------------
name <- "Wroclaw"
code <- "Breslau"
lat	<- 51.11
lon	<- 17.04
alt	<- 117
source <- "PALAEO-RA"
link   <- ""

raw <- read_excel("/scratch3/PALAEO-RA/daily_data/original/Breslau/Kanold_T3_PL_Breslau_1717-1726_subdaily.xls",
                  na="-", skip=7,
                  col_names=c("Year", "Month", "Day", "Time", "dd.orig", "p", "ta", "ta_string", "drachm", "scrup", "gran", "meta.rr", "notes")
)
raw[raw=="NA"] <- NA
raw[raw==""]   <- NA

raw <- raw %>%
  fill(Year, Month, Day, .direction="down")
head(raw)

df <- raw %>%
  mutate(
    Hour = case_when(
      Time=="F" ~ 8L,
      Time=="N" ~ 12L,
      Time=="A" ~ 17L,
      Time=="5" ~ 17L,
      .default = as.numeric(Time)
    ),
    Minute=0L,
    dd = dd2deg(dd_normalize(dd.orig)),
    
    ta = ifelse(!is.na(ta_string), paste0(ta, ta_string), ta),
    
    rr = na_if(str_c(
      if_else(is.na(drachm), "", str_c(drachm, "drachm")),
      if_else(is.na(scrup),  "", str_c(scrup,  "scrup")),
      if_else(is.na(gran),   "", str_c(gran,   "gran"))
    ),""),
    
    meta.time = case_when(
      Time=="F" ~ "orig.time=morning",
      Time=="N" ~ "orig.time=afternoon",
      Time=="A" ~ "orig.time=evening",
      TRUE ~ NA_character_
    ),
    
    meta.rr = ifelse(
      !is.na(notes) & str_detect(notes, regex("nieder", ignore_case = TRUE)),
      notes,
      NA_character_
    ),
    
    meta.dd = paste0("orig.dd=",dd.orig),
    meta.dd = if_else(
      !is.na(notes) & str_detect(notes, regex("wind", ignore_case = TRUE)),
      paste0(meta.dd, " | ", notes),
      meta.dd
    ),
    meta.dd = if_else(
      !is.na(meta.time),
      paste0(meta.time, " | ", meta.dd),
      meta.dd
    )
    
  ) %>% select(Year, Month, Day, Hour, Minute, p, ta, dd, rr, meta.dd, meta.rr, meta.time)

head(df)

base_cols <- df[c("Year","Month","Day","Hour","Minute")]

for (var in c("p", "ta")) {
  dat   <- cbind(base_cols, setNames(df[var], var))
  write_sef_f(
    dat,
    outfile = outfile.name(name, var, dat, TRUE),
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
    meta    = df$meta.time,
    metaHead = "Observer=J.Kanold",
    time_offset = time.offset(lon),
    keep_na = FALSE
  )
}

# dd
var <-"dd"
dat <- cbind(base_cols, setNames(df[var], var))
write_sef_f(
  dat,
  outfile = outfile.name(name, var, dat, TRUE),
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
  meta    = df$meta.dd,
  metaHead = "Observer=J.Kanold",
  time_offset = time.offset(lon),
  keep_na = TRUE
)

# rr
var <-"rr"
dat <- cbind(base_cols, setNames(df[var], var),meta= df$meta.rr) %>%
  mutate(Hour=24L,Minute=24L)

head(dat)
write_sef_f(
  dat,
  outfile = outfile.name(name, var, dat, TRUE),
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
  meta    = dat$meta,
  metaHead = "Observer=J.Kanold | obs.time=midnight | 1 Drachm. was 3.727 g after the Nürnberg standard | 1 Drachm. = 3 Scrup. = 60 Gran.",
  keep_na = FALSE
)


# 1717-1718 (ta, p, rr) ---------------------------------------------------

raw <- read_excel("/scratch3/PALAEO-RA/daily_data/original/Breslau/Kanold_T3_PL_Breslau_1717-1718_subdaily.xlsx",
                  na="-", skip=8,
                  col_names=c("Year", "Month", "Day", "p1", "p2", "p3", "ta1", "ta2", "ta3", "drachm", "scrup", "gran", "notes")
                 )
head(raw)

df.rr <- raw %>%
  mutate(
    Hour=NA_integer_,
    Minute=NA_integer_
  ) %>% select(Year, Month, Day, Hour, Minute, drachm, scrup, gran, notes)

df.rr[df.rr=="NA"] <- NA

df.rr <- df.rr %>%
  mutate(rr = str_c(
    if_else(is.na(drachm), "", str_c(drachm, "drachm")),
    if_else(is.na(scrup),  "", str_c(scrup,  "scrup")),
    if_else(is.na(gran),   "", str_c(gran,   "gran"))
  )) %>% select(Year, Month, Day, Hour, Minute,rr, meta=notes)
df.rr[df.rr==""] <- NA

head(df.rr)

var<-"rr"
write_sef_f(
  as.data.frame(df.rr),
  outfile = outfile.name(name, var, df.rr, subdaily=FALSE),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Breslau',
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "day",
  units   = "unknown",
  meta    = df.rr$meta,
  metaHead="The opening of the rain gauge was “almost 2/3 of ¼ of ell” in diameter, i.e. ca. 9.5 cm? | 1 Drachm. was 3.727 g after the Nürnberg standard | 1 Drachm. = 3 Scrup. = 60 Gran.",
  keep_na = FALSE
)


df <- raw %>%
  pivot_longer(
    cols = c(p1, p2, p3, ta1, ta2, ta3),
    names_to = c(".value", "tod"),
    names_pattern = "^([a-z]+)([123])$",
  ) %>%
  mutate(
    Hour = case_when(
      tod==1 ~ 7,
      tod==2 ~ 12,
      tod==3 ~ 20
    ),
    Minute = 0L,
  ) %>% select(Year, Month, Day, Hour, Minute, p, ta)

head(df)

for (var in c("p", "ta")){
  dat  <- df[c(1:5, which(names(df)==var))]
  write_sef_f(
    as.data.frame(dat),
    outfile = outfile.name(name, var, dat, TRUE),
    outpath = '/scratch3/PALAEO-RA/daily_data/final/Breslau',
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
    time_offset = time.offset(lon),
    keep_na = FALSE
  )
}






