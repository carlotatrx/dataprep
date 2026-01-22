rm(list=ls())
library(dataresqc)
library(readxl)
library(tidyr)
library(dplyr)
library(lubridate)
library(tibble)
library(stringr)
source('/scratch2/ccorbella/code/dataprep/helpfun.R')

outdir <- "/scratch3/PALAEO-RA/daily_data/final/Raciborz"


# new data from Rajmund email ---------------------------------------------

name <- "Przybylak_Raciborz"
code <- "Raciborz"
lat <- round(50+22/60,4)
lon <- round(18+11/60,4)
alt	<- 206.8
source <- "Programm der Königlichen Gymnasien zu Ratibor für die Zeit von Ostern 1879 bis dahin 1880, Ratibor 1880; Die meteorologischen Verhältnisse von Ratibor, auf Grund der Monatstabellen der meteorologischen"
link <- ""

file <- "/scratch3/PALAEO-RA/daily_data/original/Raciborz/Raciborz_Temp_1848-1879.xlsx"

raw <- read_excel(file,
                  sheet=1,
                  skip =2,
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
    meta = paste0("orig=", Value, "R"),
    Value = round(as.numeric(Value)*1.25,1),
    
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
  link    = "",
  nam     = name,
  var     = var,
  stat    = "mean",
  period  = "day",
  units   = units(var),
  meta    = df2$meta,
  keep_na = TRUE
)
