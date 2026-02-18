rm(list=ls())
library(dplyr)
library(readxl)
library(purrr)
library(tidyr)
library(dataresqc)
source('/scratch2/ccorbella/code/dataprep/helpfun.R')

outdir <- "/scratch3/PALAEO-RA/daily_data/final/Avignon"


lat <- 43.949
lon <- 4.811
alt <- 51
name <- "Yuri_Avignon"
code <- "Avignon"
source <- "Guérin, J. (1839). Observations météorologiques faites à Avignon; suivies d'un tableau monographique des taches du soleil. Imprimerie de Jacquet et Joudou. Archives nationales, Observations climatologiques et pluviométriques des postes climatologiques bénévoles, classés par département : Seine à Yonne (Cote: 19820184/73/2)."
link   <- ""
raw <- read_excel("/scratch3/PALAEO-RA/daily_data/original/Avignon/GUERIN TABLE  YEAR 1816.xls",
                  sheet=2,na="--",
                  range = "A6:L1000",
                  col_names = c("mo_name", "Day", "ta6", "ta14", "no1", "ta10", "ta12", "no2", "porig6", "porig10", "porig12", "porig14"))

head(raw)

raw <- raw %>%
  fill(mo_name) %>%                               # 1. fill down
  mutate(
    Month = recode(mo_name,                       # 2. month → number
                   "JANV" = 1,  "FEV" = 2,  "MARS" = 3,
                   "AVRIL" = 4, "MAI" = 5,  "JUIN" = 6,
                   "JUIL" = 7,  "AOÛT" = 8, "SEPT" = 9,
                   "OCT" = 10,  "NOV" = 11, "DEC" = 12
    )
  )

num_cols <- c("ta6","ta14","ta10","ta12",
              "porig6","porig10","porig12","porig14")

raw <- raw %>%
  mutate(across(all_of(num_cols),
                ~ round(as.double(.), 2)))


df <- raw %>%
  pivot_longer(
    cols = num_cols, 
    names_to = c(".value", "Hour"), 
    names_pattern = "(porig|ta)(\\d+)"
  ) %>% mutate(
    Year=1816L,
    Hour = as.integer(Hour),
    Minute = 0L,

    p = round(convert_pressure(porig, lat=lat, alt=alt, atb=ta),2),

    meta.p = paste0(meta_time(Hour, Minute), " | orig=",porig,"mm | atb=", ta, "C"),
    meta.ta = paste0(meta_time(Hour, Minute)),
  ) %>% select(Year, Month, Day, Hour, Minute, p, ta, meta.ta, meta.p)

# remove empty rows at the bottom
df <- df %>%
  filter(!(is.na(ta) & is.na(p)))

# remove rows without days (monthly /yearly means in the excel)
df <- df %>%
  filter(!is.na(Day))


head(df)

# save
var <-"p"
write_sef_f(
  as.data.frame(df[,c("Year","Month", "Day", "Hour", "Minute", "p")]),
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
  metaHead = "PTC=Y | PGC=Y",
  meta    = df$meta.p,
  time_offset = time.offset(lon),
  keep_na = TRUE
)

var <-"ta"
write_sef_f(
  as.data.frame(df[,c("Year","Month", "Day", "Hour", "Minute", "ta")]),
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
  meta    = df$meta.ta,
  time_offset = time.offset(lon),
  keep_na = TRUE
)

