rm(list=ls())
source('/scratch2/ccorbella/code/dataprep/helpfun.R')
library(dplyr)
library(dataresqc)
library(readxl)
library(dplyr)
library(purrr)
library(stringr)
library(lubridate)

outdir <- '/scratch3/PALAEO-RA/daily_data/final/Valencia'
indir <- '/scratch3/PALAEO-RA/daily_data/original/Valencia/'

files <- list.files(indir, pattern = "\\.xlsx?$|\\.xls$", full.names = TRUE)

library(janitor)

std_one <- function(df) {
  df <- df %>% janitor::clean_names()
  
  # helpers
  pick1 <- function(...) {
    x <- c(...)
    x[x %in% names(df)][1]
  }
  
  day   <- pick1("dia", "d_a", "day")
  month <- pick1("mes", "month")
  year  <- pick1("ano", "anho", "a_o", "year")
  hour  <- pick1("hora", "hour")
  
  wind_txt <- pick1("vientos")          # text wind (if exists)
  dd_txt   <- pick1("vien_let")         # orig.dd
  w_num    <- pick1("vien_num")         # orig.w
  sky      <- pick1("est_del_ciel")     # sky text
  
  df %>%
    transmute(
      Year   = as.integer(.data[[year]]),
      Month  = as.integer(.data[[month]]),
      Day    = as.integer(.data[[day]]),
      Hour   = as.integer(.data[[hour]]),
      Minute = 0L,
      
      Term   = suppressWarnings(as.numeric(.data[["term"]])),
      Bar_p  = suppressWarnings(as.numeric(.data[["bar_p"]])),
      Bar_l  = suppressWarnings(as.numeric(.data[["bar_l"]])),
      Hygro1 = suppressWarnings(as.numeric(.data[["higro1"]])),
      SecHum = suppressWarnings(as.numeric(.data[["sec_hum"]])),
      
      # wind fields you asked for
      orig_dd = coalesce(
        if (!is.na(dd_txt)) as.character(.data[[dd_txt]]) else NA_character_,
        if (!is.na(wind_txt)) as.character(.data[[wind_txt]]) else NA_character_
      ),
      orig_w  = if (!is.na(w_num)) suppressWarnings(as.numeric(.data[[w_num]])) else NA_real_,
      
      Sky = if (!is.na(sky)) as.character(.data[[sky]]) else NA_character_
    )
}


read_one <- function(f) {
  df <- read_excel(f, .name_repair = "minimal")  # keeps accents like "Día"

  df <- std_one(df) %>%
    mutate(
      file = basename(f),
      date = lubridate::make_date(Year, Month, Day)
    ) %>%
    select(file, date, everything())
  df
}

all_df <- map_dfr(files, read_one)

# drop when everything is NA
all_df <- all_df %>%
  filter(!(is.na(Year) & is.na(Month) & is.na(Day)))



# Create a list (equivalent to a dictionary in R)
meta <- list(
  ID = "Valencia",
  Name = "DominguezCastro_Valencia",
  lat = 39.47,
  lon = -0.38,
  alt = 25,
  Source = "Dominguez_Castro_2014",
  Link = "https://docta.ucm.es/entities/publication/b26dac98-5ffe-482c-9419-2f94168cc7eb"
)

df.ta <- all_df %>%
  mutate(
    meta=meta_time(Hour, Minute)
  ) %>% select(Year, Month, Day, Hour, Minute, Term, meta)

df.p <- all_df %>%
  mutate(
    p = round(convert_pressure(as.numeric(Bar_p) + as.numeric(Bar_l)/12, f=27.07, lat=meta[["lat"]], alt=meta[["alt"]],
                                   atb=Term),2),
    meta=paste0(meta_time(Hour, Minute), " | orig=", Bar_p, "in", Bar_l, "l | atb=", Term, "C")
  ) %>% select(Year, Month, Day, Hour, Minute, p, meta)

# clean winds for Valenica
dd_preclean <- function(x) {
  y <- toupper(trimws(x))
  y <- gsub("^M\\.?\\s*[0-9]+\\.?$", "M", y)
  # drop numbers, dots, stray punctuation
  y <- gsub("[0-9\\.]", "", y)
  y <- gsub("-", "", y)
  
  recode(y,
         # Valencian / Spanish names
         "L"      = "E",   "LEV"    = "E",   "LEVECHE" = "SE",
         "P"      = "W",   "PON"    = "W",
         "M"      = "S",   "MIGJORN"= "S",
         "T"      = "N",   "TRAMUNTANA" = "N",
         "TRE"    = "N",
         
         # French / mixed spellings
         "OUESTE" = "W",
         "OUEST"  = "W",
         "NOR"    = "N",
         "NORT"   = "N",
         "SUD"    = "S",
         "EST"    = "E",
         "ESTNOR" = "ENE",
         "SUDEST"="SE",
         
         # compound textual forms
         "ESTNORD" = "ENE",
         "NOREST"  = "NE",
         "NORDESTE" = "NE",
         "ESTSUD"  = "ESE",
         "EES" = "ESE",
         "SUDES"   = "SE",
         "NORNOR"  = "NNW",
         "SUDSUD"  = "SSW",
         "OUESTSUD"= "WSW",
         "ESTE" ="E",
         "NOROUE"  = "NW",
         "ESTSUDE" = "SE",
         
         # weird abbreviations
         "ON"  = "WNW",
         "NON" = "NNW",
         "ONO" = "WNW",
         "OSO" = "WSW",
         "NOO" = "NW",
         "OO"  = "W",
         "EE"  = "E",
         "OS" = "SW",
         "EN" = "NE",
         "SSS" = "S",
         "SP" = "SW",
         "MNE" = "E",
         "SS" ="S",
         "ME" = "SE",
         
         .default = y,
         .missing = NA_character_
  )
}


df.dd <- all_df %>%
  mutate(
    orig_dd_clean = dd_preclean(orig_dd),
    dd.norm = dd_normalize(orig_dd_clean),
    dd = dd2deg(dd.norm),
    meta=paste0(meta_time(Hour, Minute), " | orig=", orig_dd)
  ) %>% select(Year, Month, Day, Hour, Minute, dd, meta)

# check which winds didn't get converted
# (need to select orig_dd and orig_dd_clean columns beforehand)
df.dd %>%
  filter(is.na(dd), !is.na(orig_dd)) %>%
  select(Year, Month, Day, dd, orig_dd, orig_dd_clean, meta) 


# ta
var <- "ta"
write_sef_f(Data=as.data.frame(df.ta),
            outfile=outfile.name(meta[["Name"]], var, df.ta, subdaily=TRUE),
            outpath=outdir,
            cod=meta[["ID"]],
            variable=var,
            nam=meta[["Name"]],
            lat=meta[["lat"]],
            lon=meta[["lon"]], alt=meta[["alt"]], sou=meta[["Source"]],
            
            link=meta[["Link"]], units=units(var), stat="point",
            time_offset = time.offset(meta[["lon"]]),
            meta=df.ta$meta, keep_na = T)

# p
var<-"p"
write_sef_f(Data=as.data.frame(df.p),
            outfile=outfile.name(meta[["Name"]], var, df.p, subdaily=TRUE),
            outpath=outdir,
            cod=meta[["ID"]],
            variable=var,
            nam=meta[["Name"]],
            lat=meta[["lat"]],
            lon=meta[["lon"]], alt=meta[["alt"]], sou=meta[["Source"]],
            units=units(var),
            link=meta[["Link"]], stat="point",
            time_offset=time.offset(meta[["lon"]]),
            meta=df.p$meta, keep_na =TRUE)

# dd
var<-'dd'
write_sef_f(Data=as.data.frame(df.dd),
            outfile=outfile.name(meta[["Name"]], var, df.dd, subdaily=TRUE),
            outpath=outdir,
            cod=meta[["ID"]],
            variable=var,
            nam=meta[["Name"]],
            lat=meta[["lat"]],
            lon=meta[["lon"]], alt=meta[["alt"]], sou=meta[["Source"]],
            units=units(var),
            link=meta[["Link"]], stat="point",
            time_offset=time.offset(meta[["lon"]]),
            meta=df.dd$meta, keep_na =TRUE)

