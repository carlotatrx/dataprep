rm(list=ls())
library(dataresqc)
library(readxl)
library(tidyr)
library(dplyr)
library(lubridate)
source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')

outdir <- "/scratch3/PALAEO-RA/daily_data/final/Cartagena_Southern_Spain"
name <- "Cartagena"
code <- "Car1804-1807"
lat	<- round(37+36/60,4)
lon	<- round(-(0+58/60),4)
alt	<- 67
source <- "Diario de Cartagena. 1804-1807. Archivo Municipal de Murcia, 1-G-13/14, (http://archivodemurcia.es)."
link   <- "https://repositorio.ual.es/handle/10835/6806"

file <- list.files("/scratch3/PALAEO-RA/daily_data/original/Southern_Spain", full.names=TRUE)
raw1 <- read_excel(file, sheet="Car1804-1806", skip=8)
raw2 <- read_excel(file, sheet="Car1807", skip=8)

# need to convert to character for consistency
raw2$`Tmor(F)` <- as.character(raw2$`Tmor(F)`)
raw2$`Tmid(F)` <- as.character(raw2$`Tmid(F)`)
raw2$`Tnig(F)` <- as.character(raw2$`Tnig(F)`)

raw <- bind_rows(raw1, raw2) %>%
  arrange(Year, Month, Day)

Tcols <- c("Tmor", "Tmid", "Tnig")

# specify when R and when F
df <- raw %>%
  rename(
    Tmor = "Tmor(F)",
    Tmid = "Tmid(F)",
    Tnig = "Tnig(F)",
    prmor = "pmor(EI)",
    prmid = "pmid(EI)",
    prnig = "pnig(EI)"
  ) %>%
  mutate(across(all_of(Tcols), \(x) gsub(",", ".", x))) %>%               # fix decimal comma
  mutate(across(all_of(Tcols), \(x) trimws(x))) %>%
  mutate(across(all_of(Tcols), \(x) {
    ifelse(grepl("\\*$", x),
           as.numeric(sub("\\*$", "", x)) * 2.25 + 32,   # R → °F conversion
           as.numeric(x))
  }, .names = "F_{.col}"))

vars_period <- grep("(mor|mid|nig)$", names(df), value = TRUE)

df_long <- df %>%
  pivot_longer(
    cols = all_of(vars_period),
    names_to = "varperiod",
   # names_pattern = "^(T|p|W|WF|Q|F_T)(mor|mid|nig)$",
    values_to = "value",
    values_transform = as.character
  ) %>%
  mutate(
    period = case_when(
      grepl("mor$", varperiod) ~ "mor",
      grepl("mid$", varperiod) ~ "mid", 
      grepl("nig$", varperiod) ~ "nig"
    ),
    Hour = case_when(
      period == "mor" ~ 8L,
      period == "mid" ~ 12L,
      period == "nig" ~ 21L
    ),
    Minute = 0L,
    variable = gsub("(mor|mid|nig)$", "", varperiod)
  ) %>%
  select(-varperiod) %>%
  pivot_wider(
    names_from = variable,
    values_from = value
  )

df_final <- df_long %>% 
  mutate(
    obs.time = case_when(
      period == "mor" ~ "morning",
      period == "mid" ~ "midday",
      period == "nig" ~ "night",
      TRUE   ~  NA_character_
    ),
    ta = round((as.numeric(F_T) -32)/1.8, 1), # convert to Celsius
    p  = round(convert_pressure(as.numeric(pr), f=25.4, lat=lat, alt=alt, atb=ta), 1),

    dd.norm = dd_normalize(W),
    dd = dd2deg(dd.norm),

    meta.dd = paste0("obs.time=", obs.time, " | orig.dd=", W, " | force=", WF),
    meta.ta = paste0("obs.time=", obs.time, " | orig.ta=", 
                     gsub("\\*","", T), # remove the *
                     ifelse(Year==1804 & Month==8, "R", "C")),
    meta.p  = paste0("obs.time=", obs.time, " | orig.p=", pr, "English in | atb=", ta, "C"),
  ) %>%
  select(Year, Month, Day, Hour, Minute, ta, p, dd, meta.ta, meta.p, meta.dd)


head(df_final)

meta_map  <- list(ta = df_final$meta.ta, p = df_final$meta.p, dd = df_final$meta.dd)
base_cols <- df_final[c("Year","Month","Day","Hour","Minute")]
vars <- c("ta", "p", "dd")

for (var in vars) {
  meta  <- meta_map[[var]]
  dat   <- cbind(base_cols, setNames(df_final[var], var))
  write_sef_f(
    dat,
    outfile = outfile.name(name, var, df_final, subdaily=TRUE),
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
    metaHead= ifelse(var=="p", "PTC=Y | PGC=Y", ""),
    meta    = meta,
    keep_na = TRUE
  )
}
