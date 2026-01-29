# Excel (multi-sheets, one per year) -> SEF (ta and p), 3 obs/day
# - fills down Year (and Month name)
# - pivots by time-of-day columns (Termín: cols 5–7)
# - ta:  °C = (°R + desetiny/10) * 1.25
# - p:   mmHg = (couly + body/12 + carky/120) * 27.07
#        then p_hPa = convert_pressure(mmHg, atb = ta_atb_C)   # needs your function
# - always keeps originals in Meta: orig.time=..., orig=...

rm(list = ls())

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(dataresqc)

source("/scratch2/ccorbella/code/dataprep/helpfun.R")

# ---- user inputs ----
infile <- "/scratch3/PALAEO-RA/daily_data/original/Zitenice/Žitenice (Kreybich)_original od PZ.xls"
outdir <- "/scratch3/PALAEO-RA/daily_data/final/Zitenice/"

code <- "Zitenice"
name <- "Brug_Zitenice"
lat <- 50.55
lon <- 14.16
alt <- 223
source <- "Yuri Brugnara"
link <- ""


# ---- helpers ----
fill_down <- function(x) {
  x %>%
    mutate(across(everything(), ~ na_if(.x, ""))) %>%
    tidyr::fill(everything(), .direction = "down")
}

# Termín values -> Hour (keep orig.time)
tod_to_hour <- function(tod) {
  case_when(
    tod == 7 ~ 7L,
    tod == 2 ~ 14L,
    tod == 9 ~ 21L,
    TRUE ~ as.integer(tod)
  )
}

# ---- parse ONE sheet ----
#  column layout (by position):
#  1 MonthName, 2 Year, 3 MonthNum, 4 Day,
#  5–7 Termín (tod1,tod2,tod3),
#  8–13 Thermometer: (R1,d1,R2,d2,R3,d3),
#  14–22 Barometer:  (c1,b1,k1,c2,b2,k2,c3,b3,k3),
#  23–28 Ther.am Bar.: (R1,d1,R2,d2,R3,d3)
#  29 Witterung, 30 Poznámky  (if present; ignored)
read_year_sheet <- function(sheet_name) {
  raw <- read_excel(infile, sheet = sheet_name, skip=2, col_names = FALSE, .name_repair = "minimal")
  names(raw) <- paste0("V", seq_along(raw))
  
  df <- raw %>%
    mutate(across(everything(), ~ str_squish(as.character(.x)))) %>%
    mutate(
      MonthName = V1,
      Year = suppressWarnings(as.integer(V2)),
      Month = suppressWarnings(as.integer(V3)),
      Day = suppressWarnings(as.integer(V4))
    ) %>%
    fill(MonthName, .direction = "down") %>%
    filter(!is.na(Year), !is.na(Month), !is.na(Day))
  
  # Termín columns are V5,V6,V7 (tod slots)
  tod_long <- df %>%
    transmute(Year, Month, Day,
              tod1 = suppressWarnings(as.integer(V5)),
              tod2 = suppressWarnings(as.integer(V6)),
              tod3 = suppressWarnings(as.integer(V7))) %>%
    pivot_longer(starts_with("tod"), names_to = "tod_slot", values_to = "tod_orig") %>%
    mutate(
      tod_orig = as.integer(tod_orig),
      Hour = case_when(
        tod_slot == "tod1" & is.na(tod_orig) ~ 7L,
        sheet_name == "1790" & tod_slot == "tod1" ~ 7L,
        sheet_name == "1790" & tod_slot == "tod2" ~ 14L,
        sheet_name == "1790" & tod_slot == "tod3" ~ 21L,
        tod_orig == 7 ~ 7L,
        tod_orig == 2 ~ 14L,
        tod_orig == 9 ~ 21L,
        TRUE ~ tod_orig
      ),
      Minute = 0L,
      meta_time = case_when(
        sheet_name == "1790" & tod_slot == "tod1" ~ "orig.time=NA",
        sheet_name == "1790" & tod_slot == "tod2" ~ 'orig.time="1-2"',
        sheet_name == "1790" & tod_slot == "tod3" ~ 'orig.time="10"',
        TRUE ~ paste0("orig.time=", tod_orig)
      )
    )
  
  # Thermometer: V8..V13 = (R1,d1,R2,d2,R3,d3)
  ta_long <- df %>%
    transmute(Year, Month, Day,
              R1 = suppressWarnings(as.integer(V8)),  d1 = suppressWarnings(as.integer(V9)),
              R2 = suppressWarnings(as.integer(V10)), d2 = suppressWarnings(as.integer(V11)),
              R3 = suppressWarnings(as.integer(V12)), d3 = suppressWarnings(as.integer(V13))) %>%
    pivot_longer(cols = c(R1,d1,R2,d2,R3,d3),
                 names_to = c("part","idx"),
                 names_pattern = "([Rd])([123])",
                 values_to = "val") %>%
    pivot_wider(names_from = part, values_from = val) %>%
    mutate(
      tod_slot = paste0("tod", idx),
      ta = round((R + d/10) * 1.25,2),
      meta_ta = paste0("orig=", R, ".", d, "R")
    ) %>%
    select(Year, Month, Day, tod_slot, ta, meta_ta)
  
  ta <- tod_long %>%
    left_join(ta_long, by = c("Year","Month","Day","tod_slot")) %>%
    transmute(
      Year, Month, Day, Hour, Minute,
      Value = ta,
      Meta = ifelse(is.na(ta), meta_time, paste(meta_time, meta_ta, sep = " | "))
    )
  
  
  # Barometer: V14..V22 = (c1,b1,k1,c2,b2,k2,c3,b3,k3)
  p_long <- df %>%
    transmute(Year, Month, Day,
              c1 = suppressWarnings(as.integer(V14)), b1 = suppressWarnings(as.integer(V15)), k1 = suppressWarnings(as.integer(V16)),
              c2 = suppressWarnings(as.integer(V17)), b2 = suppressWarnings(as.integer(V18)), k2 = suppressWarnings(as.integer(V19)),
              c3 = suppressWarnings(as.integer(V20)), b3 = suppressWarnings(as.integer(V21)), k3 = suppressWarnings(as.integer(V22))) %>%
    pivot_longer(cols = c(c1,b1,k1,c2,b2,k2,c3,b3,k3),
                 names_to = c("part","idx"),
                 names_pattern = "([cbk])([123])",
                 values_to = "val") %>%
    pivot_wider(names_from = part, values_from = val) %>%
    mutate(
      tod_slot = paste0("tod", idx),
      mmHg = (c + b/12 + k/120) * 27.07,
      meta_orig_p = paste0("orig=", c, ".", b, ".", k)
    ) %>%
    select(Year, Month, Day, tod_slot, mmHg, meta_orig_p)
  
  # Thermometer at barometer: V23..V28 = (R1,d1,R2,d2,R3,d3)
  atb_long <- df %>%
    transmute(Year, Month, Day,
              R1 = suppressWarnings(as.integer(V23)), d1 = suppressWarnings(as.integer(V24)),
              R2 = suppressWarnings(as.integer(V25)), d2 = suppressWarnings(as.integer(V26)),
              R3 = suppressWarnings(as.integer(V27)), d3 = suppressWarnings(as.integer(V28))) %>%
    pivot_longer(cols = c(R1,d1,R2,d2,R3,d3),
                 names_to = c("part","idx"),
                 names_pattern = "([Rd])([123])",
                 values_to = "val") %>%
    pivot_wider(names_from = part, values_from = val) %>%
    mutate(
      tod_slot = paste0("tod", idx),
      atb_C = (R + d/10) * 1.25,
      meta_atb = paste0("atb=",R,".",d,"R" )
    ) %>%
    select(Year, Month, Day, tod_slot, atb_C, meta_atb)
  
  p <- tod_long %>%
    left_join(p_long,  by = c("Year","Month","Day","tod_slot")) %>%
    left_join(atb_long, by = c("Year","Month","Day","tod_slot")) %>%
    mutate(
      Value = ifelse(is.na(mmHg) | is.na(atb_C), NA_real_, round(convert_pressure(mmHg, lat=lat, alt=alt, atb = atb_C),2)),
      Meta  = ifelse(is.na(Value), meta_time, paste(meta_time, meta_orig_p, meta_atb, sep = " | "))
    ) %>%
    select(Year, Month, Day, Hour, Minute, Value, Meta)
  
  list(ta=ta, p=p)
}


# ---- main: read all sheets ----
sheets <- excel_sheets(infile)
sheets <- sheets[3:length(sheets)-1]


df_list <- map(sheets, ~read_year_sheet(.x))

df_ta <- bind_rows(map(df_list,"ta")) %>%
  arrange(Year, Month, Day, Hour, Minute)

df_p  <- bind_rows(map(df_list,"p")) %>%
  arrange(Year, Month, Day, Hour, Minute)

dfs <- list(ta=df_ta, p=df_p)

indir <- '/scratch3/PALAEO-RA/daily_data/final/'

### quick one
dirname <- "Zitenice/"

for (var in names(dfs)) {
  df <- dfs[[var]]
  filename <- outfile.name(name, var, df, TRUE)
  write_sef_f(
    as.data.frame(df[,c("Year","Month", "Day", "Hour", "Minute", "Value"),]),
    outfile = filename,
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
    meta    = df$Meta,
    time_offset = time.offset(lon),
    keep_na = FALSE
  )
  
  # do the QC's already
  

  filename_full <- paste0(filename, ".tsv")
  qc(glue(indir,dirname,filename_full), outpath=glue(indir, dirname))
  
  qcfilename <- paste0("qc_Zitenice_",var,"_subdaily.txt")
  write_flags_f(infile=glue(indir,dirname,filename_full),
                qcfile=glue(indir,dirname, qcfilename),
                outpath=glue(indir,dirname),
                match=FALSE)
}
  



