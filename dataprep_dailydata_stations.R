library(dataresqc)
library(readxl)
library(tidyr)
library(dplyr)
library(glue)

source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')

################################################################################
linein2mmHg <- function(inch, l) {
  val = inch + l/12
  return(round(val * 27.07, 2))
}

R2C <- function(val) {
  return(round(val * 1.25, 1))
}

C2R <- function(val) {
  return(round(val / 1.25), 2)
}

# Augsburg ----------------------------------------------------------------

filedir <- '/scratch3/PALAEO-RA/daily_data/original/Augsburg/'
files   <- list.files(filedir, pattern="^Europe_T2_DE_Augsburg_\\d{4}_subdaily.xls$", full.names = T)
augsburg_alt <- 500
augsburg_lon <- 10.8833
augsburg_lat <- 48.35

# for file in filedir:

augsburg <- function(filename) {
  df <- read_excel(filename, skip = 5)
  colnames(df) <- c("Year", "Month", "Day", "p_in7", "p_l7", "p_in14", "p_l14", "p_in21", "p_l21", "tbar_R7",
                    "tbar_R14", "tbar_R21", "ta_R7", "ta_R14", "ta_R21", "notes")
  
  df$notes <- NULL # drop column
  
  # replace all commas with points
  df[] <- lapply(df, function (x) {
    if (is.character(x) || is.factor(x)){
      x <- as.character(x)
      x <- gsub(",", ".", x)
      x    
    } else {
      x
    }
  })
  
  # replace all "-" with NA
  df[] <- lapply(df, function(x){
    if (is.character(x)) {
      x[x == "-"] <- NA
    }
    x
  })
  
  # coerce all cols to num
  df[] <- lapply(df, function (x) {
    if (is.character(x)) suppressWarnings(as.numeric(x)) else x
  })
  
  # unit conversions
  df.all <- df %>%
    mutate(
      p7  = linein2mmHg(p_in7, p_l7),
      p14 = linein2mmHg(p_in14, p_l14),
      p21 = linein2mmHg(p_in21, p_l21),
      
      tbar7  = R2C(tbar_R7),
      tbar14 = R2C(tbar_R14),
      tbar21 = R2C(tbar_R21),
      
      ta7  = R2C(ta_R7),
      ta14 = R2C(ta_R14),
      ta21 = R2C(ta_R21)
    )
  
  # pivot ta to long table
  df.ta <- df.all %>%
    select(Year, Month, Day, ta7, ta14, ta21, ta_R7, ta_R14, ta_R21) %>%
    pivot_longer(
      cols = c(ta7, ta14, ta21, ta_R7, ta_R14, ta_R21), # only pivot ta vars
      names_to = c(".value", "Hour"),
      names_pattern = "(ta|ta_R)(\\d+)",
      values_to = "Value"
    ) %>%
    mutate(
      Hour = as.integer(Hour),
      Minute = 0,
      Meta = ifelse(!is.na(ta_R), paste0("orig_ta=",ta_R,"R"), NA_character_)
    ) %>%
    select(Year, Month, Day, Hour, Minute, ta, Meta) %>%
    rename(Value = ta)
  
  # Pivot p to long format
  df.p <- df.all %>%
    select(Year, Month, Day,
           p7, p14, p21,                # converted pressures (mmHg)
           p_in7, p_in14, p_in21,       # inches
           p_l7, p_l14, p_l21,          # lines
           tbar_R7, tbar_R14, tbar_R21  # barometer ta in R
    ) %>%
    pivot_longer(
      cols = c(p7, p14, p21,
               p_in7, p_in14, p_in21,
               p_l7, p_l14, p_l21,
               tbar_R7, tbar_R14, tbar_R21), # only cols to pivot
      names_to = c("var", "Hour"),
      names_pattern = "([a-zA-Z_]+)(\\d+)"
    ) %>%
    pivot_wider(
      names_from = var,
      values_from = value
    ) %>%
    mutate(
      Hour = as.integer(Hour),
      Minute = 0,
      value = round(convert_pressure(p, lat = augsburg_lat, alt = augsburg_alt, atb=tbar_R),2),
      Meta = ifelse(!is.na(p_in) & !is.na(p_l) & !is.na(tbar_R),
                    paste0("orig_p=", p_in, ".", p_l, "Pin | ",
                           "atb=", tbar_R, "R"),
                    NA_character_)
    ) %>%
    select(Year, Month, Day, Hour, Minute, value, Meta)
  
  return(list(ta=df.ta, p=df.p))
}

augsburg.all <- purrr:::map(files, augsburg)

df.ta.augsburg <- bind_rows(purrr:::map(augsburg.all, "ta"))
df.p.augsburg  <- bind_rows(purrr:::map(augsburg.all, "p"))

write_sef_f(
  Data=as.data.frame(df.ta.augsburg), outfile="Augsburg_ta_subdaily.tsv",
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod="Augsburg",
  variable="ta",
  nam="Augsburg",
  lat=augsburg_lat, lon=augsburg_lon, alt=augsburg_alt,
  meta=df.ta.augsburg$Meta, metaHead = "Observer=Stark",
  units="C", sou="PALAEO-RA", stat="point",keep_na = F
)

write_sef_f(
  Data=as.data.frame(df.p.augsburg), outfile="Augsburg_p_subdaily.tsv",
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod="Augsburg",
  variable="p",
  nam="Augsburg",
  lat=augsburg_lat, lon=augsburg_lon, alt=augsburg_alt,
  meta=df.p.augsburg$Meta, metaHead = "Observer=Stark | PGC=Y | PTC=Y",
  units="hPa", sou="PALAEO-RA", stat="point", keep_na = F
)

## check
for (var in c('ta','p')) {
  # Run check_sef and print results
  infile <- paste0("/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/Augsburg_",var,"_subdaily.tsv")
  check_result <- check_sef(infile)
  print(check_result)
  
  # Run qc and save results
  qc(infile, outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/sef_tests")
  
  # Plot decimals and save as PNG
  plot_decimals(infile,
                outfile = paste0("/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/sef_tests/Augsburg_",var,"_subdaily_decimals"))
  
  write_flags_f(infile=infile,
                qcfile=glue('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/sef_tests/qc_Augsburg_{var}_subdaily.txt'),
                outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
                match=F)
}


# Boehringen --------------------------------------------------------------

# qc fir wubd direction
df.dd <- read_sef('/scratch3/PALAEO-RA/DataRescue/Projects/Europe/5_QCed/Boehringen/PALAEO-RA_Europe_Boehringen_17630101-17820512_dd.tsv')
