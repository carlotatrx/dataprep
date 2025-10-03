rm(list = ls())
library(lubridate)
library(dplyr)
library(readr)
library(tidyr)
library(dataresqc)
source("/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R")

files <- list.files("/scratch3/PALAEO-RA/daily_data/original/Extremadura",
                    pattern = ".tab$", full.names = TRUE)

link <- "https://doi.pangaea.de/10.1594/PANGAEA.925805"
source <- "Vaquero, José Manuel; Bravo-Paredes, Nieves; Obregón, María Angeles; Sánchez-Carrasco, Víctor Manuel; Valente, Maria Antónia; Trigo, Ricardo M; Domínguez-Castro, Fernando; Montero-Martín, Javier; Vaquero-Martínez, Javier; Antón, Manuel; García, José Agustín; Gallego, María Cruz (2020): Early meteorological records from BAALAN1 (1866). PANGAEA."

time.offset <- function(lon) {as.numeric(lon)*12/180}

outfile.name <- function(name, var, df, subdaily=T) {
  subdaily.str <- ifelse(subdaily, "subdaily", "daily")
  paste0(name,"_",get_date_range(df), "_", var,"_", subdaily.str)
}

units <- function(var) {
  case_when(
    var == "ta" ~ "C",
    var == "p"  ~ "hPa",
    var == "dd" ~ "deg",
    var == "rr" ~ "mm",
    T ~ "unknown"
  )
}

# wind conversion
directions <- c("N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
                "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW")

dd_normalize <- function(x) {
  y <- toupper(trimws(x))
  y <- recode(y,
              "NWE"  = "WNW",
              "NWW"  = "WNW",
              "SWW"  = "WSW",
              "WENW" = "WNW",
              "ENW"  = "ENE",
              "NEE"  = "ENE",
              "NN0"  = "NNW",   # OCR '0' -> 'W'
              "NNO"  = "NNW",   # 'O' (Oeste) -> W
              "SWS"  = "SSW",
              "ES"   = "SE",
              .default = y,
              .missing = NA_character_
              )
  ifelse(y %in% directions, y, NA_character_)
}
dd2deg <- function(x) {
  ifelse(is.na(x), NA_real_, 22.5 * (match(toupper(x), directions) - 1))
}
map2dd <- function(x) {
  x2 <- x %>% trimws()
  x2[x2==""] <- NA_character_
  x2.low <- tolower(x2)  # per si de cas
  out <- x2
  out[!is.na(x2.low) & x2.low=="gallego"] <- "NW"   # gallego: Dicho del viento: Procedente del noroeste, de la parte de Galicia.(RAE)
  out[!is.na(x2.low) & x2.low=="solano"]  <- "W"    # solano: Viento que sopla de donde nace el sol.(RAE)
  toupper(out)
}



# BAZAFR1 -----------------------------------------------------------------
code <- "BAZAFR1"
lat  <- 38.433333
lon  <- -6.429167
name <- "Extremadura-Zafra"
f <- files[grepl(code, files)] # select file from list

raw <- suppressWarnings(read_tsv(
  file      = f,
  skip      = 13,
  col_types = cols(.default = col_guess()),
  trim_ws   = T,
  guess_max = 100000
)) %>%
  rename(
    Date = 1,
    ta1 = 2,
    ta2 = 3,
    ta3 = 4,
  ) %>% mutate(
    Year   = year(Date),
    Month  = month(Date),
    Day    = day(Date)
  )
head(raw)

df.ta <- raw %>%
  pivot_longer(
    cols = c(ta1, ta2, ta3),
    names_to = "tod",
    names_prefix = "ta",
    values_to = "ta.orig"
  ) %>% mutate(
    # assign seasonal times:
    # Oct–Mar: 07,12,17 ; Apr–Sep: 06,12,18
    Hour = case_when(
      Month %in% c(10,11,12,1,2,3) & tod=="1" ~ 7L,
      Month %in% c(10,11,12,1,2,3) & tod=="2" ~ 12L,
      Month %in% c(10,11,12,1,2,3) & tod=="3" ~ 17L,
      Month %in% c(4,5,6,7,8,9)    & tod=="1" ~ 6L,
      Month %in% c(4,5,6,7,8,9)    & tod=="2" ~ 12L,
      Month %in% c(4,5,6,7,8,9)    & tod=="3" ~ 18L
    ),
    Minute = 0L,
    meta   = paste0("obs.#=", tod, " | time=inferred")
  ) %>% select(Year, Month, Day, Hour, Minute, Value=ta.orig, meta)

head(df.ta)

var <- "ta"
write_sef_f(
  as.data.frame(df.ta),
  outfile = outfile.name(name, var, df.ta, TRUE),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = NA,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  period  = "0",
  stat    = "point",
  units   = units(var),
  meta    = df.ta$meta,
  keep_na = FALSE
)

# BAVVEN1 -----------------------------------------------------------------
code <- "BAVVEN1"
lat  <- 38.265278
lon  <- -6.474444
name <- "Extremadura-Valencia_del_Ventoso"

f <- files[grepl(code, files)] # select file from list

raw <- suppressWarnings(read_tsv(
  file      = f,
  skip      = 25,
  col_types = cols(.default = col_guess()),
  trim_ws   = T,
  guess_max = 100000
)) %>%
  rename(
    Date        = 1,
    p1      = 2,
    p2 = 3,
    tmax = 6,
    tmin = 7,
    tmean = 8,
    dd1 = 10,
    dd2 = 11,
    w1 = 12,
    w2 = 13,
    prcp = 16
  ) %>% mutate(
    Year   = as.integer(year(Date)),
    Month  = as.integer(month(Date)),
    Day    = day(Date)
  )


df.p <- raw %>%
  pivot_longer(
    cols = c("p1","p2"),
    names_to = "tod",
    values_to = "p.orig"
  ) %>% mutate(
    Hour = if_else(tod=="p1", 9L, 15L),
    Minute= 0L,
    meta = paste0("orig.time=", sprintf("%02d", Hour), ":", sprintf("%02d", Minute))
  ) %>% select(Year, Month, Day, Hour, Minute, Value=p.orig, meta)

df.dd <- raw %>%
  pivot_longer(
    cols = c(dd1, dd2, w1, w2),
    names_to = c(".value", "tod"),
    names_pattern = "([a-z]+)([12])"
  ) %>%
  mutate(
    dd = dd2deg(dd),
    Hour = if_else(tod==1L, 9L, 15L),
    Minute = 0L,
    meta = paste0(
      "orig.dd=", dd,
      ifelse(!is.na(w) & nzchar(trimws(as.character(w))),
             paste0(" | force=", w), ""), 
      ifelse(Hour==9L, " | orig.time=morning", " | orig.time=afternoon")
    )
  ) %>% select(Year, Month, Day, Hour, Minute, Value=dd, meta)


# save dd & p
for (var in c("dd", "p")){
  dat <- if (var=="dd") df.dd else df.p
  write_sef_f(
    as.data.frame(dat),
    outfile = outfile.name(name, var, dat, TRUE),
    outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
    cod     = code,
    lat     = lat,
    lon     = lon,
    alt     = NA,
    sou     = source,
    link    = link,
    nam     = name,
    var     = var,
    stat    = "point",
    units   = ifelse(var=="dd", units(var), "arbitrary units"),
    meta    = dat$meta,
    time_offset = ifelse(var=="dd", 0, time.offset(lon)),
    keep_na = FALSE
  )
}

# save tmin & tmax
df.ta <- raw %>%
  transmute(
    Year, Month, Day,
    Hour = 24L,
    Minute = 0L,
    tmin, tmax
  )
base_cols <- df.ta[c("Year","Month","Day","Hour","Minute")]

for (var in c("tmin", "tmax")) {
  dat   <- cbind(base_cols, setNames(df.ta[var], var))
  write_sef_f(
    dat,
    outfile = outfile.name(name, var, dat, FALSE),
    outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
    cod     = code,
    lat     = lat,
    lon     = lon,
    alt     = NA,
    sou     = source,
    link    = link,
    nam     = name,
    var     = var,
    period  = "day",
    stat    = ifelse(var=="tmax", "max", "min"),
    units   = "C",
    keep_na = FALSE
  )
}

# save rr
df.rr <- raw %>%
  transmute(
    Year, Month, Day,
    Hour = 24L,
    Minute = 0L,
    prcp
  )
var <- "rr"
write_sef_f(
  as.data.frame(df.rr),
  outfile = outfile.name(name, var, df.rr, FALSE),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = NA,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  period  = "day",
  stat    = "sum",
  units   = units(var),
  keep_na = FALSE
)


# BAVALD1 -----------------------------------------------------------------
code <- "BAVALD1"
lat  <- 37.871581
lon  <- -6.522527
name <- "Extremadura-Valdesevilla"

f <- files[grepl(code, files)] # select file from list

raw <- read.delim(
  file          = f,
  skip          = 15,
  sep           = "\t",
  fileEncoding  = "UTF-8",
  check.names   = FALSE,    # keep original column names
  stringsAsFactors = FALSE,
  quote         = ""        # columns contain parentheses/commas but no quotes
) %>%
  rename(
    Date        = 1,
    p.orig      = 2,
    ta          = 3,
    dd.orig     = 4,
    w.orig      = 5
  ) %>% mutate(
    Year   = as.integer(year(Date)),
    Month  = as.integer(month(Date)),
    Day    = day(Date),
    Hour   = NA_integer_,
    Minute = NA_integer_
  )

df <- raw %>%
  mutate(
    p = round(convert_pressure(as.numeric(p.orig), lat=lat, atb=as.numeric(ta)),2),
    dd = dd2deg(dd.orig),
    meta.p = paste0("orig.p=", p.orig, "mmHg | atb=", ta),
    meta.dd = paste0("orig.dd=", dd.orig, " | force=", w.orig)
  )

meta_map  <- list(ta = "", p = df$meta.p, dd = df$meta.dd)
base_cols <- df[c("Year","Month","Day","Hour","Minute")]

for (var in c("ta", "p", "dd")) {
  meta  <- meta_map[[var]]
  dat   <- cbind(base_cols, setNames(df[var], var))
  write_sef_f(
    dat,
    outfile = outfile.name(name, var, dat, FALSE),
    outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
    cod     = code,
    lat     = lat,
    lon     = lon,
    alt     = NA,
    sou     = source,
    link    = link,
    nam     = name,
    var     = var,
    stat    = "point",
    units   = units(var),
    metaHead = "orig.coords=38.700000,-6.829167",
    meta    = meta,
    keep_na = FALSE
  )
}

# BACVAC1 -----------------------------------------------------------------
code <- "BACVAC1"
lat  <- 38.086667
lon  <- -6.421111
name <- "Extremadura-Cabeza_la_Vaca"

f <- files[grepl(code, files)] # select file from list

raw <- suppressWarnings(read_tsv(
  file      = f,
  skip      = 23,
  col_types = cols(.default = col_guess()),
  trim_ws   = T,
  guess_max = 100000
)) %>%
  rename(
    Date         = 1,
    p1.orig      = 2,
    p2.orig      = 3,
    tmax1        = 4,
    tmax2        = 5,
    tmin         = 6,
    ddorig1      = 7,
    ddorig2      = 8,
    worig1       = 9,
    worig2       = 10,
    prcp         = 13,
    evaporation  = 14
  ) %>% mutate(
    Year   = as.integer(year(Date)),
    Month  = as.integer(month(Date)),
    Day    = day(Date)
  )

df.p <- raw %>%
  pivot_longer(
    cols = c(p1.orig, p2.orig),
    names_to = "tod",
    values_to = "p.orig"
  ) %>% 
  mutate(
    Hour = case_when(
      tod=="p1.orig" ~ 9L,
      tod=="p2.orig" ~ 15L
    ),
    Minute = 0L,
    
    p = round(convert_pressure(p.orig, lat = lat,),2),
    
    meta   = paste0("orig.time=", sprintf("%02d", Hour), ":", sprintf("%02d", Minute), " | orig.p=", p.orig, "mmHg")
  ) %>% select(Year, Month, Day, Hour, Minute, Value=p, meta)

head(df.p)


df.ta <- raw %>%
  mutate(
    tmax2_is = !is.na(tmax2),
    tmax1_is = !is.na(tmax1),
    tmax2_C  = ifelse(tmax2_is, tmax2 * 1.25, NA_real_),  # convert R to C for ta other
    
    tmax = coalesce(tmax1, tmax2_C),  # prefer the first one TxTxTx
    
    meta.max = case_when(
      tmax1_is & tmax2_is ~ paste0("alternative.orig.ta=", tmax2, "R"),
      !tmax1_is & tmax2_is ~ paste0("orig.ta=", tmax2, "R"),
      TRUE ~ ""
    ),
    
    Hour   = 24L,
    Minute = 0L,
    
    # qc: flaig if Tmax <= Tmin
    flag = !is.na(tmax2) & !is.na(tmin) & (tmax <= tmin),
    meta.max = ifelse(flag, paste0(meta.max, " | qc= Tmax<=Tmin"), meta.max),
    meta.min = paste0("qc= Tmax<=Tmin")
  ) %>% select(Year, Month, Day, Hour, Minute, tmax, tmin, meta.max, meta.min)

head(df.ta)

df.dd <- raw %>%
  pivot_longer(
    cols = c(ddorig1, ddorig2, worig1, worig2),
    names_to = c(".value", "tod"),
    names_pattern = "([a-z]+)([12])"
  ) %>%
  mutate(
    dd = dd2deg(ddorig),
    Hour = if_else(tod==1L, 9L, 15L),
    Minute = 0L,
    meta = paste0(
      "orig.dd=", ddorig,
      ifelse(!is.na(worig) & nzchar(trimws(as.character(worig))),
             paste0(" | force=", worig), ""), 
      ifelse(Hour==9L, " | orig.time=morning", " | orig.time=afternoon")
    )
  ) %>% select(Year, Month, Day, Hour, Minute, Value=dd, meta)

head(df.dd)

df.rr <- raw %>%
  mutate(
    Hour = 24L,
    Minute = 0L
  ) %>% select(Year, Month, Day, Hour, Minute, Value=prcp)

# save dd
var <- "dd"
write_sef_f(
  as.data.frame(df.dd[, c("Year", "Month", "Day", "Hour", "Minute", "Value")]),
  outfile = outfile.name(name, var, df.dd, TRUE),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = NA,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "point",
  units   = units(var),
  meta    = df.dd$meta,
  time_offset = time.offset(lon),
  keep_na = FALSE
)

# save rr
var <- 'rr'
write_sef_f(
  as.data.frame(df.rr[, c("Year", "Month", "Day", "Hour", "Minute", "Value")]),
  outfile = outfile.name(name, var, df.rr, FALSE),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = NA,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "sum",
  period  = "day",
  units   = units(var),
  keep_na = TRUE
)

# save p
var <- "p"
write_sef_f(
  as.data.frame(df.p[, c("Year", "Month", "Day", "Hour", "Minute", "Value")]),
  outfile = outfile.name(name, var, df.p, TRUE),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = NA,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "point",
  units   = "hPa",
  meta    = df.p$meta,
  time_offset = time.offset(lon),
  metaHead = "PTC=N | PGC=N",
  keep_na = FALSE
)

# save tmax and tmin
var <- "tmax"
write_sef_f(
  as.data.frame(df.ta[, c("Year", "Month", "Day", "Hour", "Minute", "tmax")]),
  outfile = outfile.name(name, var, df.ta, FALSE),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = NA,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "max",
  period  = "day",
  units   = "C",
  meta    = df.ta$meta.max,
  keep_na = FALSE
)

var <- "tmin"
write_sef_f(
  as.data.frame(df.ta[, c("Year", "Month", "Day", "Hour", "Minute", "tmin")]),
  outfile = outfile.name(name, var, df.ta, FALSE),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = NA,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "min",
  period  = "day",
  units   = "C",
  meta    = df.ta$meta.min,
  keep_na = FALSE
)

# BABADA2 -----------------------------------------------------------------
code <- "BABADA2"
lat  <- 38.900000
lon  <- -6.962500
name <- "Extremadura-Badajoz"

f <- files[grepl(code, files)] # select file from list

raw <- suppressWarnings(read_tsv(
  file      = f,
  skip      = 15,
  col_types = cols(.default = col_guess()),
  trim_ws   = T,
  guess_max = 100000
)) %>%
  rename(
    Date         = 1,
    p.orig       = 2,
    ta.orig      = 3,
    dd.orig      = 4,
    dd.descr     = 5
  )

df <- raw %>%
  mutate(
    Year   = year(Date),
    Month  = month(Date),
    Day    = day(Date),
    Hour   = 24L,
    Minute = 0L,
    p.orig  = suppressWarnings(as.numeric(p.orig)),
    ta.orig = suppressWarnings(as.numeric(ta.orig)),
    p  = round(convert_pressure(p.orig, lat=lat, atb=ta.orig),2),
    dd = dd2deg(dd.orig),
    meta.p  = paste0("orig.p=", p.orig, "mmHg | atb=", ta.orig, "C"),
    meta.dd = paste0("orig.dd=", dd.orig, " | force=", dd.descr)
  ) %>% select(Year, Month, Day, Hour, Minute, p, dd, ta.orig, meta.p, meta.dd)

head(df)

value_map <- c(ta = "ta.orig", p = "p", dd = "dd")
meta_map  <- list(ta = "", p = df$meta.p, dd = df$meta.dd)
base_cols <- df[c("Year","Month","Day","Hour","Minute")]


for (var in names(value_map)) {
  value <- value_map[[var]]
  meta  <- meta_map[[var]]
  dat   <- cbind(base_cols, setNames(df[value], value))
  write_sef_f(
    dat,
    outfile = outfile.name(name, var, df, F),
    outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
    cod     = code,
    lat     = lat,
    lon     = lon,
    alt     = NA,
    sou     = source,
    link    = link,
    nam     = name,
    var     = var,
    stat    = "mean",
    period  = "day",
    units   = units(var),
    metaHead = ifelse(var=="p", "PTC=Y | PGC=N", ""),
    meta    = meta,
    keep_na = F
  )
}

# BABADA6 -----------------------------------------------------------------
code <- "BABADA6"
f <- files[grepl(code, files)] # select file from list

raw <- read.delim(
  file          = f,
  skip          = 25,
  sep           = "\t",
  fileEncoding  = "UTF-8",
  check.names   = FALSE,    # keep original column names
  stringsAsFactors = FALSE,
  quote         = ""        # columns contain parentheses/commas but no quotes
) %>%
  rename(
    Date = 1,
    ta1 = 2,
    ta2 = 3,
    ta3 = 4,
    z1 = 5,
    l1 = 6,
    z2 = 7,
    l2 = 8,
    z3 = 9,
    l3 = 10,
    dd.orig.1 = 11,
    dd.orig.2 = 13,
    dd.orig.3 = 15
  ) %>% mutate(
    Date   = ym(Date), # careful here because Date has format "YYYY-MM"
    Year   = year(Date),
    Month  = month(Date)
  ) %>%
  group_by(Year, Month) %>%
  mutate(Day = row_number()) %>%
  ungroup

head(raw)

df.ta <- raw %>%
  pivot_longer(
    cols = c(ta1, ta2, ta3),
    names_to = "tod",
    names_prefix = "ta",
    values_to = "ta.orig"
  ) %>% mutate(
    # assign seasonal times:
    # Oct–Mar: 07,12,17 ; Apr–Sep: 06,12,18
    Hour = case_when(
      Month %in% c(10,11,12,1,2,3) & tod=="1" ~ 7L,
      Month %in% c(10,11,12,1,2,3) & tod=="2" ~ 12L,
      Month %in% c(10,11,12,1,2,3) & tod=="3" ~ 17L,
      Month %in% c(4,5,6,7,8,9)    & tod=="1" ~ 6L,
      Month %in% c(4,5,6,7,8,9)    & tod=="2" ~ 12L,
      Month %in% c(4,5,6,7,8,9)    & tod=="3" ~ 18L
    ),
    Minute = 0L,
    meta   = paste0("obs.#=", tod, " | time=inferred | Date=inferred")
  ) %>% select(Year, Month, Day, Hour, Minute, Value=ta.orig, meta)

head(df.ta)


df.p <- raw %>%
  pivot_longer(
    cols = c(z1:l3, ta1:ta3),
    names_to = c(".value", "tod"),
    names_pattern = "([a-z]+)([123])"
  ) %>% 
  mutate(
    p.orig = z + l/12,
    
    # use English inch 1in=25.4mHg
    p = round(convert_pressure(p.orig, f = 25.04, lat = lat, atb = ta),2),

    # assign seasonal times:
    # Oct–Mar: 07,12,17 ; Apr–Sep: 06,12,18
    tod = parse_number(tod),
    Hour = case_when(
      Month %in% c(10,11,12,1,2,3) & tod=="1" ~ 7L,
      Month %in% c(10,11,12,1,2,3) & tod=="2" ~ 12L,
      Month %in% c(10,11,12,1,2,3) & tod=="3" ~ 17L,
      Month %in% c(4,5,6,7,8,9)    & tod=="1" ~ 6L,
      Month %in% c(4,5,6,7,8,9)    & tod=="2" ~ 12L,
      Month %in% c(4,5,6,7,8,9)    & tod=="3" ~ 18L
    ),
    Minute = 0L,
    meta   = paste0("obs.#=", tod, " | time=inferred | Date=inferred | orig.p=",
                    z, "in", l, "l | atb=", ta)
  ) %>% select(Year, Month, Day, Hour, Minute, Value=p, meta, p.orig)

head(df.p)


df.dd <- raw %>%
  pivot_longer(
    cols = c(dd.orig.1, dd.orig.2, dd.orig.3),
    names_to = "tod",
    names_prefix = "dd\\.orig\\.",
    values_to = "dd_txt"
  ) %>%
  mutate(
    tod = as.integer(tod),
    dd_txt = as.character(dd_txt),
    dd_norm = dd_normalize(dd_txt),
    dd = dd2deg(dd_norm),
    Hour = case_when(
      Month %in% c(10,11,12,1,2,3) & tod=="1" ~ 7L,
      Month %in% c(10,11,12,1,2,3) & tod=="2" ~ 12L,
      Month %in% c(10,11,12,1,2,3) & tod=="3" ~ 17L,
      Month %in% c(4,5,6,7,8,9)    & tod=="1" ~ 6L,
      Month %in% c(4,5,6,7,8,9)    & tod=="2" ~ 12L,
      Month %in% c(4,5,6,7,8,9)    & tod=="3" ~ 18L
    ),
    Minute = 0L,
    meta   = paste0("obs.#=", tod, " | time=inferred | Date=inferred | orig.dd=",dd_txt)
    
  ) %>% select(Year, Month, Day, Hour, Minute, Value=dd, meta)

head(df.dd)


dfs <- list(
  ta = df.ta,
  p  = df.p,
  dd = df.dd
)

for (var in names(dfs)) {
  dat   <- select(dfs[[var]],  Year, Month, Day, Hour, Minute, Value)
  write_sef_f(
    as.data.frame(dat),
    outfile = outfile.name(name, var, dfs[[var]], TRUE),
    outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
    cod     = code,
    lat     = lat,
    lon     = lon,
    alt     = NA,
    sou     = source,
    link    = link,
    nam     = name,
    var     = var,
    stat    = "point",
    period  = 0,
    units   = units(var),
    metaHead = ifelse(var=="p", "PTC=Y | PGC=N | assumed English inch (25.4mmHg)", ""),
    meta    = dfs[[var]]$meta,
    keep_na = FALSE
  )
}

# BAALAN1 & 2 -----------------------------------------------------------------
code <- "BAALAN2" # "BAALAN1" or "BALAAN2"
lat  <- 38.784444
lon  <- -6.243611
name <- "Extremadura-Alange"

f <- files[grepl(code, files)] # select file from list

raw <- suppressWarnings(read_tsv(
  file      = f,
  skip      = 13,
  col_types = cols(.default = col_guess()),
  trim_ws   = T,
  guess_max = 100000
  )
) %>%
  rename(
    Date         = 1,
    p.orig       = 2,
    ta.orig      = 3,
    rh           = 4
  )

df <- raw %>%
  mutate(
    Year   = year(Date),
    Month  = month(Date),
    Day    = day(Date),
    Hour   = hour(Date),
    Minute = minute(Date),
    p      = round(convert_pressure(p.orig, lat=lat, atb=ta.orig),2),
    meta.time = paste0("orig.time=", sprintf("%02d", Hour),":", sprintf("%02d",Minute)),
    meta.p = paste0("orig.p=",p.orig,"mmHg | ", meta.time)
  )
head(df)

# save_ta
var <- "ta"
write_sef_f(
  as.data.frame(df[, c("Year", "Month", "Day", "Hour", "Minute", "ta.orig")]),
  outfile = outfile.name(name, var, df, T),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = NA,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "point",
  units   = units(var),
  meta    = df$meta.time,
  metaHead = "Instrument=Thermometer",
  time_offset = time.offset(lon),
  keep_na = F
)

# save p
var <- "p"
write_sef_f(
  as.data.frame(df[, c("Year", "Month", "Day", "Hour", "Minute", "p")]),
  outfile = outfile.name(name, var, df, T),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = NA,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "point",
  units   = units(var),
  meta    = df$meta.p,
  metaHead = "PGC=N | PTC=Y",
  time_offset = time.offset(lon),
  keep_na = F
)

# save rh
var <- "rh"
units.rh <- ifelse(code=="BAALAN1", "relative units", "%")
write_sef_f(
  as.data.frame(df[, c("Year", "Month", "Day", "Hour", "Minute", "rh")]),
  outfile = outfile.name(name, var, df, T),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = NA,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "point",
  units   = units.rh,
  meta    = df$meta.time,
  metaHead = "Instrument=Saussure hygrometer",
  time_offset = time.offset(lon),
  keep_na = F
)


# BAALAN3 -----------------------------------------------------------------
code <- "BAALAN3"

f <- files[grepl(code, files)] # select file from list

raw <- suppressWarnings(read_tsv(
  file      = f,
  skip      = 13,
  col_types = cols(.default = col_guess()),
  trim_ws   = T,
  guess_max = 100000
)) %>%
  rename(
    Date         = 1,
    ta.orig      = 2,
    p.orig       = 3
  )


head(raw)
df <- raw %>%
  mutate(
    Year   = year(Date),
    Month  = month(Date),
    Day    = day(Date),
    Hour   = NA,
    Minute = NA,
    ta     = ta.orig*1.25,
    p.num  = suppressWarnings(as.numeric(p.orig)),
    p      = round(convert_pressure(p.num, lat=lat, atb=ta),2),
    meta.p = paste0("orig.p=",p.orig,"mmHg"),
    meta.ta = paste0("orig.ta=", ta.orig,"R")
  )
head(df)

# save_ta
var <- "ta"
write_sef_f(
  as.data.frame(df[, c("Year", "Month", "Day", "Hour", "Minute", "ta")]),
  outfile = outfile.name(name, var, df, F),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = NA,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "point",
  units   = units(var),
  meta    = df$meta.ta,
  metaHead = "Instrument=Reaumur thermometer",
  time_offset = time.offset(lon),
  keep_na = F
)

# save_p
var <- "p"
write_sef_f(
  as.data.frame(df[, c("Year", "Month", "Day", "Hour", "Minute", "p")]),
  outfile = outfile.name(name, var, df, F),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = NA,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "point",
  units   = units(var),
  meta    = df$meta.p,
  metaHead = "PTC=Y | PGC=N",
  time_offset = time.offset(lon),
  keep_na = F
)


# BAALAN4 -----------------------------------------------------------------
code <- "BAALAN4"

f <- files[grepl(code, files)] # select file from list

raw <- read.delim(
  file          = f,
  skip          = 20,
  sep           = "\t",
  fileEncoding  = "UTF-8",
  check.names   = FALSE,    # keep original column names
  stringsAsFactors = FALSE,
  quote         = ""        # columns contain parentheses/commas but no quotes
  ) %>%
  rename(
    Date         = 1,
    ta.09        = 2,
    ta.14        = 3,
    dd.descr.09  = 4,
    dd.descr.14  = 5,
    dd.descr.mean = 6,
    dd.09        = 7,
    dd.14        = 8,
    dd.mean      = 9,
    p.orig.09    = 10,
    p.orig.14    = 11
  ) %>%
  mutate(
    Year   = year(Date),
    Month  = month(Date),
    Day    = day(Date)
  )

head(raw)

############ wind #############

df.dd <- raw %>%
  pivot_longer(
    cols = c(dd.descr.09, dd.descr.14, dd.descr.mean),
    names_to = "tod",
    names_prefix = "dd\\.descr\\.",
    values_to = "dd.orig"
  ) %>% mutate(
    Hour = case_when(
      tod=="09" ~ 9L,
      tod=="14" ~ 14L,
      tod=="mean" ~ 24L
    ),
    Minute = 0,
    dd.clean = map2dd(dd.orig),
    Value  = dd2deg(dd.clean),
    meta   = ifelse(is.na(Value),"", paste0("orig.dd=", dd.orig))
  ) %>% select(Year, Month, Day, Hour, Minute, Value, meta, dd.orig)

head(df.dd)

df.dd.daily <- df.dd %>%
  filter(Hour==24L, !is.na(Value)) 

df.dd.subdaily <- df.dd %>%
  filter(Hour != 24L, !is.na(Value)) %>%
  mutate(
    meta = paste0(meta, " | orig.time", sprintf("%02d", Hour),":", sprintf("%02d",Minute))
  )

# save_dd_daily
var <- "dd"
write_sef_f(
  as.data.frame(df.dd.daily[, c("Year", "Month", "Day", "Hour", "Minute", "Value")]),
  outfile = outfile.name(name, var, df.dd.daily, F),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = NA,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "mean",
  period  = "day",
  units   = "deg",
  meta    = df.dd.daily$meta,
  keep_na = F
)

write_sef_f(
  as.data.frame(df.dd.subdaily[, c("Year", "Month", "Day", "Hour", "Minute", "Value")]),
  outfile = outfile.name(name, var, df.dd.subdaily, T),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = NA,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "point",
  units   = "deg",
  meta    = df.dd.subdaily$meta,
  time_offset = time.offset(lon),
  keep_na = F
)

############# ta #######################

df.ta <- raw %>%
  pivot_longer(
    cols = c(ta.09, ta.14),
    names_to = "tod",
    names_prefix = "ta\\.",
    values_to = "ta.orig"
  ) %>% mutate(
    Hour = case_when(
      tod=="09" ~ 9L,
      tod=="14" ~ 14L
    ),
    Minute = 0,
    Value  = ta.orig*1.25,
    meta   = paste0("orig.ta=", ta.orig,"R | orig.time=", sprintf("%02d", Hour),":", sprintf("%02d",Minute))
  ) %>% select(Year, Month, Day, Hour, Minute, Value, meta)
head(df.ta)

var <- "ta"
write_sef_f(
  as.data.frame(df.ta[, c("Year", "Month", "Day", "Hour", "Minute", "Value")]),
  outfile = outfile.name(name, var, df.ta, T),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = NA,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "point",
  units   = "C",
  meta    = df.ta$meta,
  time_offset = time.offset(lon),
  metaHead = "Instrument=Reaumur thermometer",
  keep_na = F
)


############# p #######################

df.p <- raw %>%
  pivot_longer(
    cols = c(p.orig.09, p.orig.14),
    names_to = "tod",
    names_prefix = "p\\.orig\\.",
    values_to = "p.orig"
  ) %>% mutate(
    Hour = case_when(
      tod=="09" ~ 9L,
      tod=="14" ~ 14L
    ),
    Minute = 0,
    Value  = round(convert_pressure(p.orig, lat=lat, atb=df.ta$Value),2),
    meta   = paste0("orig.p=", p.orig,"mmHg | orig.time=", sprintf("%02d", Hour),":", sprintf("%02d",Minute), " | atb=", df.ta$Value, "C")
  ) %>% select(Year, Month, Day, Hour, Minute, Value, p.orig, meta)
head(df.p)

var <- "p"
write_sef_f(
  as.data.frame(df.p[, c("Year", "Month", "Day", "Hour", "Minute", "Value")]),
  outfile = outfile.name(name, var, df.p, T),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Extremadura',
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = NA,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "point",
  units   = "hPa",
  meta    = df.p$meta,
  time_offset = time.offset(lon),
  metaHead = "PTC=Y | PGC=N",
  keep_na = F
)

