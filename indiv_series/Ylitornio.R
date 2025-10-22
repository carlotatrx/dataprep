rm(list = ls())
library(lubridate)
library(dplyr)
library(readr)
library(tidyr)
library(dataresqc)
library(readxl)
library(XLConnect)

source("/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R")

file <- list.files("/scratch3/PALAEO-RA/daily_data/original/Ylitornio",
                   full.names = TRUE)

code <-	"Ylitornio"
name <-	"Ylitornio"
lat <- 66.319266
lon	<- 23.670970
alt	<- 55
source <-	"Helama, S., Holopainen, J., Timonen, M., Ogurtsov, M. G., Lindholm, M., Meriläinen, J., Eronen, M. (2004). Comparison of living-tree and subfossil ringwidths with summer temperatures"
link <-	"Holopainen, J. (2006): Reconstructions of past climates from documentary and natural sources in Finland since the 18th century. Ph.D. dissertation. Publications of the Department of Geology D9, Yliopistopaino, Helsinki, 33 pp."
metaHead <- "Observer=Johan Portin"

wind_dir_map <- c(
  "1" = "N/NNO",
  "2" = "NO/ONO",
  "3" = "O/OSO",
  "4" = "SO/SSO",
  "5" = "S/SSW",
  "6" = "SW/WSW",
  "7" = "W/WNW",
  "8" = "NW/NNW",
  "9" = NA_character_
)

# read file ---------------------------------------------------------------

raw <- readWorksheetFromFile(file, sheet=1, startRow=7, header=T)

df <- raw %>%
  pivot_longer(
    cols = matches("Observation\\.hour|Air\\.pressure|Temperature(\\.\\d*)?$|Wind\\.Direction|Wind\\.Speed|Cloudiness"),
    names_to = c(".value", "tod"),
    names_pattern = "(.*?)(?:\\.(\\d+))?$"
  ) %>%
  rename(
    Hour = "Observation.hour"
  ) %>%
  mutate(
    Minute = if_else(is.na(Hour), NA_integer_, 0L),
    
    # p = round(convert_pressure(Air.pressure, f=27.07, lat=lat, alt=alt, atb=Temperature),1),
    ddorig = recode(as.character(Wind.Direction), !!!wind_dir_map),
    # keep only the first wind direction
    ddnorm = dd_normalize(sub("/.*", "", ddorig)),
    dd     = dd2deg(ddnorm),
    
    meta.dd = paste0(meta_time(Hour, Minute), " | orig.dd=",
                     Wind.Direction, " (", ddorig, ") | force=", Wind.Speed),
    meta.p  = paste0(meta_time(Hour, Minute), " | orig.p=", Air.pressure, "in | atb=", Temperature, "C"),
    meta.ta = meta_time(Hour, Minute),
    
  ) %>% select(Year, Month, Day, Hour, Minute,
               rr = "Precipitation",
               dd, meta.dd,
               p  = "Air.pressure", meta.p,
               ta = "Temperature", meta.ta)

head(df)


# precipitation -----------------------------------------------------------

# precipitation is daily, it works different
df.rr <- df %>%
  group_by(Year, Month, Day) %>%
  mutate(
    Hour = 24L,
    Minute = 0L,
    rr = ifelse(row_number() == 1, rr, NA_real_)
  ) %>%
  ungroup() %>% select(Year, Month, Day, Hour, Minute, rr)



# save files --------------------------------------------------------------
var <- "rr"
write_sef_f(
  as.data.frame(df.rr),
  outfile = outfile.name(name, var, df.rr, FALSE),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Ylitornio',
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  period  = "day",
  stat    = "sum",
  units   = "unknown",
  keep_na = FALSE
)



base_cols <- df[c("Year","Month","Day","Hour","Minute")]
meta.list <- list(dd=df$meta.dd, p=df$meta.p, ta=df$meta.ta, rr="")

var <- "dd"
df.dd <- df[c("Year", "Month", "Day", "Hour", "Minute", "dd", "meta.dd")]
df.dd <- df.dd[df.dd$Year >= 1802, ]
write_sef_f(
  as.data.frame(df.dd),
  outfile = outfile.name(name, var, df.dd, TRUE),
  outpath = '/scratch3/PALAEO-RA/daily_data/final/Ylitornio',
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
  meta    = df.dd$meta.dd,
  metaHead = metaHead,
  time_offset = time.offset(lon),
  keep_na = TRUE
)

for (var in c("p", "ta")) {
  
  dat <- cbind(base_cols, setNames(df[var], var))
  write_sef_f(
    dat,
    outfile = outfile.name(name, var, dat, TRUE),
    outpath = '/scratch3/PALAEO-RA/daily_data/final/Ylitornio',
    cod     = code,
    lat     = lat,
    lon     = lon,
    alt     = alt,
    sou     = source,
    link    = link,
    nam     = name,
    var     = var,
    stat    = "point",
    units   = ifelse(var=="p", "unknown mercury inches", units(var)),
    meta    = meta.list[[var]],
    metaHead = ifelse(var=="p",
                      paste0(metaHead, " | Instrumet=unknown barometer similar to Academyof Turku"),
                      metaHead),
    time_offset = time.offset(lon),
    keep_na = FALSE
  )
}






