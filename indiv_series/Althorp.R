rm(list=ls())
library(dplyr)
library(readxl)
library(purrr)
library(tidyr)
source('/scratch2/ccorbella/code/dataprep/helpfun.R')

outdir <- "/scratch3/PALAEO-RA/daily_data/final/Althorp"


lat <- 52.28
lon <- -1
alt <- 105
name <- "Yuri_Althorp"
code <- "Althorp"
source <- "Brugnara, Y., Auchmann, R., Brönnimann, S., Allan, R. J., Auer, I., Barriendos, M., Bergström, H., Bhend, J., Brázdil, R., Compo, G. P., Cornes, R. C., Dominguez-Castro, F., van Engelen, A. F. V., Filipiak, J., Holopainen, J., Jourdain, S., Kunz, M., Luterbacher, J., Maugeri, M., Mercalli, L., Moberg, A., Mock, C. J., Pichard, G., Řezníčková, L., van der Schrier, G., Slonosky, V., Ustrnul, Z., Valente, M. A., Wypych, A., and Yin, X.: A collection of sub-daily pressure and temperature observations for the early instrumental period with a focus on the `year without a summer` 1816, Clim. Past, 11, 1027–1047, 2015."
link   <- "https://doi.org/10.5194/cp-11-1027-2015"
raw <- read_excel("/scratch3/PALAEO-RA/daily_data/original/Althorp/Althorp.xls",
                  skip=10)
df <- raw %>%
  pivot_longer(
    cols = P1:TA2, 
    names_to = c(".value", "obs.num"), 
    names_pattern = "(P|TA)(\\d)"
  ) %>% mutate(
    Hour=case_when(
      obs.num==1 ~ 8L,
      obs.num==2 ~ 18L, 
    ),

    Minute = 0L,
    phPa = round(convert_pressure(p=P, f=25.04, lat=lat, alt=alt),2),
    taC = round((TA- 32) / 1.8,2),
    meta.time=case_when(
      obs.num==1 ~ "morning",
      obs.num==2 ~ "evening", 
    ),
    
    meta.p = paste0("obs.time=", meta.time, " | orig=",P,"English in | atb=", TA, "F"),
    meta.ta = paste0("obs.time=", meta.time, " | orig=",TA, "F"),
  ) %>% select(Year, Month, Day, Hour, Minute, phPa, taC, meta.ta, meta.p)

head(df)

# save
var <-"p"
write_sef_f(
  as.data.frame(df[,c("Year","Month", "Day", "Hour", "Minute", "phPa")]),
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
  keep_na = TRUE
)

var <-"ta"
write_sef_f(
  as.data.frame(df[,c("Year","Month", "Day", "Hour", "Minute", "taC")]),
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
  keep_na = TRUE
)

