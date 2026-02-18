rm(list=ls())
library(dplyr)
library(readxl)
library(purrr)
library(tidyr)
source('/scratch2/ccorbella/code/dataprep/helpfun.R')

outdir <- "/scratch3/PALAEO-RA/daily_data/final/Nuuk"


lat <- 64.1877
lon <- -51.6756
alt <- 10
name <- "Yuri_Nuuk"
code <- "Nuuk"
source <- "Brugnara, Y., Auchmann, R., Brönnimann, S., Allan, R. J., Auer, I., Barriendos, M., Bergström, H., Bhend, J., Brázdil, R., Compo, G. P., Cornes, R. C., Dominguez-Castro, F., van Engelen, A. F. V., Filipiak, J., Holopainen, J., Jourdain, S., Kunz, M., Luterbacher, J., Maugeri, M., Mercalli, L., Moberg, A., Mock, C. J., Pichard, G., Řezníčková, L., van der Schrier, G., Slonosky, V., Ustrnul, Z., Valente, M. A., Wypych, A., and Yin, X.: A collection of sub-daily pressure and temperature observations for the early instrumental period with a focus on the `year without a summer` 1816, Clim. Past, 11, 1027–1047, 2015."
link   <- "https://doi.org/10.5194/cp-11-1027-2015"
metaHead <- "coords=approx."
raw <- read_excel("/scratch3/PALAEO-RA/daily_data/original/Nuuk/Pressure_Greenland.xls",
                  skip=1,
                  col_names=c("year", "month", "day", "p1", "l1", "p2", "l2", "p3", "l3"))

df <- raw %>%
  pivot_longer(
    cols = p1:l3, 
    names_to = c(".value", "obs.num"), 
    names_pattern = "(p|l)(\\d)"
  ) %>% mutate(
    Hour=case_when(
      obs.num==1 ~ 7L,
      obs.num==2 ~ 13L,
      obs.num==3 ~ 19L 
    ),
    Minute = 0L,
    phPa = round(convert_pressure(p=(p+l/12), f=27.07, lat=lat, alt=alt),2),
    meta = paste0("obs.num=", obs.num, " | orig=",p,"in", l,"l")
  ) %>% select(year, month, day, Hour, Minute, phPa, meta)

head(df)

# save
var <-"p"
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
  metaHead = "PTC=N | PGC=Y",
  meta    = df$meta,
  keep_na = TRUE
)

