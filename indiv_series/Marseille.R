# from the file that Stefan sent, only a couple of months
rm(list=ls())
source('/scratch2/ccorbella/code/dataprep/helpfun.R')

raw <- read.delim('/scratch3/PALAEO-RA/daily_data/original/Marseille/marseille.txt',
                  fileEncoding = "UTF-16LE",
                  header=TRUE,
                  sep="\t")
raw

outdir <- '/scratch3/PALAEO-RA/daily_data/final/Marseille/'

lat	<- 43.2965
lon <- 5.36978
alt <- 44
code <- "Marseille"
name <- "Marseille"

df.ta <- raw %>% mutate(
  Hour = NA_integer_,
  Minute = NA_integer_,
  ta = round(Temperature,1)
) %>% select(Year, Month, Day, Hour, Minute, ta)


df.p <- raw %>% mutate(
  Hour = NA_integer_,
  Minute = NA_integer_,
  p = round(Pressure,1)
) %>% select(Year, Month, Day, Hour, Minute, p)

# ta
var<-"ta"
write_sef_f(
  as.data.frame(df.ta),
  outfile = outfile.name(name, var, df.ta, subdaily=FALSE),
  outpath = outdir,
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = "",
  link    = "",
  nam     = name,
  var     = var,
  stat    = "point",
  units   = units(var),
  keep_na = TRUE
)

# p
var<-"p"
write_sef_f(
  as.data.frame(df.p),
  outfile = outfile.name(name, var, df.p, subdaily=FALSE),
  outpath = outdir,
  cod     = code,
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = "",
  link    = "",
  nam     = name,
  var     = var,
  stat    = "point",
  units   = units(var),
  keep_na = TRUE
)
