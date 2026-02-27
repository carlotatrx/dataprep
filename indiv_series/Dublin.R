rm(list=ls())
library(dataresqc)
library(readxl)
library(tidyr)
library(dplyr)

source('/scratch2/ccorbella/code/dataprep/helpfun.R')

indir <- "/scratch3/PALAEO-RA/daily_data/original/"
name <- "Dublin"

outdir <- paste0("/scratch3/PALAEO-RA/daily_data/final/",name)

source <- "Mateus, C., Potito, A., & Curley, M. (2020). Reconstruction of a long‐term historical daily maximum and minimum air temperature network dataset for Ireland (1831‐1968). Geoscience Data Journal, 7(2), 102-115."
link   <- "https://www.edepositireland.ie/entities/publication/7271837d-dbe6-463f-b103-422f42a15446"



# RDO Dunsink Dublin ------------------------------------------------------


lat <- 53.387
lon <- -6.338
alt <- 64 # m (210 ft)

metaHead <- paste0(
  "Observer=Observatory Staff | Instrument=Max/Min thermometers | ",
  "Location=Royal Dublin Observatory, Dunsink; Hilltop rural exposure northwest of city | ",
  "Notes=0.1F precision"
)

metaHead


raw <- read.csv(paste0(indir, name, "/RDO Dunsink_1851.csv"), header=FALSE,
                skip=1,fileEncoding = "latin1",
                col.names=c("year", "month", "day", "maxF", "minF", "maxC", "minC", "extra"))

head(raw)
tail(raw)

df <- raw %>%
  mutate(
    hour = NA,
    minute = NA,
    obs_time_label = "9a.m.",
    meta.Tx = case_when(
      !is.na(maxF) ~ paste0("orig=", maxF, "F | obs.time=", obs_time_label),
      is.na(maxF)  ~ ""
    ),
    meta.Tn = case_when(
      !is.na(minF) ~ paste0("orig=", minF, "F | obs.time=", obs_time_label),
      is.na(minF)  ~ ""
    )
  ) %>% 
  select(year, month, day, hour, minute, Tx = maxC, Tn = minC, meta.Tx, meta.Tn) %>%
  filter(!is.na(Tx) | !is.na(Tn))

head(df)
tail(df)

var <- "Tx"
df.Tx <- df[c("year","month", "day", "hour", "minute","Tx", "meta.Tx")] %>% filter(!is.na(Tx))
write_sef_f(
  as.data.frame(df.Tx),
  outfile = outfile.name(paste0("ILMMT_",name, "-RDO-Dunsink"), var, df.Tx, FALSE),
  outpath = outdir,
  cod     = paste0("ILMMT-",name, "-RDO-Dunsink"),
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "minimum",
  period    = "day",
  units   = units(var),
  meta = df.Tx$meta.Tx,
  metaHead=metaHead,
  keep_na=TRUE
)


var <- "Tn"
df.Tn <- df[c("year","month", "day", "hour", "minute","Tn", "meta.Tn")] %>% filter(!is.na(Tn))
write_sef_f(
  as.data.frame(df.Tn),
  outfile = outfile.name(paste0("ILMMT_",name, "-RDO-Dunsink"), var, df.Tn, FALSE),
  outpath = outdir,
  cod     = paste0("ILMMT-",name,"-RDO-Dunsink"),
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "minimum",
  period    = "day",
  units   = units(var),
  meta = df.Tn$meta.Tn,
  metaHead=metaHead,
  keep_na=TRUE
)


# RCS Dublin --------------------------------------------------------------

lat <- 53.339
lon <- -6.261
alt <- 11 # m (approximate elevation of St. Stephen's Green)

metaHead <- paste0(
  "Observer=John Evans, supervised by Prof. James Apjohn | ",
  "Instrument=Max/Min thermometers | Location=Courtyard of Royal College of Surgeons, ",
  "Dublin | DataSource=Dublin Medical Press | Notes=whole degree precision"
)

metaHead


raw <- read.csv(paste0(indir, name, "/RCS Dublin_1841-1857.csv"), header=FALSE,
                skip=1,fileEncoding = "latin1",
                col.names=c("year", "month", "day", "maxF", "minF", "maxC", "minC", "extra"))

head(raw)
tail(raw)

df <- raw %>%
  mutate(
    hour = NA,
    minute = NA,
    obs_time_label = "9a.m.",
    meta.Tx = case_when(
      !is.na(maxF) ~ paste0("orig=", maxF, "F"),
      is.na(maxF)  ~ ""
    ),
    meta.Tn = case_when(
      !is.na(minF) ~ paste0("orig=", minF, "F"),
      is.na(minF)  ~ ""
    )
  ) %>% 
  select(year, month, day, hour, minute, Tx = maxC, Tn = minC, meta.Tx, meta.Tn) %>%
  filter(!is.na(Tx) | !is.na(Tn))

head(df)
tail(df)

var <- "Tx"
df.Tx <- df[c("year","month", "day", "hour", "minute","Tx", "meta.Tx")] %>% filter(!is.na(Tx))
write_sef_f(
  as.data.frame(df.Tx),
  outfile = outfile.name(paste0("ILMMT_",name, "-RCS"), var, df.Tx, FALSE),
  outpath = outdir,
  cod     = paste0("ILMMT-",name, "-RCS"),
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "minimum",
  period    = "day",
  units   = units(var),
  meta = df.Tx$meta.Tx,
  metaHead=metaHead,
  keep_na=TRUE
)


var <- "Tn"
df.Tn <- df[c("year","month", "day", "hour", "minute","Tn", "meta.Tn")] %>% filter(!is.na(Tn))
write_sef_f(
  as.data.frame(df.Tn),
  outfile = outfile.name(paste0("ILMMT_",name, "-RCS"), var, df.Tn, FALSE),
  outpath = outdir,
  cod     = paste0("ILMMT-",name,"-RCS"),
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "minimum",
  period    = "day",
  units   = units(var),
  meta = df.Tn$meta.Tn,
  metaHead=metaHead,
  keep_na=TRUE
)

# Grafton Street Dublin ---------------------------------------------------

lat <- 53.342
lon <- -6.259
alt <- 12 # m (approximate elevation of Grafton Street)
metaHead <- paste0("Observer=George Yeates (1843-1849) | Instrument=Rutherford’s self-registering thermometers | Location=Grafton Street | Note=George Yeates was an instrument maker;",
                   " whole degree precision | DataSource=Transactions of the Royal Irish Academy."
)

metaHead


raw <- read.csv(paste0(indir, name, "/Grafton Street Dublin_1843-1849.csv"),
                col.names=c("year", "month", "day", "maxF", "minF", "maxC", "minC"))
head(raw)
tail(raw)

df <- raw %>%
  mutate(
    hour = NA,
    minute = NA,
    obs_time_label = "10a.m.",
    # Meta logic including the standardized instrument source
    meta.Tx = case_when(
      !is.na(maxF) ~ paste0("orig=", maxF, "F | obs.time=", obs_time_label),
      is.na(maxF)  ~ ""
    ),
    meta.Tn = case_when(
      !is.na(minF) ~ paste0("orig=", minF, "F | obs.time=", obs_time_label),
      is.na(minF)  ~ ""
    )
  ) %>% 
  select(year, month, day, hour, minute, Tx = maxC, Tn = minC, meta.Tx, meta.Tn)

head(df)
tail(df)

var <- "Tn"
df.Tn <- df[c("year","month", "day", "hour", "minute","Tn", "meta.Tn")] %>% filter(!is.na(Tn))
write_sef_f(
  as.data.frame(df.Tn),
  outfile = outfile.name(paste0("ILMMT_",name, "-GraftonStreet_Yaetes"), var, df.Tn, FALSE),
  outpath = outdir,
  cod     = paste0("ILMMT-",name, "-GraftonStreet_Yaetes"),
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "minimum",
  period    = "day",
  units   = units(var),
  meta = df.Tn$meta.Tn,
  metaHead=metaHead,
  keep_na=TRUE
)


var <- "Tx"
df.Tx <- df[c("year","month", "day", "hour", "minute","Tx", "meta.Tx")] %>% filter(!is.na(Tx))

write_sef_f(
  as.data.frame(df.Tx),
  outfile = outfile.name(paste0("ILMMT_",name, "-GraftonStreet_Yaetes"), var, df.Tx, FALSE),
  outpath = outdir,
  cod     = paste0("ILMMT-",name, "-GraftonStreet_Yaetes"),
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "maximum",
  period    = "day",
  units   = units(var),
  metaHead=metaHead,
  meta=df.Tx$meta.Tx,
  keep_na=TRUE
)



# Dublin Commercial Buildings ---------------------------------------------

lat <- 53.344
lon <- -6.263
alt <- 10 # (approximate, based on Dame Street elevation)
metaHead <- paste0(
  "Observer=Clerks at the Commercial Buildings | ",
  "Instrument=Rutherford’s minimum thermometer | ",
  "Location=Dame Street | DataSource=",
  "Dublin Evening Mail newspaper archives | Note=whole degree precision."
)
raw <- read.csv(paste0(indir, name, "/Commercial Buildings Dublin_1849-1858.csv"),
                col.names=c("year", "month", "day", "maxF", "minF", "maxC", "minC"))
head(raw)
tail(raw)

df <- raw %>%
  mutate(
    hour = NA,
    minute = NA,
    obs_time_label = "10a.m.",
    # Meta logic including the standardized instrument source
    meta.Tx = case_when(
      !is.na(maxF) ~ paste0("orig=", maxF, "F | obs.time=", obs_time_label),
      is.na(maxF)  ~ ""
    ),
    meta.Tn = case_when(
      !is.na(minF) ~ paste0("orig=", minF, "F | obs.time=", obs_time_label),
      is.na(minF)  ~ ""
    )
  ) %>% 
  select(year, month, day, hour, minute, Tx = maxC, Tn = minC, meta.Tx, meta.Tn) %>%
  filter(!is.na(Tx) | !is.na(Tn))

head(df)
tail(df)

var <- "Tn"
write_sef_f(
  as.data.frame(df[c("year","month", "day", "hour", "minute","Tn")]),
  outfile = outfile.name(paste0("ILMMT_",name, "-DameStreet"), var, df, FALSE),
  outpath = outdir,
  cod     = paste0("ILMMT-",name, "-DameStreet"),
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "minimum",
  period    = "day",
  units   = units(var),
  meta = df$meta.Tn,
  metaHead=metaHead,
  keep_na=TRUE
)


var <- "Tx"
write_sef_f(
  as.data.frame(df[c("year","month", "day", "hour", "minute",var)]),
  outfile = outfile.name(paste0("ILMMT_",name,"-DameStreet"), var, df, FALSE),
  outpath = outdir,
  cod     = paste0("ILMMT-",name,"-DameStreet"),
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "maximum",
  period    = "day",
  units   = units(var),
  metaHead=metaHead,
  meta=df$meta.Tx,
  keep_na=TRUE
)

# Dublin Phoenix Park -----------------------------------------------------
lat <- 53.364
lon <- -6.348
alt <- 46

metaHead <- paste0(
  "Observer=Staff of the Royal Engineers, Ordnance Survey Office (1831-1900) | ",
  "Instrument=Octagonal wooden thermometer screen (1831-1860s), Stevenson screen (introduced late 19th century),",
  "Max and Min self-registering thermometers"
)

raw <- read.csv(paste0(indir, name, "/Phoenix Park Dublin_1831-1958.csv"),
                col.names=c("year", "month", "day", "maxF", "minF", "maxC", "minC"))

head(raw)
tail(raw)

df <- raw %>%
  mutate(
    hour=NA,
    minute=NA,
    meta.Tx=paste0("orig=",maxF, "F | obs.time=9a.m."),
    meta.Tn=paste0("orig=",minF, "F | obs.time=9a.m."),
  ) %>% select(year, month, day, hour, minute, Tx=maxC, Tn=minC, meta.Tx, meta.Tn) %>% 
  filter(
    # only drop where both are NA
    !is.na(Tx) | !is.na(Tn)
  )

head(df)
tail(df)

var <- "Tn"
write_sef_f(
  as.data.frame(df[c("year","month", "day", "hour", "minute","Tn")]),
  outfile = outfile.name(paste0("ILMMT_",name, "-PhoenixPark"), var, df, FALSE),
  outpath = outdir,
  cod     = paste0("ILMMT-",name, "-PhoenixPark"),
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "minimum",
  period    = "day",
  units   = units(var),
  metaHead = metaHead,
  meta = df$meta.Tn,
  keep_na=TRUE
)

var <- "Tx"
write_sef_f(
  as.data.frame(df[c("year","month", "day", "hour", "minute",var)]),
  outfile = outfile.name(paste0("ILMMT_",name, "-PhoenixPark"), var, df, FALSE),
  outpath = outdir,
  cod     = paste0("ILMMT-",name, "-PhoenixPark"),
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "maximum",
  period    = "day",
  units   = units(var),
  metaHead = metaHead,
  meta=df$meta.Tx,
  keep_na=TRUE
)



# Dublin Trinity College --------------------------------------------------
lat <- 53.344
lon <- -6.258
alt <- 8

metaHead <- paste0(
  "Observer=Humphrey Lloyd (1840-1851), Dr. Erasmus Dixon (1904-1920) | ",
  "Instrument=Magnetical Observatory (early), Stevenson screen (late) | ",
  "Location=Magnetical Observatory, Fellows' Garden (1840-1851), School of Physic Garden (1904-1920); Precision shift from 0.1F to whole degrees in 1915."
)

raw <- read.csv(paste0(indir, name, "/Trinity College Dublin_1840-1959.csv"),
                col.names=c("year", "month", "day", "maxF", "minF", "maxC", "minC"))

head(raw)
tail(raw)

df <- raw %>%
  mutate(
    hour=NA,
    minute=NA,
    meta.Tx=paste0("orig=",maxF, "F | obs.time=9a.m."),
    meta.Tn=paste0("orig=",minF, "F | obs.time=9a.m."),
  ) %>% select(year, month, day, hour, minute, Tx=maxC, Tn=minC, meta.Tx, meta.Tn)

head(df)
tail(df)

var <- "Tn"
write_sef_f(
  as.data.frame(df[c("year","month", "day", "hour", "minute","Tn")]),
  outfile = outfile.name(paste0("ILMMT_",name, "-TrinityCollege"), var, df, FALSE),
  outpath = outdir,
  cod     = paste0("ILMMT-",name, "-TrinityCollege"),
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "minimum",
  period    = "day",
  units   = units(var),
  metaHead = metaHead,
  meta = df$meta.Tn,
  keep_na=FALSE
)

var <- "Tx"
write_sef_f(
  as.data.frame(df[c("year","month", "day", "hour", "minute",var)]),
  outfile = outfile.name(paste0("ILMMT_",name, "-TrinityCollege"), var, df, FALSE),
  outpath = outdir,
  cod     = paste0("ILMMT-",name, "-TrinityCollege"),
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "maximum",
  period    = "day",
  units   = units(var),
  metaHead = metaHead,
  meta=df$meta.Tx,
  keep_na=FALSE
)

# Dublin botanic Gardens --------------------------------------------------
lat <- 53.372
lon <- -6.721
alt <- 18

raw <- read.csv(paste0(indir, name, "/Botanic Gardens Dublin_1834-1958.csv"),
                col.names=c("year", "month", "day", "maxF", "minF", "maxC", "minC"))
head(raw)

df <- raw %>%
  mutate(
    hour=NA,
    minute=NA,
    meta.Tx=ifelse(year>=1911, paste0("orig=",maxF, "F | obs.time=9p.m."),paste0("orig=",maxF, "F")),
    meta.Tn=ifelse(year>=1911, paste0("orig=",minF, "F | obs.time=9p.m."),paste0("orig=",minF, "F")),
  ) %>% select(year, month, day, hour, minute, Tx=maxC, Tn=minC, meta.Tx, meta.Tn)

head(df)
tail(df)

var <- "Tn"
write_sef_f(
  as.data.frame(df[c("year","month", "day", "hour", "minute","Tn")]),
  outfile = outfile.name(paste0("ILMMT_",name), var, df, FALSE),
  outpath = outdir,
  cod     = paste0("ILMMT-",name, "-BotanicGardens"),
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "minimum",
  period    = "day",
  units   = units(var),

  meta = df$meta.Tn
)


var <- "Tx"
write_sef_f(
  as.data.frame(df[c("year","month", "day", "hour", "minute",var)]),
  outfile = outfile.name(paste0("ILMMT_",name), var, df, FALSE),
  outpath = outdir,
  cod     = paste0("ILMMT-",name, "-BotanicGardens"),
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "maximum",
  period    = "day",
  units   = units(var),
  metaHead = paste0(
    "Observer=John Underwood (1819-1830s), David Moore (1848-1879), Frederick Moore (1879-1920), ",
    "David M'Ardle (1884) | Instrument=Barometer by James Lynch (1813), Stevenson Screen (verified 1880s)"
  ),
  meta=df$meta.Tx
)


# Dublin NLI series -------------------------------------------------------

raw <- read.csv(paste0(indir, name, "/Botanic Gardens Dublin NLI series_1882-1952.csv"),
                col.names=c("year", "month", "day", "maxF", "minF", "maxC", "minC"))
head(raw)

df <- raw %>%
  mutate(
    hour=NA,
    minute=NA,
    meta.Tx=paste0("orig=",maxF, "F | obs.time=9-10a.m."),
    meta.Tn=paste0("orig=",minF, "F | obs.time=9-10a.m."),
  ) %>% select(year, month, day, hour, minute, Tx=maxC, Tn=minC, meta.Tx, meta.Tn) %>% drop_na(Tx)

head(df)

tail(df)

var <- "Tn"
write_sef_f(
  as.data.frame(df[c("year","month", "day", "hour", "minute","Tn")]),
  outfile = outfile.name(paste0("ILMMT_",name,"-BotanicGardens_NLI"), var, df, FALSE),
  outpath = outdir,
  cod     = paste0("ILMMT-",name, "-BotanicGardens_NLI"),
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "minimum",
  period    = "day",
  units   = units(var),
  meta = df$meta.Tn,
  metaHead=metaHead <- paste0(
    "Observer=John Underwood (1819-1830s), David Moore (1848-1879), Frederick Moore (1879-1920), ",
    "David M'Ardle (1884) | Instrument=Max and Min thermometers (NLI series 1882-1920). | Note=Weekly manuscripts recorded at 9 or 10 a.m."
  ),
  keep_na=FALSE
)


var <- "Tx"
write_sef_f(
  as.data.frame(df[c("year","month", "day", "hour", "minute",var)]),
  outfile = outfile.name(paste0("ILMMT_",name,"-BotanicGardens_NLI"), var, df, FALSE),
  outpath = outdir,
  cod     = paste0("ILMMT-",name, "-BotanicGardens_NLI"),
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "maximum",
  period    = "day",
  units   = units(var),
  metaHead=metaHead <- paste0(
    "Observer=John Underwood (1819-1830s), David Moore (1848-1879), Frederick Moore (1879-1920), ",
    "David M'Ardle (1884) | Instrument=Max and Min thermometers (NLI series 1882-1920). | Note=Weekly manuscripts recorded at 9 or 10 a.m."
  ),
  meta=df$meta.Tx,
  keep_na=FALSE
)


# Dublin 1880 botanic gardens ---------------------------------------------

raw <- read.csv(paste0(indir, name, "/Botanic Gardens Dublin_1880 .csv"),
                col.names=c("year", "month", "day", "9maxF", "9minF", "9maxC", "9minC", "21maxF", "21minF", "21maxC", "21minC"))

head(raw)

df <- raw %>%
  pivot_longer(cols=starts_with("X"),
               names_to = "name",
               values_drop_na = TRUE
  ) %>%
  # 1. Extract the Hour and the Metric into two columns
  # This regex looks for the digits (Hour) and the remaining letters (Metric)
  extract(name, 
          into = c("hour", "Metric"), 
          regex = "X?(\\d+)(.*)") %>%
  # 2. Pivot the Metric back out so maxF, maxC, etc. are columns
  pivot_wider(
    names_from = Metric, 
    values_from = value
  ) %>% mutate(
    hour=as.integer(hour),
    minute=0L,
    meta.Tx=paste0("orig=",maxF, "F"),
    meta.Tn=paste0("orig=",minF, "F")
  ) %>% select(year, month, day, hour, minute, Tx=maxC, Tn=minC, meta.Tx, meta.Tn)

head(df)



metaHead_1880 <- paste0(
  "Observer=Frederick Moore (1880) | Instrument=Stevenson Screen, ",
  "verified Max/Min thermometers | Note=Transition year from RDS abstracts ",
  "to daily handwritten registers; screen location at Glasnevin."
)

var <- "Tn"
write_sef_f(
  as.data.frame(df[c("year","month", "day", "hour", "minute","Tn")]),
  outfile = outfile.name(paste0("ILMMT_",name,"-BotanicGardens"), var, df, TRUE),
  outpath = outdir,
  cod     = paste0("ILMMT-",name, "-BotanicGardens_NLI"),
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "minimum",
  period    = "point",
  units   = units(var),
  meta = df$meta.Tn,
  metaHead=metaHead_1880,
  keep_na=FALSE
)


var <- "Tx"
write_sef_f(
  as.data.frame(df[c("year","month", "day", "hour", "minute",var)]),
  outfile = outfile.name(paste0("ILMMT_",name,"-BotanicGardens"), var, df, TRUE),
  outpath = outdir,
  cod     = paste0("ILMMT-",name, "-BotanicGardens"),
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "maximum",
  period    = "point",
  units   = units(var),
  metaHead=metaHead_1880,
  meta=df$meta.Tx,
  keep_na=FALSE
)


# Dublin Fitzwilliam Square ---------------------------------------------
raw <- read.csv(paste0(indir, name, "/Fitzwilliam Square Dublin_1871-1937.csv"),
                col.names=c("year", "month", "day", "maxF", "minF", "maxC", "minC"))
head(raw)

lat <- 53.336
lon <- -6.252
alt <- 15.8

df <- raw %>%
  mutate(
    hour=NA,
    minute=NA,
    meta.Tx=paste0("orig=",maxF, "F | obs.time=9p.m."),
    meta.Tn=paste0("orig=",minF, "F | obs.time=9a.m."),
  ) %>% select(year, month, day, hour, minute, Tx=maxC, Tn=minC, meta.Tx, meta.Tn) %>% drop_na(Tx)

head(df)

tail(df)

metaHead <- paste0(
  "Observer=Sir John William Moore (1871-1920) | Instrument=Stevenson screen, ",
  "Casella Max/Min thermometers",
  "Location=40 Fitzwilliam Square West."
)

var <- "Tn"
write_sef_f(
  as.data.frame(df[c("year","month", "day", "hour", "minute","Tn")]),
  outfile = outfile.name(paste0("ILMMT_",name,"-FitzwilliamSquare_Moore"), var, df, FALSE),
  outpath = outdir,
  cod     = paste0("ILMMT-",name, "-FitzwilliamSquare"),
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "minimum",
  period    = "day",
  units   = units(var),
  meta = df$meta.Tn,
  metaHead=metaHead,
  keep_na=FALSE
)


var <- "Tx"
write_sef_f(
  as.data.frame(df[c("year","month", "day", "hour", "minute",var)]),
  outfile = outfile.name(paste0("ILMMT_",name,"-FitzwilliamSquare_Moore"), var, df, FALSE),
  outpath = outdir,
  cod     = paste0("ILMMT-",name, "-FitzwilliamSquare"),
  lat     = lat,
  lon     = lon,
  alt     = alt,
  sou     = source,
  link    = link,
  nam     = name,
  var     = var,
  stat    = "maximum",
  period    = "day",
  units   = units(var),
  metaHead=metaHead,
  meta=df$meta.Tx,
  keep_na=FALSE
)
