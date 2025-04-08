library(dataresqc)
library(lubridate)
library(dplyr)

source("/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/write_sef_f.R")

outdir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/'


# Bologna -----------------------------------------------------------------

df <- read.csv('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/csvs/Bologna_TMP2m.csv')

df.ta.bol <- data.frame(
  Year= year(df$Date),
  Month=month(df$Date),
  Day=day(df$Date),
  Hour='NA',
  Minute='NA',
  Value=df$TMP2m
)

meta <- list(
  ID = "169",
  Name = "Bologna",
  Lat = 44.500000,
  Lon = 11.346111,
  Alt = 53,
  Vbl = "ta",
  Units = "C",
  Link = "http://www.ecad.eu",
  Source = "Klein Tank, A.M.G. and Coauthors, 2002. Daily dataset of 20th-century surface air temperature and precipitation series for the European Climate Assessment. Int. J. of Climatol., 22, 1441-1453.-RA"
)

write_sef_f(Data=df.ta.bol, outfile="Bologna_ta.tsv",
            outpath=outdir,
            cod=meta[["ID"]],
            variable=meta[['Vbl']],
            nam=meta[["Name"]],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]], link=meta[["Link"]],
            sou=meta[["Source"]], units=meta[["Units"]], stat="point",keep_na = F
)

# CADIZ -----------------------------------------------------------------

df <- read.csv('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/csvs/Cadiz_PRMSL.csv')

df.p.cad <- data.frame(
  Year= year(df$Date),
  Month=month(df$Date),
  Day=day(df$Date),
  Hour=0,
  Minute=0,
  Value=df$PRMSL
)

meta <- list(
  ID = "Cadiz",
  Name = "Cadiz",
  Lat = 36.53,
  Lon = -5.705,
  Alt = 15,
  Vbl = "p",
  Units = "unknown",
  Source = "Final_Series_IMPROVE"
)

write_sef_f(Data=df.p.cad, outfile="Cadiz_p.tsv",
            outpath=outdir,
            cod=meta[["ID"]],
            variable=meta[['Vbl']],
            nam=meta[["Name"]],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]],
            sou=meta[["Source"]], units=meta[["Units"]], stat="point",keep_na = F
)


# CBT ---------------------------------------------------------------------

map_wrongCBT <- list(
  'gen-18' = 1801,
  '43132' = 1802,
  '43160' = 1803,
  '43191' = 1804,
  'mag-18' = 1805,
  'giu-18' = 1806,
  'lug-18' = 1807,
  'ago-18' = 1808,
  'set-18' = 1809,
  'ott-18' = 1810,
  '43405' = 1811,
  'dic-18' = 1812
)

x <- read_sef('/home/ccorbella/scratch2_symboliclink/files/1807_USBstick/not_to_use/CBT_ta_daily.tsv')
meta <- read_meta('/home/ccorbella/scratch2_symboliclink/files/1807_USBstick/not_to_use/CBT_ta_daily.tsv')

# Convert the Year column to character (if not already)
x$Year <- as.character(x$Year)

# Replace the incorrect values using your custom mapping
x$Year <- ifelse(x$Year %in% names(map_wrongCBT),
                 unlist(map_wrongCBT[x$Year]),
                 x$Year)

# Convert back to integer
x$Year <- as.integer(x$Year)

df <- data.frame(x)
df['Var'] <- NULL
df['Hour'] <- NA
df['Minute'] <- NA

# change order of columns
df <- df[,c(1,2,3,5,6,4)]

write_sef_f(Data=df, outfile="CBT_ta.tsv",
            outpath='/home/ccorbella/scratch2_symboliclink/files/1807_USBstick/',
            cod=meta[["id"]],
            variable=meta[['var']],
            nam=meta[["name"]],
            lat=meta[["lat"]],
            lon=meta[["lon"]], alt=meta[["alt"]], link=meta[['link']], period='day',
            sou=meta[["source"]], units=meta[["units"]], stat="day",keep_na = F
            )

# Central England Temperature (CET) -----------------------------------------------------------------

df <- read.csv('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/csvs/CET_TMP2m.csv')

df.ta.cet <- data.frame(
  Year= year(df$Date),
  Month=month(df$Date),
  Day=day(df$Date),
  Hour='NA',
  Minute='NA',
  Value=df$TMP2m
)

meta <- list(
  ID = "HadCET",
  Name = "Central England Temperature",
  Lat = 52.5,
  Lon = -1.9,
  Alt = 44,
  Vbl = "ta",
  Units = "C",
  Link = "https://www.metoffice.gov.uk/hadobs/hadcet/data/download.html",
  Source = "Parker, D.E., T.P. Legg, and C.K. Folland. 1992. A new daily Central England Temperature Series, 1772-1991. Int. J. Clim., Vol 12, pp 317-342"
)

write_sef_f(Data=df.ta.cet, outfile="CET_ta.tsv",
            outpath=outdir,
            cod=meta[["ID"]],
            variable=meta[['Vbl']],
            nam=meta[["Name"]],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]], link=meta[['Link']],
            sou=meta[["Source"]], units=meta[["Units"]], stat="point",keep_na = F
)


# London ------------------------------------------------------------------
path1 <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_orig/London_11_17870101-18221231_mslp.tsv'
path2 <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_orig/London_12_18230101-18411231_mslp.tsv'

sef1 <- read_sef(path1, all = TRUE)
sef2 <- read_sef(path2, all = TRUE)

meta1 <- read_meta(path1)
meta2 <- read_meta(path2)

# check that other relevant metadata are the same
meta_check_fields <- c("var", "lat", "lon", "alt", "source", "link", "units")
for (field in meta_check_fields) {
  if (!identical(meta1[[field]], meta2[[field]])) {
    warning(paste("Metadata field", field, "differs between files. Using value from file1."))
  }
}

combined_sef <- rbind(sef1, sef2)
combined_sef$Var <- NULL     # drop first col
combined_sef$Period <- NULL  # the period col is written w the write_sef fnct
combined_sef <- combined_sef[order(combined_sef$Year, combined_sef$Month, combined_sef$Day,
                                   combined_sef$Hour, combined_sef$Minute), ]

write_sef_f(combined_sef, 
            outpath = outdir, 
            outfile = "London_p_subdaily.tsv", 
            variable = "p", 
            cod = "11_and_12",                         # combined ID
            nam = meta1["name"],                       # keep name from file 1
            lat = meta1["lat"], lon = meta1["lon"], alt = meta1["alt"], 
            sou = meta1["source"], link = meta1["link"], 
            units = meta1["units"], stat = "point", 
            period = "NA", 
            metaHead = "altitude of 29.6 from 1823-01-01 to 1841-15-31 | var=mslp",
            meta = combined_sef$Meta,
            keep_na = F)


# Milan -------------------------------------------------------------------

df <- read.csv('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/csvs/Milan_PRMSL.csv')

df.p.mil <- data.frame(
  Year= year(df$Date),
  Month=month(df$Date),
  Day=day(df$Date),
  Hour='NA',
  Minute='NA',
  Value=df$PRMSL
)

meta <- list(
  ID = "IMPROVE_Milan",
  Name = "Milan",
  Lat = 45.47,
  Lon = 9.19,
  Alt = 150,
  Vbl = "p",
  Units = "unknown",
  Source = "PALAEO-RA"
)

write_sef_f(Data=df.p.mil, outfile="Milan_p.tsv",
            outpath=outdir,
            cod=meta[["ID"]],
            variable=meta[['Vbl']],
            nam=meta[["Name"]],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]],
            sou=meta[["Source"]], units=meta[["Units"]], stat="point",keep_na = F
)


# Padova -----------------------------------------------------------------

# Read metadata from the reference file
meta_ref <- read_meta('/home/ccorbella/scratch2_symboliclink/files/1807_USBstick/not_to_use/Padova_p_daily.tsv')

df <- read.csv('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/csvs/Padova_PRMSL.csv')

df.p.pad <- data.frame(
  Year= year(df$Date),
  Month=month(df$Date),
  Day=day(df$Date),
  Hour='NA',
  Minute='NA',
  Value=df$PRMSL
)

df <- read.csv('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/csvs/Padova_TMP2m.csv')

df.ta.pad <- data.frame(
  Year= year(df$Date),
  Month=month(df$Date),
  Day=day(df$Date),
  Hour='NA',
  Minute='NA',
  Value=df$TMP2m
)

# Write df to SEF using metadata from reference file
write_sef_f(df.p.pad,
            outpath = outdir,
            outfile = "Padova_p.tsv",
            variable = "p",
            cod = meta_ref["id"],
            nam = meta_ref["name"],
            lat = meta_ref["lat"],
            lon = meta_ref["lon"],
            alt = meta_ref["alt"],
            sou = meta_ref["source"],
            link = meta_ref["link"],
            units = "hPa",
            stat = "point",  # or "point" if relevant
            metaHead = meta_ref["meta"],
            keep_na = F)

write_sef_f(df.ta.pad,
            outpath = outdir,
            outfile = "Padova_ta.tsv",
            variable = "ta",
            cod = meta_ref["id"],
            nam = meta_ref["name"],
            lat = meta_ref["lat"],
            lon = meta_ref["lon"],
            alt = meta_ref["alt"],
            sou = meta_ref["source"],
            link = meta_ref["link"],
            units = "C",
            stat = "point", 
            metaHead = meta_ref["meta"],
            keep_na = F)


# Paris -------------------------------------------------------------------
df <- read.csv('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/csvs/Paris_TMP2m.csv')

df <- data.frame(
  Year= year(df$Date),
  Month=month(df$Date),
  Day=day(df$Date),
  Hour='NA',
  Minute='NA',
  Value=df$TMP2m
)

meta <- list(
  ID = "Paris",
  Name = "Paris",
  Lat = 48.86,
  Lon = 2.34,
  Alt = 42,
  Vbl = "ta",
  Units = "C",
  Link = "paris_Daily_Updated_meteo_2024_127_33_supp.xlsx from Stefan",
  Source = "Rousseau D., 2024. La série des températures journalières à Paris de 1658 à 2023"
)

write_sef_f(Data=df, outfile="Paris_ta.tsv",
            outpath=outdir,
            cod=meta[["ID"]],
            variable=meta[['Vbl']],
            nam=meta[["Name"]],
            lat=meta[["Lat"]], link=meta[['link']],
            lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]], 
            metaHead = "latlon from Peter's Paris_ta_noon.tsv",
            units=meta['Units'], stat="point",keep_na = F
)

## for pressure, change the nameof the file so it can later be read
file.rename('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/Paris_mslp.tsv',
            '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/Paris_p.tsv')

# Stockholm -----------------------------------------------------------------

df <- read.delim('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_orig/stockholm_SLP_1756_2012_hPa_hom.txt',
                 header=FALSE, fill=T, sep="", col.names= c('Year','Month','Day','PRMSL','dunno','dunno2'))

df.p.Sto <- data.frame(
  Year= df$Year,
  Month=df$Month,
  Day=df$Day,
  Hour='NA',
  Minute='NA',
  Value=df$PRMSL
)

# keep only data for period of homogenization and 3 obs/day
df <- df %>% filter(Year>=1785, Year<=1860)

meta <- list(
  ID = "IMPROVE_Stockholm",
  Name = "Stockholm Old Astronomical Observatory",
  Lat = 59.35,
  Lon = 18.05,
  Alt = 44,
  Source = "Moberg A, Bergström H, Ruiz Krigsman J, Svanered O. 2002: Daily air temperature and pressure series for Stockholm (1756-1998). Climatic Change 53: 171-212",
  Link="http://people.su.se/~amobe/stockholm/stockholm-historical-weather-observations-ver-1.0"
)

write_sef_f(Data=df.p.Sto, outfile="Stockholm_p.tsv",
            meta="hom.daymean | orig.obs=3",
            outpath=outdir,
            cod=meta[["ID"]],
            variable="p",
            nam=meta[["Name"]], link=meta[['Link']],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]],
            metaHead = "Homogenized Sea Level Pressure 1756-2012, i.e. observed air pressure reduced to 0 degC, normal gravity and sea level. Additionally 	corrected for known and supposed biased barometers and also homogenized against reference stations",
            units="hPa", stat="mean",keep_na = F
)

# TORINO -----------------------------------------------------------------

df <- read.csv('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_orig/TorinoPR792_953.csv',sep='\t')

df.p.Tor <- data.frame(
  Year=df$a,
  Month=df$m,
  Day=df$g,
  Hour=df$ora,
  Minute='NA',
  Value=df$Pressione.ridotta.al.l.d.m...hPa.
)

df.p.Tor <- df.p.Tor %>% filter(Year>=1803,Year<=1866)

meta <- list(
  ID = "TOR",
  Name = "Torino",
  Lat = 45.06813258066838,
  Lon = 7.68408326971596,
  Alt = 254,
  Source = "Yuri-stations-01/Turin/; Contatti: Società Meteorologica Italiana, info@nimbus.it",
  Link="http://people.su.se/~amobe/stockholm/stockholm-historical-weather-observations-ver-1.0"
)

write_sef_f(Data=df.p.Sto, outfile="Torino_p.tsv",
            outpath=outdir,
            cod=meta[["ID"]],
            variable="p",
            nam=meta[["Name"]], link=meta[['Link']],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]], metaHead = "quota barometro (m)=280.90",
            units="hPa", stat="point",keep_na = F
)


# UPPSALA -----------------------------------------------------------------

df <- read.csv('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/csvs/Uppsala_PRMSL.csv')

df.p.Upp <- data.frame(
  Year= year(df$Date),
  Month=month(df$Date),
  Day=day(df$Date),
  Hour='NA',
  Minute='NA',
  Value=df$PRMSL
)

meta <- list(
  ID = "IMPROVE_Uppsala",
  Name = "Uppsala",
  Lat = 59.8605,
  Lon = 17.6407,
  Alt = 15,
  Source = "PALEO-RA"
)

write_sef_f(Data=df.p.Upp, outfile="Uppsala_p.tsv",
            outpath=outdir,
            cod=meta[["ID"]],
            variable="p",
            nam=meta[["Name"]],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]], metaHead = "from Final_Series_Improve",
            units="unknown", stat="point",keep_na = F
)


# VALÈNCIA ----------------------------------------------------------------
df <- read.csv('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_orig/Valencia/Valencia_concatenated.csv',
               header=T)

# Create a list (equivalent to a dictionary in R)
meta <- list(
  ID = "Dom_Valencia",
  Name = "Valencia",
  Lat = 39.47,
  Lon = -0.38,
  Alt = 25,
  Source = "Dominguez_Castro_2014",
  Link = "https://docta.ucm.es/entities/publication/b26dac98-5ffe-482c-9419-2f94168cc7eb"
)

df.ta.Valencia <- data.frame(
  Year=df$Year,
  Month=df$Month,
  Day=df$Day,
  Hour=df$Hour,
  Minute=0,
  Value=as.numeric(df$Term) * 1.25
)

df.p.Valencia <- data.frame(
  Year=df$Year,
  Month=df$Month,
  Day=df$Day,
  Hour=df$Hour,
  Minute=0,
  Value=(as.numeric(df$Bar.p.) + as.numeric(df$Bar.l.) /12 ) * 27.07 * 1.3332239
)

# drop pressures beyond the 5 sigma range
p.mean <- mean(df.p.Valencia$Value, na.rm = T)
p.std  <- sd(df.p.Valencia$Value, na.rm = T)

df.p.Valencia <- subset(df.p.Valencia, Value < p.mean + 5*p.std)
df.p.Valencia <- subset(df.p.Valencia, Value > p.mean - 5*p.std)

write_sef_f(Data=df.ta.Valencia,outfile="Valencia_ta_subdaily.tsv",
            outpath=outdir,
            cod=meta[["ID"]],
            variable="ta",
            nam=meta[["Name"]],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]], metaHead = "orig_ta=Reaumur",
            link=meta[["Link"]], units="C", stat="point",
            meta="orig_ta=Reaumur", keep_na = F)

write_sef_f(Data=df.p.Valencia,outfile="Valencia_p_subdaily.tsv",
            outpath=outdir,
            cod=meta[["ID"]],
            variable="p",
            nam=meta[["Name"]],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]], metaHead = "orig_p=Paris inch",
            link=meta[["Link"]], units="unknown", stat="point",
            meta="orig_p=Paris inch | qc=±5σ", keep_na = F)

# Vienna ---------------------------------------------------------------
df <- read.csv('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_orig/Wien Sternwarte_5905_TAG_18150101_18171231.csv',
               header=T, na.strings = c("", "NA"))
df <- df %>% na.omit()

## Anbei die Daten 1815 – 1817 der Station Wien Sternwarte, 3 mal tägliche Messungen Lufttemperatur und Luftdruck. 
## Die Daten sind in den vorgesehenen Spalten von 7,14 und 19 Uhr, soweit es dokumentiert ist wurde aber
## um um 8,15 und 22 Uhr gemessen.

df$datum <- ymd(df$datum) # Convert date column

# Reshape temperature values to long format (t7 = 7h, t14 = 14h, t19 = 19h)
df.ta <- df %>%
  select(datum, t7, t14, t19) %>%
  pivot_longer(cols = starts_with("t"), names_to = "HourTag", values_to = "Value") %>%
  mutate(
    Year = year(datum),
    Month = month(datum),
    Day = day(datum),
    Hour = case_when(
      HourTag == "t7" ~ 8,
      HourTag == "t14" ~ 15,
      HourTag == "t19" ~ 22
    ),
    Value = as.numeric(Value)/10,
  )

df.ta <-data.frame(
  Year=df.ta$Year,
  Month=df.ta$Month,
  Day=df.ta$Day,
  Hour=df.ta$Hour,
  Minute='0',
  Value=df.ta$Value
)

meta <- list(
  ID = "Wien Sternwarte_5905",
  Name = "Vienna",
  Lat = "48.20909573207442",
  Lon = "16.37758981033942",
  Alt = "198.5",
  Source = "Yuri Brugnara",
  meta = "observer=Tiesnecker & team | obs=in 32m tall tower, 171+32m"
)
  

write_sef_f(Data=df.ta, outfile="Vienna_ta_subdaily.tsv",
            outpath=outdir,
            cod=meta[["ID"]],
            variable="ta",
            nam=meta[["Name"]],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]], metaHead = meta[['meta']],
            units="C", stat="point", keep_na = F
)


df.p <- df %>%
  select(datum, druck07, druck14, druck19) %>%
  pivot_longer(cols = starts_with("druck"), names_to = "HourTag", values_to = "Value") %>%
  mutate(
    Year = year(datum),
    Month = month(datum),
    Day = day(datum),
    Hour = case_when(
      HourTag == "druck07" ~ 8,
      HourTag == "druck14" ~ 15,
      HourTag == "druck19" ~ 22
    ),
    Value = as.numeric(Value)/1000,
  )

df.p <-data.frame(
  Year=df.p$Year,
  Month=df.p$Month,
  Day=df.p$Day,
  Hour=df.p$Hour,
  Minute='0',
  Value=df.p$Value
)

write_sef_f(Data=df.ta, outfile="Vienna_p_subdaily.tsv",
            outpath=outdir,
            cod=meta[["ID"]],
            variable="p",
            nam=meta[["Name"]],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]], metaHead = meta[['meta']],
            units="hPa", stat="point", keep_na = F
)

# Ylitornio ---------------------------------------------------------------
df <- read.csv('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/csvs/Ylitornio_all.csv',
               header=T)

df.ta.Yli <- data.frame(
  Year=df$Year,
  Month=df$Month,
  Day=df$Day,
  Hour='NA',
  Minute='NA',
  Value=df$TMP2m
)

df.p.Yli <- data.frame(
  Year=df$Year,
  Month=df$Month,
  Day=df$Day,
  Hour='NA',
  Minute='NA',
  Value=df$PRMSL
)

meta <- list(
  ID = "Ylitornio",
  Name = "Ylitornio",
  Lat = "66.319266",
  Lon = "23.670970",
  Alt = "55",
  Source = "Helama, S., Holopainen, J., Timonen, M., Ogurtsov, M. G., Lindholm, M., Meriläinen, J., Eronen, M. (2004). Comparison of living-tree and subfossil ringwidths with summer temperatures from 18th, 19th and 20th centuries in Northern Finland. Dendrochronologia 21/3, 147 - 154.",
  meta = "observer=Johan Portin | orig_ta=C | orig_p=unknown(similar to Turku Academy measurements) | latlon extracted from Peter's Yli"
)

write_sef_f(Data=df.ta.Yli, outfile="Ylitornio_ta.tsv",
            outpath=outdir,
            cod=meta[["ID"]],
            variable="ta",
            nam=meta[["Name"]],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]], metaHead = meta[['meta']],
            units="C", stat="point",
            meta="orig_ta=C", keep_na = F
)

write_sef_f(Data=df.p.Yli, outfile="Ylitornio_p.tsv",
            outpath=outdir,
            cod=meta[["ID"]],
            variable="p",
            nam=meta[["Name"]],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]], metaHead = meta[['meta']],
            units="unknown", stat="point", keep_na = F
)


# Zwanenburg --------------------------------------------------------------
library(strucchange)

f <- '/home/ccorbella/scratch2_symboliclink/files/1807_USBstick/KNMI_KNMI-42_Zwanenburg_17380101-18601231_p_daily.tsv'
x <- read_sef(f, all=T)
meta <- read_meta(f)

x <- x[!is.na(x$Value), ]

# Create a daily time series vector
dates <- as.Date(paste(x$Year, x$Month, x$Day, sep = "-"))
pressure <- x$Value

# Combine into a data frame
df <- data.frame(date = dates, pressure = pressure)

# Run breakpoint detection (intercept-only model)
bp_model <- breakpoints(pressure ~ 1, breaks=5)

# View summary and breakpoints
summary(bp_model)

# Plot results
plot(bp_model, main = "Breakpoints in Daily Pressure Time Series")
lines(df$pressure, col = "blue")  # optional: actual series overlay

# Optional: plot segmented regression fit
plot(df$date, pressure, type = "l", col = "gray",
     main = "Segmented Mean Pressure with Detected Breaks",
     xlab = "Date", ylab = "Pressure")
lines(fitted(bp_model, breaks = bp_model$breakpoints), col = "red", lwd = 2)