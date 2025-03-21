library(dataresqc)
library(lubridate)

source("/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/write_sef_f.R")

outdir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/'


# Bologna -----------------------------------------------------------------

df <- read.csv('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/Bologna_TMP2m.csv')

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

df <- read.csv('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/Cadiz_PRMSL.csv')

df.p.cad <- data.frame(
  Year= year(df$Date),
  Month=month(df$Date),
  Day=day(df$Date),
  Hour='NA',
  Minute='NA',
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

write_sef_f(Data=df.ta.bol, outfile="Cadiz_p.tsv",
            outpath=outdir,
            cod=meta[["ID"]],
            variable=meta[['Vbl']],
            nam=meta[["Name"]],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]],
            sou=meta[["Source"]], units=meta[["Units"]], stat="point",keep_na = F
)

# Central England Temperature (CET) -----------------------------------------------------------------

df <- read.csv('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/CET_TMP2m.csv')

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

write_sef_f(Data=df.ta.bol, outfile="CET_ta.tsv",
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
combined_sef <- combined_sef[order(combined_sef$Year, combined_sef$Month, combined_sef$Day,
                                   combined_sef$Hour, combined_sef$Minute), ]

write_sef(combined_sef, 
          outpath = outdir, 
          outfile = "London_p_subdaily.tsv", 
          variable = meta1["var"], 
          cod = "11_and_12",                         # combined ID
          nam = meta1["name"],                       # keep name from file 1
          lat = meta1["lat"], lon = meta1["lon"], alt = meta1["alt"], 
          sou = meta1["source"], link = meta1["link"], 
          units = meta1["units"], stat = "point", 
          period = "NA", 
          metaHead = "altitude of 29.6 from 1823-01-01 to 1841-15-31",
          meta = combined_sef$Meta,
          keep_na = F)


# Milan -------------------------------------------------------------------

df <- read.csv('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/Milan_PRMSL.csv')

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

write_sef_f(Data=df.ta.bol, outfile="Milan_p.tsv",
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

df <- read.csv('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/Padova_PRMSL.csv')

df.p.pad <- data.frame(
  Year= year(df$Date),
  Month=month(df$Date),
  Day=day(df$Date),
  Hour='NA',
  Minute='NA',
  Value=df$PRMSL
)

df <- read.csv('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/Padova_TMP2m.csv')

df.ta.pad <- data.frame(
  Year= year(df$Date),
  Month=month(df$Date),
  Day=day(df$Date),
  Hour='NA',
  Minute='NA',
  Value=df$TMP2m
)

# Write df to SEF using metadata from reference file
write_sef(df.p.pad,
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

write_sef(df.ta.pad,
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
          meta = df$Meta,  # use per-row meta if df includes it
          keep_na = F)

# UPPSALA -----------------------------------------------------------------

df <- read.csv('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/Uppsala_PRMSL.csv')

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
  Minute='NA',
  Value=as.numeric(df$Term) * 1.25
)

df.p.Valencia <- data.frame(
  Year=df$Year,
  Month=df$Month,
  Day=df$Day,
  Hour=df$Hour,
  Minute='NA',
  Value=(as.numeric(df$Bar.p.) * 23.22 + as.numeric(df$Bar.l.) * 1.935 ) * 1.3322
)


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
          lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]], metaHead = "orig_p=Castillian_inch_&_lines",
          link=meta[["Link"]], units="unknown", stat="point",
          meta="orig_p=Castillian_inch_&_lines", keep_na = F)


# Ylitornio ---------------------------------------------------------------
df <- read.csv('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/Ylitornio_all.csv',
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
  Lat = "unknown",
  Lon = "unknown",
  Alt = "unknown",
  Source = "Helama, S., Holopainen, J., Timonen, M., Ogurtsov, M. G., Lindholm, M., Meriläinen, J., Eronen, M. (2004). Comparison of living-tree and subfossil ringwidths with summer temperatures from 18th, 19th and 20th centuries in Northern Finland. Dendrochronologia 21/3, 147 - 154.",
  meta = "observer=Johan Portin | orig_ta=C | orig_p=unknown(similar to Turku Academy measurements)"
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

