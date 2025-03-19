library(dataresqc)
source("/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/write_sef_f.R")

outdir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/'


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
  Minute='0',
  Value=as.numeric(df$Term) * 1.25
)

df.p.Valencia <- data.frame(
  Year=df$Year,
  Month=df$Month,
  Day=df$Day,
  Hour=df$Hour,
  Minute='0',
  Value=(as.numeric(df$Bar.p.) * 23.22 + as.numeric(df$Bar.l.) * 1.935 ) * 1.3322
)


write_sef(Data=df.ta.Valencia,
          outpath=outdir,
          cod=meta[["ID"]],
          variable="ta",
          nam=meta[["Name"]],
          lat=meta[["Lat"]],
          lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]], metaHead = "orig_ta=Reaumur",
          link=meta[["Link"]], units="C", stat="point",
          meta="orig_ta=Reaumur", keep_na = F)

write_sef(Data=df.p.Valencia,
          outpath=outdir,
          cod=meta[["ID"]],
          variable="p",
          nam=meta[["Name"]],
          lat=meta[["Lat"]],
          lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]], metaHead = "orig_p=Castillian_inch_&_lines",
          link=meta[["Link"]], units="C", stat="point",
          meta="orig_p=Castillian_inch_&_lines", keep_na = F)


# Ylitornio ---------------------------------------------------------------
df <- read.csv('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/Ylitornio_all.csv',
               header=T)

df.ta.Yli <- data.frame(
  Year=df$Year,
  Month=df$Month,
  Day=df$Day,
  Hour='9',
  Minute='0',
  Value=df$TMP2m
)

df.p.Yli <- data.frame(
  Year=df$Year,
  Month=df$Month,
  Day=df$Day,
  Hour='9',
  Minute='0',
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

write_sef_f(Data=df.ta.Yli, outfile="Yli_ta.tsv",
          outpath=outdir,
          cod=meta[["ID"]],
          variable="ta",
          nam=meta[["Name"]],
          lat=meta[["Lat"]],
          lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]], metaHead = meta[['meta']],
          units="C", stat="point",
          meta="orig_ta=C", keep_na = F
)

write_sef_f(Data=df.ta.Yli, outfile="Yli_p.tsv",
            outpath=outdir,
            cod=meta[["ID"]],
            variable="p",
            nam=meta[["Name"]],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]], metaHead = meta[['meta']],
            units="unknown", stat="point", keep_na = F
)

