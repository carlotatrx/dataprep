## Only daily temperature formatted

library(dataresqc)

lat <- 47.75
lon <- 7.34
ele <- 240

## SUB-DAILY

source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')
sef_test_path <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/sef_tests/'
outpath_preprocessed <- "/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed"

## Read data
template <- read.table("/scratch3/PALAEO-RA/DataRescue/Incoming/DigiHom/Mulhouse/mulhouse_1794_1812_DC.txt", header=TRUE, fill=TRUE, quote="",
                       sep="\t", stringsAsFactors=FALSE, na.strings="-9999", fileEncoding="latin1")
                       
## Initialize data frame
out <- data.frame(year = rep(template[,2],each=3), 
                  month = rep(template[,3],each=3), 
                  day = rep(template[,4],each=3), 
                  hour = rep(NA,3*nrow(template)), 
                  minute = 0)
                  
## Temperature
ta.C <- cbind(template$Temp_1stObs_Cels, template$Temp_2ndObs_Cels, template$Temp_3rdObs_Cels)
ta.orig <- cbind(template$Temp_1stObs_Orig, template$Temp_2ndObs_Orig, template$Temp_3rdObs_Orig)
ta.orig.flatten <- c(t(ta.orig))  # flatten column-wise
tod <- rep(c("matin", "apres midi", "soir"), nrow(template))

out$ta <- round(c(t(ta.C)), 1)
out$meta.ta <- paste0("orig_time=", tod,
                   " | orig_ta=", ta.orig.flatten, "R")
out$hour <- rep(c(8, 14, 20))


## Pressure
p.zoll <- cbind(template$Baro_1stObs_Zoll, template$Baro_2ndObs_Zoll, template$Baro_3rdObs_Zoll)
p.line <- cbind(template$Baro_1stObs_Linien, template$Baro_2ndObs_Linien, template$Baro_3rdObs_Linien)
p.zoll.flatten <- c(t(p.zoll))
p.line.flatten <- c(t(p.line))


temp <- cbind(template$Baro_1stObs_Zoll_Decimal, template$Baro_2ndObs_Zoll_Decimal, template$Baro_3rdObs_Zoll_Decimal)
out$p <- convert_pressure(c(t(temp)), 27.07, lat, ele)
out$p <- round(out$p, 1)
out$meta.p <- paste0("orig.time=", tod,
                     " | orig_p=", p.zoll.flatten, ".", p.line.flatten, "in")

write_sef_f(
  Data=out[,c("year","month","day","hour","minute","ta")], outfile="Mulhouse_ta_subdaily.tsv",
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod="Mulhouse",
  variable="ta",
  nam="Mulhouse",
  lat=lat, lon=lon, alt=ele,
  meta=out$meta.ta,
  units="C", sou="DigiHom", stat="point",keep_na = F
)

write_sef_f(
  Data=out[,c("year","month","day","hour","minute","p")], outfile="Mulhouse_p_subdaily.tsv",
  outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
  cod="Mulhouse",
  variable="p",
  nam="Mulhouse",
  lat=lat, lon=lon, alt=ele,
  meta=out$meta.p, metaHead = "PGC=Y | PTC=N",
  units="hPa", sou="DigiHom", stat="point",keep_na = F
)
