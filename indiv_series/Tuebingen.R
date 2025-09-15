rm(list=ls())
library(dataresqc)
library(XLConnect)
library(dplyr)
library(tidyr)
source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')

name <- "Tuebingen"
cod <- "Europe_Tuebingen_1"

file <- paste0('/scratch3/PALAEO-RA/DataRescue/Projects/Europe/3_Revised/',name,'/Europe_T3_DE_',name,'_1691-1694_subdaily.xlsx')

lat	<- 48.5216364
lon	<- 9.0576448
alt	<- NA
metaHead <- "Observer=Camerarius"


fill_val <- function(x){
  
  ind <- which(!is.na(x))
  #ind <- switch((ind[length(ind)]==length(x))+1, c(ind),ind)
  for(i in 1:(length(ind))){
    if(i==length(ind)){
      x[ind[i]:length(x)] <- x[ind[i]]
    } else {
      x[ind[i]:(ind[i+1]-1)] <- x[ind[i]]
    }
  }
  return(x)
}
raw <- readWorksheetFromFile(file, sheet=1, startRow=8, header=FALSE, colTypes="character")



meta_frac <- function(g, b1, b2) {
  if (is.na(b1) | is.na(b2)) {
    as.character(g)  # just return the integer part
  } else {
    paste0(g, ".", b1, "/", b2)
  }
}

p_frac <- function(g, b1, b2) {
  if (is.na(b1) | is.na(b2)) {
    as.integer(g)  # just return the integer part
  } else {
    round(as.numeric(g) + as.numeric(b1)/as.numeric(b2)/10 ,2)
  }
}

Year <- as.integer(fill_val(raw[,1]))
Month <- as.integer(fill_val(raw[,2]))
Day <- as.integer(fill_val(raw[,3]))
hours <- as.integer(substr(raw[,4],12,13))
minutes <- as.integer(substr(raw[,4],15,16))
meta_time <- raw[,4]
j <- which(nchar(meta_time)>1)
meta_time[j] <- substr(meta_time[j], 12, 16)
meta_time <- paste0("orig.time=", meta_time)
p.raw <- as.integer(raw[,8])
p.frac1 <- as.integer(raw[,9])
p.frac2 <- as.integer(raw[,10])
meta_p <- paste0('orig_p=', mapply(meta_frac, p.raw, p.frac1, p.frac2))
p <- mapply(p_frac, p.raw, p.frac1, p.frac2)

time.offset <- as.numeric(lon)*12/180

## Write SEF files
write_sef_f(data.frame(Year, Month, Day, hours, minutes,p),
            outfile=paste0(name, "_p_subdaily.tsv"),
            outpath="/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/",
            variable = "p",
            cod = cod,
            nam = name,
            lat = lat,
            lon = lon,
            alt = alt,
            sou = "PALAEO-RA",
            units = "unknown",
            stat = "point",
            meta = paste0(meta_time, " | ", meta_p),
            metaHead = metaHead,
            period = 0,
            time_offset = time.offset,
            keep_na = F)
