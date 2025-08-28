library(dataresqc)
source("/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R")

indir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/'
outdir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/'

# Ukraine

station.names <- c('Dnipro','Kamyanets','Kharkiv','Kherson','Kyiv','Lugansk','Poltava','Odesa')

files <- as.vector(outer(station.names, c('ta', 'p'),
                   function(station, var) paste0(station,"_",var,"_subdaily_qc.tsv")))

lapply(files, function(file) write_sef_time_offset(paste0(indir,file), outdir))
       