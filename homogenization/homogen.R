### Homogenize all pressure data

# 1. create input for homogenization files
# 2. do homogenization using RHclim

library(maps)
library(RColorBrewer)
library(dataresqc)
library(ncdf4)
library(terra)

rm(list=ls())

source("/scratch3/nimfeld/wear/swiss_arm/code/data_prep/RHtests-master/RHtests-master/V4_files/RHtestsV4_20190301.r")
source("/scratch3/nimfeld/wear/swiss_arm/code/data_prep/scripts/000_helpfun.R")

file.name <- '/home/ccorbella/scratch2_symboliclink/files/1807_USBstick/KNMI_KNMI-42_Zwanenburg_17380101-18601231_p_daily.tsv'
x <- read_sef(file.name)
x$Date <- as.Date(with(x, paste(Year, Month, Day, sep = "-")), format = "%Y-%m-%d")
meta <- read_meta(file.name)

wd <- "~/scratch2_symboliclink/code/KF_assimilation/dataprep/homogenization/"

## calculate monthly averages for pressure data from meteoswiss
# Create a complete monthly time sequence
all_months <- expand.grid(Year = min(x$Year):max(x$Year), Month = 1:12)
all_months <- all_months[order(all_months$Year, all_months$Month), ]

# Calculate monthly means where data coverage is sufficient
valid_means <- aggregate(x$Value, by = list(Year = x$Year, Month = x$Month), mean_Xpercna, na = 0.2)

# Merge the complete calendar with your valid means
mondat <- merge(all_months, valid_means, by = c("Year", "Month"), all.x = TRUE)



datf <- data.frame(year = mondat$Year,
                   mon = mondat$Month,
                   day = 0,
                   value= mondat$x)

write.table(datf,file=paste0(wd,"hom_input/",meta['name'],".dat"),quote=F,row.names = F,col.names = F)

######################################################
###### create input reference series
######################################################

modera_file <- "/mnt/climstor/ERC_PALAEO/ModE-RA/outdata/lowres_20mem_Set_1420-3_1850-1/abs/ensstat/by_var/mon/ModE-RA_lowres_20mem_Set_1420-3_1850-1_ensmean_slp_abs_1421-2008_mon.nc"
modera_slp <- terra::rast(modera_file,"slp")
modera_time <- time(modera_slp)

mlat <- as.numeric(meta['lat'])
mlon <- as.numeric(meta['lon'])

stat <- cellFromXY(modera_slp, cbind(mlon,mlat))
slpch <- extract(modera_slp, y = stat)

tind <- modera_time > "1762-12-31"
datf <- data.frame(year=lubridate::year(modera_time),
                   mon= lubridate::month(modera_time),
                   day=rep(00,length(modera_time)),
                   value= unlist(slpch)/100)[tind,]

write.table(datf,file=paste0(wd,"hom_input/ModERA_slp.dat"),quote=F,row.names = F,col.names = F)


###############################################
############## Homogenization all stations
##############################################

plot <- F
hmon <- 6

## read the hom input files
files <- list.files(path=paste0(wd,"/hom_input"))
files <- files[!grepl("ModERA",files)]
id <- gsub(".dat","",files)

#files_ta_p <- files[grepl("ta.dat|p.dat",files)] 

## check homogenity of reference series using FindU
if (!"ModERA_slp" %in% list.dirs(paste0(wd,"/hom_output"),full.names=FALSE)) dir.create(path = paste0(wd,"/hom_output/ModERA_slp/"))
FindU("hom_input/ModERA_slp.dat",MissingValueCode = "NA", p.lev=0.95,Iadj=10000,Mq=10,Ny4a=0,
      output = paste0("hom_output/ModERA_slp/","ModERA_slp"))

selref <- c()
ff <- 1 # now only 1 file
print(id[ff])
vals <- read.table(paste0(wd,"/hom_input/",id[ff],".dat"))[,4]
init <-gsub("_ta.*|_p.*","",id[ff])
ref <- gsub(init,"",id[ff])
selref[ff] <- gsub("_.*","",substr(ref,2, nchar(ref)))

if (sum(!is.na(vals))>12*hmon){
  if (!id[ff] %in% list.dirs(paste0(wd,"/hom_output"),full.names=FALSE)) {dir.create(path = paste0(wd,"/hom_output/",id[ff],"/"))}
  FindU.wRef(Bseries=paste0(wd,"/hom_input/",files[ff]),MissingValueCode = "NA", p.lev=0.95,Iadj=10000,Mq=10,Ny4a=0,
             output = paste0(wd,"/hom_output/",id[ff],"/",id[ff]),
             Rseries=paste0(wd,"/hom_input/ModERA_slp",".dat"))
  }


# ## additional change points
# for(ff in 1:(length(files))){
#   #for(ff in 71){
#   print(id[ff])
#   
#   if (sum(!is.na(vals))>12*hmon){
#     FindUD.wRef(Bseries=paste0(wd,"/hom_input/",files[ff]),MissingValueCode = "NA", p.lev=0.95,Iadj=10000,Mq=10,Ny4a=0,
#                 InCs=paste0(wd,"/hom_output/",id[ff],"/",id[ff],"_1Cs.txt"),
#               output = paste0(wd,"/hom_output/",id[ff],"/",id[ff]),
#               Rseries=paste0(wd,"/hom_input/ModERA_slp",".dat"))
# }
# }

#########################################################################################################
## only keep the significant change points and re-estimate the significance and magnitude of change points
#########################################################################################################

for (ff in 1:(length(files))){
  print(id[ff])  
  vals <- read.table(paste0(wd,"/hom_input/",id[ff],".dat"))[,4]
  
  if (sum(!is.na(vals))>12*hmon){
    
    StepSize.wRef(Bseries=paste0(wd,"/hom_input/",files[ff]),
                  Rseries=paste0(wd,"/hom_input/ModERA_slp.dat"),
                  InCs = paste0(wd,"/hom_output/",id[ff],"/",id[ff],"_mCs.txt"),
                  MissingValueCode = "NA", p.lev=0.95,Iadj=10000,Mq=10,Ny4a=0,output = paste0(wd,"/hom_output/",id[ff],"/",id[ff]))
  }
}


#########################################################################################################
## create final homogenized with adjustements estimated above 
#########################################################################################################
sefdata <- data.frame(Year=x$Year,
                      Month=x$Month,
                      Day=x$Day)

print(id[ff])  
vals <- read.table(paste0(wd,"/hom_input/",id[ff],".dat"))[,4]

if (sum(!is.na(vals))>12*hmon){
  
  # orig qced data
  ts <- x$Value
  
  # read hom vals
  Fdat <- read.table(paste0(wd,"/hom_output/",id[ff],"/",id[ff],"_F.dat"))
  qmt <- 5 #ifelse(all(is.na(Fdat[,11])),5,11) -> don't use the qm based estimates, because they are estimated based on the EKF reference
  ## set Fdat values to NA, if they are followed by a large number of NAs...
  #fdathelp <-  ifelse(!is.na(Fdat[,qmt]),1,0)
  
  homadj_monthly <- Fdat[,qmt] - Fdat[,3] 
  # create daily adjustments using a spline?
  mimo <- Fdat[,2]+15
  mohelp <- as.Date(paste0(substr(mimo,1,4),"-",substr(mimo,5,6),"-",substr(mimo,7,8)))
  lmo <- lubridate::days_in_month(as.Date(paste0(substr(Fdat[nrow(Fdat),2],1,4),"-",substr(Fdat[nrow(Fdat),2],5,6),"-01")))
  dhelp <- seq(as.Date(paste0(substr(Fdat[1,2],1,4),"-",substr(Fdat[1,2],5,6),"-01")),
               as.Date(paste0(substr(Fdat[nrow(Fdat),2],1,4),"-",substr(Fdat[nrow(Fdat),2],5,6),"-",lmo)),
               by="day")
  xout <- which(dhelp %in% x$Date)
  which(!x$Date %in% dhelp[xout])
  stopifnot(length(xout)==length(sefdata$Day))
  
  ## don't use spline. because it makes weird step to next period...
  homadj_daily <- homadj_monthly[match(substr(dhelp[xout],1,7),substr(mohelp,1,7))]
  # add hom adjust
  sefdata[,'Value'] <- ts + homadj_daily
  #pipe <- ifelse(sefdata[,"Meta"]=="","","|")
  sefdata[,"Meta"] <- paste0("homadj=",round(homadj_daily,2))
  
  # check for complete series
  sefdata <- check_dates(sefdata)
  
  write_sef(Data = sefdata,
            outpath = wd,
            variable = 'p',
            cod = meta["id"],
            nam = meta["name"],
            lat = meta["lat"],
            lon = meta["lon"],
            alt = meta["alt"],
            link = meta["link"],
            sou = meta["source"],
            units = meta["units"],
            stat = "mean",
            metaHead = meta["meta"],
            meta=sefdata$Meta,
            period = "day",
            time_offset =0,
            keep_na = TRUE,
            note="hom")

}



#save(TOThom,file = "/scratch3/nimfeld/wear/swiss_grids/data/analogue_dates/analogue_stationdata_spatpress_wtnew_pressqced_presshom.RData")
#load("/scratch3/nimfeld/wear/swiss_grids/data/analogue_dates/analogue_stationdata_spatpress_wtnew_pressqced_presshom.RData")

### ex-change the values of TOT sp/spm with homogenized value
colnames(TOThom)[grepl("spm|sp",colnames(TOThom))]

TOThom$NSG_sp <- TOThom$ALT_p - TOThom$LUG_p
TOThom$EWG_sp <- TOThom$GVE_p - TOThom$SMA_p
TOThom$LUG_sp  <- TOThom[,"LUG_p"] - TOThom[,"ALT_p"] 
TOThom["MIL_spm"]  <- TOThom[,"MIL_p"] - TOThom[,"SMA_p"] 
TOThom["TOR_spm"]  <- TOThom[,"TOR_p"] - TOThom[,"SMA_p"] 

stns <- c("ALT_spm","BAS_spm","BER_spm","CHU_spm","GSB_spm","DAV_spm","CDF_spm","LUG_spm","NEU_spm","SHA_spm","SAE_spm","SMA_spm","BUS_spm",
          "GOT_spm","GVE_spm","HOH_spm","KAR_spm","LUZ_spm","STG_spm")
stns_p <- paste0(substr(stns,1,3),"_p")

TOThom[,stns] <- TOThom[,stns_p] - TOThom[,"MIL_p"]

save(TOThom,file = "/scratch3/nimfeld/wear/swiss_grids/data/analogue_dates/analogue_stationdata_spatpress_wtnew_pressqced_presshom.RData")


#########################################################################################################
## create plots of homogenized and real data
#########################################################################################################
if(plot){
  
  for(ff in 1:length(files)){
    
    print(id[ff])  
    png(paste0("plots/",id[ff],"_hom.png"), width = 4000, height = 2000, res=250)
    par(mfrow=c(2,1), mar=c(2,2,1,1))
    tind <- which(!is.na(TOTqc[,id[ff]]))[1]:nrow(TOTqc)
    plot(TOTqc$date[tind], TOTqc[tind,id[ff]], type="l", xlab="",ylab="hPa")
    points(TOTqc$date[tind], TOThom[tind,id[ff]], type="l", col="cornflowerblue")
    mtext(id[ff], line = -1)
    legend("bottom", c("qc","hom"),col=c("black","cornflowerblue"), pch=20, bty="n", ncol = 2)
    plot(TOTqc$date[tind], TOTqc[tind,id[ff]] - TOThom[tind,id[ff]], type="l", xlab="", ylab="Diff")
    dev.off()
    
  }
  
}
