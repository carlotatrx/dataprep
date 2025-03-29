## Calculate air temperature and pressure daily means using ERA5-Land (or, in alternative, a MeteoCH station)
## to adjust for the daily cycle (different daily cycle for each month)
## Observations with missing time are ignored otherwise.
## For pressure, I simply take the average of the available observations when one or
## more times are unknown, with no correction for daily cycle.

rm(list=ls())
library(dataresqc)
library(geosphere)
library(ncdf4)
library(dplyr)
library(tidyr)
source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/write_sef_f.R')

era5_t <- "/home/ccorbella/scratch2_symboliclink/files/ERA5/ERA5_1950-1958_t2m.nc"
era5_p <- "/home/ccorbella/scratch2_symboliclink/files/ERA5/ERA5_1950-1958_prmsl.nc"

outpath <- "/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/"


# temperature -------------------------------------------------------------

## Temperature
# List all the temperature files that need daily adjustment in the current directory
ta_sub_files <- list.files('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed',
                           pattern="ta_subdaily.*\\.tsv$", full.names=T, recursive=F)

nct <- nc_open(era5_t)
lons <- ncvar_get(nct, "lon")
lons[lons>180] <- lons[lons>180] - 360
lats <- ncvar_get(nct, "lat")
# create coordinate df with all possible grid points
coords <- data.frame(lon=rep(lons,each=length(lats)), lat=rep(lats,times=length(lons)))

for (f in ta_sub_files) {
  #  x <- read.table(file = f, skip = 12, header = TRUE, fill = TRUE,
  #                  sep = "\t", stringsAsFactors = FALSE, quote = "")
  x <- read_sef(f, all=TRUE)
  x <- x[which(!is.na(x$Value)), ]
  ko <- grep("qc", x$Meta)
  if (length(ko) > 0) x <- x[-ko, ] # Remove flagged values
  meta <- read_meta(f)
  
  x$Time <- x$Hour + x$Minute/60   # convert time to decimal hours
  # aggregate available times per day
  ntimes <- aggregate(x$Time, list(x$Year,x$Month,x$Day), function(x) sum(!is.na(x)))
  ntimes <- ntimes[order(ntimes[,1],ntimes[,2],ntimes[,3]), ]
  
  out <- data.frame(Year=ntimes$Group.1,
                    Month=ntimes$Group.2,
                    Day=ntimes$Group.3,
                    Hour=rep(24, nrow(ntimes)),
                    Minute=rep(0, nrow(ntimes)),
                    Value=rep(NA, nrow(ntimes)))
  
  if (max(ntimes$x) > 0) {    
    x$Value[which(is.na(x$Time))] <- NA
    ## Select closest grid point
    xlon <- as.numeric(meta["lon"])
    xlat <- as.numeric(meta["lat"])
    i_grid <- which.min(distGeo(c(xlon,xlat), coords))
    i_lon <- which(lons == coords$lon[i_grid])
    i_lat <- which(lats == coords$lat[i_grid])
    
    ## Get daily cycle
    dc <- ncvar_get(nct, "t2m", c(i_lon, i_lat,1 ), c(1,1,-1))
    
    ## Convert to anomalies from the daily mean
    for (m in 1:12) dc[((m-1)*24+1):(m*24)] <- 
      dc[((m-1)*24+1):(m*24)] - mean(dc[((m-1)*24+1):(m*24)]) # indices corresponding to that month - monthly mean
    dc <- append(dc, dc[1]) # Add first value to end to ensure smooth interpolation in seasonal models
    # now dc contains hourly ta anomalies for each month
    
    ## Calculate daily means
    xmean <- aggregate(x$Value, list(x$Year,x$Month,x$Day), mean, na.rm=TRUE)
    xmean <- xmean[order(xmean$Group.1, xmean$Group.2, xmean$Group.3), ]
    # now xmean contains raw daily mean ta for each day
    
    ## Apply correction for daily cycle
    for (i_day in 1:nrow(xmean)) {                            # for each day in xmean
      xtimes <- x$Time[which(x$Year==xmean$Group.1[i_day] &   # identify obs times for each day
                             x$Month==xmean$Group.2[i_day] & 
                             x$Day==xmean$Group.3[i_day])]
      xtimes <- xtimes[!is.na(xtimes)]
      if (length(xtimes) > 0) {                               # as long as there are obs
        m <- as.integer(xmean$Group.2[i_day])                 # extract month number for current day
        dc_day <- dc[((m-1)*24+1):(m*24)]                     # extract 24-hourly vals for month m
        dc_day <- append(dc_day, dc_day[1])                   # append last value for circualrity
        corrections <- dc_day[as.integer(xtimes)+1] * (1-xtimes+as.integer(xtimes)) + # contribution from lower hour
                       dc_day[as.integer(xtimes)+2] * (xtimes-as.integer(xtimes))     # contribution from upper hour
        
        # this weighted average with (1-xtimes+as.integer(xtimes)) for the lower hour and (xtimes-as.integer(xtimes))
        # is only relevant if our hours are not integers. Otherwise, the lower hour will always have 1 as the weight.
        out$Value[i_day] <- round(xmean$x[i_day] - mean(corrections), 1)
        # finally, take raw daily mean (xmean$x[i_day]) for the day, subract the average of the interpolated
        # daily corrections, and round result to 1 decimal place.
      }
    }
  }
  # Write daily file
  # define meta column
  if (nchar(x$Meta[1]) > 0) {
    meta_col <- paste(x$Meta[1], "daily_cycle=T", sep=" | ")
  } else {
    meta_col <- "daily_cycle=T"
  }
  # Update meta["metaHead"]
  if (nchar(meta["meta"]) > 0) {
    meta["meta"] <- paste(meta["meta"], "daily means adjusted to ERA5", sep=" | ")
  } else {
    meta["meta"] <- "daily means adjusted to ERA5"
  }  
  
  #  out <- out[which(!is.na(out$Value)), ]
  if (nrow(out) > 0) {
    write_sef_f(out, outpath, 
                outfile=sub("subdaily", "daily", basename(f)), # replace "subdaily" with "daily" in the original filename
                meta["var"], meta["id"], meta["name"], 
                meta["lat"], meta["lon"], meta["alt"], meta["source"], meta["link"], 
                meta["units"], "mean", period="day", note="daily", keep_na=F,
                meta=meta_col,
                metaHead=meta["meta"])
  }
  # move original subdaily files to not_to_use/ directory
  file.rename(f, paste0(outpath,"not_to_use/", basename(f)))
}
nc_close(nct)


# pressure ----------------------------------------------------------------
## Pressure
p_sub_files <- list.files('/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed',
                          pattern="p_subdaily.*\\.tsv$", full.names=T, recursive=F)
ncp <- nc_open(era5_p)
lons <- ncvar_get(ncp, "longitude")
lons[lons>180] <- lons[lons>180] - 360
lats <- ncvar_get(ncp, "latitude")
coords <- data.frame(lon=rep(lons,each=length(lats)), lat=rep(lats,times=length(lons)))

for (f in p_sub_files) {
  x <- read_sef(f, all=TRUE)
  x <- x[which(!is.na(x$Value)), ]
  ko <- grep("qc", x$Meta)
  if (length(ko) > 0) x <- x[-ko, ] # Remove flagged values
  meta <- read_meta(f)
  
  x$Time <- x$Hour + x$Minute/60
  ntimes <- aggregate(x$Time, list(x$Year,x$Month,x$Day), function(x) sum(!is.na(x)))
  ntimes <- ntimes[order(ntimes[,1],ntimes[,2],ntimes[,3]), ]
  out <- data.frame(Year=ntimes$Group.1,
                    Month=ntimes$Group.2,
                    Day=ntimes$Group.3,
                    Hour=rep(24, nrow(ntimes)),
                    Minute=rep(0, nrow(ntimes)),
                    Value=rep(NA, nrow(ntimes)))
  if (max(ntimes$x) > 0) {
    ## Select closest grid point
    xlon <- as.numeric(meta["lon"])
    xlat <- as.numeric(meta["lat"])
    i_grid <- which.min(distGeo(c(xlon,xlat), coords))
    i_lon <- which(lons == coords$lon[i_grid])
    i_lat <- which(lats == coords$lat[i_grid])
    ## Get daily cycle
    dc <- ncvar_get(ncp, "msl", c(i_lon,i_lat,1), c(1,1,-1)) / 100
    ## Convert to anomalies from the daily mean
    for (m in 1:12) dc[((m-1)*24+1):(m*24)] <- 
      dc[((m-1)*24+1):(m*24)] - mean(dc[((m-1)*24+1):(m*24)]) 
    dc <- append(dc, dc[1])
    ## Calculate daily means
    xmean <- aggregate(x$Value, list(x$Year,x$Month,x$Day), mean, na.rm=TRUE)
    xmean <- xmean[order(xmean$Group.1, xmean$Group.2, xmean$Group.3), ]
    
    ## Apply correction for daily cycle (if all times are known)
    for (i_day in 1:nrow(xmean)) {
      xtimes <- x$Time[which(x$Year==xmean$Group.1[i_day] & 
                               x$Month==xmean$Group.2[i_day] & 
                               x$Day==xmean$Group.3[i_day])]
      if (sum(is.na(xtimes)) == 0) {
        m <- as.integer(xmean$Group.2[i_day])
        dc_day <- dc[((m-1)*24+1):(m*24)]
        dc_day <- append(dc_day, dc_day[1])
        corrections <- dc_day[as.integer(xtimes)+1] * (1-xtimes+as.integer(xtimes)) + 
          dc_day[as.integer(xtimes)+2] * (xtimes-as.integer(xtimes))
        out$Value[i_day] <- round(xmean$x[i_day] - mean(corrections), 1)
      } else {
        out$Value[i_day] <- round(xmean$x[i_day], 1)
      }
    }
  }
  
  # Write daily file
  if (nchar(x$Meta[1]) > 0) {
    meta_col <- paste(x$Meta[1], "daily_cycle=T", sep=" | ")
  } else {
    meta_col <- "daily_cycle=T"
  }

  # Update meta["metaHead"]
  if (nchar(meta["meta"]) > 0) {
    meta["meta"] <- paste(meta["meta"], "daily means adjusted to ERA5", sep=" | ")
  } else {
    meta["meta"] <- "daily means adjusted to ERA5"
  }  
  
  #  out <- out[which(!is.na(out$Value)), ]
  if (nrow(out) > 0) {
    write_sef_f(out, outpath, 
                outfile=sub("subdaily", "daily", basename(f)), # replace "subdaily" with "daily" in the original filename
                meta["var"], meta["id"], meta["name"], 
                meta["lat"], meta["lon"], meta["alt"], meta["source"], meta["link"], 
                meta["units"], "mean", period="day", note="daily", keep_na=F,
                meta=meta_col,
                metaHead=meta["meta"])
  }
  # move original subdaily files to not_to_use/ directory
  file.rename(f, paste0(outpath,"not_to_use/", basename(f)))
}
nc_close(ncp)
