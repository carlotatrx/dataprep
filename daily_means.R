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

era5_t <- "../../reanalysis-era5-land_2m_temperature_daily_cycle_2001-2018.nc"
era5_p <- "../../reanalysis-era5-land_surface_pressure_daily_cycle_2001-2018.nc"

SEF_path <- "/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/Ukraine"
outpath <- "/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/"


## Temperature
files <- list.files(SEF_path, pattern=".tsv", full.names=TRUE)
nct <- nc_open("/home/ccorbella/scratch2_symboliclink/files/station_timeseries_orig/Ukraine/an_sfc_ERA5_1950-01-01.nc")
lons <- ncvar_get(nct, "lon")
lons[lons>180] <- lons[lons>180] - 360
lats <- ncvar_get(nct, "lat")
# create coordinate df with all possible grid points
coords <- data.frame(lon=rep(lons,each=length(lats)), lat=rep(lats,times=length(lons)))

for (f in files) {
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
  #  out <- out[which(!is.na(out$Value)), ]
  if (nrow(out) > 0) {
    write_sef(out, outpath, meta["var"], meta["id"], meta["name"], 
              meta["lat"], meta["lon"], meta["alt"], meta["source"], meta["link"], 
              meta["units"], "mean", period="day", meta["meta"], note="daily", keep_na=TRUE,
              outfile=meta["name"])
  }
}
  

# Temperature dailies
file <- list.files(SEF_path, pattern="ta_daily", full.names=TRUE)
lons <- ncvar_get(nct, "longitude")
lons[lons>180] <- lons[lons>180] - 360
lats <- ncvar_get(nct, "latitude")
coords <- data.frame(lon=rep(lons,each=length(lats)), lat=rep(lats,times=length(lons)))
x <- read_sef(file, all=TRUE)
x <- x[which(!is.na(x$Value)), ]
ko <- grep("qc", x$Meta)
if (length(ko) > 0) x <- x[-ko, ] # Remove flagged values
meta <- read_meta(file)
#get time from other subdaily file
time_x <- read_sef(list.files(SEF_path, pattern="ta.tsv|ta_qc", full.names=TRUE), all=TRUE)
time_x$Time <- time_x$Hour + time_x$Minute/60
ntimes <- data.frame(x[,c("Year","Month","Day")], 3)
time <- data.frame(rep(x$Year,each=3),rep(x$Month,each=3),rep(x$Day,each=3))
if (1784 %in% time_x$Year) {
  time <- data.frame(time, time_x$Time[which(1784==time_x$Year)][1:3])
} else if (1789 %in% time_x$Year) {
  time <- data.frame(time, time_x$Time[which(1789==time_x$Year)][1:3])
} else if (1786 %in% time_x$Year) {
  time <- data.frame(time, time_x$Time[which(1786==time_x$Year)][1:3])
}
colnames(time) <- c("Year","Month","Day","Time")
#initialize dataframe
out <- data.frame(Year=ntimes$Year,
                  Month=ntimes$Month,
                  Day=ntimes$Day,
                  Hour=rep(24, nrow(ntimes)),
                  Minute=rep(0, nrow(ntimes)),
                  Value=rep(NA, nrow(ntimes)))
#get daily values (averages in SMP tables)
xmean <- x[,c("Year","Month","Day","Value")]
## Select closest grid point
xlon <- as.numeric(meta["lon"])
xlat <- as.numeric(meta["lat"])
i_grid <- which.min(distGeo(c(xlon,xlat), coords))
i_lon <- which(lons == coords$lon[i_grid])
i_lat <- which(lats == coords$lat[i_grid])
## Get daily cycle
dc <- ncvar_get(nct, "t2m", c(i_lon,i_lat,1), c(1,1,-1))
## Convert to anomalies from the daily mean
for (m in 1:12) dc[((m-1)*24+1):(m*24)] <- 
  dc[((m-1)*24+1):(m*24)] - mean(dc[((m-1)*24+1):(m*24)]) 
dc <- append(dc, dc[1])
## Apply correction for daily cycle
for (i_day in 1:nrow(xmean)) {
  xtimes <- time$Time[which(time$Year==xmean$Year[i_day] & 
                              time$Month==xmean$Month[i_day] & 
                              time$Day==xmean$Day[i_day])]
  if (length(xtimes) > 0) {
    m <- as.integer(xmean$Month[i_day])
    dc_day <- dc[((m-1)*24+1):(m*24)]
    dc_day <- append(dc_day, dc_day[1])
    corrections <- dc_day[as.integer(xtimes)+1] * (1-xtimes+as.integer(xtimes)) + 
      dc_day[as.integer(xtimes)+2] * (xtimes-as.integer(xtimes))
    out$Value[i_day] <- round(xmean$Value[i_day] - mean(corrections), 1)
  }
}
## Write daily file
if (nrow(out) > 0) {
  write_sef(out, outpath, meta["var"], meta["id"], meta["name"], 
            meta["lat"], meta["lon"], meta["alt"], meta["source"], meta["link"], 
            meta["units"], "mean", period="day", meta["meta"], note="daily", keep_na=TRUE)
}
nc_close(nct)


## Pressure
  files <- list.files(SEF_path, pattern="p.tsv|p_qc", full.names=TRUE)
  ncp <- nc_open(era5_p)
  lons <- ncvar_get(ncp, "longitude")
  lons[lons>180] <- lons[lons>180] - 360
  lats <- ncvar_get(ncp, "latitude")
  coords <- data.frame(lon=rep(lons,each=length(lats)), lat=rep(lats,times=length(lons)))
  for (f in files) {
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
      dc <- ncvar_get(ncp, "sp", c(i_lon,i_lat,1), c(1,1,-1)) / 100
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
    ## Write daily file
    #  out <- out[which(!is.na(out$Value)), ]
    #  if (nrow(out) > 0) {
    write_sef(out, outpath, meta["var"], meta["id"], meta["name"], 
              meta["lat"], meta["lon"], meta["alt"], meta["source"], meta["link"], 
              meta["units"], "mean", period="day", meta["meta"], note="daily", keep_na=TRUE)
    #  }
  }
  nc_close(ncp)