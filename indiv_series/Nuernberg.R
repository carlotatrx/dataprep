## SCRIPT EXTRACTED FROM PROJECTS/SCRIPTS/

## NOTES:
## Temperature (unknown Fahrenheit's thermometer) was not formatted
## Monthly pressure data were not formatted as they do not fit with the daily data
## and would thus introduce an inhomogeneity (only 2 months are missing from the daily)

library(dataresqc)
library(XLConnect)

## Create station code
name_folder <- strsplit(getwd(), "/")[[1]][8]
name_group <- strsplit(getwd(), "/")[[1]][6]
station_code <- paste(name_group, name_folder, sep = "_")

## Create output directory (if necessary)
if (!name_folder %in% list.dirs("../../4_Formatted/",full.names=FALSE)) {
  dir.create(paste0("../../4_Formatted/", name_folder))
}

## Read inventory
inventory <- readWorksheetFromFile("/scratch3/PALAEO-RA/DataRescue/Inventory_PALAEO-RA_digi.xls",
                                   sheet=1, colTypes="character")

## Define function to fill missing values with the value on the previous row
fill_variable <- function(timeseries, strings) {
  for (i in 2:length(timeseries)) {
    if (timeseries[i] %in% strings) {
      if(!timeseries[i-1] %in% strings) {
        timeseries[i] <- timeseries[i-1]
      }
    }
  }
  return(timeseries)
}


## MONTHLY
## Read data
# template <- readWorksheetFromFile("../../3_Revised/Nuernberg/Europe_T3_DE_Nuernberg_1734-1735_monthly-pres.xlsx",
#                                   sheet=1, startRow=8, header=FALSE, colTypes="character")
# metadata <- inventory[which(inventory$Code == station_code & inventory$Resolution=="monthly"), ]
# template <- template[which(!is.na(template[,2])), ]
# template[,1] <- fill_variable(template[,1], NA)
# 
# ## Arrange output
# pres <- convert_pressure(as.numeric(template[,4]), 25.4,
#                          as.numeric(metadata$Latitude), as.numeric(metadata$Altitude))
# out <- data.frame(Year = as.integer(template[,1]),
#                   Month = as.integer(template[,2]),
#                   Day = rep("", nrow(template)),
#                   Hour = rep("", nrow(template)),
#                   Minute = rep("", nrow(template)),
#                   Value = round(pres, 1),
#                   Meta = paste0("orig=", template[,4], "in"),
#                   stringsAsFactors = FALSE)
# 
# ## Write SEF files
# write_sef(out[, 1:6],
#           outpath = paste0("../../4_Formatted/", name_folder),
#           variable = "p",
#           cod = metadata$Code,
#           nam = metadata$Station,
#           lat = metadata$Latitude,
#           lon = metadata$Longitude,
#           alt = metadata$Altitude,
#           sou = "PALAEO-RA",
#           link = metadata$C3S.Link,
#           units = "hPa",
#           stat = "mean",
#           meta = out$Meta,
#           metaHead = paste0("Observer=", metadata$Observer, " | PTC=N | PGC=Y"),
#           period = "month",
#           note = "monthly",
#           keep_na = TRUE)
# 
# file.rename(list.files(paste0("../../4_Formatted/", name_folder), full.names=TRUE), 
#             gsub(" ", "", list.files(paste0("../../4_Formatted/", name_folder), full.names=TRUE)))


## SUBDAILY
## Loop over files
files <- list.files(paste0("../../3_Revised/", name_folder), pattern="subdaily", full.names=TRUE)
out <- list()
for (i in 1:length(files)) {
  
  ## Read data
  template <- readWorksheetFromFile(files[i], sheet=1, startRow=8, header=FALSE, colTypes="character")
  metadata <- inventory[which(inventory$Code == station_code & inventory$Resolution=="subdaily"), ]
  
  ## Pressure
  j <- ifelse(i==1, 7, 4)
  template[,j] <- fill_variable(template[,j], NA)
  pres <- as.numeric(template[,j]) + as.numeric(template[,j+1])/10 + as.numeric(template[,j+2])/100
  pres <- convert_pressure(pres, 25.4, as.numeric(metadata$Latitude), as.numeric(metadata$Altitude))
  
  ## Wind direction
  directions <- c("N", "NNO", "NO", "ONO", "O", "OSO", "SO", "SSO", "S",
                  "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW")
  wdir <- 22.5 * (match(toupper(template[,16]), directions) - 1)
  
  ## Create output data frame
  out[[i]] <- data.frame(Year = as.integer(template[,1]),
                         Month = as.integer(template[,2]),
                         Day = as.integer(template[,3]),
                         Hour = rep(NA, nrow(template)),
                         Minute = rep(NA, nrow(template)),
                         p = round(pres, 1),
                         p_meta = paste0("orig=", template[,j], ".", template[,j+1], template[,j+2], "in"),
                         dd = round(wdir, 0),
                         dd_meta = paste0("orig=", template[,16]),
                         stringsAsFactors = FALSE)
}

## Merge data frames
out <- do.call(rbind, out)

## Write SEF files
write_sef(out[, 1:6],
          outpath = paste0("../../4_Formatted/", name_folder),
          variable = "p",
          cod = metadata$Code,
          nam = metadata$Station,
          lat = metadata$Latitude,
          lon = metadata$Longitude,
          alt = metadata$Altitude,
          sou = "PALAEO-RA",
          link = metadata$C3S.Link,
          units = "hPa",
          stat = "point",
          meta = out$p_meta,
          metaHead = paste0("Observer=", metadata$Observer, " | PTC=N | PGC=Y"),
          period = 0,
          keep_na = TRUE)

write_sef(out[, c(1:5,8)],
          outpath = paste0("../../4_Formatted/", name_folder),
          variable = "dd",
          cod = metadata$Code,
          nam = metadata$Station,
          lat = metadata$Latitude,
          lon = metadata$Longitude,
          alt = metadata$Altitude,
          sou = "PALAEO-RA",
          link = metadata$C3S.Link,
          units = "degree",
          stat = "point",
          meta = out$dd_meta,
          metaHead = paste0("Observer=", metadata$Observer),
          period = 0,
          keep_na = TRUE)

## Check SEF files
for (f in list.files(paste0("../../4_Formatted/", name_folder), full.names=TRUE)) check_sef(f)