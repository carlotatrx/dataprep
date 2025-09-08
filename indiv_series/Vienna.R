library(dataresqc)
library(XLConnect)

#### functions for sef files...
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


files <- list.files("/scratch3/PALAEO-RA/daily_data/original/Vienna", full.names=T)

i <- 1

template <- readWorksheetFromFile(files[i], sheet=ifelse(i==1,3,1), startRow=8, 
                                  header=FALSE, colTypes="character")

Year <- as.integer(fill_val(template[,1]))
Month <- as.integer(fill_val(template[,2]))
Day <- fill_val(as.integer(template[,3]))

if (i == 1) {
  hours <- rep(NA, nrow(template))
  minutes <- hours
  meta_time <- ""
}

head(template)
