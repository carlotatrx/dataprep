library(dataresqc)
library(readxl)
library(tidyr)
library(dplyr)

source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')

################################################################################
linein2mmHg <- function(inch, l) {
  val = inch + l/12
  return(round(val * 27.07, 2))
}

R2C <- function(val) {
  return(round(val * 1.25, 1))
}

# Augsburg ----------------------------------------------------------------

filedir <- '/scratch3/PALAEO-RA/DataRescue/Projects/Europe/2_Digitized/Augsburg/'

# for file in filedir:

df <- read_excel(paste0(filedir,"Europe_T2_DE_Augsburg_1813_subdaily.xls"), skip = 5)
colnames(df) <- c("Year", "Month", "Day", "p7in", "p7l", "p14in", "p14l", "p21in", "p21l", "tbar7",
                  "tbar14", "tbar21", "ta7", "ta14", "ta21", "notes")

df$notes <- NULL # drop column

# replace all commas with points
df[] <- lapply(df, function (x) {
  if (is.character(x) || is.factor(x)){
    x <- as.character(x)
    x <- gsub(",", ".", x)
    x    
  } else {
    x
  }
})

# replace all "-" with NA
df[] <- lapply(df, function(x){
  if (is.character(x)) {
    x[x == "-"] <- NA
  }
  x
})

# coerce all cols to num
df[] <- lapply(df, function (x) {
  if (is.character(x)) suppressWarnings(as.numeric(x)) else x
})


df <- df %>%
  mutate(
    p7  = linein2mmHg(p7in, p7l),
    p14 = linein2mmHg(p14in, p14l),
    p21 = linein2mmHg(p21in, p21l),
  
    tbar7C  = R2C(tbar7),
    tbar14C = R2C(tbar14),
    tbar21C = R2C(tbar21),
    
    ta7C  = R2C(ta7),
    ta14C = R2C(ta14),
    ta21C = R2C(ta21)
  )
# pivot to long table

df.long <- df %>%
  pivot_longer(
    cols = -c(Year, Month, Day),
    names_to = c("var", "Hour"),
    names_pattern = "^([a-zA-Z0-9_]+?)(\\d+)$",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = var,
    values_from = value
  ) %>%
  mutate(Hour = as.numeric(Hour)) %>%
  arrange(Year, Month, Day, Hour)

df <- df[, c(1,2,3,4,8,7,5,6)] # reorder columns

head(df.long)
