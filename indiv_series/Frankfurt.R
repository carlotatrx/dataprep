rm(list=ls())
library(dataresqc)
library(readxl)
library(dplyr)
library(tidyr)
library(hms)
source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')

file1 <- '/scratch3/PALAEO-RA/DataRescue/BORIS/XLS/Europe/Frankfurt/Europe_T3_DE_Frankfurt_1749_subdaily.xlsx'
file2 <- '/scratch3/PALAEO-RA/DataRescue/BORIS/XLS/Europe/Frankfurt/Europe_T3_DE_Frankfurt_1750-1752_subdaily.xlsx'
file3 <- '/scratch3/PALAEO-RA/DataRescue/BORIS/XLS/Europe/Frankfurt/Europe_T3_DE_Frankfurt_1753-1755_subdaily.xlsx'

raw1 <- read_excel(file1, sheet = 1, skip  = 6)
raw2 <- read_excel(file2, sheet = 1, skip  = 6)
raw3 <- read_excel(file3, sheet = 1, skip = 6)

# Keep only the first 10 columns and rename
df1 <- raw1 %>%
  select(1:9) %>%
  rename(
    Year = 1, Month = 2, Day = 3, Hour = 4, p = 5,
    ta = 6, ta.corrected = 7, direction = 8, force=9
  ) %>%
  # Fill down Year/Month/Day
  mutate(across(where(is.character), ~ na_if(trimws(.x), ""))) %>%
  fill(Year, Month, Day, .direction = "down") %>%
  mutate(
    Year = as.integer(Year),
    Month = as.integer(Month),
    Day = as.integer(Day),
    Hour = round(suppressWarnings(as.numeric(Hour))* 86400/ 3600, 0), # convert stupidity
    p = suppressWarnings(as.numeric(p)),
    p.zoll = as.integer(p%/%10), # first 2 digits
    p.in = as.integer(round((p%%10) * 10)), # everything after decimal * 10
    ta = suppressWarnings(as.numeric(ta)),
    ta.corrected = suppressWarnings(as.numeric(ta.corrected)),
    force = suppressWarnings(as.numeric(force))
  )
head(df1)

# Keep only the first 10 columns and rename
df2 <- raw2 %>%
  select(1:10) %>%
  rename(
    Year = 1, Month = 2, Day = 3, Hour = 4, direction = 5, force=6,
    ta = 7, p.zoll = 8, p.in = 9, hygrom = 10 
  ) %>%
  # Fill down Year/Month/Day
  mutate(across(where(is.character), ~ na_if(trimws(.x), ""))) %>%
  fill(Year, Month, Day, .direction = "down") %>%
  mutate(
    Year = as.integer(Year),
    Month = as.integer(Month),
    Day = as.integer(Day),
    Hour = suppressWarnings(as.integer(Hour)),
    p.zoll = suppressWarnings(as.numeric(p.zoll)),
    p.in = suppressWarnings(as.numeric(p.in)),
    ta = suppressWarnings(as.numeric(ta)),
    force = suppressWarnings(as.numeric(force)),
    direction = ifelse(direction=="-", NA, direction)
  )
head(df2)

df3 <- raw3 %>%
  select(1:11) %>%
  rename(
    Year = 1, Month = 2, Day = 3, Hour = 4, direction = 5, force=6,
    ta = 7, p.zoll = 8, p.in = 9, p.zoll.min = 10, p.in.min = 11
  ) %>%
  # Fill down Year/Month/Day
  mutate(across(where(is.character), ~ na_if(trimws(.x), ""))) %>%
  fill(Year, Month, Day, p.zoll, p.zoll.min, .direction = "down") %>%
  mutate(
    Year = as.integer(Year),
    Month = as.integer(Month),
    Day = as.integer(Day),
    Hour = suppressWarnings(as.integer(Hour)),
    p.zoll = suppressWarnings(as.numeric(p.zoll)),
    p.in = suppressWarnings(as.numeric(p.in)),
    p.zoll.min = suppressWarnings(as.numeric(p.zoll.min)),
    p.in.min = suppressWarnings(as.numeric(p.in.min)),
    ta = suppressWarnings(as.numeric(ta)),
    force = suppressWarnings(as.numeric(force)),
    direction = ifelse(direction=="-", NA, direction)
  )
head(df3)



