library(dataresqc)
library(dplyr)
library(tidyr)
library(readxl)

source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/write_sef_f.R')

indir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_orig/Ukraine/'
outdir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/'


read_meta_nonofficial <- function (file = file.choose(), parameter = NULL) {
  header <- read.table(file = file, quote = "", comment.char = "", 
                       sep = "\t", nrows = 12, stringsAsFactors = FALSE, fill = TRUE)
  pars <- c("version", "id", "name", "lat", "lon", "alt", "source", 
            "link", "var", "stat", "units", "meta")
  if (is.null(parameter)) {
    out <- header[, 2]
    names(out) <- pars
  }
  else {
    out <- header[match(parameter, pars), 2]
    names(out) <- parameter
  }
  return(out)
}

# Dnipro -----------------------------------------------------------------
infile <- 'Dnipro.tsv'
df <- read.delim(paste0(indir, infile), header=T, sep='\t', stringsAsFactors = F,
                 skip=12)

meta <- read_meta_nonofficial(paste0(indir,infile))

# Create temperature dataframe
df.ta.Dnipro <- df %>%
  separate(Hour, into = c("Hour", "Minute"), sep = ":", convert = TRUE) %>%
  mutate(value = T * 1.25, Minute='NA') %>%
  select(Year, Month, Day, Hour, Minute, value)

# df.p.Dnipro <- df %>%
#   separate(Hour, into = c("Hour", "Minute"), sep = ":", convert = TRUE) %>%
#   mutate(
#   value = case_when(
#       Year >= 1833 & Year <= 1838 ~ P * 33.8639 ,        # inches → hPa
#       Year %in% c(1839:1842, 1850) ~ P * 1.27/(750.06)*1000,       
#       # R.s.l. → hPa ccording to (Shostin, 1975;Lamb, 1986) this unit can be converted into millimetres
#       # based on the relation 1 R.s.l. = 1.27 mm. That is, 1,000 hPa  = 1,000 mbar = 750.06 mmHg = 590.60 R.s.l
#       TRUE ~ NA_real_                                   # fallback if unexpected year
#     )
#   ) %>%
#   select(Year, Month, Day, Hour, Minute, value)

df.p.Dnipro <- df %>%
  separate(Hour, into = c("Hour", "Minute"), sep = ":", convert = TRUE) %>%
  mutate(
    Lmm = case_when(
      Year >= 1833 & Year <= 1838 ~ P * 25.399999704976,        # inches → mmHg
      Year %in% c(1839:1842, 1850) ~ P * 1.27,       
      # R.s.l. → hPa ccording to (Shostin, 1975;Lamb, 1986) this unit can be converted into millimetres
      # based on the relation 1 R.s.l. = 1.27 mm. That is, 1,000 hPa  = 1,000 mbar = 750.06 mmHg = 590.60 R.s.l
      TRUE ~ NA_real_                                   # fallback if unexpected year
    ),
    unit_label = case_when(
      Year >= 1833 & Year <= 1838 ~ "in",
      Year %in% c(1839:1842, 1850) ~ "R.s.l.",
      TRUE ~ "unknown"
    ),
    ta_out = T * 1.25,
    value = round(convert_pressure(p=Lmm, f=1, lat=as.numeric(meta[['lat']]), alt=as.numeric(meta[['alt']]), atb=ta_out),2),
    p_diff = round(value - Lmm*1.3332239, 2),
    meta = paste("orig_p=", unit_label, " | orig_ta_out=", round(ta_out, 1), " | Δp=", p_diff, "hPa", sep = "")
  ) %>%
  select(Year, Month, Day, Hour, Minute, value, meta)

write_sef_f(Data=df.ta.Dnipro,
            outpath=outdir, outfile="Dnipro_ta_subdaily.tsv",
            cod=meta[["id"]],
            variable='ta',
            nam=meta[["name"]],
            lat=meta[["lat"]],
            lon=meta[["lon"]], alt=meta[["alt"]], sou=meta[["source"]], metaHead = 'orig_ta=Reaumur',
            link=meta[["link"]], units='C', stat="point",
            meta="", keep_na = F)

write_sef_f(Data=df.p.Dnipro,
            outpath=outdir, outfile="Dnipro_p_subdaily.tsv",
            cod=meta[["id"]],
            variable='ta',
            nam=meta[["name"]],
            lat=meta[["lat"]],
            lon=meta[["lon"]], alt=meta[["alt"]], sou=meta[["source"]], metaHead = 'PGC=Y | PTC=Y',
            link=meta[["link"]], units='hPa', stat="point",
            meta=df.p.Dnipro$meta, keep_na = F)

# Kamyanets-Podilskyi -----------------------------------------------------------------

infile <- 'Kamyanets-Podilskyi.tsv'
df <- read.delim(paste0(indir, infile), header=T, sep='\t', stringsAsFactors = F,
                 skip=12)

meta <- read_meta_nonofficial(paste0(indir,infile))
lat_Kamyanets <- as.numeric(meta[['lat']])
alt_Kamyanets <- as.numeric(meta[['alt']])

# Create temperature dataframe
df.ta.Kamyanets <- df %>%
  separate(Time, into = c("Hour", "Minute"), sep = ":", convert = TRUE) %>%
  mutate(value = T * 1.25, Minute='NA') %>%
  select(Year, Month, Day, Hour, Minute, value)

df.p.Kamyanets <- df %>%
  separate(Time, into = c("Hour", "Minute"), sep = ":", convert = TRUE) %>%
  mutate(
    ta_out = T * 1.25,
    meta_ta = paste0("origa_ta_outside=", round(ta_out, 1)),
    Lmm = P * 25.4,     # inHg → mmHg
    
    value = round(convert_pressure(p=Lmm, f=1, lat=lat_Kamyanets, alt=alt_Kamyanets, atb=ta_out),2),
    p_diff = round(value - Lmm*1.3332239, 2),
    
    meta = paste("orig_p=in | ", meta_ta, " | Δp=", p_diff, "hPa", sep = "")
  ) %>%
  select(Year, Month, Day, Hour, Minute, value, meta)

write_sef_f(Data=df.ta.Kamyanets,
            outpath=outdir, outfile="Kamyanets_ta_subdaily.tsv",
            cod=meta[["id"]],
            variable='ta',
            nam=meta[["name"]],
            lat=meta[["lat"]],
            lon=meta[["lon"]], alt=meta[["alt"]], sou=meta[["source"]], metaHead = 'orig_ta=Reaumur',
            link=meta[["link"]], units='C', stat="point",
            meta="", keep_na = F)

write_sef_f(Data=df.p.Kamyanets,
            outpath=outdir, outfile="Kamyanets_p_subdaily.tsv",
            cod=meta[["id"]],
            variable='p',
            nam=meta[["name"]],
            lat=meta[["lat"]],
            lon=meta[["lon"]], alt=meta[["alt"]], sou=meta[["source"]], metaHead = 'PGC=Y | PTC=Y',
            link=meta[["link"]], units='hPa', stat="point",
            meta=df.p.Kamyanets$meta, keep_na = F)

# Kharkiv_University -----------------------------------------------------------------
infile <- 'Kharkiv_University.tsv'
df <- read.delim(paste0(indir, infile), header=T, sep='\t', stringsAsFactors = F,
                 skip=12)

meta <- read_meta_nonofficial(paste0(indir,infile))

# Create temperature dataframe
df.ta.Kharkiv <- df %>%
  separate(Hour, into = c("Hour", "Minute"), sep = ":", convert = TRUE) %>%
  mutate(value = ifelse(T==-999.9, NA_real_,T * 1.25), Minute='NA') %>%
  select(Year, Month, Day, Hour, Minute, value)

df.p.Kharkiv <- df %>%
  separate(Hour, into = c("Hour", "Minute"), sep = ":", convert = TRUE) %>%
  mutate(Lmm = ifelse(P==-999.9, NA_real_,P * 1.27),
         ta_out = T * 1.25,
         value = round(convert_pressure(p=Lmm, f=1, lat=as.numeric(meta['lat']), alt=as.numeric(meta['alt']), atb=ta_out),
                       2),
         p_diff = round(value - Lmm*1.3332239, 2),
         meta = paste("orig_p=in | orig_ta_out=", round(ta_out,2), " | Δp=", p_diff, "hPa", sep = "")
  ) %>%     # R.s.l. → hPa
  select(Year, Month, Day, Hour, Minute, value, meta)

write_sef_f(Data=df.ta.Kharkiv,
            outpath=outdir, outfile="Kharkiv_ta_subdaily.tsv",
            cod=meta[["id"]],
            variable='ta',
            nam=meta[["name"]],
            lat=meta[["lat"]],
            lon=meta[["lon"]], alt=meta[["alt"]], sou=meta[["source"]], metaHead = 'orig_ta=Reaumur',
            link=meta[["link"]], units='C', stat="point",
            meta="", keep_na = F)

write_sef_f(Data=df.p.Kharkiv,
            outpath=outdir, outfile="Kharkiv_p_subdaily.tsv",
            cod=meta[["id"]],
            variable='ta',
            nam=meta[["name"]],
            lat=meta[["lat"]],
            lon=meta[["lon"]], alt=meta[["alt"]], sou=meta[["source"]], metaHead = 'orig_p=R.s.l.',
            link=meta[["link"]], units='hPa', stat="point",
            meta=df.p.Kharkiv$meta, keep_na = F)



# Kherson -----------------------------------------------------------------

infile <- 'Kherson_pre1825.tsv'
df <- read.delim(paste0(indir, infile), header=T, sep='\t', stringsAsFactors = F,
                 skip=12)

# convert Reaumur to C, in to mmHg
df$T_on_N <- ifelse(df$T_on_N == -999.9, NA, df$T_on_N * 1.25)
df$T_on_S <- ifelse(df$T_on_S == -999.9, NA, df$T_on_S * 1.25)
df$P <-      ifelse(df$P == -999.9, NA, df$P * 25.4)

# df <- df %>% filter(P < 3000, P > 950)

# extract hour and minute
df$Minute <- as.integer(sub(".*:","", df$Hour))
df$Hour <- as.integer(sub(":.*", "", df$Hour))

### add 12 days to each date
# df$date <- as.Date(sprintf("%04d-%02d-%02d", df$Year, df$Month, df$Day))
# df$date <- df$date + 12 # add 12 days

# Update the Year, Month, and Day columns based on the new date
# df$Year <- as.integer(format(df$date, "%Y"))
# df$Month <- as.integer(format(df$date, "%m"))
# df$Day <- as.integer(format(df$date, "%d"))

# read meta from lines
meta_lines <- strsplit(readLines(paste0(indir, infile), n=12), '\t')

meta <- list()
for (line in meta_lines) {
  if (length(line) >= 2 ** nzchar(line[2])) {
    key <- line[1]
    value <- line[2]
    meta[[key]] <- value
  } else { # for lines wo a value (e.g. SEF, ID)
    meta[[line[1]]] <- NA
  }
}

# create two different files for T on N and T on S
df_T_on_N <- df[, c("Year","Month","Day","Hour","Minute","T_on_N")]
colnames(df_T_on_N)[colnames(df_T_on_N)=="T_on_N"] <- "Value"
df_T_on_S <- df[, c("Year","Month","Day","Hour","Minute","T_on_S")]
colnames(df_T_on_S)[colnames(df_T_on_S)=="T_on_S"] <- "Value"

# generate TSV for T_on_N
meta[['Vbl']] <- 'ta'
meta[['Units']] <- 'C'

write_sef_f(Data=df_T_on_S,
            outpath=outdir, outfile="Kherson_ta_subdaily_TonS.tsv",
            cod=meta[["ID"]],
            variable=meta[["Vbl"]],
            nam=meta[["Name"]],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]],
            metaHead = "T_on_N | orig_ta=Reaumur | Hours correspond to original day periods morning, midday, evening",
            link=meta[["Link"]], units=meta[["Units"]], stat="point", keep_na = F)
write_sef_f(Data=df_T_on_N,
            outpath=outdir, outfile="Kherson_ta_subdaily_TonN.tsv",
            cod=meta[["ID"]],
            variable=meta[["Vbl"]],
            nam=meta[["Name"]],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]], 
            metaHead = "T_on_N | orig_ta=Reaumur | Hours correspond to original day periods morning, midday, evening",
            link=meta[["Link"]], units=meta[["Units"]], stat="point", keep_na = F)

# create pressure file
df.p.Kherson <- df %>%
  mutate(value = round(convert_pressure(p=P, f=1, lat=as.numeric(meta$Lat), alt=as.numeric(meta$Alt), atb=T_on_N), 2),
         p_diff = round(value - P*1.3332239,2),
         meta = paste("orig_p=in | orig_ta_onN=", round(T_on_N,2), " | Δp=", p_diff, "hPa", sep = "")
  ) %>%     # R.s.l. → hPa
  select(Year, Month, Day, Hour, Minute, value, meta)

write_sef_f(Data=df.p.Kherson,
            outpath=outdir, outfile="Kherson_p_subdaily",
            cod=meta[["ID"]],
            variable='p',
            nam=meta[["Name"]],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]], metaHead = "PGC=Y | PTC=Y",
            link=meta[["Link"]], units='hPa', stat="point",
            meta=df.p.Kherson$meta, keep_na = F)

# Kyiv --------------------------------------------------------------------

infile <- 'Kyiv_ta_pre1837.tsv'
df <- read.delim(paste0(indir, infile), header=T, sep='\t', stringsAsFactors = F,
                 skip=12)

# convert Reamur to Kelvin
df$Value <- ifelse(df$Value == -999.9, NA, df$Value * 1.25)

# extract hour and minute
df$Minute <- as.integer(sub(".*:","", df$Hour))
df$Hour <- as.integer(sub(":.*", "", df$Hour))

# rearrange columns
df <- df[, c("Year","Month","Day","Hour","Minute","Value")]

# read meta for SEF
meta_lines <- strsplit(readLines(paste0(indir, infile), n=12), '\t')

meta <- list()
for (line in meta_lines) {
  if (length(line) >= 2 ** nzchar(line[2])) {
    key <- line[1]
    value <- line[2]
    meta[[key]] <- value
  } else { # for lines wo a value (e.g. SEF, ID)
    meta[[line[1]]] <- NA
  }
}
meta[['Vbl']] <- 'ta'
meta[['Units']] <- 'C'
meta[['metaHead']] <- 'Hours correspond to original day periods morning, midday, evening'

write_sef_f(Data=df,
            outpath=outdir, outfile='Kyiv_ta_subdailyNO.tsv',
            cod=meta[["ID"]],
            variable=meta[["Vbl"]],
            nam=meta[["Name"]],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]],
            metaHead = meta[['metaHead']],
            link=meta[["Link"]], units=meta[["Units"]], stat="point",
            meta="orig_ta=Reaumur", keep_na = F)

# Lugansk --------------------------------------------------------------------

alt_Lugansk=49
lat_Lugansk=48.565556

filepath <- paste0(indir, 'Lugansk.xlsx')
sheets <- setdiff(excel_sheets(paste0(indir,'Lugansk.xlsx')),"Meta")

df <- bind_rows(
  lapply(sheets, function(sheet) {
    df <- tryCatch(
        {
          suppressMessages(suppressWarnings(
            read_excel(paste0(indir,'Lugansk.xlsx'), sheet = sheet, range = cell_cols(1:7), .name_repair = "minimal")
          ))
        },
        warning = function(w) {
          message("⚠️  Warning in sheet: ", sheet, " → ", conditionMessage(w))
          suppressMessages(suppressWarnings(read_excel(file_path, sheet = sheet, range = cell_cols(1:7), .name_repair = "minimal")))
        }
    )
    time_cols <- grep("^Time", names(df))
    if (length(time_cols) > 1) {
      # remove first "time" column (morning/midday/evening) if duplicated
      df <- df[, -time_cols[1]]
    }
    
    # Ensure required columns exist
    if (!"P, in" %in% names(df)) {
      df$`P, in` <- NA_real_
    }
    if (!"P, R.s.l." %in% names(df)) {
      df$`P, R.s.l.` <- NA_real_
    }
    if (!"Pt, R" %in% names(df)) {
      df$`Pt, R` <- NA_real_
    }
    df <- df[, intersect(names(df), c("Year", "Month", "Day", "Time", "T, R", "P, in", "P, R.s.l.", "Pt, R"))]
    df <- df %>%
      mutate(across(c("T, R", "P, in", "P, R.s.l.", "Pt, R"), ~na_if(., -999.9))) # -999.9 to NA
  
    df <- df %>%
      mutate(
        Hour = as.integer(sub(":.*", "", Time)),
        Minute = ifelse(grepl(":", Time), as.integer(sub(".*:", "", Time)), NA_integer_)
      )
    # rearrange and drop time col
    df <- df %>%
      select(Year, Month, Day, Hour, Minute, everything(), -Time)
    df
  })
)

df.ta.Lugansk <- df %>%
  mutate(value = `T, R` * 1.25) %>%
  select(Year, Month, Day, Hour, Minute, value)

df.ta.Lugansk <- as.data.frame(df.ta.Lugansk)

df.p.Lugansk <- df %>%
  # Check if there are any rows where both "P, in" and "P, R.s.l." are not NA (where we have obs in both units)
  filter(!is.na(`P, in`) | !is.na(`P, R.s.l.`)) %>%
  mutate(
    ta_bar = `Pt, R` * 1.25, # barometer Temperature in C
    ta_out = `T, R`  * 1.25, # outside temperature in C
    ta_used = case_when(     # use barometer Ta when available, outside when not
      !is.na(ta_bar) ~ ta_bar,
      TRUE ~ ta_out
    ),
    meta_ta = case_when(    # prepare meta column
      !is.na(ta_bar) ~ paste0("orig_ta_bar=", round(ta_bar,1)),
      !is.na(ta_out) ~ paste0("orig_ta_outside=", round(ta_out, 1)),
      TRUE ~ "orig_ta_unknown"
    ),
    # gamma = 1.82e-4, # thermal expansion coeff of mercury at 0°C
    
    # original pressure in mmHg
    Lmm = case_when(
      !is.na(`P, in`) ~ `P, in` * 25.4, # convert inHg → mmHg
      !is.na(`P, R.s.l.`) ~ (`P, R.s.l.` * 1.27 / 750.06 * 1000) / 33.8639
    ),
    
    # corrected P in hPa
    value = convert_pressure(p=Lmm, f=1, lat=lat_Lugansk, alt=alt_Lugansk, atb=ta_used),
    uncorrected = convert_pressure(p = Lmm, f = 1, lat = 48.57, alt = 90, atb = rep(0, n())),  # temp = 0 for reference

    # difference in pressure:
    p_diff = round(value - uncorrected, 2),
    
    meta_orig = case_when(
      !is.na(`P, in`) ~ "orig_p=in",
      !is.na(`P, R.s.l.`) ~ "orig_p=R.s.l."
    ),
    meta= paste(meta_orig, " | ", meta_ta, " | Δp=", p_diff,"hPa",sep="")
  ) %>%
  select(Year, Month, Day, Hour, Minute, value, meta)

df.p.Lugansk <- as.data.frame(df.p.Lugansk)

write_sef_f(Data=df.ta.Lugansk,
            outpath=outdir, outfile='Lugansk_ta_subdaily.tsv',
            cod='Lugansk',
            variable='ta',
            nam='Lugansk',
            lat=lat_Lugansk,
            lon='39.2275', alt=alt_Lugansk, sou="Ukrainian early (pre-1850) historical weather observations, Skrynk et al.",
            link="https://doi.org/10.15407/uhmi.report.01", units="C", stat="point",
            metaHead="relocation(1843)",
            meta="orig_ta=Reaumur", keep_na = F)

write_sef_f(Data=df.p.Lugansk,
            outpath=outdir, outfile='Lugansk_p_subdaily.tsv',
            cod='Lugansk',
            variable='p',
            nam='Lugansk',
            lat='48.565556',
            lon='39.2275', alt='59', sou="Ukrainian early (pre-1850) historical weather observations, Skrynk et al.",
            link="https://doi.org/10.15407/uhmi.report.01", units="hPa", stat="point",
            metaHead="relocation(1843) | PGC=Y | PTC=Y",
            meta=df.p.Lugansk$meta, keep_na = F)

# Odesa --------------------------------------------------------------------
lat_Odesa <- NA
lon_Odesa <- NA 
alt_Odesa <- NA
sheets <- setdiff(excel_sheets(paste0(indir,'Odesa.xlsx')),"Meta")

df <- bind_rows(
  lapply(sheets, function(sheet) {
    df <- read_excel(paste0(indir,'Odesa.xlsx'), sheet = sheet, range = cell_cols(1:6))
    time_cols <- grep("^Time", names(df))
    if (length(time_cols) > 1) {
      # remove first "time" column (morning/midday/evening) if duplicated
      df <- df[, -time_cols[1]]
    }
    df 
  })
)

df <- df %>%
  separate(Time, into = c("Hour", "Minute"), sep = ":", convert = TRUE) %>%
  mutate(Minute=NA_integer_)

df.ta.Odesa <- df %>%
  mutate(value = `T, R` * 1.25) %>%
  select(Year, Month, Day, Hour, Minute, value)

df.ta.Odesa <- as.data.frame(df.ta.Odesa)

df.p.Odesa <- df %>%
  # Check if there are any rows where both "P, in" and "P, R.s.l." are not NA (where we have obs in both units)
  filter(!is.na(`P, in`) | !is.na(`P, R.s.l.`)) %>%
  mutate(
    ta_out = `T, R` * 1.25, # outside temperature in °C
    meta_ta = paste0("orig_ta_outside=", round(ta_out,1)),
    
    Lmm = case_when(
      !is.na(`P, in`) ~ `P, in`  * 25.4,                                    # inHg to mmHg
      !is.na(`P, R.s.l.`) ~ (`P, R.s.l.` * 1.27 / 750.06 * 1000) / 33.8639  # hPa to mmHg
    ),
    
    # Pressure corrected to 0°C + gravity → hPa
    value = convert_pressure(p = Lmm, f = 1, lat = lat_Odesa, alt = alt_Odesa, atb = ta_out),
    
    # Uncorrected pressure for delta_p (using atb = 0)
    uncorrected = convert_pressure(p = Lmm, f = 1, lat = lat_Odesa, alt = alt_Odesa, atb = rep(0, n())),
    p_diff = round(value - uncorrected, 2),
    
    meta_orig = case_when(
      !is.na(`P, in`) ~ "orig_p=in",
      !is.na(`P, R.s.l.`) ~ "orig_p=R.s.l."
    ),
    
    meta = paste(meta_orig, " | ", meta_ta, " | Δp=", p_diff, "hPa", sep = "")
  ) %>%
  select(Year, Month, Day, Hour, Minute, value, meta)

df.p.Odesa <- as.data.frame(df.p.Odesa)

write_sef_f(Data=df.ta.Odesa,
            outpath=outdir, outfile='Odesa_ta_subdaily.tsv',
            cod='Odesa',
            variable='ta',
            nam='Odesa',
            lat=lat_Odesa,
            lon=lon_Odesa, alt=alt_Odesa, sou="Ukrainian early (pre-1850) historical weather observations, Skrynk et al.",
            link="https://doi.org/10.15407/ uhmi.report.01", units="C", stat="point",
            meta="orig_ta=Reaumur", keep_na = F)

write_sef_f(Data=df.p.Odesa,
            outpath=outdir, outfile='Odesa_p_subdaily.tsv',
            cod='Odesa',
            variable='p',
            nam='Odesa',
            lat=lat_Odesa,
            lon=lon_Odesa, alt=alt_Odesa, sou="Ukrainian early (pre-1850) historical weather observations, Skrynk et al.",
            link="https://doi.org/10.15407/ uhmi.report.01", units="hPa", stat="point",
            metaHead="PGC=Y | PTC=Y",
            meta=df.p.Odesa$meta, keep_na = F)


# Poltava --------------------------------------------------------------------
lat_Poltava <- 49.609444
alt_Poltava <- 160
sheets <- setdiff(excel_sheets(paste0(indir,'Poltava.xlsx')),"Meta")

df <- bind_rows(
  lapply(sheets, function(sheet) {
    df <- read_excel(paste0(indir,'Poltava.xlsx'), sheet = sheet, range = cell_cols(1:6))
    time_cols <- grep("^Time", names(df))
    if (length(time_cols) > 1) {
      # remove first "time" column (morning/midday/evening) if duplicated
      df <- df[, -time_cols[1]]
    }
    df 
  })
)

df <- df %>%
  separate(Time, into = c("Hour", "Minute"), sep = ":", convert = TRUE) %>%
  mutate(Minute=NA_integer_)
  
df.ta.Poltava <- df %>%
  mutate(value = `T, R` * 1.25) %>%
  select(Year, Month, Day, Hour, Minute, value)

df.ta.Poltava <- as.data.frame(df.ta.Poltava)

df.p.Poltava <- df %>%
  # Check if there are any rows where both "P, in" and "P, R.s.l." are not NA (where we have obs in both units)
  filter(!is.na(`P, in`) | !is.na(`P, R.s.l.`)) %>%
  mutate(
    ta_out = `T, R` * 1.25, # outside temperature in °C
    meta_ta = paste0("orig_ta_outside=", round(ta_out,1)),
    
    Lmm = case_when(
      !is.na(`P, in`) ~ `P, in`  * 25.4,                                    # inHg to mmHg
      !is.na(`P, R.s.l.`) ~ (`P, R.s.l.` * 1.27 / 750.06 * 1000) / 33.8639  # hPa to mmHg
    ),
    
    # Pressure corrected to 0°C + gravity → hPa
    value = convert_pressure(p = Lmm, f = 1, lat = lat_Poltava, alt = alt_Poltava, atb = ta_out),
    
    # Uncorrected pressure for delta_p (using atb = 0)
    uncorrected = convert_pressure(p = Lmm, f = 1, lat = lat_Poltava, alt = alt_Poltava, atb = rep(0, n())),
    p_diff = round(value - uncorrected, 2),
    
    meta_orig = case_when(
      !is.na(`P, in`) ~ "orig_p=in",
      !is.na(`P, R.s.l.`) ~ "orig_p=R.s.l."
    ),
    
    meta = paste(meta_orig, " | ", meta_ta, " | Δp=", p_diff, "hPa", sep = "")
  ) %>%
  select(Year, Month, Day, Hour, Minute, value, meta)

df.p.Poltava <- as.data.frame(df.p.Poltava)

write_sef_f(Data=df.ta.Poltava,
            outpath=outdir, outfile='Poltava_ta_subdaily.tsv',
            cod='Poltava',
            variable='ta',
            nam='Poltava',
            lat=lat_Poltava,
            lon='34.544722', alt=alt_Poltava, sou="Ukrainian early (pre-1850) historical weather observations, Skrynk et al.",
            link="https://doi.org/10.15407/ uhmi.report.01", units="C", stat="point",
            meta="orig_ta=Reaumur", keep_na = F)

write_sef_f(Data=df.p.Poltava,
            outpath=outdir, outfile='Poltava_p_subdaily.tsv',
            cod='Poltava',
            variable='p',
            nam='Poltava',
            lat='49.609444',
            lon='34.544722', alt='160', sou="Ukrainian early (pre-1850) historical weather observations, Skrynk et al.",
            link="https://doi.org/10.15407/ uhmi.report.01", units="hPa", stat="point",
            metaHead="PGC=Y | PTC=Y",
            meta=df.p.Poltava$meta, keep_na = F)
