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
  mutate(value = T * 1.25, Minute='NA', meta = paste0('orig_ta=',T,"R | orig_time=",Time)) %>%
  select(Year, Month, Day, Hour, Minute, value, meta)

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
    meta = paste0('orig_p=', P, unit_label, ' | orig_time=', Time)
  ) %>%
  select(Year, Month, Day, Hour, Minute, value, meta)

write_sef_f(Data=df.ta.Dnipro,
            outpath=outdir, outfile="Dnipro_ta_subdaily.tsv",
            cod=meta[["id"]],
            variable='ta',
            nam=meta[["name"]],
            lat=meta[["lat"]],
            lon=meta[["lon"]], alt=meta[["alt"]], sou=meta[["source"]],
            link=meta[["link"]], units='C', stat="point",
            meta=df.ta.Dnipro$meta, keep_na = F)

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
  mutate(value = T * 1.25, Minute='NA', meta=paste0('orig_ta=',T,'R')) %>%
  separate(Time, into = c("Hour", "Minute"), sep = ":", convert = TRUE) %>%
  select(Year, Month, Day, Hour, Minute, value, meta)

df.p.Kamyanets <- df %>%
  mutate(
    ta_out = T * 1.25,
    meta_ta = paste0("origa_ta_outside=", round(ta_out, 1)),
    Lmm = P * 25.4,     # inHg → mmHg
    value = round(convert_pressure(p=Lmm, f=1, lat=lat_Kamyanets, alt=alt_Kamyanets, atb=ta_out),2),
    p_diff = round(value - Lmm*1.3332239, 2),
    meta = paste0('orig_p=', P, 'in')
    ) %>%
  separate(Time, into = c("Hour", "Minute"), sep = ":", convert = TRUE) %>%
  select(Year, Month, Day, Hour, Minute, value, meta)

write_sef_f(Data=df.ta.Kamyanets,
            outpath=outdir, outfile="Kamyanets_ta_subdaily.tsv",
            cod=meta[["id"]],
            variable='ta',
            nam=meta[["name"]],
            lat=meta[["lat"]],
            lon=meta[["lon"]], alt=meta[["alt"]], sou=meta[["source"]],
            link=meta[["link"]], units='C', stat="point",
            meta=df.ta.Kamyanets$meta, keep_na = F)

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
# ''' there are two files with Kharkiv. See page 65 of Ukranian data paper. Only temperature data are different
#     for the first months of 1843 while all the values of air pressure and air temperature for the remaining months
#     are almost the same.'''

infile <- 'Kharkiv_University.tsv'
dfA <- read.delim(paste0(indir, infile), header=T, sep='\t', stringsAsFactors = F,
                 skip=12) %>%
       separate(Hour, into = c("Hour", "Minute"), sep = ":", convert = TRUE)

filepath <- paste0(indir, 'Kharkiv.xlsx')
sheets <- setdiff(excel_sheets(filepath),"Meta")
   
dfB <- bind_rows(
  lapply(sheets, function(sheet) {
    dfB <- read_excel(filepath,
                     sheet = sheet,
                     range = cell_cols(1:6),
                     col_types = c("numeric", "numeric", "numeric", "text", "numeric", "numeric"))
  })) %>%
  mutate(
    Time = case_when(
      Time == "0.625" ~ "15:00",
      Time == "0.375" ~ "09:00",
      Time == "0.875" ~ "21:00",
      TRUE ~ Time
    )
  ) %>%
  separate(Time, into = c("Hour", "Minute"), sep = ":", convert = TRUE)

meta <- read_meta_nonofficial(paste0(indir,infile))

merge.Kharkiv <- function(dfA, dfB, var = c("ta", "p")) {
  var <- match.arg(var)
  
  # Full join on date-time columns
  merged <- full_join(dfA %>% rename(value_A = value, meta_A = meta),
                      dfB %>% rename(value_B = value, meta_B = meta),
                      by = c("Year", "Month", "Day", "Hour", "Minute"))
  
  # Decide final value and updated meta
  result <- merged %>%
    mutate(
      value = ifelse(!is.na(value_A), value_A, value_B),
      meta = case_when(
        !is.na(value_A) & !is.na(value_B) ~ paste0(meta_A, " | altern_", meta_B),
        !is.na(value_A) ~ meta_A,
        TRUE ~ meta_B
      )
    ) %>%
    select(Year, Month, Day, Hour, Minute, value, meta)
  
  return(result)
}

# Create temperature dataframe
dfA.ta.Kharkiv <- dfA %>%
  mutate(value = ifelse(T==-999.9, NA_real_,T * 1.25), meta=paste0('orig_ta=',T,'R')) %>%
  select(Year, Month, Day, Hour, Minute, value, meta)

dfB.ta.Kharkiv <- dfB %>%
  rename('T' = "T, R") %>%
  mutate(value = ifelse(T==-999.9, NA_real_,T * 1.25), meta=paste0('orig_ta=',T,'R')) %>%
  select(Year, Month, Day, Hour, Minute, value, meta)

df.ta.Kharkiv <- merge.Kharkiv(dfA.ta.Kharkiv, dfB.ta.Kharkiv, var='ta')


dfA.p.Kharkiv <- dfA %>%
  mutate(Lmm = ifelse(P==-999.9, NA_real_,P * 1.27),
         ta_out = T * 1.25,
         value = round(convert_pressure(p=Lmm, f=1, lat=as.numeric(meta['lat']), alt=as.numeric(meta['alt']), atb=ta_out),
                       2),
         meta = paste0("orig_p=",P,"R.s.l.")
  ) %>%     # R.s.l. → hPa
  select(Year, Month, Day, Hour, Minute, value, meta)

dfB.p.Kharkiv <- dfB %>%
  rename('T'='T, R', 'P'='P, R.s.l.') %>%
  mutate(Lmm = ifelse(P==-999.9, NA_real_,P * 1.27),
         ta_out = T * 1.25,
         value = round(convert_pressure(p=Lmm, f=1, lat=as.numeric(meta['lat']), alt=as.numeric(meta['alt']), atb=ta_out),
                       2),
         meta = paste0("orig_p=",P,"R.s.l.")
  ) %>%     # R.s.l. → hPa
  select(Year, Month, Day, Hour, Minute, value, meta)

df.p.Kharkiv <- merge.Kharkiv(dfA.p.Kharkiv, dfB.p.Kharkiv, var='p')


write_sef_f(Data=df.ta.Kharkiv,
            outpath=outdir, outfile="Kharkiv_ta_subdaily.tsv",
            cod=meta[["id"]],
            variable='ta',
            nam=meta[["name"]],
            lat=meta[["lat"]],
            lon=meta[["lon"]], alt=meta[["alt"]], sou=meta[["source"]],
            link=meta[["link"]], units='C', stat="point",
            meta=df.ta.Kharkiv$meta, keep_na = F)

write_sef_f(Data=df.p.Kharkiv,
            outpath=outdir, outfile="Kharkiv_p_subdaily.tsv",
            cod=meta[["id"]],
            variable='ta',
            nam=meta[["name"]],
            lat=meta[["lat"]],
            lon=meta[["lon"]], alt=meta[["alt"]], sou=meta[["source"]], metaHead = 'PGC=Y | PTC=Y',
            link=meta[["link"]], units='hPa', stat="point",
            meta=df.p.Kharkiv$meta, keep_na = F)



# Kherson -----------------------------------------------------------------

alt_Kherson=47
lat_Kherson=46.73833

filepath <- paste0(indir, 'Kherson.xlsx')
sheets <- setdiff(excel_sheets(paste0(indir,'Kherson.xlsx')),"Meta")

df <- bind_rows(
  lapply(sheets, function(sheet) {
    df <- tryCatch(
      {
        suppressMessages(suppressWarnings(
          read_excel(paste0(indir,'Kherson.xlsx'), sheet = sheet, range = cell_cols(1:8), .name_repair = "minimal")
        ))
      },
      warning = function(w) {
        message("⚠️  Warning in sheet: ", sheet, " → ", conditionMessage(w))
        suppressMessages(suppressWarnings(read_excel(file_path, sheet = sheet, range = cell_cols(1:8), .name_repair = "minimal")))
      }
    )
    #df <- df %>% mutate(meta=NA_character_)
    #df$meta <- NA_character_
    
    time_cols <- grep("^Time", names(df))
    if (length(time_cols) > 1) {
      # save first "time" column (morning/midday/evening) if duplicated
      df <- df[, -time_cols[1]]
    }
    
    # Ensure required columns exist
    if (!"T, on N, R" %in% names(df)) {
      df$`T, on N, R` <- NA_real_
    }
    if (!"T, on S, R" %in% names(df)) {
      df$`T, on S, R` <- NA_real_
    }
    if (!"P, in" %in% names(df)) {
      df$`P, in` <- NA_real_
    }
    if (!"T, R" %in% names(df)) {
      df$`T, R` <- NA_real_
    }
    if (!"Tm/Tn, R" %in% names(df)) {
      df$`Tm/Tn, R` <- NA_real_
    }
    df <- df[, intersect(names(df), c("Year", "Month", "Day", "Time", "T, on N, R", "T, on S, R", "T, R", "P, in", "Tm/Tn, R"))]
    df <- df %>%
      mutate(across(c("T, on N, R", "T, on S, R", "T, R", "P, in", "Tm/Tn, R"), ~na_if(., -999.9))) # -999.9 to NA
    
    df <- df %>%
      mutate(
        Hour = as.integer(sub(":.*", "", Time)),
        Minute = ifelse(grepl(":", Time), as.integer(sub(".*:", "", Time)), NA_integer_),
        meta = case_when(
          Hour == 6 ~ "orig_time=morning",
          Hour == 14 ~ "orig_time=midday",
          Hour == 22 ~ "orig_time=evening",
          Hour == 10 ~ "orig_time=morning"
        )
      )

    df <- df %>%
      select(Year, Month, Day, Hour, Minute, meta, everything(), -Time)
  })
)

df.ta.Kherson <- df %>%
  
  mutate(
    ta_onN = `T, on N, R` * 1.25,
    ta_onS = `T, on S, R` * 1.25,
    ta     = `T, R` * 1.25,
    value = case_when(
      !is.na(ta) ~ ta,
      !is.na(ta_onN) ~ ta_onN,
      !is.na(ta_onS) ~ ta_onS,
      TRUE ~ NA_real_
    ),
    meta_ta = case_when( # if N and S walls available, use N wall for value and S wall as alternative
      !is.na(ta_onN) ~ paste0(meta, " | orig_ta=",ta_onN,"R | altern_orig_ta=",ta_onS,"R"),
      TRUE ~ paste0(meta, " | orig_ta=",ta,"R")
    )
  ) %>%
  select(Year, Month, Day, Hour, Minute, value, meta_ta)

df.ta.Kherson <- as.data.frame(df.ta.Kherson)

df.p.Kherson <- df %>%
  mutate(
    ta_onN = `T, on N, R` * 1.25,
    ta_onS = `T, on S, R` * 1.25,
    ta = `T, R` * 1.25,
    ta_used = case_when(
      !is.na(ta) ~ ta,
      !is.na(ta_onN) ~ ta_onN,
      !is.na(ta_onS) ~ ta_onS,
    ),
    meta_ta = case_when(    # prepare meta column
      !is.na(ta) ~ paste0("orig_ta=", round(ta_used,1)),
      !is.na(ta_onN) ~ paste0("orig_ta='TonN',", round(ta_used, 1)),
      !is.na(ta_onS) ~ paste0("orig_ta='TonS',", round(ta_used, 1)),
      TRUE ~ "orig_ta_unknown"
    ),
    
    # corrected P in hPa
    value = round(convert_pressure(p=`P, in` * 25.4, f=1, lat=lat_Kherson, alt=alt_Kherson, atb=ta_used),2),
    
    # difference in pressure:
    p_diff = round(value - `P, in` * 25.4 * 1.33322368, 2),

    meta= paste(meta_ta, " | Δp=", p_diff,"hPa",sep="")
  ) %>%
  select(Year, Month, Day, Hour, Minute, value, meta)

df.p.Kherson <- as.data.frame(df.p.Kherson)

# read meta for SEF
meta_lines <- strsplit(readLines(paste0(indir, 'Kherson_pre1825.tsv'), n=12), '\t')

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

write_sef_f(Data=df.ta.Kherson,
            outpath=outdir, outfile='Kherson_ta_subdaily.tsv',
            cod=meta[["ID"]],
            variable='ta',
            nam=meta[["Name"]],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]],
            meta = df.ta.Kherson$meta_ta,
            link=meta[["Link"]], units='C', stat="point", keep_na = F)

write_sef_f(Data=df.p.Kherson,
            outpath=outdir, outfile='Kherson_p_subdaily.tsv',
            cod=meta[["ID"]],
            variable='p',
            nam=meta[["Name"]],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]],
            metaHead = paste0(meta[['metaHead']], " | orig_p=R.s.l. | PGC=Y | PTC=Y"),
            link=meta[["Link"]], units='hPa', stat="point",
            meta=df.p.Kherson$meta, keep_na = F)

# Kyiv --------------------------------------------------------------------

alt_Kyiv = 166
lat_Kyiv = 50.39222

filepath <- paste0(indir, 'Kyiv.xlsx')
sheets <- setdiff(excel_sheets(paste0(indir,'Kyiv.xlsx')),"Meta")

df <- bind_rows(
  lapply(sheets, function(sheet) {
    df <- tryCatch(
      {
        suppressMessages(suppressWarnings(
          read_excel(paste0(indir,'Kyiv.xlsx'), sheet = sheet, range = cell_cols(1:7), .name_repair = "minimal")
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
    
    if (!"P, R.s.l." %in% names(df)) { # if sheet doesn't have pressure, then leave it empty
      df$`P, R.s.l.` <- NA_real_
    }

    df <- df[, intersect(names(df), c("Year", "Month", "Day", "Time", "T, R", "P, R.s.l."))]
    df <- df %>%
      mutate(across(c("T, R", "P, R.s.l."), ~na_if(., -999.9))) # -999.9 to NA
    
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


df.ta.Kyiv <- df %>%
  mutate(value = `T, R` * 1.25) %>%
  select(Year, Month, Day, Hour, Minute, value)

df.ta.Kyiv <- as.data.frame(df.ta.Kyiv)

df.p.Kyiv <- df %>%
  mutate(
    ta = `T, R`  * 1.25,    # outside temperature in C
    meta_ta = if_else(!is.na(ta), paste0("orig_ta=", round(ta, 1)), "orig_ta=NA"),
    
    # original pressure in mmHg
    Lmm = `P, R.s.l.` * 1.27, # R.s.l. -> mmHg ## / 750.06 * 1000) / 33.8639

    # corrected P in hPa
    value = round(convert_pressure(p=Lmm, f=1, lat=lat_Kyiv, alt=alt_Kyiv, atb=ta),2),
    
    # difference in pressure:
    p_diff = round(value - Lmm * 1.33322368, 2),
    

    meta= paste(meta_ta, " | Δp=", p_diff,"hPa",sep="")
  ) %>%
  select(Year, Month, Day, Hour, Minute, value, meta)

df.p.Kyiv <- as.data.frame(df.p.Kyiv)

# read meta for SEF
meta_lines <- strsplit(readLines(paste0(indir, 'Kyiv_all_post1837.tsv'), n=12), '\t')

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

meta[['metaHead']] <- "Hours correspond to 'morning/midday/evening'"

write_sef_f(Data=df.ta.Kyiv,
            outpath=outdir, outfile='Kyiv_ta_subdaily.tsv',
            cod=meta[["ID"]],
            variable='ta',
            nam=meta[["Name"]],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]],
            metaHead = paste(meta[['metaHead']],"| orig_ta=Reaumur"),
            link=meta[["Link"]], units='C', stat="point", keep_na = F)

write_sef_f(Data=df.p.Kyiv,
            outpath=outdir, outfile='Kyiv_p_subdaily.tsv',
            cod=meta[["ID"]],
            variable='p',
            nam=meta[["Name"]],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]],
            metaHead = paste(meta[['metaHead']],"| orig_p=R.s.l. | PGC=Y | PTC=Y"),
            link=meta[["Link"]], units='hPa', stat="point",
            meta=df.p.Kyiv$meta, keep_na = F)

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
      !is.na(ta_out) ~ paste0("orig_ta_out=", round(ta_out, 1)),
      TRUE ~ "orig_ta_unknown"
    ),
    # gamma = 1.82e-4, # thermal expansion coeff of mercury at 0°C
    
    # original pressure in mmHg
    Lmm = case_when(
      !is.na(`P, in`) ~ `P, in` * 25.4, # convert inHg → mmHg
      !is.na(`P, R.s.l.`) ~ (`P, R.s.l.` * 1.27) # R.s.l. -> mmHg ## / 750.06 * 1000) / 33.8639
    ),
    
    # corrected P in hPa
    value = round(convert_pressure(p=Lmm, f=1, lat=lat_Lugansk, alt=alt_Lugansk, atb=ta_used),2),

    # difference in pressure:
    p_diff = round(value - Lmm * 1.33322368, 2),
    
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
            metaHead="relocation(1843) | orig_ta=Reaumur | Hours correspond to 'morning/midday/evening'", keep_na = F)

write_sef_f(Data=df.p.Lugansk,
            outpath=outdir, outfile='Lugansk_p_subdaily.tsv',
            cod='Lugansk',
            variable='p',
            nam='Lugansk',
            lat='48.565556',
            lon='39.2275', alt='59', sou="Ukrainian early (pre-1850) historical weather observations, Skrynk et al.",
            link="https://doi.org/10.15407/uhmi.report.01", units="hPa", stat="point",
            metaHead="relocation(1843) | PGC=Y | PTC=Y | Hours correspond to 'morning/midday/evening'",
            meta=df.p.Lugansk$meta, keep_na = F)

# Odesa --------------------------------------------------------------------
lat_Odesa <- 46.440833
lon_Odesa <- 30.770278 
alt_Odesa <- 42
sheets <- setdiff(excel_sheets(paste0(indir,'Odesa.xlsx')),"Meta")

df <- bind_rows(
  lapply(sheets, function(sheet) {
    df <- read_excel(paste0(indir,'Odesa.xlsx'), 
                     sheet = sheet,
                     range = cell_cols(1:6),
                     col_types = c("numeric", "numeric", "numeric", "text", "numeric", "numeric"))
  })
)

df <- df %>%
  mutate(
    Time = case_when(
      Time == "0.625" ~ "15:00",
      Time == "0.66666666666666663" ~ "16:00",
      Time == "0.75" ~ "18:00",
      TRUE ~ Time
    ),
    Day = ifelse(Month == 6 & Day == 31, 30, Day)
  )


df <- df %>%
  separate(Time, into = c("Hour", "Minute"), sep = ":", convert = TRUE) %>%
  mutate(Minute=NA_integer_)

df.ta.Odesa <- df %>%
  mutate(value = ifelse(`T, R`==-999.9, NA_real_, `T, R` * 1.25)) %>%
  select(Year, Month, Day, Hour, Minute, value)

df.ta.Odesa <- as.data.frame(df.ta.Odesa)

df.p.Odesa <- df %>%
  mutate(
    ta_out = ifelse(`T, R`==-999.9, NA_real_, `T, R` * 1.25), # outside temperature in °C
    meta_ta = paste0("orig_ta_out=", round(ta_out,1)),
    
    P = ifelse(`P, in`==-999.9, NA_real_, `P, in`*25.4),
    # Pressure corrected to 0°C + gravity → hPa
    value = convert_pressure(p = P, f = 1, lat = lat_Odesa, alt = alt_Odesa, atb = ta_out),
    p_diff = round(value - P, 2),
    
    meta = paste(meta_ta, " | Δp=", p_diff, "hPa", sep = "")
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
            link="https://doi.org/10.15407/uhmi.report.01", units="C", stat="point",
            metaHead="orig_ta=Reaumur", keep_na = F)

write_sef_f(Data=df.p.Odesa,
            outpath=outdir, outfile='Odesa_p_subdaily.tsv',
            cod='Odesa',
            variable='p',
            nam='Odesa',
            lat=lat_Odesa,
            lon=lon_Odesa, alt=alt_Odesa, sou="Ukrainian early (pre-1850) historical weather observations, Skrynk et al.",
            link="https://doi.org/10.15407/uhmi.report.01", units="hPa", stat="point",
            metaHead="orig_p=in | PGC=Y | PTC=Y",
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
    meta_ta = paste0("orig_ta_out=", round(ta_out,1)),
    
    Lmm = if_else(!is.na(`P, in`), `P, in` * 25.4, `P, R.s.l.` * 1.27), #  / 750.06 * 1000) / 33.8639  # hPa to mmHg
    
    # Pressure corrected to 0°C + gravity → hPa
    value = round(convert_pressure(p = Lmm, f = 1, lat = lat_Poltava, alt = alt_Poltava, atb = ta_out),2),
    p_diff = round(value - Lmm * 1.333, 2),

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
            link="https://doi.org/10.15407/uhmi.report.01", units="C", stat="point",
            metaHead="Hours correspond to 'at sunrise/midday/at sunset' | orig_ta=Reaumur", keep_na = F)

write_sef_f(Data=df.p.Poltava,
            outpath=outdir, outfile='Poltava_p_subdaily.tsv',
            cod='Poltava',
            variable='p',
            nam='Poltava',
            lat=lat_Poltava,
            lon='34.544722', alt='160', sou="Ukrainian early (pre-1850) historical weather observations, Skrynk et al.",
            link="https://doi.org/10.15407/uhmi.report.01", units="hPa", stat="point",
            metaHead="Hours correspond to 'at sunrise/midday/at sunset' | PGC=Y | PTC=Y",
            meta=df.p.Poltava$meta, keep_na = F)
