library(dataresqc)
library(dplyr)
library(tidyr)
library(readxl)

source('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/helpfun.R')

indir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_orig/Ukraine/'
outdir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/'



# Dnipro -----------------------------------------------------------------
infile <- 'Dnipro.tsv'
df <- read.delim(paste0(indir, infile), header=T, sep='\t', stringsAsFactors = F,
                 skip=12)

meta <- read_meta_nonofficial(paste0(indir,infile))

# Create temperature dataframe
df.ta.Dnipro <- df %>%
  separate(Hour, into = c("Hour", "Minute"), sep = ":", convert = TRUE) %>%
  mutate(value = ifelse(T==-999.9, NA_real_,T * 1.25), meta = paste0('orig_ta=',T,"R | orig_time=",Time)) %>%
  select(Year, Month, Day, Hour, Minute, value, meta)

df.p.Dnipro <- df %>%
  separate(Hour, into = c("Hour", "Minute"), sep = ":", convert = TRUE) %>%
  mutate(
    P = ifelse(P==-999.9, NA_real_, P),
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
    ta_out = ifelse(T==-999.9, NA_real_, T * 1.25),
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
            link="https://doi.org/10.15407/uhmi.report.01", units='C', stat="point",
            meta=df.ta.Dnipro$meta, keep_na = F)

write_sef_f(Data=df.p.Dnipro,
            outpath=outdir, outfile="Dnipro_p_subdaily.tsv",
            cod=meta[["id"]],
            variable='p',
            nam=meta[["name"]],
            lat=meta[["lat"]],
            lon=meta[["lon"]], alt=meta[["alt"]], sou=meta[["source"]], metaHead = 'PGC=Y | PTC=Y',
            link="https://doi.org/10.15407/uhmi.report.01", units='hPa', stat="point",
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
  mutate(value = ifelse(T==-999.9, NA_real_, T * 1.25), meta=paste0('orig_ta=',T,'R')) %>%
  separate(Time, into = c("Hour", "Minute"), sep = ":", convert = TRUE) %>%
  select(Year, Month, Day, Hour, Minute, value, meta)

df.p.Kamyanets <- df %>%
  mutate(
    ta_out = ifelse(T==-999.9, NA_real_, T * 1.25),
    meta_ta = paste0("origa_ta_outside=", round(ta_out, 1)),
    P = ifelse(P==-999.9, NA_real_, P),
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
            lon=meta[["lon"]], alt=meta[["alt"]], link="https://doi.org/10.15407/uhmi.report.01",
            sou=meta[["source"]], units='C', stat="point",
            meta=df.ta.Kamyanets$meta, keep_na = F)

write_sef_f(Data=df.p.Kamyanets,
            outpath=outdir, outfile="Kamyanets_p_subdaily.tsv",
            cod=meta[["id"]],
            variable='p',
            nam=meta[["name"]],
            lat=meta[["lat"]],
            lon=meta[["lon"]], alt=meta[["alt"]], link="https://doi.org/10.15407/uhmi.report.01", metaHead = 'PGC=Y | PTC=Y',
            sou=meta[["source"]], units='hPa', stat="point",
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
        !is.na(value_B) ~ paste0("altern_", meta_B),
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

# sort df
df.ta.Kharkiv <- df.ta.Kharkiv %>%
  arrange(Year, Month, Day, Hour, Minute)

dfA.p.Kharkiv <- dfA %>%
  mutate(Lmm = ifelse(P==-999.9, NA_real_,P * 1.27),
         ta_out = ifelse(T==-999.9, NA_real_, T * 1.25),
         value = round(convert_pressure(p=Lmm, f=1, lat=as.numeric(meta['lat']), alt=as.numeric(meta['alt']), atb=ta_out),
                       2),
         meta = paste0("orig_p=",P,"R.s.l.")
  ) %>%     # R.s.l. → hPa
  select(Year, Month, Day, Hour, Minute, value, meta)

dfB.p.Kharkiv <- dfB %>%
  rename('T'='T, R', 'P'='P, R.s.l.') %>%
  mutate(Lmm = ifelse(P==-999.9, NA_real_,P * 1.27),
         ta_out = ifelse(T==-999.9, NA_real_, T * 1.25),
         value = round(convert_pressure(p=Lmm, f=1, lat=as.numeric(meta['lat']), alt=as.numeric(meta['alt']), atb=ta_out),
                       2),
         meta = paste0("orig_p=",P,"R.s.l.")
  ) %>%     # R.s.l. → hPa
  select(Year, Month, Day, Hour, Minute, value, meta)

df.p.Kharkiv <- merge.Kharkiv(dfA.p.Kharkiv, dfB.p.Kharkiv, var='p')

df.p.Kharkiv <- df.p.Kharkiv %>%
  arrange(Year, Month, Day, Hour, Minute)

write_sef_f(Data=df.ta.Kharkiv,
            outpath=outdir, outfile="Kharkiv_ta_subdaily.tsv",
            cod=ifelse(meta[['id']]=="", 'Kharkiv', meta[['id']]),
            variable='ta',
            nam=meta[["name"]],
            lat=meta[["lat"]],
            lon=meta[["lon"]], alt=meta[["alt"]], sou=meta[["source"]],
            link=meta[["link"]], units='C', stat="point",
            meta=df.ta.Kharkiv$meta, keep_na = F)

write_sef_f(Data=df.p.Kharkiv,
            outpath=outdir, outfile="Kharkiv_p_subdaily.tsv",
            cod=ifelse(meta[['id']]=="", 'Kharkiv', meta[['id']]),
            variable='p',
            nam=meta[["name"]],
            lat=meta[["lat"]],
            lon=meta[["lon"]], alt=meta[["alt"]], sou=meta[["source"]], metaHead = 'PGC=Y | PTC=Y',
            link=meta[["link"]], units='hPa', stat="point",
            meta=df.p.Kharkiv$meta, keep_na = F)



# Kherson -----------------------------------------------------------------

alt_Kherson=47
lat_Kherson=46.73833

filepath <- paste0(indir, 'Kherson_forR.xlsx')
sheets <- setdiff(excel_sheets(filepath),"Meta")

df <- bind_rows(
  lapply(sheets, function(sheet) {
    df <- tryCatch(
      {
        suppressMessages(suppressWarnings(
          read_excel(filepath, sheet = sheet, range = cell_cols(1:8), .name_repair = "minimal")
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
    # meta_ta = case_when(    # prepare meta column
    #   !is.na(ta) ~ paste0("orig_ta=", round(ta_used,1)),
    #   !is.na(ta_onN) ~ paste0("orig_ta='TonN',", round(ta_used, 1)),
    #   !is.na(ta_onS) ~ paste0("orig_ta='TonS',", round(ta_used, 1)),
    #   TRUE ~ "orig_ta_unknown"
    # ),
    # 
    # corrected P in hPa
    value = round(convert_pressure(p=`P, in` * 25.4, f=1, lat=lat_Kherson, alt=alt_Kherson, atb=ta_used),2),
    

    meta_p= paste(meta, " | orig_p=",`P, in`,"in", sep="")
  ) %>%
  select(Year, Month, Day, Hour, Minute, value, meta_p)

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
            lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[['Source']],
            meta = df.ta.Kherson$meta_ta,
            link=meta[["Link"]], units='C', stat="point", keep_na = F)

write_sef_f(Data=df.p.Kherson,
            outpath=outdir, outfile='Kherson_p_subdaily.tsv',
            cod=meta[["ID"]],
            variable='p',
            nam=meta[["Name"]],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[['Source']],
            metaHead = "PGC=Y | PTC=Y",
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
        Minute = ifelse(grepl(":", Time), as.integer(sub(".*:", "", Time)), NA_integer_),
        meta = case_when(
          Hour == 6 ~ "orig_time=morning",
          Hour == 14 ~ "orig_time=midday",
          Hour == 22 ~ "orig_time=evening",
          Hour == 10 ~ "orig_time=morning"
          )
        )
        
    # rearrange and drop time col
    df <- df %>%
      select(Year, Month, Day, Hour, Minute, everything(), -Time)
    df
  })
)


df.ta.Kyiv <- df %>%
  mutate(value = `T, R` * 1.25, meta_ta=paste0(meta, " | orig_ta=", `T, R`,"R")) %>%
  select(Year, Month, Day, Hour, Minute, value, meta_ta)

df.ta.Kyiv <- as.data.frame(df.ta.Kyiv)

df.p.Kyiv <- df %>%
  mutate(
    ta = `T, R`  * 1.25,    # outside temperature in C
    
    # original pressure in mmHg
    Lmm = `P, R.s.l.` * 1.27, # R.s.l. -> mmHg ## / 750.06 * 1000) / 33.8639

    # corrected P in hPa
    value = round(convert_pressure(p=Lmm, f=1, lat=lat_Kyiv, alt=alt_Kyiv, atb=ta),2),

    meta_p= paste(meta, " | orig_p=", `P, R.s.l.`,"R.s.l.",sep="")
  ) %>%
  select(Year, Month, Day, Hour, Minute, value, meta_p)

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

write_sef_f(Data=df.ta.Kyiv,
            outpath=outdir, outfile='Kyiv_ta_subdaily.tsv',
            cod=meta[["ID"]],
            variable='ta',
            nam=meta[["Name"]],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]],
            meta=df.ta.Kyiv$meta,
            link=meta[["Link"]], units='C', stat="point", keep_na = F)

write_sef_f(Data=df.p.Kyiv,
            outpath=outdir, outfile='Kyiv_p_subdaily.tsv',
            cod=meta[["ID"]],
            variable='p',
            nam=meta[["Name"]],
            lat=meta[["Lat"]],
            lon=meta[["Lon"]], alt=meta[["Alt"]], sou=meta[["Source"]],
            metaHead = "PGC=Y | PTC=Y",
            link=meta[["Link"]], units='hPa', stat="point",
            meta=df.p.Kyiv$meta, keep_na = F)

# Lugansk --------------------------------------------------------------------

alt_Lugansk=49
lat_Lugansk=48.565556

filepath <- paste0(indir, 'Lugansk_forR.xlsx')
sheets <- setdiff(excel_sheets(filepath),"Meta")

df <- bind_rows(
  lapply(sheets, function(sheet) {
    df <- tryCatch(
        {
          suppressMessages(suppressWarnings(
            read_excel(filepath, sheet = sheet, range = cell_cols(1:7), .name_repair = "minimal")
          ))
        },
        warning = function(w) {
          message("⚠️  Warning in sheet: ", sheet, " → ", conditionMessage(w))
          suppressMessages(suppressWarnings(read_excel(file_path, sheet = sheet, range = cell_cols(1:7), .name_repair = "minimal")))
        }
    )
    
    # Ensure required columns exist
    if (!"TimeA" %in% names(df)) {
      df$TimeA <- NA_character_
    }
    if (!"P, in" %in% names(df)) {
      df$`P, in` <- NA_real_
    }
    if (!"P, R.s.l." %in% names(df)) {
      df$`P, R.s.l.` <- NA_real_
    }
    if (!"Pt, R" %in% names(df)) {
      df$`Pt, R` <- NA_real_
    }
    df <- df[, intersect(names(df), c("Year", "Month", "Day", "TimeA", "Time", "T, R", "P, in", "P, R.s.l.", "Pt, R"))]
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
  mutate(value = ifelse(`T, R`==-999.9, NA_real_, `T, R` * 1.25),
         meta = if_else(
           is.na(TimeA),
           paste0('orig_ta=',`T, R`,"R"),
           paste0('orig_time=', TimeA," | orig_ta=",`T, R`,"R")
           )) %>%
  select(Year, Month, Day, Hour, Minute, value, meta)

df.ta.Lugansk <- as.data.frame(df.ta.Lugansk)

df.p.Lugansk <- df %>%
  
  # Flag special case: data from 1837-12-09 onwards are actually in inches
  mutate(
    is_inch_hack = (Year == 1837 & Month == 12 & Day >= 9)
  ) %>%
  
  # Check if there are any rows where both "P, in" and "P, R.s.l." are not NA (where we have obs in both units)
  filter(!is.na(`P, in`) | !is.na(`P, R.s.l.`) | is_inch_hack) %>%
  mutate(
    ta_bar = `Pt, R` * 1.25, # barometer Temperature in C
    ta_out = `T, R`  * 1.25, # outside temperature in C
    ta_used = case_when(     # use barometer Ta when available, outside when not
      !is.na(ta_bar) ~ ta_bar,
      TRUE ~ ta_out
    ),

    # original pressure in mmHg
    Lmm = case_when(
      !is.na(`P, in`) ~ `P, in` * 25.4, # convert inHg → mmHg
      is_inch_hack ~ `P, R.s.l.` * 25.4, # for the vals in 1837 that are recorded in inches
      !is.na(`P, R.s.l.`) ~ (`P, R.s.l.` * 1.27) # R.s.l. -> mmHg ## / 750.06 * 1000) / 33.8639
    ),
    
    # corrected P in hPa
    value = round(convert_pressure(p=Lmm, f=1, lat=lat_Lugansk, alt=alt_Lugansk, atb=ta_used),2),
    
    meta_orig = case_when(
      !is.na(`P, in`) ~ paste0("orig_p=", `P, in`,"in"),
      is_inch_hack ~ paste0("orig_p=", `P, in`, "in"),
      !is.na(`P, R.s.l.`) ~ paste0("orig_p=", `P, R.s.l.`, "R.s.l.")
    ),
    
    meta = if_else(
      is.na(TimeA),
      meta_orig,
      paste0('orig_time=', TimeA," | ", meta_orig)
    )
    
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
            meta=df.ta.Lugansk$meta,
            metaHead="station relocation in 1843", keep_na = F)

write_sef_f(Data=df.p.Lugansk,
            outpath=outdir, outfile='Lugansk_p_subdaily.tsv',
            cod='Lugansk',
            variable='p',
            nam='Lugansk',
            lat='48.565556',
            lon='39.2275', alt='59', sou="Ukrainian early (pre-1850) historical weather observations, Skrynk et al.",
            link="https://doi.org/10.15407/uhmi.report.01", units="hPa", stat="point",
            metaHead="station relocation in 1843 | PGC=Y | PTC=Y",
            meta=df.p.Lugansk$meta, keep_na = F)

# Odessa --------------------------------------------------------------------
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
  mutate(value = ifelse(`T, R`==-999.9, NA_real_, `T, R` * 1.25),
         meta = paste0("orig_ta=", `T, R`, "R")) %>%
  select(Year, Month, Day, Hour, Minute, value, meta)

df.ta.Odesa <- as.data.frame(df.ta.Odesa)

df.p.Odesa <- df %>%
  mutate(
    ta_out = ifelse(`T, R`==-999.9, NA_real_, `T, R` * 1.25), # outside temperature in °C

    P = ifelse(`P, in`==-999.9, NA_real_, `P, in`*25.4),
    # Pressure corrected to 0°C + gravity → hPa
    value = round(convert_pressure(p = P, f = 1, lat = lat_Odesa, alt = alt_Odesa, atb = ta_out),2),

    meta = paste("orig_p=", `P, in`, "in", sep = "")
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
            meta=df.ta.Odesa$meta, keep_na = F)

write_sef_f(Data=df.p.Odesa,
            outpath=outdir, outfile='Odesa_p_subdaily.tsv',
            cod='Odesa',
            variable='p',
            nam='Odesa',
            lat=lat_Odesa,
            lon=lon_Odesa, alt=alt_Odesa, sou="Ukrainian early (pre-1850) historical weather observations, Skrynk et al.",
            link="https://doi.org/10.15407/uhmi.report.01", units="hPa", stat="point",
            metaHead="PGC=Y | PTC=Y",
            meta=df.p.Odesa$meta, keep_na = F)


# Poltava --------------------------------------------------------------------
lat_Poltava <- 49.609444
alt_Poltava <- 160
indir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_orig/Ukraine/'
outdir <- '/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/'

sheets <- setdiff(excel_sheets(paste0(indir,'Poltava_forR.xlsx')),"Meta")

df <- bind_rows(
  lapply(sheets, function(sheet) {
    df <- read_excel(paste0(indir,'Poltava_forR.xlsx'), sheet = sheet)

    # Ensure required columns exist
    if (!"TimeA" %in% names(df)) {
      df$TimeA <- NA_character_
    }
    if (!"P, in" %in% names(df)) {
      df$`P, in` <- NA_real_
    }
    if (!"P, R.s.l." %in% names(df)) {
      df$`P, R.s.l.` <- NA_real_
    }
    df <- df[, intersect(names(df), c("Year", "Month", "Day", "TimeA", "Time", "T, R", "P, in", "P, R.s.l."))]
    df <- df %>%
      mutate(across(c("T, R", "P, in", "P, R.s.l.", ), ~na_if(., -999.9))) # -999.9 to NA
    
    return(df)
  })
)

df.ta.Poltava <- df %>%
  separate(Time, into = c("Hour", "Minute"), sep = ":", convert = TRUE) %>%

  mutate(value = round(`T, R` * 1.25, 2),
         meta = if_else(
           is.na(TimeA),
           paste0('orig_ta=',`T, R`,"R"),
           paste0('orig_time=', TimeA," | orig_ta=",`T, R`,"R")
         )) %>%
  select(Year, Month, Day, Hour, Minute, value, meta)

df.ta.Poltava <- as.data.frame(df.ta.Poltava)

df.p.Poltava <- df %>%
  separate(Time, into = c("Hour", "Minute"), sep = ":", convert = TRUE) %>%

  # Check if there are any rows where both "P, in" and "P, R.s.l." are not NA (where we have obs in both units)
  # filter(!is.na(`P, in`) | !is.na(`P, R.s.l.`)) %>%
  mutate(
    ta_out = `T, R` * 1.25, # outside temperature in °C

    Lmm = if_else(!is.na(`P, in`), `P, in` * 25.4, `P, R.s.l.` * 1.27), #  / 750.06 * 1000) / 33.8639  # hPa to mmHg
    
    # Pressure corrected to 0°C + gravity → hPa
    value = round(convert_pressure(p = Lmm, f = 1, lat = lat_Poltava, alt = alt_Poltava, atb = ta_out),2),

    meta_orig = case_when(
      !is.na(`P, in`) ~ paste0("orig_p=", `P, in`,"in"),
      !is.na(`P, R.s.l.`) ~ paste0("orig_p=", `P, R.s.l.`, "R.s.l.")
    ),
    
    meta = if_else(
      is.na(TimeA),
      meta_orig,
      paste0('orig_time=', TimeA," | ", meta_orig)
    )
  ) %>%
  select(Year, Month, Day, Hour, Minute, value, meta)

df.p.Poltava <- as.data.frame(df.p.Poltava)

write_sef_f(Data=df.ta.Poltava,
            outpath=outdir, outfile='Poltava_ta_subdaily.tsv',
            cod='Poltava',
            variable='ta',
            nam='Poltava',
            lat=lat_Poltava,
            lon='34.544722', alt=alt_Poltava,
            link="https://doi.org/10.15407/uhmi.report.01", units="C", stat="point",
            meta=df.ta.Poltava$meta, keep_na = F)

write_sef_f(Data=df.p.Poltava,
            outpath=outdir, outfile='Poltava_p_subdaily.tsv',
            cod='Poltava',
            variable='p',
            nam='Poltava',
            lat=lat_Poltava,
            lon='34.544722', alt='160', sou="Ukrainian early (pre-1850) historical weather observations, Skrynk et al.",
            link="https://doi.org/10.15407/uhmi.report.01", units="hPa", stat="point",
            metaHead="PGC=Y | PTC=Y",
            meta=df.p.Poltava$meta, keep_na = F)
