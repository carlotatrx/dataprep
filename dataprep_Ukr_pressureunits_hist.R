library(readr); library(dplyr); library(stringr); library(purrr); library(tidyr)

'''
  analyze historical Ukrainian pressure in SEF format to verify 
  the original measurement units used before metrication.
  1. Read subdaily SEF pressure files from Ukraine (needs to be specified).
  2. Extract the original pressure readings from `Meta` column.
  3. Convert modern hPa values back to mmHg and computes the implied mm per inch or
     mm per R.s.l. for each observation.
  4. Compare these implied conversion factors against known historical definitions:
       - British/Russian inch = 25.40 mm
       - Prussian inch        = 26.15 mm
       - Rhineland inch       = 26.34 mm
       - Russian semi-line    = 1.27 mm
  5. Plot the distribution of implied mm per inch in a nice histogram.
  
'''

dir <- "/home/ccorbella/scratch2_symboliclink/files/station_timeseries_preprocessed/"

# helper
hPa_to_mmHg <- function(h) h / 1.333224

# read one SEF (auto-skip metadata, read the tabular part)
read_sef_quick <- function(fp){
  lines <- read_lines(fp)
  head_row <- which(str_detect(lines, "^Year\\tMonth\\tDay\\tHour"))
  stopifnot(length(head_row)==1)
  read_tsv(fp, skip = head_row-1, show_col_types = FALSE) %>%
    mutate(file = basename(fp))
}

# load all pressure SEFs
names <- c('Kamyanets','Kharkiv','Dnipro','Kyiv','Odesa','Lugansk','Kherson')
files <- paste0(dir, names, '_p_subdaily_qc.tsv')

dat <- map_dfr(files, read_sef_quick)

# parse station name, units etc. from the header (optional but nice)
get_header_field <- function(fp, key){
  ln <- read_lines(fp, n_max = 60)
  hit <- ln[str_detect(ln, paste0("^", key, "\\t"))]
  if(!length(hit)) return(NA_character_)
  str_trim(str_remove(hit, paste0("^", key, "\\t")))
}
meta <- tibble(file = basename(files)) %>%
  mutate(Name  = map_chr(files, get_header_field, "Name"),
         Units = map_chr(files, get_header_field, "Units"))

dat <- dat %>% left_join(meta, by="file")

# extract original readings from "Meta" column
# supports "orig_p=29.9in", "orig_p=590Rsl", "orig_p=750mm", etc.
parsed <- dat %>%
  mutate(
    orig_p_str = str_extract(Meta, "(?<=orig_p=)[^|]+"),
    orig_num   = as.numeric(str_extract(orig_p_str, "[0-9]+\\.?[0-9]*")),
    orig_unit  = case_when(
      str_detect(orig_p_str, "in")  ~ "in",
      str_detect(orig_p_str, "R.s.l|R.s.l.") ~ "rsl",
      str_detect(orig_p_str, "mm")  ~ "mm",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(orig_unit), !is.na(orig_num))

# compute implied conversions
res <- parsed %>%
  transmute(
    Station = Name,
    Date = as.Date(sprintf("%04d-%02d-%02d", Year, Month, Day)),
    Value_hPa = Value,
    mmHg = hPa_to_mmHg(Value),
    orig_num, orig_unit
  ) %>%
  mutate(
    # implied mm per "inch" OR per R.s.l.
    mm_per_inch = ifelse(orig_unit=="in",  mmHg / orig_num, NA_real_),
    mm_per_rsl  = ifelse(orig_unit=="rsl", mmHg / orig_num, NA_real_)
  )

# --- SUMMARY: inches ---------------------------------------------------------
inch_stats <- res %>% filter(!is.na(mm_per_inch)) %>%
  group_by(Station) %>%
  summarize(
    n = n(),
    median_mm_per_inch = median(mm_per_inch, na.rm=TRUE),
    mean_mm_per_inch   = mean(mm_per_inch, na.rm=TRUE),
    sd_mm_per_inch     = sd(mm_per_inch, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    diff_vs_25_4   = median_mm_per_inch - 25.4,
    diff_vs_26_15  = median_mm_per_inch - 26.15,
    diff_vs_26_34  = median_mm_per_inch - 26.34
  )

print(inch_stats)

# --- SUMMARY: R.s.l. ---------------------------------------------------------
rsl_stats <- res %>% filter(!is.na(mm_per_rsl)) %>%
  group_by(Station) %>%
  summarize(
    n = n(),
    median_mm_per_rsl = median(mm_per_rsl, na.rm=TRUE),
    mean_mm_per_rsl   = mean(mm_per_rsl, na.rm=TRUE),
    sd_mm_per_rsl     = sd(mm_per_rsl, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  mutate(diff_vs_1_27 = median_mm_per_rsl - 1.27)

print(rsl_stats)

overall <- res %>% summarize(
  n_in   = sum(!is.na(mm_per_inch)),
  med_in = median(mm_per_inch, na.rm=TRUE),
  n_rsl  = sum(!is.na(mm_per_rsl)),
  med_rsl= median(mm_per_rsl,  na.rm=TRUE)
)
print(overall)

cat("\nBenchmarks:\n  British/Russian inch: 25.40 mm\n  Prussian inch: 26.15 mm\n  Rhineland inch: 26.34 mm\n  R.s.l.: 1.27 mm\n")

# safety check
stopifnot(abs(median(res$mm_per_inch, na.rm=TRUE) - 25.4) < 0.1)
stopifnot(abs(median(res$mm_per_rsl,  na.rm=TRUE) - 1.27) < 0.01)


# plot --------------------------------------------------------------------

library(ggplot2)

# Data for reference lines
inch_lines <- data.frame(
  x   = c(25.40, 26.15, 26.34),
  lbl = c("British/Russian inch", "Prussian inch", "Rhineland inch"),
  lty = c("solid", "dashed", "dashed"),
  y   = c(780,540,600)
)

# save plot
pdf('/home/ccorbella/scratch2_symboliclink/code/KF_assimilation/dataprep/image/Ukraine_calendar_shifts/inches_to_mmHg.pdf',
    width=8,height=7)


# Base histogram
p <- ggplot(res %>% filter(!is.na(mm_per_inch)),
            aes(mm_per_inch, fill=Station)) +
  geom_histogram(bins=60, alpha=0.5, position="identity") +
  labs(x = "Implied mm per inch", y = "Number of Observations")

# Add lines + labels + bigger text
p +
  geom_vline(data = inch_lines,
             aes(xintercept = x, linetype = lty),
             linewidth = 0.8, show.legend = FALSE) +
  scale_linetype_identity() +
  geom_text(data = inch_lines,
            aes(x = x - 0.01, y = y, label = lbl),
            angle = 90, vjust = -.2, hjust = 1, size = 5,
            inherit.aes = FALSE) +
  theme_minimal(base_size = 18) +
  theme(
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 14),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 14)
  )
dev.off()
