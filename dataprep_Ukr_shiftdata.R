library(dataresqc)
library(dplyr)
library(lubridate)
library(readr)

# Define input and output files
input_file <- "Dnipro_ta_subdaily.tsv"
output_file <- "Dnipro_ta_subdaily_corrected.tsv"

# Read SEF file
sef <- read_sef(input_file)

# Convert to data frame for easier manipulation
df <- sef$data

# Create full date column
df$full_date <- as.Date(with(df, sprintf("%04d-%02d-%02d", Year, Month, Day)))

# Define cutoff date
cutoff_date <- as.Date("1838-12-15")

# Apply date shift where needed
df <- df %>%
  mutate(
    apply_shift = full_date > cutoff_date,
    full_date = if_else(apply_shift, full_date + 12, full_date),
    Meta = if_else(apply_shift,
                   if_else(is.na(Meta) | Meta == "", "Δdate=+12d", paste(Meta, "Δdate=+12d", sep = ";")),
                   Meta)
  )

# Update Year, Month, Day
df$Year <- year(df$full_date)
df$Month <- month(df$full_date)
df$Day <- day(df$full_date)

# Drop helper columns
df$full_date <- NULL
df$apply_shift <- NULL

# Replace the data in the sef object
sef$data <- df

# Add metaHead flag
sef$metaHead$calendar_correction <- "Y"

# Write updated SEF file
write_sef(sef, output_file)
