library(dataresqc)
library(ggplot2)
library(patchwork)

# Load homogenized SEF file
final_file <- "~/scratch2_symboliclink/code/KF_assimilation/dataprep/homogenization/KNMI_KNMI-42_Zwanenburg_17380101-18601231_p_hom.tsv"
x_final <- read_sef(final_file)

# Build Date column
x_final$Date <- as.Date(with(x_final, sprintf("%04d-%02d-%02d", Year, Month, Day)))

# Plot
ggplot(x_final, aes(x = Date, y = Value)) +
  geom_line(color = "blue") +
  labs(title = "Homogenized Pressure Time Series",
       x = "Date", y = "Pressure (hPa)") +
  theme_minimal()


# Load original SEF file
original_file <- "/home/ccorbella/scratch2_symboliclink/files/1807_USBstick/KNMI_KNMI-42_Zwanenburg_17380101-18601231_p_daily.tsv"
x_orig <- read_sef(original_file)

# Build Date column
x_orig$Date <- as.Date(with(x_orig, sprintf("%04d-%02d-%02d", Year, Month, Day)))

# Plot
ggplot(x_orig, aes(x = Date, y = Value)) +
  geom_line(color = "darkred") +
  labs(title = "Original Pressure Time Series",
       x = "Date", y = "Pressure (hPa)") +
  theme_minimal()

# Create Date column
x_orig$Date <- as.Date(with(x_orig, sprintf("%04d-%02d-%02d", Year, Month, Day)))
x_final$Date <- as.Date(with(x_final, sprintf("%04d-%02d-%02d", Year, Month, Day)))


# plots -------------------------------------------------------------------


# Filter from 1780 onward
x_orig <- subset(x_orig, Date >= as.Date("1780-01-01"))
x_final <- subset(x_final, Date >= as.Date("1780-01-01"))

# Merge for difference line
x_diff <- merge(x_orig[, c("Date", "Value")],
                x_final[, c("Date", "Value")],
                by = "Date", suffixes = c("_orig", "_hom"))
x_diff$diff <- x_diff$Value_orig - x_diff$Value_hom

# Common y-limits
ymin <- min(c(x_orig$Value, x_final$Value), na.rm = TRUE)
ymax <- max(c(x_orig$Value, x_final$Value), na.rm = TRUE)

# Top plot: Original + Correction line
p1 <- ggplot() +
  geom_line(data = x_orig, aes(x = Date, y = Value, group=1), na.rm=F, color = "dodgerblue4") +  # top plot
  geom_line(data = x_diff, aes(x = Date, y = 1007 + diff, group=1), na.rm=F, color = "magenta", size=1) +
  # annotate("text", x = as.Date("1840-01-01"), y = ymax - 5, label = "Applied Correction", color = "magenta", hjust = 0, size = 4.2, fontface = "italic") +
  
  annotate("segment", x = as.Date("1833-01-01"), xend = as.Date("1836-01-01"),
           y = ymax - 5, yend = ymax - 5, color = "magenta", size = 1) +
  annotate("text", x = as.Date("1837-01-01"), y = ymax - 5, 
           label = "Correction Applied", color = "magenta", hjust = 0, size = 4.2, fontface = "italic") +
  annotate("text", x = as.Date("1780-06-01"), y = ymax - 2, 
           label = "Original", fontface = "bold", hjust = 0, size = 5) +
  labs(y = "Pressure [hPa]", x = NULL) +
  ylim(ymin, ymax) +
  theme_minimal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

  
# Bottom plot: Homogenized
p2 <- ggplot(x_final, aes(x = Date, y = Value)) +
  geom_line(color = "dodgerblue", alpha = 1) +  # Blue tone 2
  annotate("text", x = as.Date("1780-06-01"), y = ymax - 2, 
           label = "Homogenized", fontface = "bold", hjust = 0, size = 5) +

  labs(y = "Pressure [hPa]", x = NULL) +
  ylim(ymin, ymax) +
  theme_minimal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

# Combine with global title
zwa.plot <- (p1 / p2) + 
  plot_annotation(
    title = "Pressure Observations at Zwanenburg (Netherlands)",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.margin = margin(t = 5, r = 10, b = 10, l = 10)  # t = top margin in pts
    )
  )

#  plot_annotation(title = "Pressure Observations at Zwanenburg ", theme = theme(plot.title = element_text(hjust = 0.5, size=14)))

ggsave("/scratch2/ccorbella/code/KF_assimilation/dataprep/homogenization/zwanenburg_hom_plot2.png", zwa.plot, dpi = 600, width = 10, height = 6)

