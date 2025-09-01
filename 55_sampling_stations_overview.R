# Author: dlu @ veelab
# Version: 2025-08-29

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)
library(googlesheets4)
library(tidyr)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



# load data ---------------------------------------------------------------

sheet_url <- "https://docs.google.com/spreadsheets/d/113hSsqFV73bfdHTs5WoFFQU1DvnAcV5kPYmtanhROsY/edit#gid=1190755110"
sheet_name <- "GV, PLV, VPH presence in each sample"
data <- read_sheet(sheet_url, sheet = sheet_name)

# Data transformation for plotting
data_long <- data %>%
  rename(
    Station = `Sampling Station`,
    VPH = `Number of VPHs`,
    PLV = `Number of PLVs`,
    NCV = NCVs
  ) %>%
  pivot_longer(
    cols = c(VPH, PLV, NCV),
    names_to = "VirusType",
    values_to = "Count"
  )

# make names shorter
data_long$Station <- str_remove(data_long$Station, pattern = "\\_.*$")


custom_colors <- c(
  "VPH" = "goldenrod1",
  "PLV" = "hotpink",
  "NCV" = "steelblue"
)

# Create the ggplot2 plot
# We use geom_point to create the circles, with size and color aesthetics mapped
# to the 'Count' and 'VirusType' columns.
presence_plot <- ggplot(data_long, aes(x = Station, y = VirusType)) +
  geom_point(
    aes(size = ifelse(Count == 0, NA, Count), fill = VirusType), # hide 0 counts
    color = "black",
    shape = 21,
    stroke = 1.5,
    na.rm = TRUE
  ) +
  scale_size_continuous(range = c(2, 10)) +
  scale_fill_manual(values = custom_colors, guide = "none") + # remove VirusType legend
  labs(
    x = "Sampling Station",
    y = NULL, # remove y-axis label
    size = "Number of Viruses"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

presence_plot


ggsave(plot = presence_plot, file = "final/NCV_VPH_PLV_per_station.png", width = 7, height = 2.5)
ggsave(plot = presence_plot, file = "final/NCV_VPH_PLV_per_station.pdf", width = 7, height = 2.5)
ggsave(plot = presence_plot, file = "final/NCV_VPH_PLV_per_station.svg", width = 7, height = 2.5)

