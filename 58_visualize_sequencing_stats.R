# Author: dlu @ veelab
# Version: 2025-10-23

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)
library(scales)
library(dplyr)
library(patchwork)
library(lubridate)
library(cowplot)
library(data.table)
library(tidyverse)


# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



# load data ---------------------------------------------------------------

data <- fread("helperfiles/supp_data_2_from_singleton_2021.csv")


# visualize proportions ---------------------------------------------------


# 1. Reshape the data from wide to long format
# This is crucial for grouping the bars in ggplot2
data_long <- data %>%
  # Select the identifier and the two proportion columns
  select(WWTP, proportion_illumina_trimmed_reads_mapped_to_assembly, proportion_of_reads_mapped_NP) %>%
  
  # Pivot the proportion columns into two new columns: Read_Type (the former column name)
  # and Proportion_Mapped (the value)
  pivot_longer(
    cols = starts_with("proportion"),
    names_to = "Read_Type",
    values_to = "Proportion_Mapped"
  ) %>%
  
  # Optionally, clean up the 'Read_Type' names for better plot readability
  mutate(
    Read_Type = recode(
      Read_Type,
      "proportion_illumina_trimmed_reads_mapped_to_assembly" = "Illumina Trimmed Reads",
      "proportion_of_reads_mapped_NP" = "Nanopore Trimmed Reads"
    )
  )

# 2. Create the grouped bar plot
ggplot(data_long, aes(x = WWTP, y = Proportion_Mapped, fill = Read_Type)) +
  geom_col(position = "dodge", colour = "black", alpha = 0.8) +
  
  # Add informative labels and title
  labs(
    title = "Proportion of Trimmed Reads Mapped",
    x = "Wastewater Treatment Plant (WWTP)",
    y = "Proportion of Reads Mapped",
    fill = "Read Type"
  ) +
  ylim(c(0, 1)) +  
  # Format the y-axis to display as a percentage (optional, but helpful for proportions)
  scale_y_continuous(labels = scales::percent) +
  theme_cowplot() +

  
  # Adjust the theme to rotate the x-axis labels
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom" 
  )


ggsave(plot = last_plot(), height = 5, width = 8, 
       filename = "final/mapping_percentages.png")


# visualize sequencing depth ----------------------------------------------

GBP_SCALE <- 1e9

illumina_data_plot <- bp_data %>%
  mutate(
    Illumina_Removed_bp = total_Illumina_data_bp - total_Illumina_data_trimmed_bp
  ) %>%
  select(
    WWTP, 
    Remaining_bp = total_Illumina_data_trimmed_bp, 
    Removed_bp = Illumina_Removed_bp
  ) %>%
  pivot_longer(
    cols = c(Remaining_bp, Removed_bp),
    names_to = "Status",
    values_to = "BP_Count"
  ) %>%
  mutate(BP_Gbp = BP_Count / GBP_SCALE) %>%
  mutate(
    Status = factor(Status, levels = c("Removed_bp", "Remaining_bp"))
  )

plot_illumina <- ggplot(illumina_data_plot, aes(x = WWTP, y = BP_Gbp, fill = Status)) +
  geom_col(position = "stack", color = "black", linewidth = 0.2) +
  scale_fill_manual(
    values = c("Remaining_bp" = "#4c72b0", "Removed_bp" = "#ff9900"),
    labels = c("Remaining_bp" = "Remaining (Trimmed)", "Removed_bp" = "Removed by Trimming")
  ) +
  labs(
    title = "Illumina Reads: Trimmed vs. Removed BP",
    x = "",
    y = "Base Pairs (Gbp)",
    fill = "Read Status"
  ) +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

print(plot_illumina)

nanopore_data_plot <- bp_data %>%
  mutate(
    NP_Removed_bp = total_nanopore_data_bp - totalbp_NP_data_trimmed_bp
  ) %>%
  select(
    WWTP, 
    Remaining_bp = totalbp_NP_data_trimmed_bp, 
    Removed_bp = NP_Removed_bp
  ) %>%
  pivot_longer(
    cols = c(Remaining_bp, Removed_bp),
    names_to = "Status",
    values_to = "BP_Count"
  ) %>%
  mutate(BP_Gbp = BP_Count / GBP_SCALE) %>%
  mutate(
    Status = factor(Status, levels = c("Removed_bp", "Remaining_bp"))
  )

plot_nanopore <- ggplot(nanopore_data_plot, aes(x = WWTP, y = BP_Gbp, fill = Status)) +
  geom_col(position = "stack", color = "black", linewidth = 0.2) +
  scale_fill_manual(
    values = c("Remaining_bp" = "#4c72b0", "Removed_bp" = "#ff9900"),
    labels = c("Remaining_bp" = "Remaining (Trimmed)", "Removed_bp" = "Removed by Trimming")
  ) +
  labs(
    title = "Nanopore Reads: Trimmed vs. Removed BP",
    x = "",
    y = "Base Pairs (Gbp)",
    fill = "Read Status"
  ) +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

print(plot_nanopore)


ggsave(plot = plot_illumina, height = 5, width = 8, 
       filename = "final/gbp_sequencing_stats_ill.png")

ggsave(plot = plot_nanopore, height = 5, width = 8, 
       filename = "final/gbp_sequencing_stats_np.png")
