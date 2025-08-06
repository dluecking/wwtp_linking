# Author: dlu @ veelab
# Version: 2025-04-22

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)
library(googlesheets4)
library(cowplot)
library(patchwork)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# read info ---------------------------------------------------------------


# URL of the Google Sheet
sheet_url <- "https://docs.google.com/spreadsheets/d/1QLNiqSt0XOS4xVPAeZAppwVjjjPIKdEE6w6f2_Qm55c"

# Read the specified sheet and convert to data.table
gv_data <- read_sheet(sheet_url, sheet = "Final GVs overview") %>% 
  filter(personal_assessment_order != "remove")


# for ViBioM poster -------------------------------------------------------

gv_data <- as.data.frame(gv_data)
gv_data <- gv_data %>% 
  mutate(gc = as.numeric(gv_data$gc))


p1 <- ggplot(gv_data, aes(x = personal_assessment_order, y = length, fill = personal_assessment_order)) +
  geom_boxplot(alpha = 0.9) +
  xlab("") +
  ylab("Length [bp]") +
  theme_cowplot() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 9),
        legend.position = "none") +
  scale_fill_manual(
    name = "Category",
    values = c(
      "Yaravirales" = "#1b9e77",
      "Asfuvirales" = "#d95f02",
      "Imitervirales" = "#66a61e",
      "Pandoravirales" = "#e7298a",
      "Pimascovirales" = "#7570b3",
      "unknown" = "grey"
    )
  ) +
  scale_y_continuous(labels = scales::scientific)


p2 <- ggplot(gv_data, aes(x = personal_assessment_order, y = gc, fill = personal_assessment_order)) +
  geom_boxplot(alpha = 0.9) +
  xlab("") +
  ylab("GC [%]") +
  theme_cowplot() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
  scale_fill_manual(
    name = "Category",
    values = c(
      "Yaravirales" = "#1b9e77",
      "Asfuvirales" = "#d95f02",
      "Imitervirales" = "#66a61e",
      "Pandoravirales" = "#e7298a",
      "Pimascovirales" = "#7570b3",
      "unknown" = "grey"
    )
  )

p_composite <- p1 / p2
p_composite

ggsave(plot = p_composite, file = "final/gv_stats.pdf", height = 5, width = 4)
ggsave(plot = p_composite, file = "final/gv_stats.svg", height = 5, width = 4)
ggsave(plot = p_composite, file = "final/gv_stats.png", height = 5, width = 4)








# # old ---------------------------------------------------------------------
# 
# 
# # plot 1 completeness -----------------------------------------------------
# 
# # Summarize the counts per taxonomic order and completeness status
# # Make sure we're working with a data.table
# setDT(gv_data)
# 
# # Count completeness per taxonomic order
# completeness_summary <- gv_data[, .N, by = .(personal_assessment_order, completeness)]
# 
# # Plot
# ggplot(completeness_summary, aes(x = personal_assessment_order, y = N, fill = completeness)) +
#   geom_bar(stat = "identity", position = "fill", color = "black") +  # proportions instead of raw count
#   theme_minimal() +
#   labs(x = "Taxonomic Order", y = "Proportion of Genomes", fill = "Completeness") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   theme_cowplot() +
#   coord_flip() +
#   xlab("") +
#   scale_y_continuous(labels = scales::percent) +  # show y-axis as percentage
#   scale_fill_manual(values = c(
#     "complete" = "#1b7837",           # dark green
#     "likely complete" = "#a6dba0",    # light green
#     "incomplete" = "#636363",         # dark grey
#     "likely incomplete" = "#bdbdbd",  # light grey
#     "ask alejandro" = "#636363"       # grey again
#   ))
# 
# ggsave(file = "GV_completeness_plot.png", plot = last_plot(), height = 4, width = 8)
# 
# # plot 2: length, GC and coding density -----------------------------------
# 
# gv_data <- as.data.frame(gv_data)
# gv_data <- gv_data %>% 
#   mutate(gc = as.numeric(gv_data$gc))
# 
# p1 <- ggplot(gv_data, aes(x = personal_assessment_order, y = gc)) +
#   geom_boxplot() +
#   coord_flip() +
#   xlab("") +
#   ylab("GC [%]") +
#   theme_cowplot()
# 
# p2 <- ggplot(gv_data, aes(x = personal_assessment_order, y = length)) +
#   geom_boxplot() +
#   coord_flip() +
#   xlab("") +
#   ylab("Length [bp]") +
#   theme_cowplot() +
#   theme(
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank()
#   )
# 
# p3 <- p1 + p2
# ggsave(file = "GV_gc_and_length_plot.png", plot = p3, height = 4, width = 8)
