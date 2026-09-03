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
sheet_url <- "https://docs.google.com/spreadsheets/d/159KaIRGjUGnN8uIJVsG42HIvWQ8MG2dysxNbNqw0VBg/edit?gid=697074736#gid=697074736"

# Read the specified sheet and convert to data.table
gv_data <- read_sheet(sheet_url, sheet = "Table S1")



# viz ---------------------------------------------------------------------

# colors
order_colors <- c(
  '"Yaravirales"' = "#0B3D5C",
  "Asfuvirales" = "#F2DACEFF",
  "Imitervirales" = "#b4934b",
  '"Pandoravirales"' = "#CCFDFFFF",
  "Pimascovirales" = "#00d6ed",
  "unknown" = "grey80"
)

gv_data <- as.data.frame(gv_data)


p1 <- ggplot(gv_data, aes(x = order, y = length, fill = order)) +
  geom_boxplot(alpha = 0.9) +
  xlab("") +
  ylab("Length [Mbp]") +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        legend.position = "none") +
  scale_fill_manual(
    name = "Category",
    values = order_colors
  ) +
  scale_y_continuous(
    labels = scales::label_number(scale = 1/1e6),
    breaks = c(100e3, 500e3, 1e6, 1500e3)
  ) +
  coord_flip()


p2 <- ggplot(gv_data, aes(x = order, y = `%GC`, fill = order)) +
  geom_boxplot(alpha = 0.9) +
  xlab("") +
  ylab("GC [%]") +
  theme_cowplot() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10)) +
  scale_fill_manual(
    name = "Category",
    values = order_colors
  ) +
  coord_flip()

# add ORFan information for each GV ---------------------------------------

annot_files <- paste0("intermediate/gv_annotations/", gv_data$shortname, "_polished_annotation.csv")
annot_df_list <- lapply(annot_files, fread)
names(annot_df_list) <- gv_data$shortname

gv_data$annotated_ORFs <- NA

for(i in 1:nrow(gv_data)){
  gv_data$annotated_ORFs[i] <- nrow(annot_df_list[[gv_data$shortname[i]]])
}

gv_data$ORFan_perc <- round((gv_data$ORFs - gv_data$annotated_ORFs) / gv_data$ORFs * 100, 2)

p3 <- ggplot(gv_data, aes(x = order, y = ORFan_perc, fill = order)) +
  geom_boxplot(alpha = 0.9) +
  xlab("") +
  ylab("ORFans [%]") +
  theme_cowplot() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        legend.position = "none") +
  scale_fill_manual(
    name = "Category",
    values = order_colors
  ) +
  ylim(c(33, 100)) +
  coord_flip()




p_composite <- p2 + p1 + p3 
p_composite

ggsave(plot = p_composite, file = "final/gv_stats.pdf", height = 3, width = 8)
ggsave(plot = p_composite, file = "final/gv_stats.svg", height = 3, width = 8)
ggsave(plot = p_composite, file = "final/gv_stats.png", height = 3, width = 8)



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
# completeness_summary <- gv_data[, .N, by = .(order, completeness)]
# 
# # Plot
# ggplot(completeness_summary, aes(x = order, y = N, fill = completeness)) +
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
# p1 <- ggplot(gv_data, aes(x = order, y = gc)) +
#   geom_boxplot() +
#   coord_flip() +
#   xlab("") +
#   ylab("GC [%]") +
#   theme_cowplot()
# 
# p2 <- ggplot(gv_data, aes(x = order, y = length)) +
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
