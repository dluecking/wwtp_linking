# Author: dlu @ veelab
# Version: 2025-06-16

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)
library(tidyr)
library(cowplot)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# load sample data --------------------------------------------------------

sample_names <- fread("helperfiles/WWTP_Denmark_SRA_PRJNA629478.txt", header = F)
names(sample_names) <- c("sample", "ill", "ont")

sample_names$sample_short <- str_remove(sample_names$sample, "\\_.*$")



# load occurance data -----------------------------------------------------

coverage_dir <- "intermediate/occurance/mapping/ont"
coverage_files <- list.files(coverage_dir, pattern = "coverage.txt", full.names = TRUE)

list_of_dfs <- list()

if (length(coverage_files) == 0) {
  message(paste0("No 'coverage.txt' files found in '", coverage_dir, "'."))
} else {
  for (file_path in coverage_files) {
    file_name <- basename(file_path)
    current_ont_id_from_file <- str_remove(file_name, pattern = "\\_mapped.*$")
    
    match_index <- match(current_ont_id_from_file, sample_names$ont)
    
    if (is.na(match_index) || length(match_index) == 0) {
      warning(paste0("No matching 'ont' ID found in 'sample_names' for file: ", file_name, ". Skipping."))
      next
    }
    
    sample_short_name <- sample_names$sample_short[match_index]
    
    if (length(sample_short_name) > 1) {
      warning(paste0("Multiple matching 'sample_short' names for '", current_ont_id_from_file, "'. Using the first: '", sample_short_name[1], "'."))
      sample_short_name <- sample_short_name[1]
    }
    
    tmp_df <- fread(file_path) %>%
      select(`#rname`, meandepth) %>%
      rename(contig_id = `#rname`)
    
    names(tmp_df)[names(tmp_df) == "meandepth"] <- sample_short_name
    tmp_df$contig_id <- as.character(tmp_df$contig_id)
    
    list_of_dfs[[sample_short_name]] <- tmp_df
  }
}

if (length(list_of_dfs) > 0) {
  occurance_df <- list_of_dfs[[1]]
  
  if (length(list_of_dfs) > 1) {
    for (i in 2:length(list_of_dfs)) {
      occurance_df <- merge(occurance_df, list_of_dfs[[i]], by = "contig_id", all.x = TRUE)
    }
  }
} else {
  occurance_df <- data.table("contig_id" = as.character())
}
rm(list_of_dfs, tmp_df, sample_names, current_ont_id_from_file, file_name, file_path, i, match_index, coverage_dir, coverage_files, sample_short_name)


# only keep connected contigs ---------------------------------------------

lcs_to_keep <- fread("intermediate/network/list_of_connected_lcs.txt", header = F)
names(lcs_to_keep) <- c("contig_id")

# Filter the dataframe
occurance_filtered <- occurance_df %>%
  filter(
    # Condition 1: Keep rows where col1 does NOT end with "_lc"
    !grepl("_lc$", contig_id) |
      # Condition 2: OR (if col1 DOES end with "_lc") keep it if it's found in contigs_to_keep
      (grepl("_lc$", contig_id) & contig_id %in% lcs_to_keep$contig_id)
  )




# calculate spearman ------------------------------------------------------

numeric_matrix <- occurance_filtered %>%
  tibble::column_to_rownames("contig_id") %>%
  as.matrix()

# Transpose the matrix so contigs are columns and samples are rows
transposed_matrix <- t(numeric_matrix)

print("Dimensions of matrix for correlation (samples x contigs):")
print(dim(transposed_matrix))

# calc spearman
cor_matrix <- cor(transposed_matrix, method = "spearman") 

# how does it look?
print(head(cor_matrix[, 1:4]))

# remove duplicate information
cor_matrix[upper.tri(cor_matrix, diag = TRUE)] <- NA

# prepare for from to network
edge_list <- as.data.frame(as.table(cor_matrix)) %>%
  na.omit() %>%
  rename(from = Var1, to = Var2, spearman_ont = Freq)

edge_list$absolut_spearman <- abs(edge_list$spearman_ont)
edge_list <- edge_list %>% 
  mutate(correlation = case_when(
    spearman_ont >= 0 ~ "positive",
    spearman_ont < 0 ~ "negative"
  ))


hist(x = edge_list$spearman_ont,
     freq = FALSE, # Plot densities instead of frequencies
     main = "Histogram with Density Curve",
     xlab = "Values",
     col = "lightblue",
     border = "black")
lines(density(edge_list$spearman_ont), col = "red", lwd = 2)

qqnorm(edge_list$spearman_ont,
       main = "Normal Q-Q Plot",
       xlab = "Theoretical Quantiles",
       ylab = "Sample Quantiles")
qqline(edge_list$spearman_ont, col = "blue", lwd = 2)





# save edgelist -----------------------------------------------------------

fwrite(edge_list %>% filter(absolut_spearman >= 0.425), "intermediate/network/occurance_ont.csv")


# # quick exploration of what different cutoffs mean ------------------------
# 
# cutoffs <- c(0.5, 0.6, 0.7, 0.8, 0.9)
# 
# df <- data.table()
# 
# for(CUTOFF in cutoffs){
#   a <- edge_list %>% 
#     filter(absolut_spearman >= CUTOFF)
#   all_contigs <- data.table(contig_id = unique(c(as.character(a$from), as.character(a$to))),
#                             type = "")
#   all_contigs <- all_contigs %>%
#     mutate(
#       type = case_when(
#         str_ends(contig_id, "vph") ~ "vph",
#         str_ends(contig_id, "lc")  ~ "lc",
#         str_ends(contig_id, "plv") ~ "plv",
#         TRUE ~ "gv"
#       )
#     )
#   
#   tmp_df <- as.data.table(table(all_contigs$type))
#   tmp_df$cutoff <- CUTOFF
#   
#   df <- rbind(df, tmp_df)
#   
# }
# 
# names(df) <- c("type", "n", "cutoff")
# 
# df$total_of_that_type <- 0 
# df$total_of_that_type[df$type == "vph"] <- 15
# df$total_of_that_type[df$type == "lc"] <- 4001
# df$total_of_that_type[df$type == "plv"] <- 14
# df$total_of_that_type[df$type == "gv"] <- 63
# 
# df$percent <- df$n / df$total_of_that_type * 100
# 
# ggplot(df, aes(x = cutoff, fill = type, y = percent)) +
#   geom_bar(stat = "identity", color = "black", alpha = 0.8, position = "dodge") +
#   theme_cowplot() +
#   ggtitle("[%] of contigs included in network, per spearman cutoff")
