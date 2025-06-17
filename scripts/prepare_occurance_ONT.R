# Author: dlu @ veelab
# Version: 2025-06-16

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)
library(tidyr)
library(Hmisc)

# set working directory
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("/lisc/scratch/dome/luecking/projects/wwtp_linking")


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


# calculate spearman ------------------------------------------------------

numeric_matrix <- occurance_df %>%
  tibble::column_to_rownames("contig_id") %>%
  as.matrix()


# new try wih hmisc
transposed_matrix <- t(numeric_matrix)
print("Dimensions of matrix for correlation (samples x contigs):")
print(dim(transposed_matrix))

# make sure its matrix
transposed_matrix_mat <- as.matrix(transposed_matrix)
# calc spearman correlation
cor_results <- rcorr(transposed_matrix_mat, type = "spearman")

cor_matrix <- cor_results$r # Correlation matrix
# cor_pvalues <- cor_results$P # P-value matrix

# remove duplicate information
cor_matrix[upper.tri(cor_matrix, diag = TRUE)] <- NA
print(head(cor_matrix[, 1:4]))


# # Transpose the matrix so contigs are columns and samples are rows
# transposed_matrix <- t(numeric_matrix)
# 
# print("Dimensions of matrix for correlation (samples x contigs):")
# print(dim(transposed_matrix))
# 
# # calc spearman
# cor_matrix <- cor(transposed_matrix, method = "spearman", use = "pairwise.complete.obs") # Or "complete.obs"
# # how does it look?
# print(head(cor_matrix[, 1:4]))

# remove duplicate information
cor_matrix[upper.tri(cor_matrix, diag = TRUE)] <- NA

# prepare for from to network
edge_list <- as.data.frame(as.table(cor_matrix)) %>%
  na.omit() %>%
  rename(from = Var1, to = Var2, spearman_ont = Freq)

fwrite(edge_list, "intermediate/network/occurance_ont.csv")