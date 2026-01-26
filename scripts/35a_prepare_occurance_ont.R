# Author: dlu @ veelab
# Version: 2025-06-16

# start time
start_time <- Sys.time()

# Packages
library(dplyr)
library(data.table)
library(stringr)
library(tidyr)
library(foreach)
library(doParallel)


cat("[INFO] loaded libraries!\n")

# set working directory
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("/lisc/data/scratch/dome/willemsen/luecking/projects/wwtp_linking")
# setwd("/run/user/1000/gvfs/sftp:host=login01.lisc.univie.ac.at,user=luecking/lisc/home/user/luecking/luecking_scratch/projects/wwtp_linking")

# cores
n_cores <- 16
registerDoParallel(cores = n_cores)

# test mode?
TEST_MODE <- FALSE
if(TEST_MODE){
  cat("[DEBUG] TEST MODE IS ACTIVATED, ONLY 1k CONTIGS SELECTED!!!")
}


cat("[INFO] set cores to 16!\n")
cat("[TIME] Elapsed:", format(Sys.time() - start_time), "\n")

# functions ---------------------------------------------------------------

# Manual RMT threshold calculation
calculate_rmt_threshold <- function(cor_matrix, n_samples) {
  
  # 1. Calculate eigenvalues of correlation matrix
  eigenvalues <- eigen(cor_matrix, only.values = TRUE)$values
  
  # 2. Marchenko-Pastur parameters
  n_features <- nrow(cor_matrix)
  Q <- n_features / n_samples  # Ratio
  
  # 3. Theoretical maximum eigenvalue for random matrix
  # Assuming noise variance σ² = 1 (correlation matrix is standardized)
  lambda_max <- (1 + sqrt(Q))^2
  
  # 4. Find empirical threshold
  # Count eigenvalues above theoretical maximum
  n_significant <- sum(eigenvalues > lambda_max)
  
  # 5. Convert eigenvalue threshold to correlation threshold
  # This is approximate - you can also use the eigenvalue directly
  # For now, use the eigenvalue corresponding to the cutoff
  
  if(n_significant > 0) {
    threshold_eigenvalue <- sort(eigenvalues, decreasing = TRUE)[n_significant]
    
    # Rough conversion (this is approximate)
    # Better: use the eigenvalue to filter the correlation matrix directly
    threshold_corr <- sqrt(threshold_eigenvalue / n_samples)
  } else {
    threshold_corr <- 0.5  # Fallback
  }
  
  return(list(
    threshold = threshold_corr,
    lambda_max_theory = lambda_max,
    n_significant_eigenvalues = n_significant,
    eigenvalues = eigenvalues
  ))
}

# Parallel correlation
calculate_cor_parallel <- function(matrix, method = "spearman", n_cores = n_cores) {
  
  n_contigs <- ncol(matrix)
  
  # Calculate in parallel chunks
  cor_matrix <- foreach(i = 1:n_contigs, .combine = rbind, 
                        .packages = c("stats")) %dopar% {
                          cor(matrix[, i], matrix, method = method)
                        }
  
  rownames(cor_matrix) <- colnames(matrix)
  colnames(cor_matrix) <- colnames(matrix)
  
  return(cor_matrix)
}

clr_transform <- function(x) {
  log_x <- log(x + 1)
  log_x - mean(log_x)  # Center by geometric mean
}


# load sample data --------------------------------------------------------

sample_names <- fread("helperfiles/WWTP_Denmark_SRA_PRJNA629478.txt", header = F)
names(sample_names) <- c("sample", "ill", "ont")

sample_names$sample_short <- str_remove(sample_names$sample, "\\_.*$")

cat("[INFO] loaded helperfile and extracted sample names!\n")


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

cat("[INFO] loaded coverage files!\n")


################################################################################
# NOTE
# WE DID THIS BEFORE; BUT NOW WANT TO CIRCUMVENT THIS BY CALCULATING BETTER!

# # only keep connected contigs ---------------------------------------------
# 
# lcs_to_keep <- fread("intermediate/network/list_of_connected_lcs.txt", header = F)
# names(lcs_to_keep) <- c("contig_id")
# 
# # Filter the dataframe
# occurance_filtered <- occurance_df %>%
#   filter(
#     # Condition 1: Keep rows where col1 does NOT end with "_lc"
#     !grepl("_lc$", contig_id) |
#       # Condition 2: OR (if col1 DOES end with "_lc") keep it if it's found in contigs_to_keep
#       (grepl("_lc$", contig_id) & contig_id %in% lcs_to_keep$contig_id)
#   )
################################################################################


# remove contigs for which we have a weak signal --------------------------
# the idea is: if you are only present in very few samples, you create noise

X_threshold <- 1.0 # Replace with your desired minimum value (X)
Y_columns <- 3     # Replace with your desired minimum number of columns (Y)

# Assuming 'occurance_filtered' is your data frame

occurance_filtered <- occurance_df %>%
  rowwise() %>% # Process row by row
  mutate(
    # Count how many non-contig_id columns have a value > X_threshold
    count_above_X = sum(c_across(-contig_id) > X_threshold, na.rm = TRUE)
  ) %>%
  filter(count_above_X >= Y_columns) %>% # Keep rows where the count is at least Y_columns
  select(-count_above_X) # Remove the temporary count_above_X column if not needed

cat("[INFO] filtered occurance DF!\n")
cat("[TIME] Elapsed:", format(Sys.time() - start_time), "\n")


# test case ---------------------------------------------------------------

if(TEST_MODE) {
  cat("[TEST] Running in TEST MODE with subset of data\n")
  set.seed(42)
  test_contigs <- sample(occurance_filtered$contig_id, size = 1000)
  occurance_filtered <- occurance_filtered %>% 
    filter(contig_id %in% test_contigs)
  
  cat("[TEST] Using only", nrow(occurance_filtered), "contigs\n")
}


# calculate spearman ------------------------------------------------------

numeric_matrix <- occurance_filtered %>%
  tibble::column_to_rownames("contig_id") %>%
  as.matrix()

# clr transform
clr_matrix <- t(apply(numeric_matrix, 1, clr_transform))

# remove emtpy or inf values:
cat("[DEBUG] Checking CLR matrix...\n")
cat("  - Dimensions:", dim(clr_matrix), "\n")
cat("  - NA values:", sum(is.na(clr_matrix)), "\n")
cat("  - Inf values:", sum(is.infinite(clr_matrix)), "\n")

# Remove rows (contigs) that are all NA or have Inf
valid_rows <- apply(clr_matrix, 1, function(x) {
  !all(is.na(x)) && !any(is.infinite(x)) && var(x, na.rm = TRUE) > 0
})

cat("[INFO] Removing", sum(!valid_rows), "problematic contigs\n")
clr_matrix <- clr_matrix[valid_rows, ]

# Transpose the matrix so contigs are columns and samples are rows
transposed_matrix <- t(clr_matrix)

cat("[INFO] transformed into CLR matrix!\n")
cat("[TIME] Elapsed:", format(Sys.time() - start_time), "\n")

# calc spearman
cor_matrix <- calculate_cor_parallel(transposed_matrix, 
                                     method = "spearman", 
                                     n_cores = n_cores)

cat("[INFO] calculated cor matrix!\n")
cat("[TIME] Elapsed:", format(Sys.time() - start_time), "\n")

# calculate RMT value
rmt_result <- calculate_rmt_threshold(cor_matrix, n_samples = ncol(transposed_matrix))
cat("[INFO] calculated RMT threshold!\n")
cat("[TIME] Elapsed:", format(Sys.time() - start_time), "\n")

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

# filter edgle_list based on RMT value
cat("[INFO] RMT threshold:", rmt_result$threshold, "\n")
cat("[INFO] Number of significant eigenvalues:", rmt_result$n_significant_eigenvalues, "\n")
cat(paste0("[INFO] Size of edge_list BEFORE filtereing above RMT threshold: ", nrow(edge_list), "\n"))

# Filter edges
edge_list_rmt <- edge_list %>% 
  filter(absolut_spearman >= rmt_result$threshold)
cat(paste0("[INFO] Size of edge_list BEFORE filtereing above RMT threshold: ", nrow(edge_list_rmt)), "\n")

# save that to a file -----------------------------------------------------

fwrite(edge_list_rmt, "intermediate/network/occurance_ont.csv")
cat("[INFO] final edgelist saved to: intermediate/network/occurance_ont.csv")
cat("[TIME] Elapsed:", format(Sys.time() - start_time), "\n")