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
cat("[TIME] Elapsed:", format(Sys.time() - start_time), "\n")


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


# prepare for network -----------------------------------------------------
# this part runs into an integer overflow, so we adjusted by running blocks:


# # remove duplicate information
# cor_matrix[upper.tri(cor_matrix, diag = TRUE)] <- NA
# 
# # prepare for from to network
# edge_list <- as.data.frame(as.table(cor_matrix)) %>%
#   na.omit() %>%
#   rename(from = Var1, to = Var2, spearman_ont = Freq)

# prepare for network -----------------------------------------------------

cat("[INFO] Converting correlation matrix to edge list (memory-efficient approach)...\n")

# Get dimensions
n_contigs <- nrow(cor_matrix)
contig_names <- rownames(cor_matrix)

# Pre-filter: only keep edges above threshold to save memory
SPEARMAN_CUTOFF <- 0.60

# Process in chunks to avoid memory issues
chunk_size <- 1000
edge_list_chunks <- list()

for(start_idx in seq(1, n_contigs, by = chunk_size)) {
  end_idx <- min(start_idx + chunk_size - 1, n_contigs)
  
  cat(sprintf("[INFO] Processing contigs %d to %d of %d...\n", 
              start_idx, end_idx, n_contigs))
  
  # Extract chunk of correlation matrix (rows = current chunk, cols = all contigs up to current row)
  chunk_mat <- cor_matrix[start_idx:end_idx, , drop = FALSE]
  
  # For each row in chunk, check all columns that come BEFORE this contig (lower triangle only)
  chunk_edges <- foreach(i = start_idx:end_idx, .combine = rbind) %do% {
    row_idx <- i - start_idx + 1
    
    # Only look at columns before this row (lower triangle to avoid duplicates)
    # This ensures we capture ALL pairwise comparisons exactly once
    if(i > 1) {
      cols_to_check <- 1:(i-1)  # All contigs with index < current contig
      correlations <- chunk_mat[row_idx, cols_to_check]
      abs_cors <- abs(correlations)
      
      # Pre-filter by threshold
      keep <- abs_cors >= SPEARMAN_CUTOFF
      
      if(any(keep)) {
        data.frame(
          from = contig_names[i],
          to = contig_names[cols_to_check[keep]],
          spearman_ont = correlations[keep],
          stringsAsFactors = FALSE
        )
      } else {
        NULL
      }
    } else {
      NULL
    }
  }
  
  if(!is.null(chunk_edges) && nrow(chunk_edges) > 0) {
    edge_list_chunks[[length(edge_list_chunks) + 1]] <- chunk_edges
  }
  
  # Clean up
  rm(chunk_mat, chunk_edges)
  gc()
}

# Combine all chunks
edge_list <- rbindlist(edge_list_chunks)

cat("[INFO] Created edge list with", nrow(edge_list), "edges\n")
cat("[TIME] Elapsed:", format(Sys.time() - start_time), "\n")


edge_list$absolut_spearman <- abs(edge_list$spearman_ont)
edge_list <- edge_list %>% 
  mutate(correlation = case_when(
    spearman_ont >= 0 ~ "positive",
    spearman_ont < 0 ~ "negative"
  ))


# info for later ----------------------------------------------------------

thresholds <- seq(0.1, 0.7, by = 0.05)

for(t in thresholds) {
  n <- sum(edge_list$absolut_spearman >= t)
  cat(sprintf("[DEBUG] Threshold %.2f: %10d edges (%.2f%% of max)\n", 
              t, n, n / (51767 * 51766 / 2) * 100))
}

SPEARMAN_CUTOFF <- 0.60

# filter edgle_list based on RMT value
cat("[INFO] SEPARMAN CUTOFF selected:", SPEARMAN_CUTOFF, "\n")
cat(paste0("[INFO] Size of edge_list BEFORE filtereing above threshold: ", nrow(edge_list), "\n"))

# Filter edges
edge_list_rmt <- edge_list %>% 
  filter(absolut_spearman >= SPEARMAN_CUTOFF)
cat(paste0("[INFO] Size of edge_list AFTER filtereing above threshold: ", nrow(edge_list_rmt)), "\n")

# save that to a file -----------------------------------------------------

fwrite(edge_list_rmt, paste0("intermediate/network/occurance_ont_", Y_columns,".csv"))
cat("[INFO] final edgelist saved to: ", paste0("intermediate/network/occurance_ont_", Y_columns,".csv"))
cat("[TIME] Elapsed:", format(Sys.time() - start_time), "\n")
