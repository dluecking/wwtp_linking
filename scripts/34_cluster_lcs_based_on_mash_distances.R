# Author: dlu @ veelab
# Version: 2026-02-11

# Packages
library(argparse)
library(data.table)
library(dplyr)
library(igraph)

# set working directory
# setwd("/lisc/data/scratch/dome/willemsen/luecking/projects/wwtp_linking")

# Parse arguments
parser <- ArgumentParser(description = "Cluster LC contigs based on Mash distances")

parser$add_argument("--threshold", type = "double", default = 0.03,
                    help = "Distance threshold for clustering")
parser$add_argument("--mash_distances_filtered", type = "character", required = TRUE,
                    help = "Path to filtered Mash distances file")
parser$add_argument("--all_lc_contig_ids", type = "character", required = TRUE,
                    help = "Path to file with all LC contig IDs")
parser$add_argument("--membership_output", type = "character", required = TRUE,
                    help = "Output path for cluster membership file")

opt <- parser$parse_args()



# # debugging ---------------------------------------------------------------
setwd("/run/user/1000/gvfs/sftp:host=login01.lisc.univie.ac.at,user=luecking/lisc/home/user/luecking/luecking_scratch/projects/wwtp_linking/")

opt$threshold <- 0.03
opt$mash_distances_filtered <- "intermediate/mash/lc_contigs_only_mash_distances_filtered.txt"
opt$mash_distances_filtered <- "intermediate/mash/mash_distances_filtered_005.txt"
opt$all_lc_contig_ids <- "intermediate/network/all_contig_ids.txt"
opt$membership_output <- "intermediate/mash/lc_contigs_only_cluster_membership_005.csv"


# clustering --------------------------------------------------------------

cat("=== LC Contig Clustering ===\n")
cat("Threshold:", opt$threshold, "\n\n")

# Read data# Read data
cat("Loading data...\n")
distances <- fread(opt$mash_distances_filtered, 
                   col.names = c("contig1", "contig2", "distance", "p_value", "shared_hashes"))
all_contigs <- fread(opt$all_lc_contig_ids, header = FALSE)$V1  # Should be ALL, not just LC

# Separate LC from non-LC
lc_contigs <- all_contigs[grepl("_lc$", all_contigs)]
non_lc_contigs <- all_contigs[!grepl("_lc$", all_contigs)]

cat("Similar pairs:", nrow(distances), "\n")
cat("Total contigs:", length(all_contigs), "\n")
cat("  - LC contigs:", length(lc_contigs), "\n")
cat("  - Non-LC contigs (VPH/PLV/GV):", length(non_lc_contigs), "\n\n")

# Build graph and find clusters (LCs only)
cat("Finding clusters for LC contigs...\n")
g <- graph_from_data_frame(distances %>% select(contig1, contig2), directed = FALSE)
clusters <- components(g)

# Create membership for clustered LCs
membership <- data.frame(
  contig_id = V(g)$name,
  cluster_id = clusters$membership
) %>%
  mutate(cluster_id = sprintf("LC_CLUSTER_%05d", cluster_id))

# Add LC singletons
lc_singletons <- setdiff(lc_contigs, membership$contig_id)
if (length(lc_singletons) > 0) {
  max_id <- max(clusters$membership)
  singleton_df <- data.frame(
    contig_id = lc_singletons,
    cluster_id = sprintf("LC_CLUSTER_%05d", max_id + 1:length(lc_singletons))
  )
  membership <- bind_rows(membership, singleton_df)
}

# Add cluster size for LCs
membership <- membership %>%
  group_by(cluster_id) %>%
  mutate(cluster_size = n()) %>%
  ungroup()

# Add non-LC contigs (each is its own "cluster")
if (length(non_lc_contigs) > 0) {
  non_lc_df <- data.frame(
    contig_id = non_lc_contigs,
    cluster_id = non_lc_contigs,  # Use their own ID as cluster ID
    cluster_size = 1
  )
  membership <- bind_rows(membership, non_lc_df)
}

# Final arrangement
membership <- membership %>% arrange(cluster_id)

# save membership
fwrite(membership, opt$membership_output, col.names = TRUE)

# Stats (LC clusters only)
cat("\n=== LC Clustering Results ===\n")
cat("Total LC contigs:", length(lc_contigs), "\n")
cat("LC clusters:", n_distinct(membership$cluster_id[grepl("LC_CLUSTER", membership$cluster_id)]), "\n")
cat("LC singletons:", sum(membership$cluster_size == 1 & grepl("LC_CLUSTER", membership$cluster_id)), "\n")
cat("Largest LC cluster:", max(membership$cluster_size[grepl("LC_CLUSTER", membership$cluster_id)]), "contigs\n")

cat("\n=== Final Membership ===\n")
cat("Total contigs:", nrow(membership), "\n")
cat("  - LC clusters:", sum(grepl("LC_CLUSTER", membership$cluster_id)), "\n")
cat("  - Non-LC (individual):", sum(!grepl("LC_CLUSTER", membership$cluster_id)), "\n\n")