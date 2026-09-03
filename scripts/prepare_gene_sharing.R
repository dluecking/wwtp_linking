# Author: dlu @ veelab
# Version: 2025-06-17

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)

# set working directory
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# setwd("/lisc/data/scratch/dome/willemsen/luecking/projects/wwtp_linking")


# what fraction of genes needs to be shared so I keep it?
CUTOFF <- 0.05 # this is 5% right now

# load blastp output ------------------------------------------------------

blast_out <- fread("intermediate/blastp/all_vs_all_blastp.tsv")
# blast_out <- fread("../intermediate/blastp/all_vs_all_blastp_subset.tsv")
names(blast_out) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp")


# add contig ids ----------------------------------------------------------

blast_out$q_contig_id <- str_remove(blast_out$qseqid, "\\_\\d*$")
blast_out$s_contig_id <- str_remove(blast_out$sseqid, "\\_\\d*$")


# filtering ---------------------------------------------------------------

blast_out <- blast_out %>% 
  filter(q_contig_id != s_contig_id) %>% 
  filter(pident >= 80, qcovhsp >= 80)



# summarize and save ------------------------------------------------------

summary_table <- blast_out %>% 
  group_by(q_contig_id, s_contig_id) %>% 
  summarise(genes_shared = n(), .groups = "drop") %>% 
  rename(from = q_contig_id, to = s_contig_id)


# count genes for each contig ---------------------------------------------

protein_file <- "intermediate/proteins/all_proteins.faa"
protein_accessions <- system2("grep", args = c("'>'", protein_file), stdout = TRUE)
protein_accessions <- as.data.table(protein_accessions)

# get to contig_id
protein_accessions$contig_id <- str_remove(str_remove(str_remove(protein_accessions$protein_accessions, "\\s.*$"), "\\_\\d*$"), ">")
protein_numbers <- as.data.table(table(protein_accessions$contig_id))


# calculate the gene sharing value ----------------------------------------

summary_table$genes_per_contig <- protein_numbers$N[match(summary_table$from, protein_numbers$V1)]
summary_table$gene_sharing <- summary_table$genes_shared / summary_table$genes_per_contig

# fwrite(summary_table %>% filter(gene_sharing >= CUTOFF) %>% select(from, to, gene_sharing), file = "intermediate/network/gene_sharing.csv")



# quick exploration: ------------------------------------------------------
# does contig length correlate with number of genes shared?
# 
# plot_data <- summary_table %>%
#   group_by(from) %>%
#   summarise(
#     total_genes_shared = sum(genes_shared),
#     genes_per_contig = first(genes_per_contig) # Assuming this value is the same for all rows of a 'from' contig
#   )
# 
# ggplot(plot_data, aes(x = genes_per_contig, y = total_genes_shared)) +
#   geom_point(alpha = 0.5, color = "steelblue") +
#   # Add linear model with confidence intervals
#   geom_smooth(method = "lm", formula = y ~ x, color = "firebrick", fill = "lightgray") +
#   # Add R-squared and equation to the plot
#   annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, 
#            label = paste("R^2 == ", round(summary(lm(total_genes_shared ~ genes_per_contig, data = plot_data))$r.squared, 3)),
#            parse = TRUE) +
#   labs(
#     title = "Correlation: Genes vs. Total Shared Genes",
#     x = "Number of Genes on Contig",
#     y = "Total Number of Shared Genes (Summed)",
#     subtitle = "Linear regression with 95% confidence interval,\n10M subset, excluded 1 outlier."
#   ) +
#   theme_minimal()
# 
# ggsave(plot = last_plot(), filename = "../final/gene_sharing_correlation_genes_vs_genes_shared.png", height = 5, width = 5)



