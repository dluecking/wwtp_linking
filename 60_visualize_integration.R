# Author: dlu @ veelab
# Version: 2025-08-25

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)
library(forcats)
library(cowplot)



# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# load data and retrieve mapped_to ----------------------------------------

blast_out <- fread("intermediate/integration/blast_results/all_overhangs_blastn.tsv")
names(blast_out) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                      "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                      "qlen", "slen")

blast_out$read <- str_remove(blast_out$qseqid, "\\_to\\_.*$")
blast_out$mapped_to <- str_remove(str_remove(blast_out$qseqid, "^.*\\_to\\_"), "\\_overhang.*$")

# filter out self hits
blast_out <- blast_out %>% 
  filter(mapped_to != sseqid)

blast_out$qcov <- ( blast_out$qend - blast_out$qstart ) /blast_out$qlen * 100

blast_out8080 <- blast_out %>% filter(pident >= 80, qcov >= 80)



# data wrangling ----------------------------------------------------------

# GV_OF_INTEREST <- "Vibo_2_3_4"
# GV_OF_INTEREST <- "Hjor_1"
# GV_OF_INTEREST <- "Bjer_2_3"
# GV_OF_INTEREST <- "OdNE_1"
# GV_OF_INTEREST <- "Hjor_2_3_4_6"

for(GV_OF_INTEREST in c("Vibo_2_3_4", "Hjor_1", "Bjer_2_3", "OdNE_1", "Hjor_2_3_4_6")){
  print(GV_OF_INTEREST)
  df <- blast_out8080 %>% 
    filter(mapped_to == GV_OF_INTEREST) 
  
  df <- df %>%
    mutate(read_type = case_when(
      str_detect(qseqid, "_m$") ~ "middle",
      TRUE ~ "boundary"
    ))
  
  
  df <- df %>% 
    select(qseqid, sseqid, mapped_to, sstart, send, slen, read_type)
  
  # how many top matching contigs do we want to keep?
  N = 5
  
  
  # keep only the 8 most frequent sseqids
  top <- df %>%
    dplyr::count(sseqid) %>%                 # count occurrences
    arrange(desc(n)) %>%              # sort manually
    slice_head(n = N) %>%             # take top 8
    pull(sseqid)
  
  df_top <- df %>%
    filter(sseqid %in% top)
  
  # make qseqid a factor for plotting order
  df_top <- df_top %>%
    mutate(qseqid = factor(qseqid))
  
  # summarize the number of reads per sseqid
  label_data <- df_top %>%
    group_by(sseqid) %>%
    summarise(n_reads = n())
  
  # plot mappings per sseqid
  ggplot(df_top) +
    geom_rect(aes(xmin = sstart, xmax = send,
                  ymin = as.numeric(qseqid) - 0.4,
                  ymax = as.numeric(qseqid) + 0.4,
                  fill = read_type),
              alpha = 0.7) +
    scale_fill_manual(values = c("boundary" = "red", "middle" = "black")) +
    facet_wrap(~ fct_infreq(sseqid), scales = "free_y", ncol = 1) +
    scale_x_continuous(name = "Position on contig") +
    scale_y_continuous(name = "Mapped reads", breaks = NULL) +
    ggtitle(paste0("Integration of ", df_top$mapped_to[1])) +
    theme_bw() +
    theme(
      strip.text = element_text(face = "bold", size = 6),
      panel.spacing = unit(0.5, "lines")
    ) +
    # add labels to the top right corner of each facet
    geom_text(
      data = label_data,
      aes(
        x = Inf, y = Inf,
        label = paste0("n=", n_reads)
      ),
      hjust = 1.1, vjust = 1.5,   # position slightly inside the corner
      size = 3.2,
      color = "black",
      inherit.aes = FALSE
    )
  
  ggsave(file = paste0("final/integration_plots/", GV_OF_INTEREST, "integration_plot.png"), 
         plot = last_plot(),
         height = 6, width = 6)
}
