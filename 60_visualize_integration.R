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
library(Rsamtools)
library(patchwork)
library(tidyr)


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

df <- blast_out8080 %>% 
  filter(sseqid == "OdNE_1") 

df <- df %>%
  mutate(read_type = case_when(
    str_detect(qseqid, "_m$") ~ "middle",
    TRUE ~ "boundary"
  ))


df <- df %>%
  select(qseqid, sseqid, mapped_to, sstart, send, slen, read_type)

# input -------------------------------------------------------------------

bam_path    <- "intermediate/integration/overhang_mappings/OdNE/overhanging_Aved_tig00303955-10-54120_vph.bam"
vph_contig  <- unique(df$mapped_to)          # e.g. "Aved_tig00303955-10-54120_vph"
host_contig <- unique(df$sseqid)             # e.g. "OdNE_1"



# parse bam ---------------------------------------------------------------

bam  <- BamFile(bam_path)
open(bam)

param <- ScanBamParam(
  what  = scanBamWhat(),
  which = GRanges(vph_contig, IRanges(1, 1e7))
)

raw <- scanBam(bam, param = param)[[1]]
close(bam)

# parse soft clips from CIGAR, compute ref_end, assign junction cluster
parse_clips <- function(cigar, pos) {
  sc5 <- as.integer(str_extract(cigar, "^(\\d+)(?=S)", group = 1))
  sc3 <- as.integer(str_extract(cigar, "(\\d+)(?=S$)", group = 1))
  sc5[is.na(sc5)] <- 0L
  sc3[is.na(sc3)] <- 0L
  
  # compute ref_end from CIGAR (M, D, N consume reference)
  ref_consumed <- sapply(cigar, function(cig) {
    m <- gregexpr("(\\d+)[MDN=X]", cig, perl = TRUE)
    lens <- regmatches(cig, m)[[1]]
    sum(as.integer(str_extract(lens, "\\d+")))
  })
  
  list(sc5 = sc5, sc3 = sc3, ref_end = pos + ref_consumed - 1L)
}

clips <- parse_clips(raw$cigar, raw$pos)

reads_raw <- tibble(
  read_id   = raw$qname,
  ref_start = raw$pos,
  ref_end   = clips$ref_end,
  sc5       = clips$sc5,
  sc3       = clips$sc3,
  strand    = ifelse(bitwAnd(raw$flag, 16L) > 0, "-", "+")
) %>%
  filter(!is.na(ref_start)) %>%
  mutate(
    junction = case_when(
      ref_start <= 20                    ~ "J1",
      ref_start >= 3740 & ref_start <= 3800 ~ "J2",
      TRUE                               ~ "other"
    ),
    clip_type = case_when(
      sc5 >= 100 & sc3 >= 100 ~ "both",
      sc5 >= 100              ~ "5only",
      sc3 >= 100              ~ "3only",
      TRUE                    ~ "none"
    ),
    read_group = paste0(junction, "_", clip_type),
    row        = row_number()
  )


# positions ---------------------------------------------------------------

vph_len  <- 5422
host_hits <- df %>%
  mutate(
    xmin      = pmin(sstart, send),
    xmax      = pmax(sstart, send),
    row       = row_number()
  )

host_xmin <- min(host_hits$xmin) - 2000
host_xmax <- max(host_hits$xmax) + 2000


# colors ------------------------------------------------------------------
pal <- c(
  J1_5only  = "#F0A882",   # orange — J1, 5' clip exits left
  J2_both   = "#D46B3A",   # green  — J2, spans both sides
  J2_3only  = "#D46B3A",   # blue   — J2, 3' clip exits right
  J2_5only  = "#D46B3A",   # light green
  other_none = "#D46B3A"
)

BOX_H <- 0.8


# plot 1, vph reads -------------------------------------------------------
# junction marker positions
j1_pos <- 1
j2_pos <- median(reads_raw$ref_start[reads_raw$junction == "J2"])

p_vph <- ggplot() +
  
  # virophage backbone
  annotate("rect", xmin = 1, xmax = vph_len, ymin = -BOX_H, ymax = BOX_H,
           fill = "goldenrod1", color = "black", linewidth = 0.6) +
  
  # junction lines
  geom_vline(xintercept = c(j1_pos, j2_pos),
             color = "#c0392b", linewidth = 0.8, linetype = "dashed") +
  
  # 5' soft clip tails (dashed, left of ref_start)
  geom_rect(
    data = reads_raw %>% filter(sc5 >= 100),
    aes(xmin = ref_start - sc5 * 0.15,   # scaled for display
        xmax = ref_start,
        ymin = row + 0.5,
        ymax = row + 0.5 + BOX_H,
        fill = read_group),
    color = "black", linewidth = 0.25, linetype = "dashed", alpha = 0.5
  ) +
  
  # mapped body
  geom_rect(
    data = reads_raw,
    aes(xmin = ref_start, xmax = ref_end,
        ymin = row + 0.5, ymax = row + 0.5 + BOX_H,
        fill = read_group),
    color = "black", linewidth = 0.25
  ) +
  
  # 3' soft clip tails (dashed, right of ref_end)
  geom_rect(
    data = reads_raw %>% filter(sc3 >= 100),
    aes(xmin = ref_end,
        xmax = ref_end + sc3 * 0.15,
        ymin = row + 0.5,
        ymax = row + 0.5 + BOX_H,
        fill = read_group),
    color = "black", linewidth = 0.25, linetype = "dashed", alpha = 0.5
  ) +
  
  # contig label
  # annotate("text", x = 1, y = max(reads_raw$row) + 3.5,
  #          label = vph_contig, size = 3, hjust = 0,
  #          color = "black", fontface = "bold") +
  
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_continuous(labels = scales::comma,
                     breaks = scales::pretty_breaks(6)) +
  coord_cartesian(xlim = c(-200, vph_len + 200)) +
  labs(x = "Position on virophage (bp)", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid   = element_blank(),
        plot.title   = element_text(face = "bold", size = 12))


# plot 2: blast hits on NCV -----------------------------------------------

# insertion locus = min/max of all hits
ins_xmin <- min(host_hits$xmin)
ins_xmax <- max(host_hits$xmax)
pivot    <- median(c(host_hits$xmin, host_hits$xmax))  # ~943967

p_host <- ggplot() +
  
  # NCV backbone
  annotate("rect", xmin = host_xmin, xmax = host_xmax, ymin = -BOX_H, ymax = BOX_H,
           fill = "steelblue", color = "black", linewidth = 0.6) +

  # pivot line
  geom_vline(xintercept = pivot,
             color = "#c0392b", linewidth = 0.8, linetype = "dashed") +
  
  # BLAST hit rectangles
  geom_rect(data = host_hits,
            aes(xmin = xmin, xmax = xmax,
                ymin = row * (BOX_H + 0.1) + 0.3,
                ymax = row * (BOX_H + 0.1) + 0.3 + BOX_H,
                fill = read_type),
            color = "black", linewidth = 0.25) +
  
  scale_fill_manual(values = c("boundary" = "#e8956d", "middle" = "#D46B3A"),
                    name = NULL) +
  
  # pivot label
  annotate("text", x = pivot, y = max(host_hits$row) * (BOX_H + 0.1) + 1.2,
           label = scales::comma(round(pivot)),
           size = 3, hjust = -0.1, color = "#c0392b", fontface = "bold") +
  
  # contig label
  annotate("text", x = host_xmin,
           y = max(host_hits$row) * (BOX_H + 0.1) + 1.5,
           label = host_contig, size = 3, hjust = 0,
           color = "black", fontface = "bold") +
  
  scale_x_continuous(labels = scales::comma,
                     breaks = scales::pretty_breaks(5)) +
  labs(x = sprintf("Position on %s (bp)", host_contig), y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.y   = element_blank(),
        axis.ticks.y  = element_blank(),
        panel.grid    = element_blank(),
        legend.position = "none")


# combine plots -----------------------------------------------------------

p_vph / p_host +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(
    title    = sprintf("%s  →  %s", vph_contig, host_contig),
    subtitle = sprintf("%d reads · J1 (pos %d) + J2 (pos %d) · BLAST pivot ~%s bp",
                       nrow(reads_raw), j1_pos, round(j2_pos),
                       scales::comma(round(pivot))),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10, color = "grey50")
    )
  )





# old

# old ---------------------------------------------------------------------
# 
# 
# # data wrangling ----------------------------------------------------------
# 
# # GV_OF_INTEREST <- "Vibo_2_3_4"
# # GV_OF_INTEREST <- "Hjor_1"
# # GV_OF_INTEREST <- "Bjer_2_3"
# GV_OF_INTEREST <- "OdNE_1"
# # GV_OF_INTEREST <- "Hjor_2_3_4_6"
# 
# for(GV_OF_INTEREST in c("Vibo_2_3_4", "Hjor_1", "Bjer_2_3", "OdNE_1", "Hjor_2_3_4_6")){
#   print(GV_OF_INTEREST)
#   df <- blast_out8080 %>% 
#     filter(sseqid == GV_OF_INTEREST) 
#   
#   df <- df %>%
#     mutate(read_type = case_when(
#       str_detect(qseqid, "_m$") ~ "middle",
#       TRUE ~ "boundary"
#     ))
#   
#   
#   df <- df %>% 
#     select(qseqid, sseqid, mapped_to, sstart, send, slen, read_type)
#   
#   # how many top matching contigs do we want to keep?
#   N = 1
#   
#   
#   # keep only the 8 most frequent sseqids
#   top <- df %>%
#     dplyr::count(sseqid) %>%                 # count occurrences
#     arrange(desc(n)) %>%              # sort manually
#     slice_head(n = N) %>%             # take top 8
#     pull(sseqid)
#   
#   df_top <- df %>%
#     filter(sseqid %in% top)
#   
#   # make qseqid a factor for plotting order
#   df_top <- df_top %>%
#     mutate(qseqid = factor(qseqid))
#   
#   # summarize the number of reads per sseqid
#   label_data <- df_top %>%
#     group_by(sseqid) %>%
#     summarise(n_reads = n())
#   
#   # plot mappings per sseqid
#   ggplot(df_top) +
#     geom_rect(aes(xmin = sstart, xmax = send,
#                   ymin = as.numeric(qseqid) - 0.4,
#                   ymax = as.numeric(qseqid) + 0.4,
#                   fill = read_type),
#               alpha = 0.7) +
#     scale_fill_manual(values = c("boundary" = "red", "middle" = "black")) +
#     facet_wrap(~ fct_infreq(sseqid), scales = "free_y", ncol = 1) +
#     scale_x_continuous(name = "Position on contig") +
#     scale_y_continuous(name = "Mapped reads", breaks = NULL) +
#     ggtitle(paste0("Integration of ", df_top$mapped_to[1])) +
#     theme_bw() +
#     theme(
#       strip.text = element_text(face = "bold", size = 6),
#       panel.spacing = unit(0.5, "lines")
#     ) +
#     # add labels to the top right corner of each facet
#     geom_text(
#       data = label_data,
#       aes(
#         x = Inf, y = Inf,
#         label = paste0("n=", n_reads)
#       ),
#       hjust = 1.1, vjust = 1.5,   # position slightly inside the corner
#       size = 3.2,
#       color = "black",
#       inherit.aes = FALSE
#     )
#   
#   ggsave(file = paste0("final/integration_plots/", GV_OF_INTEREST, "integration_plot.png"), 
#          plot = last_plot(),
#          height = 3, width = 6)
# }
# 
# 
# 
# 
# # quick exploratory thing -------------------------------------------------
# # we need to vis the starting positions of vibo_middle reads on vibo
# pos <- fread("testing/vibo_middle_read_positions_on_vibo.txt")
# 
# ggplot(pos, aes(x = V2)) +
#   geom_rug(alpha = 0.4) +
#   geom_histogram(aes(y = after_stat(count)), bins = 200, alpha = 0.8, color = "black") +
#   scale_x_continuous(limits = c(0, 140000), name = "Position on contig (zoomed)") +
#   ylab("Count / rug ticks") +
#   ggtitle("Read mapping positions on Vibo_2_3_4") +
#   theme_bw()
# 
# ggsave(file = paste0("final/integration_plots/vibo_middle_read_positions_on_vibo.png"), 
#        plot = last_plot(),
#        height = 3, width = 6)
