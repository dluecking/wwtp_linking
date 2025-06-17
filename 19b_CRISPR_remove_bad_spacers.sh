#!/bin/bash

module load conda
conda activate crispr-env

prinseq-lite.pl \
    -fasta intermediate/minced/minced_results_spacers.fa \
    -lc_method dust \
    -lc_threshold 30 \
    -out_good intermediate/minced/minced_results_spacers_filtered_lc \
    -out_bad intermediate/minced/minced_results_spacers_low_complexity