#!/bin/bash
module load conda
conda activate spacerextractor

TARGET_FILE="intermediate/contigs/combined_gv_plv_vph.fna"
TARGET_DB="intermediate/CRISPR/map_spacers_to_targets/target_db"
SPACER_FILE="intermediate/minced/minced_results_spacers_filtered_lc.fasta"


# outdir and working directory
mkdir -p intermediate/CRISPR/map_spacers_to_targets

# create combined target file
cat data/gv_contigs/* intermediate/contigs/vph/* intermediate/contigs/plv/* > $TARGET_FILE

# create target database
spacerextractor create_target_db -i $TARGET_FILE -d $TARGET_DB -t 8 --replace_spaces

# map spacers to targets
spacerextractor map_to_target -i $SPACER_FILE -d $TARGET_DB -o intermediate/CRISPR/map_spacers_to_targets -t 8 

