#!/bin/bash

module load Conda
conda activate /lisc/data/scratch/dome/willemsen/luecking/conda/conda_envs/R

Rscript scripts/CRISPR_create_list_of_contigs_with_array.R

echo "done"
