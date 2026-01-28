#!/bin/bash
#SBATCH --job-name=prep_gene_sharing
#SBATCH --output=log/36_parse_lc_tax_info_%A_%a.out
#SBATCH --error=log/36_parse_lc_tax_info_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16    
#SBATCH --time=08:00:00     
#SBATCH --mem=48G             

module load Conda
conda activate /lisc/data/scratch/dome/willemsen/luecking/conda/conda_envs/R

Rscript scripts/36_parse_lc_tax_info.R 

echo "done!"
