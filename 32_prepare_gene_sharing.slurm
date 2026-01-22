#!/bin/bash
#SBATCH --job-name=prep_gene_sharing
#SBATCH --output=log/32_prep_gene_sharing_%A_%a.out
#SBATCH --error=log/32_prep_gene_sharing_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12    
#SBATCH --time=05:00:00     
#SBATCH --mem=256G             

module load Conda
conda activate /lisc/data/scratch/dome/willemsen/luecking/conda/conda_envs/R

Rscript scripts/prepare_gene_sharing.R 

echo "done!"

