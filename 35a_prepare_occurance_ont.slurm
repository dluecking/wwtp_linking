#!/bin/bash
#SBATCH --job-name=prep_gene_sharing
#SBATCH --output=log/35a_prep_occurance_ont_%A_%a.out
#SBATCH --error=log/35a_prep_occurance_ont_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16    
#SBATCH --time=12:00:00     
#SBATCH --mem=256G             

module load Conda
conda activate /lisc/data/scratch/dome/willemsen/luecking/conda/conda_envs/R_parallel

Rscript scripts/35a_prepare_occurance_ont.R 

echo "done!"

