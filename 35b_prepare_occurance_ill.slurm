#!/bin/bash
#SBATCH --job-name=prep_gene_sharing
#SBATCH --output=log/35b_prep_occurance_ill_%A_%a.out
#SBATCH --error=log/35b_prep_occurance_ill_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16    
#SBATCH --time=04:00:00     
#SBATCH --mem=128G             

module load Conda
conda activate /lisc/data/scratch/dome/willemsen/luecking/conda/conda_envs/R_parallel

Rscript scripts/35b_prepare_occurance_ill.R 

echo "done!"

