#!/bin/bash
#SBATCH --job-name=prep_occurance_ont
#SBATCH --output=log/prep_occurance_ont_%A_%a.out
#SBATCH --error=log/prep_occurance_ont_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12    
#SBATCH --time=01:00:00     
#SBATCH --mem=128G             

module load conda
conda activate R_new

Rscript scripts/prepare_occurance_ONT.R 

echo "done!"

