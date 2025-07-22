#!/bin/bash
#SBATCH --job-name=iqtree2
#SBATCH --output=log/iqtree2_output_%A_%a.out
#SBATCH --error=log/iqtree2_error_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --mem=16G


module load iqtree

iqtree2 -T 8 -s intermediate/treebuilding/all_public_and_my_own_MCP_proteins_cleaned.aln -B 1000
