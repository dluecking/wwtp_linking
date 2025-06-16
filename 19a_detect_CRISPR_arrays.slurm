#!/bin/bash
#SBATCH --job-name=run_minced
#SBATCH --output=log/minced_job_%j.out
#SBATCH --error=log/minced_job_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=03:00:00

module load conda
conda activate minced 

# Define input and output directories
INPUT_CONTIGS="intermediate/contigs/all_contigs.fna"
OUTPUT_DIR="intermediate/minced"
OUTPUT_BASE="minced_results" # Base name for output files

# Create the output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Run minced with options for comprehensive GFF output
minced -spacers -gffFull \
       "${INPUT_CONTIGS}" \
       "${OUTPUT_DIR}/${OUTPUT_BASE}.txt" \
       "${OUTPUT_DIR}/${OUTPUT_BASE}.gff"