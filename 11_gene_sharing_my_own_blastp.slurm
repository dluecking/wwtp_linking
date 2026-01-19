#!/bin/bash
#SBATCH --job-name=own_blastp
#SBATCH --output=log/11_own_blastp_%A_%a.out # OUTPUT THAT WOULD PRINT TO SCREEN GOES HERE
#SBATCH --error=log/11_own_blastp_%A_%a.err # LOG FILE GOES HERE
#SBATCH --ntasks=1 # FOR PARRALELIZATION, CAN SAFELY KEEP AT 1
#SBATCH --cpus-per-task=12 # NUMBER OF CPUS YOU NEED (good number: 4-16)
#SBATCH --time=3-00:00:00  # TIME YOU NEED DAY:HOURS:MINUTES:SECONDS (MORE -> it finishes probably, LESS -> you get priority on SLURM)
#SBATCH --mem=256G         # AMOUNT OF RAM YOU NEED (good number: 8GB for small jobs, 64GB for medium jobs, 256+GB for assemblies)


## --- Load necessary modules ---
module load DIAMOND/2.1.9-GCC-13.2.0

## --- Define File Paths ---
ALL_PROTEINS="intermediate/proteins/all_proteins.faa"
DB_DIR="intermediate/blastp"
DB_PATH="${DB_DIR}/all_proteins.dmnd"
OUTPUT_FILE="${DB_DIR}/all_vs_all_blastp.tsv"


mkdir -p $DB_DIR

echo "Starting DIAMOND all-vs-all BLASTp job"
echo "Using ${SLURM_CPUS_PER_TASK} cores."


echo "--- Creating DIAMOND database ---"
diamond makedb --in $ALL_PROTEINS \
               -d $DB_PATH \
               --threads $SLURM_CPUS_PER_TASK


echo "--- Running all-vs-all BLASTp search ---"
diamond blastp --query $ALL_PROTEINS \
               --db $DB_PATH \
               --out $OUTPUT_FILE \
               --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp \
               --min-score 50 \
               --sensitive \
               --threads $SLURM_CPUS_PER_TASK

echo "Job finished successfully!"