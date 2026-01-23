#!/bin/bash
#SBATCH --job-name=CRISPRcasID
#SBATCH --output=log/21_crisprcasid_job_%j.out
#SBATCH --error=log/21_crisprcasid_job_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4      # Request more CPUs for potentially faster processing
#SBATCH --mem=16G              # Request more memory, especially for larger input files
#SBATCH --time=16:00:00        # Set a reasonable time limit (e.g., 4 hours)

module load Conda
conda activate /lisc/data/scratch/dome/willemsen/luecking/conda/conda_envs/crispr-env

# --- Define paths ---
CRISPRCASID_BIN="/lisc/home/user/luecking/bin/CRISPRcasIdentifier-1.1.0/CRISPRcasIdentifier.py"
INPUT_FASTA="intermediate/CRISPR/candidate_contigs.fna"
OUTPUT_FILE="intermediate/CRISPR/crisprcasid_output.txt"

# --- Run CRISPRcasIdentifier ---
echo "Starting CRISPRcasIdentifier on ${INPUT_FASTA}..."
echo "Output will be saved to ${OUTPUT_DIR}"

python "${CRISPRCASID_BIN}" \
    -f "${INPUT_FASTA}" \
    -o "${OUTPUT_FILE}" \
    -sc partial \
    -st dna 

echo "CRISPRcasIdentifier run complete."


mv output/* intermediate/CRISPR/
