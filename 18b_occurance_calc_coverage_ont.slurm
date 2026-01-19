#!/bin/bash
#SBATCH --job-name=cov_ont
#SBATCH --output=log/18b_cov_ont_%A_%a.out
#SBATCH --error=log/18b_cov_ont_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1    
#SBATCH --time=08:00:00     
#SBATCH --mem=64G             
#SBATCH --array=1-23 # up to 23

# load module
module load SAMtools/1.23-GCC-14.2.0

READLIST="helperfiles/ont_readfiles.txt"
SRR=$(sed -n "$((SLURM_ARRAY_TASK_ID ))p" "$READLIST")

# input files
BAM_SORTED="intermediate/occurance/mapping/ont/${SRR}_mapped_sorted.bam"
# output files
COV_FILE="intermediate/occurance/mapping/ont/${SRR}_mapped_coverage.txt"

# actually calculate coverage
echo "Calculating coverage for sample: $SRR"
echo "Using BAM file: $BAM_SORTED"

samtools coverage $BAM_SORTED > $COV_FILE
