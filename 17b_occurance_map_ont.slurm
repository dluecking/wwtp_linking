#!/bin/bash
#SBATCH --job-name=map_ont
#SBATCH --output=log/17b_map_ont_%A_%a.out
#SBATCH --error=log/17b_map_ont_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16    
#SBATCH --time=1-00:00:00     
#SBATCH --mem=256G             
#SBATCH --array=21%6 # up to 23

# load modules
module load SAMtools/1.23-GCC-14.2.0
module load BEDTools/2.31.1-GCC-13.3.0
module load minimap2/2.30-GCCcore-14.3.0

# inputs
REF_FILE="intermediate/contigs/all_contigs.fna"
READLIST="helperfiles/ont_readfiles.txt"

SRR=$(sed -n "$((SLURM_ARRAY_TASK_ID ))p" "$READLIST")
READFILE="data/reads/ont/${SRR}_pchop.fastq.gz"

echo "Mapping reads for sample: $SRR"
echo "Reads: $READFILE"
echo "Using reference: $REF_FILE"

# create output dir
OUTDIR="intermediate/occurance/mapping/ont"
mkdir -p "$OUTDIR"

# outputs
SAM="${OUTDIR}/${SRR}_mapped.sam"
BAM="${OUTDIR}/${SRR}_mapped.bam"
BAM_SORTED="${OUTDIR}/${SRR}_mapped_sorted.bam"
COV_FILE="${OUTDIR}/${SRR}_mapped_coverage.txt"

# map reads to ref file
minimap2 -ax map-ont -t $SLURM_CPUS_PER_TASK $REF_FILE $READFILE > $SAM
# convert to bam
samtools view -b -F 4 -@ $SLURM_CPUS_PER_TASK -T $REF_FILE $SAM > $BAM

# sort and index
samtools sort -@ $SLURM_CPUS_PER_TASK $BAM -o $BAM_SORTED
samtools index -@ $SLURM_CPUS_PER_TASK $BAM_SORTED

# calculate coverage # we will be doing this later
# samtools coverage $BAM_SORTED > $COV_FILE

# remove sam and unsorted bam file, they are too big!
rm $SAM
rm $BAM