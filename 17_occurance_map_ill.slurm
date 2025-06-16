#!/bin/bash
#SBATCH --job-name=map_ill
#SBATCH --output=log/map_ill_%A_%a.out
#SBATCH --error=log/map_ill_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16    
#SBATCH --time=1-00:00:00     
#SBATCH --mem=256G             
#SBATCH --array=21-21 # up to 23

# load module
module load bowtie2
module load samtools
module load bedtools

INDEX_FILE="intermediate/contigs/all_contigs_index"
READLIST="helperfiles/illumina_readfiles.txt"

SRR=$(sed -n "$((SLURM_ARRAY_TASK_ID ))p" "$READLIST")

R1="intermediate/occurance/reads_derep/illumina/${SRR}_1_derep.fastq.gz"
R2="intermediate/occurance/reads_derep/illumina/${SRR}_2_derep.fastq.gz"

echo "Mapping reads for sample: $SRR"
echo "R1: $R1"
echo "R2: $R2"
echo "Using index: $INDEX_FILE"

# OUTPUT FILES
SAM="intermediate/occurance/mapping/illumina/${SRR}_mapped.sam"
BAM="intermediate/occurance/mapping/illumina/${SRR}_mapped.bam"
BAM_SORTED="intermediate/occurance/mapping/illumina/${SRR}_mapped_sorted.bam"

# map reads to reference
bowtie2 --sensitive-local --no-mixed --no-unal -I 0 -X 1000 -p $SLURM_CPUS_PER_TASK -x "$INDEX_FILE" -1 "$R1" -2 "$R2" -S "$SAM"

# convert sam to bam
samtools view -b -F 4 $SAM > $BAM
   
# sort and index
samtools sort -@ $SLURM_CPUS_PER_TASK $BAM -o $BAM_SORTED
samtools index -@ $SLURM_CPUS_PER_TASK $BAM_SORTED

# delete sam and unsorted bam file, they are too big!
rm $SAM
rm $BAM

