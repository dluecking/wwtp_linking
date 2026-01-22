#!/bin/bash
#SBATCH --job-name=map_integration
#SBATCH --output=log/24_map_integration_%A_%a.out
#SBATCH --error=log/24_map_integration_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16    
#SBATCH --time=02:00:00     
#SBATCH --mem=256G             
#SBATCH --array=2-23%5 # up to 23


# load modules
module load minimap2/2.30-GCCcore-14.3.0
module load BEDTools/2.31.1-GCC-13.3.0
module load SAMtools/1.23-GCC-14.2.0


# combine vph and plv files to REF file
REF_FILE="intermediate/contigs/combined_gv_plv_vph.fna"
cat intermediate/contigs/vph/* intermediate/contigs/plv/* > "$REF_FILE"

# read info file
INFO_FILE="helperfiles/WWTP_Denmark_SRA_PRJNA629478.txt"

# Get the relevant line for the current array task
LINE_INFO=$(sed -n "$((SLURM_ARRAY_TASK_ID))p" "$INFO_FILE")

# Extract the full sample name from the first column of the line
FULL_SAMPLE_NAME=$(echo "$LINE_INFO" | awk '{print $1}')

# Extract the SRR accession from the third column of the line
SRR_ID=$(echo "$LINE_INFO" | awk '{print $3}')

# Extract the short sample name (first word before '_') from the full sample name
SHORT_SAMPLE_NAME=$(echo "$FULL_SAMPLE_NAME" | awk -F'_' '{print $1}')

# Construct the full path to the FASTQ read file using the SRR_ID
READFILE="data/reads/ont/${SRR_ID}_pchop.fastq.gz"

echo "Mapping reads for sample: $SHORT_SAMPLE_NAME (SRA: $SRR_ID)"
echo "Reads: $READFILE"
echo "Using reference: $REF_FILE"

# outputs - now using SHORT_SAMPLE_NAME for output file prefix
SAM="intermediate/integration/combined_mappings/${SHORT_SAMPLE_NAME}_mapped.sam"
BAM="intermediate/integration/combined_mappings/${SHORT_SAMPLE_NAME}_mapped.bam"
BAM_SORTED="intermediate/integration/combined_mappings/${SHORT_SAMPLE_NAME}_mapped_sorted.bam"


# map reads to ref file
minimap2 -ax map-ont -t "$SLURM_CPUS_PER_TASK" "$REF_FILE" "$READFILE" > "$SAM"
# convert to bam
samtools view -b -F 4 -@ "$SLURM_CPUS_PER_TASK" -T "$REF_FILE" "$SAM" > "$BAM"

# sort and index
samtools sort -@ "$SLURM_CPUS_PER_TASK" "$BAM" -o "$BAM_SORTED"
samtools index -@ "$SLURM_CPUS_PER_TASK" "$BAM_SORTED"

# remove sam and unsorted bam file, they are too big!
rm "$SAM"
rm "$BAM"
