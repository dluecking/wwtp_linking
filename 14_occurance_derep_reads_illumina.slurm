#!/bin/bash
#SBATCH --job-name=derep_illumina
#SBATCH --output=log/14_derep_illumina_%A_%a.out
#SBATCH --error=log/14_derep_illumina_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=256G
#SBATCH --array=2-23 # up ot 23

module load CD-HIT/4.8.1-GCC-12.3.0

INPUT_READS_DIR="data/reads/illumina/"
OUTPUT_DIR="intermediate/occurance/reads_derep/illumina/"
READLIST="helperfiles/illumina_readfiles.txt"

SRR=$(sed -n "$((SLURM_ARRAY_TASK_ID ))p" "$READLIST")

r1_gz="${INPUT_READS_DIR}/${SRR}_1.fastq.gz"
r2_gz="${INPUT_READS_DIR}/${SRR}_2.fastq.gz"

tmp_r1="$TMPDIR/${SRR}_1.fastq"
tmp_r2="$TMPDIR/${SRR}_2.fastq"
tmp_out1="$TMPDIR/${SRR}_1_derep.fastq"
tmp_out2="$TMPDIR/${SRR}_2_derep.fastq"
tmp_clstr="$TMPDIR/${SRR}_dup.clstr"

out1="${OUTPUT_DIR}/${SRR}_1_derep.fastq.gz"
out2="${OUTPUT_DIR}/${SRR}_2_derep.fastq.gz"
clstr="${OUTPUT_DIR}/${SRR}_dup.clstr"

mkdir -p "$OUTPUT_DIR"

gunzip -c "$r1_gz" > "$tmp_r1"
gunzip -c "$r2_gz" > "$tmp_r2"

cd-hit-dup -i "$tmp_r1" -i2 "$tmp_r2" -o "$tmp_out1" -o2 "$tmp_out2" -e 0

gzip "$tmp_out1"
gzip "$tmp_out2"

mv "${tmp_out1}.gz" "$out1"
mv "${tmp_out2}.gz" "$out2"
mv "$tmp_clstr" "$clstr"
