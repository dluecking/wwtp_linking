#!/bin/bash
#SBATCH --job-name=annotate_vph
#SBATCH --output=log/annotate_vph_output_%A_%a.out
#SBATCH --error=log/annotate_vph_error_%A_%a.err
#SBATCH --array=2-15%50
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --mem=8G

########## Setup ##########
protein_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" helperfiles/vph_proteins_to_annotate.txt)
contig_file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" helperfiles/vph_contigs_to_annotate.txt)

protein_basename=${protein_file%.faa}
mkdir -p intermediate/proteins/cleaned/vph

########## Clean asterisks from protein file ##########
echo "Cleaning ${protein_file} to remove asterisks"
sed 's/\*//g' intermediate/proteins/vph/${protein_file} > intermediate/proteins/cleaned/vph/${protein_file}

########## INTERPROSCAN ##########
module load interproscan

echo "Running InterProScan for ${protein_file}"
mkdir -p intermediate/annotations/interpro/vph
interproscan.sh \
  -i intermediate/proteins/cleaned/vph/${protein_file} \
  -b intermediate/annotations/interpro/vph/${protein_basename} \
  --cpu 4

########## PROKKA ##########
module load conda
conda activate prokka

echo "Running Prokka for ${contig_file}"
mkdir -p intermediate/annotations/prokka/vph
prokka \
  --hmms data/hmm/vph_markers_roux.hmm \
  --cpus 4 \
  --kingdom Viruses \
  --coverage 50 \
  --evalue 0.001 \
  --outdir intermediate/annotations/prokka/vph/${contig_file%.fna} \
  --prefix ${contig_file%.fna} \
  --rawproduct --debug --addgenes \
  intermediate/contigs/vph/${contig_file}


