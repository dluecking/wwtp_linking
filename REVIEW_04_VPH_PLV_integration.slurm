#!/bin/bash
#SBATCH --job-name=vph_plv_integration
#SBATCH --output=log/REVIEW_04_vph_plv_integration_%A.out
#SBATCH --error=log/REVIEW_04_vph_plv_integration_%A.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=4:00:00
#SBATCH --mem=32G

# REVIEW_04_VPH_PLV_integration.slurm
#
# Search all VPH/PLV contigs against ALL microbial contigs with blastn megablast,
# to look for evidence of VPH/PLV sequence integrated into microbial assemblies
# (near-identical nucleotide matches), analogous to the GV-vs-all-contigs overlap
# check in 05_blastn_GVs_vs_all_others.sh.
#
# Reuses the "all_contigs" BLAST db already built by 26_integration_blast_overhangs.sh
# (intermediate/contigs/all_contigs.fna = lc/*_filtered.fna + vph/* + plv/* +
# data/gv_contigs/*, see 15_occurance_concat_all_contigs.sh) instead of rebuilding one.
# Since that db also contains the viral sequences themselves, VPH/PLV/GV accessions are
# excluded from the search via -negative_seqidlist so every hit is against a microbial
# contig, not a self-hit or a virus-vs-virus hit.

set -e  # Exit on error

# --- Load necessary tools ---
module load BLAST+/2.17.0-gompi-2024a

# --- Configuration ---
WORK_DIR="intermediate/REVIEW_vph_plv_integration"

VPH_DIR="intermediate/contigs/vph"
PLV_DIR="intermediate/contigs/plv"
GV_DIR="data/gv_contigs"

ALL_CONTIGS_DB="intermediate/contigs/all_contigs" # built by 26_integration_blast_overhangs.sh

QUERY_FASTA="${WORK_DIR}/vph_plv_contigs.fna"
VIRAL_SEQID_LIST="${WORK_DIR}/viral_seqids.txt"
VIRAL_SEQID_LIST_BSL="${WORK_DIR}/viral_seqids.bsl"

BLAST_OUTPUT_FILE="${WORK_DIR}/vph_plv_vs_microbial_megablast.tsv"

# Minimum percent identity for megablast (matches 05_blastn_GVs_vs_all_others.sh)
MIN_PERCENT_IDENTITY=90

echo "Working directory: ${WORK_DIR}"
mkdir -p "${WORK_DIR}"

if [ ! -f "${ALL_CONTIGS_DB}.nal" ] && [ ! -f "${ALL_CONTIGS_DB}.nhr" ]; then
    echo "ERROR: Expected BLAST database ${ALL_CONTIGS_DB} not found. Run 26_integration_blast_overhangs.sh (or makeblastdb) first." >&2
    exit 1
fi

# --- Step A: Combine all VPH + PLV contigs into a single query FASTA ---
echo "Step A: Combining VPH and PLV contigs into query FASTA..."
cat "${VPH_DIR}"/*.fna "${PLV_DIR}"/*.fna > "${QUERY_FASTA}"
echo "Query FASTA written to: ${QUERY_FASTA}"

# --- Step B: Build a negative seqid list (VPH + PLV + GV) so hits can only land on microbial contigs ---
echo "Step B: Building negative seqid list of viral accessions (VPH/PLV/GV)..."
grep -h "^>" "${VPH_DIR}"/*.fna "${PLV_DIR}"/*.fna "${GV_DIR}"/* | sed -e 's/^>//' -e 's/ .*//' > "${VIRAL_SEQID_LIST}"
blastdb_aliastool -seqid_file_in "${VIRAL_SEQID_LIST}" -seqid_file_out "${VIRAL_SEQID_LIST_BSL}"

# --- Step C: Run BLASTN (megablast) of VPH/PLV contigs against ALL microbial contigs ---
echo "Step C: Running BLASTN (megablast)..."
blastn -query "${QUERY_FASTA}" \
       -db "${ALL_CONTIGS_DB}" \
       -negative_seqidlist "${VIRAL_SEQID_LIST_BSL}" \
       -out "${BLAST_OUTPUT_FILE}" \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
       -task megablast \
       -num_threads "${SLURM_CPUS_PER_TASK}" \
       -perc_identity "${MIN_PERCENT_IDENTITY}" \
       -max_target_seqs 1000

echo "BLASTN comparison complete. Results saved to: ${BLAST_OUTPUT_FILE}"
