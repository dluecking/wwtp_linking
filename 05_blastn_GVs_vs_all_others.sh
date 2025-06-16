#!/bin/bash
#SBATCH --job-name=blast_contig_overlap
#SBATCH --output=log/blast_contig_overlap_%A.out
#SBATCH --error=log/blast_contig_overlap_%A.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16 # Adjust based on available cores and desired speed
#SBATCH --time=4:00:00 # Adjust time based on data size
#SBATCH --mem=32G # Adjust memory based on database and query size. BLAST can be memory intensive.

# --- Configuration ---
GV_CONTIGS_DIR="data/gv_contigs/"
QUERY_CONTIGS_BASE_DIR="intermediate/contigs/" # Contains subdirs vph, lc, plv
BLAST_DB_DIR="intermediate/blast_db/" # Directory to store the created BLAST database
BLAST_OUTPUT_DIR="intermediate/blast_results/" # Directory for BLAST output

# Define the paths for the combined FASTA files
GV_COMBINED_FASTA="${BLAST_DB_DIR}gv_contigs_db.fna"
QUERY_COMBINED_FASTA="${BLAST_DB_DIR}query_contigs_for_blast.fna"

# --- BLAST Parameters ---
# Minimum percent identity for megablast (e.g., 95 for very similar, 90 for strong overlap)
MIN_PERCENT_IDENTITY=95
# Minimum query coverage (e.g., 80% of query length must align)
MIN_QUERY_COVERAGE=80

# --- Load necessary tools ---
module load ncbiblastplus
module load conda
conda activate bioinf # Ensure blast+ (makeblastdb, blastn) is in this environment

echo "Starting BLAST comparison for contig overlap..."
echo "GV Contigs directory: ${GV_CONTIGS_DIR}"
echo "Query Contigs base directory: ${QUERY_CONTIGS_BASE_DIR}"
echo "BLAST Database directory: ${BLAST_DB_DIR}"
echo "BLAST Results directory: ${BLAST_OUTPUT_DIR}"
echo "---"

# Create necessary directories
mkdir -p "${BLAST_DB_DIR}"
mkdir -p "${BLAST_OUTPUT_DIR}"

# --- Step A: Create database of contigs in data/gv_contigs/ ---
echo "Step A: Combining GV contigs and creating BLAST database..."

# Combine all .fna files from GV_CONTIGS_DIR into a single FASTA for the database
find "${GV_CONTIGS_DIR}" -name "*.fasta" -print0 | xargs -0 cat > "${GV_COMBINED_FASTA}"

if [ ! -s "${GV_COMBINED_FASTA}" ]; then
    echo "ERROR: No .fasta files found in ${GV_CONTIGS_DIR} or combined FASTA is empty. Exiting."
    exit 1
fi

# Create the BLAST database
makeblastdb -in "${GV_COMBINED_FASTA}" -dbtype nucl -out "${BLAST_DB_DIR}gv_contigs_db" -parse_seqids

if [ $? -ne 0 ]; then
    echo "ERROR: makeblastdb failed. Exiting."
    exit 1
fi
echo "BLAST database created: ${BLAST_DB_DIR}gv_contigs_db"
echo "---"

# --- Step B: Blast all contigs in intermediate/{plv,lc,vph}/* against this database ---
echo "Step B: Combining query contigs and running BLASTN..."

# Combine all .fna files from intermediate/contigs/vph, lc, plv into a single query FASTA
find "${QUERY_CONTIGS_BASE_DIR}" -type d -name "vph" -o -name "lc" -o -name "plv" | \
while read -r subdir; do
    find "${subdir}" -name "*.fna" -print0
done | xargs -0 cat > "${QUERY_COMBINED_FASTA}"

if [ ! -s "${QUERY_COMBINED_FASTA}" ]; then
    echo "ERROR: No .fna files found in ${QUERY_CONTIGS_BASE_DIR}{vph,lc,plv} or combined query FASTA is empty. Exiting."
    exit 1
fi

# Define the output file for BLAST results
BLAST_RESULTS_FILE="${BLAST_OUTPUT_DIR}query_vs_gv_contigs_megablast.tsv"

# Run blastn (megablast task)
blastn -query "${QUERY_COMBINED_FASTA}" \
       -db "${BLAST_DB_DIR}gv_contigs_db" \
       -out "${BLAST_RESULTS_FILE}" \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
       -task megablast \
       -num_threads "${SLURM_CPUS_PER_TASK}" \
       -perc_identity "${MIN_PERCENT_IDENTITY}" \
       -qcov_hsp_perc "${MIN_QUERY_COVERAGE}" \
       -max_target_seqs 1000 # Limit hits per query; adjust if you expect more relevant matches

if [ $? -ne 0 ]; then
    echo "ERROR: blastn failed. Check BLAST output for details."
    exit 1
fi

echo "BLASTN comparison complete. Results saved to: ${BLAST_RESULTS_FILE}"