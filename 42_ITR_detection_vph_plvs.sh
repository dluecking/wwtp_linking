#!/bin/bash


# STEP 1 retrieve ends of PLVs and VPHs
echo "Retrieving ends of PLVs and VPHs..."
module load conda
conda activate R

Rscript scripts/retrieve_VPH_PLV_ends.R intermediate/ITRs/vph_plv_ends.fasta
echo "Ends of PLVs and VPHs retrieved successfully."


# STEP 2 create blast db
echo "Creating BLAST database for VPH/PLV ends..."
module load ncbiblastplus/2.16.0

makeblastdb -in intermediate/ITRs/vph_plv_ends.fasta \
            -dbtype nucl \
            -out intermediate/ITRs/vph_plv_ends_DB \
            -title "1000 bp ends of vph/plvs BLAST DB"
echo "BLAST database created successfully."


# STEP 3 blast ends against each other
echo "Running BLAST to find ITRs in VPHs and PLVs..."

blastn -query intermediate/ITRs/vph_plv_ends.fasta \
       -db intermediate/ITRs/vph_plv_ends_DB \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs scovs" \
       -num_threads 2 \
       -evalue 1e-3 \
       -word_size 7 \
       -strand minus \
       -out "intermediate/ITRs/ends_vs_ends.out"
echo "BLAST completed successfully."


# STEP 4 filter blast results
echo "Filtering BLAST results to identify ITRs in VPHs and PLVs..."

conda activate R
Rscript scripts/filter_VPH_PLV_ITR_blast_results.R intermediate/ITRs/ends_vs_ends.out intermediate/ITRs/vph_plv_ITR_info.tsv

echo "Filtering completed successfully. ITRs information saved to intermediate/ITRs/vph_plv_ITR_info.tsv."
