#!/bin/bash
# Load the required module
module load ncbiblastplus/2.16.0

mkdir -p intermediate/non_CRISPR/blastn_out

# Loop through all contig files in the folder
for contig in data/gv_contigs/*.fasta; do
    # Extract base name without path and extension
    base_name=$(basename "$contig" .fasta)

    # Run blastn
    blastn -query "$contig" \
           -db intermediate/non_CRISPR/vph_plv_chunks_db \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs scovs" \
           -num_threads 8 \
           -evalue 1e-3 \
           -out "intermediate/non_CRISPR/blastn_out/${base_name}_vs_vph_plv_chunks.out"

    echo "BLAST for $contig completed."
done

echo "All BLAST jobs finished!"