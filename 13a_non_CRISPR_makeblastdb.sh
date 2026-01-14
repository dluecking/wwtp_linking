#!/bin/bash

module load BLAST+/2.17.0-gompi-2024a


cat intermediate/non_CRISPR/chunks/* > intermediate/non_CRISPR/combined_chunks.fa

makeblastdb -in intermediate/non_CRISPR/combined_chunks.fa \
            -dbtype nucl \
            -out intermediate/non_CRISPR/vph_plv_chunks_db \
            -title "40 bp chunks of vphs plvs BLAST DB"


echo "BLAST database created at intermediate/non_CRISPR/vph_plv_chunks_db"