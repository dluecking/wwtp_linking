#!/bin/bash

module load Bowtie2/2.5.4-GCC-13.3.0

CONTIG_FILE="intermediate/contigs/all_contigs.fna"
INDEX_FILE="intermediate/contigs/all_contigs_index"

# build reference
bowtie2-build $CONTIG_FILE $INDEX_FILE