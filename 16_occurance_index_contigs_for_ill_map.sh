#!/bin/bash

module load bowtie2

CONTIG_FILE="intermediate/contigs/all_contigs.fna"
INDEX_FILE="intermediate/contigs/all_contigs_index"

# build reference
bowtie2-build $CONTIG_FILE $INDEX_FILE