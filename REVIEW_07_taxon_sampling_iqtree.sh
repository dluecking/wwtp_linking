#!/bin/bash

module load IQ-TREE/2.4.0

marker="GVOGm0760"

NUPHYLO_dir="intermediate/GV_tree/REVIEW_taxon_sampling/nuphylo"
IQTREE_dir="intermediate/GV_tree/REVIEW_taxon_sampling/iqtree"

mkdir -p "$IQTREE_dir"

aln="${NUPHYLO_dir}/${marker}/allseqs.trimmed.aln"

iqtree2 -T 8 -s "$aln" -B 1000 --prefix "${IQTREE_dir}/${marker}" -redo

cp "${IQTREE_dir}/${marker}.treefile" "${IQTREE_dir}/${marker}.nwk"
