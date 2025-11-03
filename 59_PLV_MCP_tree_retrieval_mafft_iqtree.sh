#!/bin/bash
module load conda
module load mafft
module load iqtree

# Filter sequences containing "_plv_" in the header
grep -A 1 "_plv_" intermediate/proteins/all_public_and_my_own_MCP_proteins.faa > intermediate/treebuilding/plv_tree/PLV_MCP_proteins.faa

# Ensure the number of sequences is 14
seq_count=$(grep -c "^>" intermediate/treebuilding/plv_tree/PLV_MCP_proteins.faa)
if [ "$seq_count" -ne 14 ]; then
    echo "Error: Number of sequences is not 14. Found $seq_count."
    exit 1
fi

# Align sequences using MAFFT
mafft --auto intermediate/treebuilding/plv_tree/PLV_MCP_proteins.faa > intermediate/treebuilding/plv_tree/PLV_MCP_proteins_aligned.faa

# Build a tree using IQ-TREE2 with 1000 bootstrap replicates
iqtree2 -s intermediate/treebuilding/plv_tree/PLV_MCP_proteins_aligned.faa -T 4 --prefix intermediate/treebuilding/plv_tree/PLV_MCP_tree -B 1000