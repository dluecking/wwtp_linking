#! /bin/bash

module load conda
conda activate aster

# concat the nuphylo trees
# cat data/tree_data/*.nwk > combined.nwk
# this was done previously

# rename the different genes in the tree
# sed -E 's/([^[:space:];(),:]+)\.([^[:space:];(),:]+)\.polished_[^[:space:];(),:]+\.copy[0-9]+/\[\1_\2\]/g' combined.nwk > combined_renamed.nwk
# also done beforehand

# run the astral command
astral -i combined.nwk -o combined_astra_tree.tre


