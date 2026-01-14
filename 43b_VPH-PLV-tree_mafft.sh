#!/bin/bash

module load MAFFT/7.526-GCC-13.3.0-with-extensions

mafft --auto intermediate/proteins/all_public_and_my_own_MCP_proteins.faa > intermediate/treebuilding/all_public_and_my_own_MCP_proteins.aln
