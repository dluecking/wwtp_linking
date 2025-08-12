#!/bin/bash

ANNOMAZING_SLURM="/lisc/home/user/luecking/bin/AnnoMazing/run_annomazing.slurm"
HELPERFILE="helperfiles/list_of_gv_protein_files.txt"

mkdir -p intermediate/gv_annotations/

for file in `cat ${HELPERFILE}`; do
    echo "Processing $file"
    current_GV=$(basename $file .faa)
    sbatch ${ANNOMAZING_SLURM} intermediate/proteins/gv/${file} intermediate/gv_annotations/${current_GV}_annotation.csv
done
    


