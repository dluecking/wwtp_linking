# Author: dlu @ veelab
# Version: 2025-07-30

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



# iterate over files, rename and save to renamed --------------------------

NCLDV_OUTPUT_DIR <- "intermediate/GV_tree/ncldv_output"

for(file in list.files(NCLDV_OUTPUT_DIR, pattern = "gvs_and_yara_GVOG.*.faa$", full.names = TRUE)){
  old_fasta <- read.fasta(file)
  
  # this gets the annotation for each sequence, extracts the second name (which is Hade_1 for example)
  # then remove everything after the last "_" (could be polished, could be gene number...)
  write.fasta(sequences = getSequence(old_fasta), 
              names = str_remove(sapply(strsplit(as.character(unlist(getAnnot(old_fasta))), " "), function(x) x[2]), "\\_[^\\_]*$"),
              file.out = str_replace(file, 
                                     pattern = "gvs_and_yara_", 
                                     replacement = "gvs_and_yara_renamed_"))
}







