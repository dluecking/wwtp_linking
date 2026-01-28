# Author: dlu @ veelab
# Version: 2025-07-02

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)
library(taxize)

# set working directory
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# setwd("/lisc/data/scratch/dome/willemsen/luecking/projects/wwtp_linking")



# load all data -----------------------------------------------------------

data <- rbindlist(lapply(list.files("intermediate/lc_tax/", pattern = "_filtered.csv", full.names = TRUE), fread))

# set key to retrieve taxinfo
ENTREZ_KEY="ed38a68e4cac02507e4bc585e8913bab5a08"
Sys.setenv(ENTREZ_KEY = ENTREZ_KEY)


# get taxinfo and save to intermediate, load if exists --------------------

# add taxonomy but only if we cant load it from locally
if(!file.exists("intermediate/lc_tax/lc_tax_info_df.csv")){
  tax_info <- classification(unique(data$majority_organism), db = "ncbi")
  saveRDS(tax_info, file = "intermediate/lc_tax/lc_tax_info_df.csv")
}else{
  tax_info <- readRDS("intermediate/lc_tax/lc_tax_info_df.csv")
}
