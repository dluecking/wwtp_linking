# Author: dlu @ veelab
# Version: 2025-08-26

# Packages
library(dplyr)
library(seqinr)
library(ggplot2)
library(data.table)
library(stringr)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



# load general data -------------------------------------------------------

# URL of the Google Sheet
sheet_url <- "https://docs.google.com/spreadsheets/d/1QLNiqSt0XOS4xVPAeZAppwVjjjPIKdEE6w6f2_Qm55c"

# Read the specified sheet and convert to data.table
gv_data <- read_sheet(sheet_url, sheet = "Final GVs overview") %>% 
  filter(completeness == "complete")


