# Window Methylation Analysis ####
# Autism Cord Blood Methylation Project
# Charles Mordaunt
# 9/24/18

# Packages ####

# Functions ####

# Data ####
samples <- read.delim(file = "Samples/Merged Cord Blood WGBS Database.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
permeth <- read.delim(file = "Tables/windows_10kb_methylation_ASD_CordBlood.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
