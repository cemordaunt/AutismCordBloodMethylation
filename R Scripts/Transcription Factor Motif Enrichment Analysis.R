# Transcription Factor Motif Enrichment Analysis -----------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 7/9/19

# Packages ####
sapply(c("tidyverse"), require, character.only = TRUE)

# Functions ####
source("R Scripts/DMR Analysis Functions.R")

# Make BED Input Files ----------------------------------------------------
# Load Regions ####
chromsXY <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")
chromsX <- c(paste("chr", 1:22, sep = ""), "chrX", "chrM")
DMRs <- list(MalesDisc = loadRegions("DMRs/Discovery/Diagnosis Males 50/DMRs_Dx_Discovery50_males.csv",
                                     chroms = chromsXY, sort = TRUE, DMRid = TRUE),
             FemalesDisc = loadRegions("DMRs/Discovery/Diagnosis Females 50/DMRs_Dx_Discovery50_females.csv",
                                       chroms = chromsX, sort = TRUE, DMRid = TRUE),
             MalesRep = loadRegions("DMRs/Replication/Diagnosis Males 50/DMRs_Dx_Replication50_males.csv",
                                    chroms = chromsXY, sort = TRUE, DMRid = TRUE),
             FemalesRep = loadRegions("DMRs/Replication/Diagnosis Females 100/DMRs_Dx_Replication100_females.csv",
                                      chroms = chromsX, sort = TRUE, DMRid = TRUE))
Background <- list(MalesDisc = loadRegions("DMRs/Discovery/Diagnosis Males 50/bsseq_background_Discovery50_males.csv",
                                           chroms = chromsXY, sort = TRUE),
                   FemalesDisc = loadRegions("DMRs/Discovery/Diagnosis Females 50/bsseq_background_Discovery50_females.csv",
                                             chroms = chromsX, sort = TRUE),
                   MalesRep = loadRegions("DMRs/Replication/Diagnosis Males 50/bsseq_background_Replication50_males.csv",
                                          chroms = chromsXY, sort = TRUE),
                   FemalesRep = loadRegions("DMRs/Replication/Diagnosis Females 100/bsseq_background_Replication100_females.csv",
                                            chroms = chromsX, sort = TRUE))
HyperDMRs <- lapply(DMRs, subset, percentDifference > 0)
HypoDMRs <- lapply(DMRs, subset, percentDifference < 0)

# Write BED Files ####
mapply(writeBED, regions = DMRs, 
       file = paste("UCSC Tracks/", c("Discovery_Males", "Discovery_Females", "Replication_Males", "Replication_Females"),
                    "_Diagnosis_DMRs.bed", sep = ""))
mapply(writeBED, regions = HyperDMRs, 
       file = paste("UCSC Tracks/", c("Discovery_Males", "Discovery_Females", "Replication_Males", "Replication_Females"),
                    "_Diagnosis_Hyper_DMRs.bed", sep = ""))
mapply(writeBED, regions = HypoDMRs, 
       file = paste("UCSC Tracks/", c("Discovery_Males", "Discovery_Females", "Replication_Males", "Replication_Females"),
                    "_Diagnosis_Hypo_DMRs.bed", sep = ""))
mapply(writeBED, regions = Background, 
       file = paste("UCSC Tracks/", c("Discovery_Males", "Discovery_Females", "Replication_Males", "Replication_Females"),
                    "_background.bed", sep = ""))

# Get HOMER Enrichments ---------------------------------------------------
# HOMER Version v4.10
# See Shell Scripts/HOMER_TF_Motif_Enrichment.sh
# Parameters
# genome = hg38
# size = given (size of regions for motifs)
# cpg (Normalized for % CpG content)
# N = 100000 (Randomly pick 100K background regions for comparison, backgrounds are 200K - 400K)


