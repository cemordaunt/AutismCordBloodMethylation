# Transcription Factor Motif Enrichment Analysis -----------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 7/11/19

# Packages ####
sapply(c("tidyverse", "GenomicRanges", "LOLA"), require, character.only = TRUE)

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

# Redefine DMRs in Terms of Background ####
GR_Background <- lapply(Background, makeGRange, direction = "all")

reDMRs <- mapply(function(x, y){
        makeGRange(x, direction = "all") %>% GRangesList() %>% redefineUserSets(userUniverse = y) %>% as.data.frame()
}, x = DMRs, y = GR_Background, SIMPLIFY = FALSE)

reHyperDMRs <- mapply(function(x, y){
        makeGRange(x, direction = "hyper") %>% GRangesList() %>% redefineUserSets(userUniverse = y) %>% as.data.frame()
}, x = DMRs, y = GR_Background, SIMPLIFY = FALSE)

reHypoDMRs <- mapply(function(x, y){
        makeGRange(x, direction = "hypo") %>% GRangesList() %>% redefineUserSets(userUniverse = y) %>% as.data.frame()
}, x = DMRs, y = GR_Background, SIMPLIFY = FALSE)

# Write BED Files ####
mapply(writeBED, regions = reDMRs, 
       file = paste("UCSC Tracks/", c("Discovery_Males", "Discovery_Females", "Replication_Males", "Replication_Females"),
                    "_Diagnosis_DMRs_HOMER.bed", sep = ""))
mapply(writeBED, regions = reHyperDMRs, 
       file = paste("UCSC Tracks/", c("Discovery_Males", "Discovery_Females", "Replication_Males", "Replication_Females"),
                    "_Diagnosis_Hyper_DMRs_HOMER.bed", sep = ""))
mapply(writeBED, regions = reHypoDMRs, 
       file = paste("UCSC Tracks/", c("Discovery_Males", "Discovery_Females", "Replication_Males", "Replication_Females"),
                    "_Diagnosis_Hypo_DMRs_HOMER.bed", sep = ""))
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
# Queued on barbera 7/11

# Load Results ####
males_disc <- list(All = loadHOMER("DMRs/Discovery/Diagnosis Males 50/HOMER_All/knownResults.txt"), 
                   Hyper = loadHOMER("DMRs/Discovery/Diagnosis Males 50/HOMER_Hyper/knownResults.txt"),
                   Hypo = loadHOMER("DMRs/Discovery/Diagnosis Males 50/HOMER_Hypo/knownResults.txt"))
males_rep <- list(All = loadHOMER("DMRs/Replication/Diagnosis Males 50/HOMER_All/knownResults.txt"), 
                  Hyper = loadHOMER("DMRs/Replication/Diagnosis Males 50/HOMER_Hyper/knownResults.txt"),
                  Hypo = loadHOMER("DMRs/Replication/Diagnosis Males 50/HOMER_Hypo/knownResults.txt"))
females_disc <- list(All = loadHOMER("DMRs/Discovery/Diagnosis Females 50/HOMER_All/knownResults.txt"), 
                     Hyper = loadHOMER("DMRs/Discovery/Diagnosis Females 50/HOMER_Hyper/knownResults.txt"),
                     Hypo = loadHOMER("DMRs/Discovery/Diagnosis Females 50/HOMER_Hypo/knownResults.txt"))
females_rep <- list(All = loadHOMER("DMRs/Replication/Diagnosis Females 100/HOMER_All/knownResults.txt"), 
                    Hyper = loadHOMER("DMRs/Replication/Diagnosis Females 100/HOMER_Hyper/knownResults.txt"),
                    Hypo = loadHOMER("DMRs/Replication/Diagnosis Females 100/HOMER_Hypo/knownResults.txt"))

# Workflow ####
# Melt into data.frame
# Subset enriched motifs (cutoffs?)
# Look at distribution of all data
test <- loadHOMER("DMRs/Discovery/Diagnosis Males 50/HOMER/knownResults.txt")
hist(test$log_qvalue, breaks = 100)
quantile(test$log_qvalue, c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 1))
#      0%      25%      50%      75%      90%      95%    97.5%      99%     100% 
# 0.00000 12.88950 16.70934 19.44746 22.07320 22.83016 23.52913 26.72516 29.38300 

hist(test$Fold_Enrichment, breaks = 100)
quantile(test$Fold_Enrichment, c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 1))
#       0%      25%      50%      75%      90%      95%    97.5%      99%     100% 
# 1.049823 1.210628 1.341379 1.518910 1.696113 1.824599 1.960478 2.130395 2.661098 

table(test$Fold_Enrichment > 1.5, test$log_qvalue > 20) # near 3rd quartile of each
#       FALSE TRUE
# FALSE   227   69
# TRUE     96   22
# 22 / 414 motifs (5.3% of total)

test$Transcription_Factor[test$Fold_Enrichment > 1.5 & test$log_qvalue > 20]
#  [1] "WT1"         "ZIC3"        "TLX"         "ZNF415"      "ZNF692"      "HIF-1A"      "FXR"         "IRF:BATF"   
#  [9] "NUR77"       "EGR1"        "ELF1"        "ETS1-Distal" "NF1"         "PAX6"        "TBX20"       "ETS"        
# [17] "ZNF675"      "VDR"         "HIF2A"       "KLF3"        "TBOX:SMAD"   "ERE"      

# Subset for motifs enriched in discovery and replication, by direction

# Heatmap of replicated motifs, clustered by log pvalues










