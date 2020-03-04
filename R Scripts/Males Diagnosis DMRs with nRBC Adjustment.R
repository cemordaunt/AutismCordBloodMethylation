# Males Diagnosis DMRs with nRBC Adjustment -------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 3/4/20

# Load Packages and Setup R Session ####
.libPaths("/share/lasallelab/programs/DMRichR/R_3.6")
sapply(c("tidyverse", "openxlsx", "bsseq", "dmrseq", "DMRichR"), require, character.only = TRUE)
set.seed(5)
register(MulticoreParam(1))

# Males Discovery DMRs with nRBC Adjustment --------------------------
# Add nRBCs to pData ####
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Discovery")
bs.filtered <- readRDS("Dx_Males/Filtered_BSseq_Discovery50_males.rds")
samples <- read.csv("MARBLES EARLI WGBS Sample Merged Covariate Database with Demo and Cell Type.csv", stringsAsFactors = FALSE)
samples <- samples[,c("Sequencing_ID", "nRBC")]
pData <- pData(bs.filtered)
table(samples$Sequencing_ID[match(rownames(pData), samples$Sequencing_ID)] == rownames(pData)) # TRUE
pData$nRBC <- samples$nRBC[match(rownames(pData), samples$Sequencing_ID)]
pData(bs.filtered) <- pData

# Call DMRs ####
regions <- dmrseq(bs = bs.filtered, cutoff = 0.05, minNumRegion = 3, maxPerms = 10, testCovariate = "Diagnosis", 
                  adjustCovariate = "nRBC", matchCovariate = NULL)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "Dx_Males/CandidateRegions_Dx_nRBC_Discovery50_males.csv")
gr2csv(sigRegions, "Dx_Males/DMRs_Dx_nRBC_Discovery50_males.csv")

# Males Replication DMRs with nRBC Adjustment --------------------------
# Add nRBCs to pData ####
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Replication")
bs.filtered <- readRDS("Dx_Males/Filtered_BSseq_Replication50_males.rds")
samples <- read.csv("MARBLES EARLI WGBS Sample Merged Covariate Database with Demo and Cell Type.csv", stringsAsFactors = FALSE)
samples <- samples[,c("Sequencing_ID", "nRBC")]
pData <- pData(bs.filtered)
table(samples$Sequencing_ID[match(rownames(pData), samples$Sequencing_ID)] == rownames(pData)) # TRUE
pData$nRBC <- samples$nRBC[match(rownames(pData), samples$Sequencing_ID)]
pData(bs.filtered) <- pData

# Call DMRs ####
regions <- dmrseq(bs = bs.filtered, cutoff = 0.05, minNumRegion = 3, maxPerms = 10, testCovariate = "Diagnosis", 
                  adjustCovariate = "nRBC", matchCovariate = NULL)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "Dx_Males/CandidateRegions_Dx_nRBC_Replication50_males.csv")
gr2csv(sigRegions, "Dx_Males/DMRs_Dx_nRBC_Replication50_males.csv")
