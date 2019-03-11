# Diagnosis and Sex Block DMRs Replication ----------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 3/11/19

# Packages ####
.libPaths("/share/lasallelab/programs/DMRichR/R_3.5")
sapply(c("tidyverse", "openxlsx", "bsseq", "dmrseq", "DMRichR"), require, character.only = TRUE)

# Global Variables ####
minCpGs <- as.numeric(3)
maxPerms <- as.numeric(10)
testCovariate <- as.character("Diagnosis")
set.seed(5)
register(MulticoreParam(1))
cutoff <- as.numeric(0.01)

# Diagnosis Block DMRs All Samples ----------------------------------------
# (Running on Barbera 3/11)
bs.filtered <- readRDS("BSseq_05Group.rds")
bs.filtered <- chrSelectBSseq(bs.filtered, seqnames = c(paste("chr", 1:22, sep = ""), "chrM")) # Remove chrX, chrY
blocks <- dmrseq(bs = bs.filtered, testCovariate = testCovariate, adjustCovariate = NULL, cutoff = cutoff,
                 minNumRegion = minCpGs, bpSpan = 5e4, minInSpan = 500, maxGapSmooth = 1e6, maxGap = 5e3, 
                 maxPerms = maxPerms, block = TRUE)
blocks$percentDifference <- round(blocks$beta/pi * 100)
sigBlocks <- blocks[blocks$pval < 0.05,]
gr2csv(blocks, "CandidateBlocks_DxNoXY_Replication50.csv")
gr2csv(sigBlocks, "DifferentialBlocks_DxNoXY_Replication50.csv")
rm(bs.filtered, blocks, sigBlocks)

# Diagnosis and Sex Block DMRs All Samples ----------------------------------------
# (Running on Barbera 3/11)
bs.filtered <- readRDS("BSseq_05Group.rds")
blocks <- dmrseq(bs = bs.filtered, testCovariate = testCovariate, adjustCovariate = "Sex", cutoff = cutoff,
                 minNumRegion = minCpGs, bpSpan = 5e4, minInSpan = 500, maxGapSmooth = 1e6, maxGap = 5e3, 
                 maxPerms = maxPerms, block = TRUE)
blocks$percentDifference <- round(blocks$beta/pi * 100)
sigBlocks <- blocks[blocks$pval < 0.05,]
gr2csv(blocks, "CandidateBlocks_DxAdjSex_Replication50.csv")
gr2csv(sigBlocks, "DifferentialBlocks_DxAdjSex_Replication50.csv")
rm(bs.filtered, blocks, sigBlocks)

# Diagnosis Block DMRs Males ----------------------------------------
# (Running on Barbera 3/11)
bs.filtered <- readRDS("Filtered_BSseq_Replication50_males.rds")
blocks <- dmrseq(bs = bs.filtered, testCovariate = testCovariate, adjustCovariate = NULL, cutoff = cutoff,
                 minNumRegion = minCpGs, bpSpan = 5e4, minInSpan = 500, maxGapSmooth = 1e6, maxGap = 5e3, 
                 maxPerms = maxPerms, block = TRUE)
blocks$percentDifference <- round(blocks$beta/pi * 100)
sigBlocks <- blocks[blocks$pval < 0.05,]
gr2csv(blocks, "CandidateBlocks_Dx_Replication50_males.csv")
gr2csv(sigBlocks, "DifferentialBlocks_Dx_Replication50_males.csv")
rm(bs.filtered, blocks, sigBlocks)

# Diagnosis Block DMRs Females ----------------------------------------
# (Running on Barbera 3/11)
bs.filtered <- readRDS("Filtered_BSseq_Replication100_females.rds")
blocks <- dmrseq(bs = bs.filtered, testCovariate = testCovariate, adjustCovariate = NULL, cutoff = cutoff,
                 minNumRegion = minCpGs, bpSpan = 5e4, minInSpan = 500, maxGapSmooth = 1e6, maxGap = 5e3, 
                 maxPerms = maxPerms, block = TRUE)
blocks$percentDifference <- round(blocks$beta/pi * 100)
sigBlocks <- blocks[blocks$pval < 0.05,]
gr2csv(blocks, "CandidateBlocks_Dx_Replication50_females.csv")
gr2csv(sigBlocks, "DifferentialBlocks_Dx_Replication50_females.csv")
rm(bs.filtered, blocks, sigBlocks)






