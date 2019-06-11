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
# (Finished, Reran background without coordinate shifting)
bs.filtered <- readRDS("Dx_All/BSseq_05Group.rds")
bs.filtered <- chrSelectBSseq(bs.filtered, seqnames = c(paste("chr", 1:22, sep = ""), "chrM")) # Remove chrX, chrY
background <- getBackground(bs.filtered, minNumRegion = minCpGs, maxGap = 5000)
background <- subset(background, width >= 5000) # Minimum block size
write.table(background, file = "Dx_All/bsseq_block_background_Replication50.csv", sep = ",", quote = FALSE, row.names = FALSE)

blocks <- dmrseq(bs = bs.filtered, testCovariate = testCovariate, adjustCovariate = NULL, cutoff = cutoff,
                 minNumRegion = minCpGs, bpSpan = 5e4, minInSpan = 500, maxGapSmooth = 1e6, maxGap = 5e3, 
                 maxPerms = maxPerms, block = TRUE)
blocks$percentDifference <- round(blocks$beta/pi * 100)
sigBlocks <- blocks[blocks$pval < 0.05,]
gr2csv(blocks, "CandidateBlocks_DxNoXY_Replication50.csv")
gr2csv(sigBlocks, "DifferentialBlocks_DxNoXY_Replication50.csv")
rm(bs.filtered, blocks, sigBlocks)

# Diagnosis and Sex Block DMRs All Samples ----------------------------------------
# (Finished, Background is with coordinate shifting)
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
# (Finished, Reran background without coordinate shifting)
bs.filtered <- readRDS("Dx_Males/Filtered_BSseq_Replication50_males.rds")
background <- getBackground(bs.filtered, minNumRegion = minCpGs, maxGap = 5000)
background <- subset(background, width >= 5000) # Minimum block size
write.table(background, file = "Dx_Males/bsseq_block_background_Replication50_males.csv", sep = ",", quote = FALSE, row.names = FALSE)

blocks <- dmrseq(bs = bs.filtered, testCovariate = testCovariate, adjustCovariate = NULL, cutoff = cutoff,
                 minNumRegion = minCpGs, bpSpan = 5e4, minInSpan = 500, maxGapSmooth = 1e6, maxGap = 5e3, 
                 maxPerms = maxPerms, block = TRUE)
blocks$percentDifference <- round(blocks$beta/pi * 100)
sigBlocks <- blocks[blocks$pval < 0.05,]
gr2csv(blocks, "CandidateBlocks_Dx_Replication50_males.csv")
gr2csv(sigBlocks, "DifferentialBlocks_Dx_Replication50_males.csv")
rm(bs.filtered, blocks, sigBlocks)

# Diagnosis Block DMRs Females ----------------------------------------
# (Finished, Reran background without coordinate shifting)
bs.filtered <- readRDS("Dx_Females_100/Filtered_BSseq_Replication100_females.rds")
background <- getBackground(bs.filtered, minNumRegion = minCpGs, maxGap = 5000)
background <- subset(background, width >= 5000) # Minimum block size
write.table(background, file = "Dx_Females_100/bsseq_block_background_Replication100_females.csv", sep = ",", quote = FALSE, row.names = FALSE)

blocks <- dmrseq(bs = bs.filtered, testCovariate = testCovariate, adjustCovariate = NULL, cutoff = cutoff,
                 minNumRegion = minCpGs, bpSpan = 5e4, minInSpan = 500, maxGapSmooth = 1e6, maxGap = 5e3, 
                 maxPerms = maxPerms, block = TRUE)
blocks$percentDifference <- round(blocks$beta/pi * 100)
sigBlocks <- blocks[blocks$pval < 0.05,]
gr2csv(blocks, "CandidateBlocks_Dx_Replication100_females.csv")
gr2csv(sigBlocks, "DifferentialBlocks_Dx_Replication100_females.csv")
rm(bs.filtered, blocks, sigBlocks)






