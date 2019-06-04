# Diagnosis and Sex DMRs Discovery ####
# Autism Cord Blood Methylation
# Charles Mordaunt
# 5/31/19
# Excluded JLCM032B and JLCM050B

# cmordaunt@epigenerate:/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Discovery
# THREADS=${SLURM_NTASKS}
# MEM=$(expr ${SLURM_MEM_PER_CPU} / 1024)
# echo "Allocated threads: " $THREADS
# echo "Allocated memory: " $MEM

# Load Packages ####
.libPaths("/share/lasallelab/programs/DMRichR/R_3.6")
sapply(c("tidyverse", "openxlsx", "bsseq", "dmrseq", "DMRichR"), require, character.only = TRUE)

# Global Variables ####
genome <- as.character("hg38")
coverage <- as.numeric(1)
perGroup <- (as.numeric(50)/100)
minCpGs <- as.numeric(3)
maxPerms <- as.numeric(10)
testCovariate <- as.character("Diagnosis")
cores <- 5
set.seed(5)
register(MulticoreParam(1))

# Annotation Databases (Done) ####
sapply(c("BSgenome.Hsapiens.UCSC.hg38", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"), require, 
       character.only = TRUE)
goi <- BSgenome.Hsapiens.UCSC.hg38
TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene
annoDb <- "org.Hs.eg.db"
annoTrack <- getAnnot(genome)
saveRDS(annoTrack, "hg38_annoTrack.rds")

# Meth ~ Diagnosis, exclude chrX and Y, DMRs ####
# Load and Process Samples (Rerun without JLCM032B and JLCM050B, Complete)
bs.filtered <- processBismark(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                              meta = read.xlsx("sample_info.xlsx", colNames = TRUE) %>% mutate_if(is.character,as.factor),
                              groups = testCovariate, Cov = coverage, mc.cores = cores, per.Group = perGroup)
bs.filtered <- chrSelectBSseq(bs.filtered, seqnames = c(paste("chr", 1:22, sep = ""), "chrM")) # Remove chrX, chrY
saveRDS(bs.filtered, "Dx_All/Filtered_BSseq_Discovery50.rds")

# Background Regions (Rerun without JLCM032B and JLCM050B, Complete)
background <- getBackground(bs.filtered, minNumRegion = minCpGs, maxGap = 1000)
write.table(background, file = "Dx_All/bsseq_background_Discovery50.csv", sep = ",", quote = FALSE, row.names = FALSE)

# DMRs and Raw Methylation (Rerun without JLCM032B and JLCM050B, Complete)
regions <- dmrseq(bs = bs.filtered, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms, testCovariate = testCovariate)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "Dx_All/CandidateRegions_DxNoXY_Discovery50.csv")
gr2csv(sigRegions, "Dx_All/DMRs_DxNoXY_Discovery50.csv")
raw <- as.data.frame(getMeth(BSseq = bs.filtered, regions = sigRegions, type = "raw", what = "perRegion"))
raw <- cbind(sigRegions, raw)
write.table(raw, "Dx_All/DMR_raw_methylation_DxNoXY_Discovery50.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Smoothed Methylation (Rerun without JLCM032B and JLCM050B, Complete)
bs.filtered.bsseq <- BSmooth(bs.filtered, BPPARAM = MulticoreParam(workers = cores, progressbar = TRUE))
pData <- pData(bs.filtered.bsseq)
pData$col <- NULL
pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "#3366CC"
pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <-  "#FF3366"
pData$label <- NULL
pData$label[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "TD"
pData$label[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <-  "ASD"
pData(bs.filtered.bsseq) <- pData
saveRDS(bs.filtered.bsseq, "Dx_All/Filtered_Smoothed_BSseq_Discovery50.rds")
smoothed <- getSmooth(bsseq = bs.filtered.bsseq, regions = sigRegions)
write.table(smoothed, "Dx_All/DMR_smoothed_methylation_DxNoXY_Discovery50.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Plots (Rerun without JLCM032B and JLCM050B, Complete)
annoTrack <- readRDS("hg38_annoTrack.rds")
pdf("Dx_All/DMRs_DxNoXY_Discovery50.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate, 
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2,
          addRegions = sigRegions, annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE,
          horizLegend = TRUE)
dev.off()

pdf("Dx_All/DMRs_DxNoXY_Discovery50_noPoints.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate, 
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2,
          addRegions = sigRegions, annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE,
          horizLegend = TRUE, addPoints = FALSE)
dev.off()

# Meth ~ Diagnosis + AdjSex DMRs and Plots ####
# (Includes JLCM032B and JLCM050B)
# New R session
bs.filtered <- readRDS("Filtered_BSseq_Discovery50.rds")
regions <- dmrseq(bs = bs.filtered, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms,
                  testCovariate = testCovariate, adjustCovariate = "Sex", matchCovariate = NULL)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "CandidateRegions_DxAdjSex_Discovery50.csv")
gr2csv(sigRegions, "DMRs_DxAdjSex_Discovery50.csv")

bs.filtered.bsseq <- readRDS("Filtered_Smoothed_BSseq_Discovery50.rds")
annoTrack <- readRDS("hg38_annoTrack.rds")
smoothed <- getSmooth(bsseq = bs.filtered.bsseq, regions = sigRegions, 
                      out = "DMR_smoothed_methylation_DxAdjSex_Discovery50.txt")

pdf("DMRs_DxAdjSex_Discovery50.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate, 
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2,
          addRegions = sigRegions, annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE,
          horizLegend = TRUE)
dev.off()

pdf("DMRs_DxAdjSex_Discovery50_noPoints.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate, 
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2,
          addRegions = sigRegions, annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE,
          horizLegend = TRUE, addPoints = FALSE)
dev.off()

# Get raw DMR methylation (Done)
bs.filtered <- readRDS("Filtered_BSseq_Discovery50.rds")
sigRegions <- read.csv("DMRs_DxAdjSex_Discovery50.csv", header = TRUE, stringsAsFactors = FALSE)
raw <- as.data.frame(getMeth(BSseq = bs.filtered, regions = data.frame2GRanges(sigRegions), type = "raw", what = "perRegion"))
raw <- cbind(sigRegions, raw)
write.table(raw, "DMR_raw_methylation_DxAdjSex_Discovery50.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# DMRs with Males Only (Rerun without JLCM032B and JLCM050B) ####
# New R session
# Load and Process Samples (Rerun without JLCM032B and JLCM050B, Complete)
bs.filtered <- processBismark(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                              meta = read.xlsx("sample_info_males.xlsx", colNames = TRUE) %>% mutate_if(is.character,as.factor),
                              groups = testCovariate, Cov = coverage, mc.cores = cores, per.Group = perGroup)
saveRDS(bs.filtered, "Dx_Males/Filtered_BSseq_Discovery50_males.rds")

# Background Regions (Rerun without JLCM032B and JLCM050B, Complete)
background <- getBackground(bs.filtered, minNumRegion = minCpGs, maxGap = 1000)
write.table(background, file = "Dx_Males/bsseq_background_Discovery50_males.csv", sep = ",", quote = FALSE, row.names = FALSE)

# Meth ~ Diagnosis DMRs and Raw Methylation (Rerun without JLCM032B and JLCM050B, Complete)
regions <- dmrseq(bs = bs.filtered, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms,
                  testCovariate = testCovariate, adjustCovariate = NULL, matchCovariate = NULL)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "Dx_Males/CandidateRegions_Dx_Discovery50_males.csv")
gr2csv(sigRegions, "Dx_Males/DMRs_Dx_Discovery50_males.csv")
raw <- as.data.frame(getMeth(BSseq = bs.filtered, regions = sigRegions, type = "raw", what = "perRegion"))
raw <- cbind(sigRegions, raw)
write.table(raw, "Dx_Males/DMR_raw_methylation_Dx_Discovery50_males.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Smoothed Methylation (Rerun without JLCM032B and JLCM050B, Complete)
bs.filtered.bsseq <- BSmooth(bs.filtered, BPPARAM = MulticoreParam(workers = cores, progressbar = TRUE))
pData <- pData(bs.filtered.bsseq)
pData$col <- NULL
pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "#3366CC"
pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <-  "#FF3366"
pData$label <- NULL
pData$label[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "TD"
pData$label[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <-  "ASD"
pData(bs.filtered.bsseq) <- pData
saveRDS(bs.filtered.bsseq, "Dx_Males/Filtered_Smoothed_BSseq_Discovery50_males.rds")
smoothed <- getSmooth(bsseq = bs.filtered.bsseq, regions = sigRegions)
write.table(smoothed, "Dx_Males/DMR_smoothed_methylation_Dx_Discovery50_males.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# Plots (Rerun without JLCM032B and JLCM050B, Complete)
annoTrack <- readRDS("hg38_annoTrack.rds")
pdf("Dx_Males/DMRs_Dx_Discovery50_males_noPoints.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate,
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2, addRegions = sigRegions, 
          annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE, horizLegend = TRUE, addPoints = FALSE)
dev.off()

pdf("Dx_Males/DMRs_Dx_Discovery50_males.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate,
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2, addRegions = sigRegions, 
          annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE, horizLegend = TRUE)
dev.off()

# DMRs with Females Only ####
# New R session
# Load and Process Samples (Done)
name <- gsub( "_.*$","", list.files(path = getwd(), pattern = "*.txt.gz"))
bs.filtered <- processBismark(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                              meta = read.xlsx("sample_info_females.xlsx", colNames = TRUE) %>% mutate_if(is.character,as.factor),
                              groups = testCovariate, Cov = coverage, mc.cores = cores, per.Group = perGroup)
saveRDS(bs.filtered, "Filtered_BSseq_Discovery50_females.rds")

# Background Regions (Done)
background <- getBackground(bs.filtered, minNumRegion = minCpGs, maxGap = 1000)
write.table(background, file = "bsseq_background_Discovery50_females.csv", sep = ",", quote = FALSE, row.names = FALSE)

# Meth ~ Diagnosis DMRs (Done)
regions <- dmrseq(bs = bs.filtered, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms,
                  testCovariate = testCovariate, adjustCovariate = NULL, matchCovariate = NULL)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "CandidateRegions_Dx_Discovery50_females.csv")
gr2csv(sigRegions, "DMRs_Dx_Discovery50_females.csv")

# Smoothed Methylation and Plots (Done)
bs.filtered.bsseq <- BSmooth(bs.filtered, BPPARAM = MulticoreParam(workers = cores, progressbar = TRUE))
pData <- pData(bs.filtered.bsseq)
pData$col <- NULL
pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "#3366CC"
pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <-  "#FF3366"
pData$label <- NULL
pData$label[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "TD"
pData$label[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <-  "ASD"
pData(bs.filtered.bsseq) <- pData
saveRDS(bs.filtered.bsseq, "Filtered_Smoothed_BSseq_Discovery50_females.rds")

smoothed <- getSmooth(bsseq = bs.filtered.bsseq, regions = sigRegions, 
                      out = "DMR_smoothed_methylation_Dx_Discovery50_females.txt")

bs.filtered.bsseq <- readRDS("Filtered_Smoothed_BSseq_Discovery50_females.rds")
testCovariate <- "Diagnosis"
sigRegions <- read.csv("DMRs_Dx_Discovery50_females.csv", header = TRUE, stringsAsFactors = FALSE)
sigRegions <- data.frame2GRanges(sigRegions)
annoTrack <- readRDS("hg38_annoTrack.rds")
# pasted in plotFunctions.R
pdf("DMRs_Dx_Discovery50_females_noPoints.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate, 
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2, addRegions = sigRegions, 
          annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE, horizLegend = TRUE, 
          addPoints = FALSE)
dev.off()

pdf("DMRs_Dx_Discovery50_females.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate, 
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2, addRegions = sigRegions, 
          annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE, horizLegend = TRUE, 
          addPoints = TRUE)
dev.off()

# Get raw DMR methylation (Done)
bs.filtered <- readRDS("Filtered_BSseq_Discovery50_females.rds")
sigRegions <- read.csv("DMRs_Dx_Discovery50_females.csv", header = TRUE, stringsAsFactors = FALSE)
raw <- as.data.frame(getMeth(BSseq = bs.filtered, regions = data.frame2GRanges(sigRegions), type = "raw", what = "perRegion"))
raw <- cbind(sigRegions, raw)
write.table(raw, "DMR_raw_methylation_Dx_Discovery50_females.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# DMRs with EARLI Samples Only ####
# Includes JLCM032B
# Load and Process Samples (Done)
cores <- 5
name <- gsub( "_.*$","", list.files(path = getwd(), pattern = "*.txt.gz"))
bs.filtered <- processBismark(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                              meta = read.xlsx("sample_info_EARLI.xlsx", colNames = TRUE) %>% mutate_if(is.character,as.factor),
                              groups = testCovariate, Cov = coverage, mc.cores = cores, per.Group = perGroup)
saveRDS(bs.filtered, "Filtered_BSseq_Discovery50_EARLI.rds") # Includes chrX, chrY

# Background Regions (Done)
background <- getBackground(bs.filtered, minNumRegion = minCpGs, maxGap = 1000) # Includes chrX, chrY
write.table(background, file = "bsseq_background_Discovery50_EARLI.csv", sep = ",", quote = FALSE, row.names = FALSE)

# Meth ~ Diagnosis DMRs and Plots (Done)
bs.filtered <- chrSelectBSseq(bs.filtered, seqnames = c(paste("chr", 1:22, sep = ""), "chrM")) # Remove chrX, chrY
regions <- dmrseq(bs = bs.filtered, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms,
                  testCovariate = testCovariate, adjustCovariate = NULL, matchCovariate = NULL)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "CandidateRegions_DxNoXY_Discovery50_EARLI.csv")
gr2csv(sigRegions, "DMRs_DxNoXY_Discovery50_EARLI.csv")

# Smoothed Methylation and Plots (Done)
bs.filtered.bsseq <- BSmooth(bs.filtered, BPPARAM = MulticoreParam(workers = cores, progressbar = TRUE))
pData <- pData(bs.filtered.bsseq)
pData$col <- NULL
pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "#3366CC"
pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <-  "#FF3366"
pData$label <- NULL
pData$label[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "TD"
pData$label[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <-  "ASD"
pData(bs.filtered.bsseq) <- pData
saveRDS(bs.filtered.bsseq, "Filtered_Smoothed_BSseq_Discovery50_EARLI_noXY.rds")

annoTrack <- readRDS("hg38_annoTrack.rds")

smoothed <- getSmooth(bsseq = bs.filtered.bsseq, regions = sigRegions, 
                      out = "DMR_smoothed_methylation_DxNoXY_Discovery50_EARLI.txt")

pdf("DMRs_DxNoXY_Discovery50_EARLI_noPoints.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate,
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2, addRegions = sigRegions, 
          annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE, horizLegend = TRUE, addPoints = FALSE)
dev.off()

pdf("DMRs_DxNoXY_Discovery50_EARLI.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate,
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2, addRegions = sigRegions, 
          annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE, horizLegend = TRUE)
dev.off()

# Get raw DMR methylation (Done)
bs.filtered <- readRDS("Filtered_BSseq_Discovery50_EARLI.rds")
sigRegions <- read.csv("DMRs_DxNoXY_Discovery50_EARLI.csv", header = TRUE, stringsAsFactors = FALSE)
raw <- as.data.frame(getMeth(BSseq = bs.filtered, regions = data.frame2GRanges(sigRegions), type = "raw", what = "perRegion"))
raw <- cbind(sigRegions, raw)
write.table(raw, "DMR_raw_methylation_DxNoXY_Discovery50_EARLI.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# DMRs with Reduced Number of Males (Rerun without JLCM032B and JLCM050B) ####
# Reduce number of males (laptop)
set.seed(5)
meta <- read.xlsx("DMRs/Discovery/Diagnosis Males 50/sample_info_males.xlsx", colNames = TRUE) 
# Ctrl_TD Exp_ASD 
#      39      37 
metaRep <- read.xlsx("DMRs/Replication/Diagnosis Males 50/sample_info_males.xlsx", colNames = TRUE)
# Ctrl_TD Exp_ASD 
#      17      21
reducedTD <- sample(meta$Name[meta$Diagnosis == "Ctrl_TD"], size = table(metaRep$Diagnosis)["Ctrl_TD"]) %>% sort
reducedASD <- sample(meta$Name[meta$Diagnosis == "Exp_ASD"], size = table(metaRep$Diagnosis)["Exp_ASD"]) %>% sort
metaReduced <- subset(meta, Name %in% reducedTD | Name %in% reducedASD)
write.xlsx(metaReduced, file = "DMRs/Discovery/Diagnosis Males 50/sample_info_reduced_males.xlsx") # Excluded JLCM032B and JLCM050B from this

# Load and Process Samples (Rerun without JLCM032B and JLCM050B, Running on Barbera 6/3)
bs.filtered <- processBismark(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                              meta = read.xlsx("sample_info_reduced_males.xlsx", colNames = TRUE) %>% 
                                      mutate_if(is.character, as.factor), groups = testCovariate, Cov = coverage, 
                              mc.cores = cores, per.Group = perGroup)
saveRDS(bs.filtered, "Dx_Males/Filtered_BSseq_Discovery50_reduced_males.rds")

# Background Regions (Rerun without JLCM032B and JLCM050B, Running on Barbera 6/3)
background <- getBackground(bs.filtered, minNumRegion = minCpGs)
write.csv(background, file = "Dx_Males/bsseq_background_Discovery50_reduced_males.csv", quote = FALSE, row.names = FALSE)

# Meth ~ Diagnosis DMRs (Rerun without JLCM032B and JLCM050B, Running on Barbera 6/3)
regions <- dmrseq(bs = bs.filtered, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms, testCovariate = testCovariate)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "Dx_Males/CandidateRegions_Dx_Discovery50_reduced_males.csv")
gr2csv(sigRegions, "Dx_Males/DMRs_Dx_Discovery50_reduced_males.csv")

# DMR Comparison (Rerun without JLCM032B and JLCM050B) ####
library(ChIPpeakAnno)
source("R Scripts/DMR Analysis Functions.R")
samples <- read.xlsx("DMRs/Discovery/Diagnosis 50/sample_info.xlsx", colNames = TRUE) # JLCM032B and JLCM050B excluded
table(samples$Diagnosis, samples$Sex)
#          F  M
# Ctrl_TD 17 39
# Exp_ASD 15 35

# Load DMRs
regions <- list(All = loadRegions("DMRs/Discovery/Diagnosis 50/DMRs_DxNoXY_Discovery50.csv", 
                                  chroms = c(paste("chr", 1:22, sep = ""), "chrM"), DMRid = TRUE),
                Males = loadRegions("DMRs/Discovery/Diagnosis Males 50/DMRs_Dx_Discovery50_males.csv", 
                                    chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), DMRid = TRUE),
                Females = loadRegions("DMRs/Discovery/Diagnosis Females 50/DMRs_Dx_Discovery50_females.csv",
                                      chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), DMRid = TRUE))

# Load Background, all background is without coordinate shifting              
background <- list(All = loadRegions("DMRs/Discovery/Diagnosis 50/bsseq_background_Discovery50.csv", 
                                  chroms = c(paste("chr", 1:22, sep = ""), "chrM")),
                Males = loadRegions("DMRs/Discovery/Diagnosis Males 50/bsseq_background_Discovery50_males.csv", 
                                    chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")),
                Females = loadRegions("DMRs/Discovery/Diagnosis Females 50/bsseq_background_Discovery50_females.csv",
                                      chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM")))

# Region Stats
regionStats <- getRegionStats(regions, background = background, n = c(106, 74, 32))
write.table(regionStats, "Tables/DMR Region Stats All Chr Dx Discovery 50.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Write BED files
mapply(function(x, y){writeBED(regions = x, file = y)}, x = regions, 
       y = c("UCSC Tracks/Discovery Diagnosis DMRs.bed", "UCSC Tracks/Discovery Diagnosis Males DMRs.bed",
             "UCSC Tracks/Discovery Diagnosis Females DMRs.bed"))

# Region Overlap Venn Diagrams
GR_regions_hyper <- sapply(regions, makeGRange, direction = "hyper")
GR_regions_hypo <- sapply(regions, makeGRange, direction = "hypo")

DMRoverlapVenn(GR_regions_hyper, NameOfPeaks = c("All", "Males", "Females"), file = "Figures/Hyper DMR Overlap Dx Discovery 50.pdf",
               cat.pos = c(345, 15, 0), cat.dist = c(0.05, 0.05, 0.03), fill = c("lightblue", "lightpink", "lightgreen"))
DMRoverlapVenn(GR_regions_hypo, NameOfPeaks = c("All", "Males", "Females"), file = "Figures/Hypo DMR Overlap Dx Discovery 50.pdf",
               cat.pos = c(345, 15, 0), cat.dist = c(0.05, 0.05, 0.03), fill = c("lightblue", "lightpink", "lightgreen"))

# Male Female Overlap
DMRoverlapVenn(GR_regions_hyper[c("Males", "Females")], NameOfPeaks = c("Males", "Females"), 
               file = "Figures/Hyper DMR Overlap Dx Discovery 50 Males vs Females.pdf", rotation.degree = 180,
               cat.pos = c(180, 180), cat.dist = c(0.02, 0.02), fill = c("lightblue", "lightpink"))
DMRoverlapVenn(GR_regions_hypo[c("Males", "Females")], NameOfPeaks = c("Males", "Females"), 
               file = "Figures/Hypo DMR Overlap Dx Discovery 50 Males vs Females.pdf", rotation.degree = 180,
               cat.pos = c(180, 180), cat.dist = c(0.02, 0.02), fill = c("lightblue", "lightpink"))

HyperMF <- GenomicRanges::intersect(x = GR_regions_hyper$Males, y = GR_regions_hyper$Females) %>% as.data.frame
HypoMF <- GenomicRanges::intersect(x = GR_regions_hypo$Males, y = GR_regions_hypo$Females) %>% as.data.frame
MF <- rbind(HyperMF, HypoMF)

HyperMalesDMRsInF <- subset(regions$Males[regions$Males$percentDifference > 0,], GR_regions_hyper$Males %over% GR_regions_hyper$Females)
HypoMalesDMRsInF <- subset(regions$Males[regions$Males$percentDifference < 0,], GR_regions_hypo$Males %over% GR_regions_hypo$Females)
MalesDMRsInF <- rbind(HyperMalesDMRsInF, HypoMalesDMRsInF)
write.table(MalesDMRsInF, "Tables/Dx Discovery Males DMRs Overlapping with Females.txt", sep = "\t", quote = FALSE, row.names = FALSE)

HyperFemalesDMRsInM <- subset(regions$Females[regions$Females$percentDifference > 0,], 
                              GR_regions_hyper$Females %over% GR_regions_hyper$Males)
HypoFemalesDMRsInM <- subset(regions$Females[regions$Females$percentDifference < 0,], 
                             GR_regions_hypo$Females %over% GR_regions_hypo$Males)
FemalesDMRsInM <- rbind(HyperFemalesDMRsInM, HypoFemalesDMRsInM)
write.table(FemalesDMRsInM, "Tables/Dx Discovery Females DMRs Overlapping with Males.txt", sep = "\t", quote = FALSE, row.names = FALSE)

rm(HyperMF, HypoMF, HyperMalesDMRsInF, HypoMalesDMRsInF, HyperFemalesDMRsInM, HypoFemalesDMRsInM, venn)









