# Diagnosis and Sex DMRs Discovery ####
# Autism Cord Blood Methylation
# Charles Mordaunt
# 1/22/19

# cmordaunt@epigenerate:/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Discovery
# THREADS=${SLURM_NTASKS}
# MEM=$(expr ${SLURM_MEM_PER_CPU} / 1024)
# echo "Allocated threads: " $THREADS
# echo "Allocated memory: " $MEM

# Load Packages ####
.libPaths("/share/lasallelab/programs/DMRichR/R_3.5")
sapply(c("tidyverse", "openxlsx", "bsseq", "dmrseq", "DMRichR"), require, character.only = TRUE)

# Global Variables ####
genome <- as.character("hg38")
coverage <- as.numeric(1)
perGroup <- (as.numeric(50)/100)
minCpGs <- as.numeric(3)
maxPerms <- as.numeric(10)
testCovariate <- as.character("Diagnosis")
#adjustCovariate <- NULL
#matchCovariate <- NULL
cores <- 10
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

# DMRs with All Samples ####
# Load and Process Samples (Done)
name <- gsub( "_.*$","", list.files(path = getwd(), pattern = "*.txt.gz"))
bs.filtered <- processBismark(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                              meta = read.xlsx("sample_info.xlsx", colNames = TRUE) %>% mutate_if(is.character,as.factor),
                              groups = testCovariate, Cov = coverage, mc.cores = cores, per.Group = perGroup)
saveRDS(bs.filtered, "Filtered_BSseq_Discovery50.rds")

# Background Regions (Done)
background <- getBackground(bs.filtered, minNumRegion = minCpGs, maxGap = 1000)
write.table(background, file = "bsseq_background_Discovery50.csv", sep = ",", quote = FALSE, row.names = FALSE)

# Smoothed Methylation (Done)
bs.filtered <- readRDS("Filtered_BSseq_Discovery50.rds")
cores <- 5
register(MulticoreParam(1))
testCovariate <- "Diagnosis"
bs.filtered.bsseq <- BSmooth(bs.filtered, BPPARAM = MulticoreParam(workers = cores, progressbar = TRUE))

pData <- pData(bs.filtered.bsseq)
pData$col <- NULL
pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "#3366CC"
pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <-  "#FF3366"
pData$label <- NULL
pData$label[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "TD"
pData$label[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <-  "ASD"
pData(bs.filtered.bsseq) <- pData

saveRDS(bs.filtered.bsseq, "Filtered_Smoothed_BSseq_Discovery50.rds")

# Meth ~ Diagnosis DMRs and Plots ####
# (Done)
bs.filtered <- readRDS("Filtered_BSseq_Discovery50.rds")
regionsDx <- dmrseq(bs = bs.filtered, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms,
                    testCovariate = testCovariate, adjustCovariate = NULL, matchCovariate = NULL)
regionsDx$percentDifference <- round(regionsDx$beta/pi * 100)
sigRegionsDx <- regionsDx[regionsDx$pval < 0.05,]
gr2csv(regionsDx, "CandidateRegions_Dx_Discovery50.csv")
gr2csv(sigRegionsDx, "DMRs_Dx_Discovery50.csv")

bs.filtered.bsseq <- readRDS("Filtered_Smoothed_BSseq_Discovery50.rds")
annoTrack <- readRDS("hg38_annoTrack.rds")
# Pasted in plotFunctions.R
smoothed <- getSmooth(bsseq = bs.filtered.bsseq, regions = sigRegionsDx, 
                      out = "DMR_smoothed_methylation_Dx_Discovery50.txt")
sigRegionsDx <- read.csv("DMRs_Dx_Discovery50.csv", header = TRUE, stringsAsFactors = FALSE)
sigRegionsDx <- data.frame2GRanges(sigRegionsDx)
pdf("DMRs_Dx_Discovery50.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegionsDx, testCovariate = testCovariate, 
          extend = 5000 - (end(sigRegionsDx) - start(sigRegionsDx) + 1)/2,
          addRegions = sigRegionsDx, annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE,
          horizLegend = TRUE)
dev.off()

pdf("DMRs_Dx_Discovery50_noPoints.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegionsDx, testCovariate = testCovariate, 
          extend = 5000 - (end(sigRegionsDx) - start(sigRegionsDx) + 1)/2,
          addRegions = sigRegionsDx, annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE,
          horizLegend = TRUE, addPoints = FALSE)
dev.off()

# Get raw DMR methylation (Done)
bs.filtered <- readRDS("Filtered_BSseq_Discovery50.rds")
sigRegionsDx <- read.csv("DMRs_Dx_Discovery50.csv", header = TRUE, stringsAsFactors = FALSE)
raw <- as.data.frame(getMeth(BSseq = bs.filtered, regions = data.frame2GRanges(sigRegionsDx), type = "raw", what = "perRegion"))
raw <- cbind(sigRegionsDx, raw)
write.table(raw, "DMR_raw_methylation_Dx_Discovery50.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Meth ~ Diagnosis + AdjSex DMRs and Plots ####
# (Done)
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

# Meth ~ Diagnosis, exclude chrX and Y, DMRs and Plots ####
# (Done)
bs.filtered <- readRDS("Filtered_BSseq_Discovery50.rds")
bs.filtered <- chrSelectBSseq(bs.filtered, seqnames = c(paste("chr", 1:22, sep = ""), "chrM")) # Remove chrX, chrY
regions <- dmrseq(bs = bs.filtered, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms,
                  testCovariate = testCovariate, adjustCovariate = NULL, matchCovariate = NULL)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "CandidateRegions_DxNoXY_Discovery50.csv")
gr2csv(sigRegions, "DMRs_DxNoXY_Discovery50.csv")

bs.filtered.bsseq <- readRDS("Filtered_Smoothed_BSseq_Discovery50.rds")
annoTrack <- readRDS("hg38_annoTrack.rds")
smoothed <- getSmooth(bsseq = bs.filtered.bsseq, regions = sigRegions, 
                      out = "DMR_smoothed_methylation_DxNoXY_Discovery50.txt")

pdf("DMRs_DxNoXY_Discovery50.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate, 
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2,
          addRegions = sigRegions, annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE,
          horizLegend = TRUE)
dev.off()

pdf("DMRs_DxNoXY_Discovery50_noPoints.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate, 
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2,
          addRegions = sigRegions, annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE,
          horizLegend = TRUE, addPoints = FALSE)
dev.off()

# Get raw DMR methylation (Done)
bs.filtered <- readRDS("Filtered_BSseq_Discovery50.rds")
sigRegions <- read.csv("DMRs_DxNoXY_Discovery50.csv", header = TRUE, stringsAsFactors = FALSE)
raw <- as.data.frame(getMeth(BSseq = bs.filtered, regions = data.frame2GRanges(sigRegions), type = "raw", what = "perRegion"))
raw <- cbind(sigRegions, raw)
write.table(raw, "DMR_raw_methylation_DxNoXY_Discovery50.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# DMRs with Males Only ####
# New R session
# Load and Process Samples (Done)
name <- gsub( "_.*$","", list.files(path = getwd(), pattern = "*.txt.gz"))
bs.filtered <- processBismark(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                              meta = read.xlsx("sample_info_males.xlsx", colNames = TRUE) %>% mutate_if(is.character,as.factor),
                              groups = testCovariate, Cov = coverage, mc.cores = cores, per.Group = perGroup)
saveRDS(bs.filtered, "Filtered_BSseq_Discovery50_males.rds")

# Background Regions (Done)
background <- getBackground(bs.filtered, minNumRegion = minCpGs, maxGap = 1000)
write.table(background, file = "bsseq_background_Discovery50_males.csv", sep = ",", quote = FALSE, row.names = FALSE)

# Meth ~ Diagnosis DMRs and Plots (Done)
regions <- dmrseq(bs = bs.filtered, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms,
                  testCovariate = testCovariate, adjustCovariate = NULL, matchCovariate = NULL)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "CandidateRegions_Dx_Discovery50_males.csv")
gr2csv(sigRegions, "DMRs_Dx_Discovery50_males.csv")

# Smoothed Methylation and Plots (Done)
bs.filtered <- readRDS("Filtered_BSseq_Discovery50_males.rds")
sigRegions <- read.csv("DMRs_Dx_Discovery50_males.csv", header = TRUE, stringsAsFactors = FALSE)
sigRegions <- data.frame2GRanges(sigRegions)

# Pasted in plotFunctions.R

bs.filtered.bsseq <- BSmooth(bs.filtered, BPPARAM = MulticoreParam(workers = cores, progressbar = TRUE))
pData <- pData(bs.filtered.bsseq)
pData$col <- NULL
pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "#3366CC"
pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <-  "#FF3366"
pData$label <- NULL
pData$label[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "TD"
pData$label[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <-  "ASD"
pData(bs.filtered.bsseq) <- pData
saveRDS(bs.filtered.bsseq, "Filtered_Smoothed_BSseq_Discovery50_males.rds")

annoTrack <- readRDS("hg38_annoTrack.rds")

smoothed <- getSmooth(bsseq = bs.filtered.bsseq, regions = sigRegions, 
                      out = "DMR_smoothed_methylation_Dx_Discovery50_males.txt")

pdf("DMRs_Dx_Discovery50_males_noPoints.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate,
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2, addRegions = sigRegions, 
          annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE, horizLegend = TRUE, addPoints = FALSE)
dev.off()

pdf("DMRs_Dx_Discovery50_males.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate,
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2, addRegions = sigRegions, 
          annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE, horizLegend = TRUE)
dev.off()

# Get raw DMR methylation (Done)
bs.filtered <- readRDS("Filtered_BSseq_Discovery50_males.rds")
sigRegions <- read.csv("DMRs_Dx_Discovery50_males.csv", header = TRUE, stringsAsFactors = FALSE)
raw <- as.data.frame(getMeth(BSseq = bs.filtered, regions = data.frame2GRanges(sigRegions), type = "raw", what = "perRegion"))
raw <- cbind(sigRegions, raw)
write.table(raw, "DMR_raw_methylation_Dx_Discovery50_males.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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

# DMR Comparison ####
samples <- read.xlsx("DMRs/Discovery/Diagnosis 50/sample_info.xlsx", colNames = TRUE)
table(samples$Diagnosis, samples$Sex)
#          F  M
# Ctrl_TD 17 39
# Exp_ASD 15 37

# Load DMRs
DxAllDMRs <- read.csv("DMRs/Discovery/Diagnosis 50/DMRs_Dx_Discovery50.csv", header = TRUE, 
                      stringsAsFactors = FALSE)
DxAllCand <- read.csv("DMRs/Discovery/Diagnosis 50/CandidateRegions_Dx_Discovery50.csv", header = TRUE, 
                      stringsAsFactors = FALSE)
DxAllBack <- read.csv("DMRs/Discovery/Diagnosis 50/bsseq_background_Discovery50.csv", header = TRUE, 
                      stringsAsFactors = FALSE)
DxSexAllDMRs <- read.csv("DMRs/Discovery/Diagnosis and Sex 50/DMRs_DxAdjSex_Discovery50.csv", header = TRUE, 
                      stringsAsFactors = FALSE)
DxSexAllCand <- read.csv("DMRs/Discovery/Diagnosis and Sex 50/CandidateRegions_DxAdjSex_Discovery50.csv", header = TRUE, 
                      stringsAsFactors = FALSE)
DxSexAllBack <- read.csv("DMRs/Discovery/Diagnosis 50/bsseq_background_Discovery50.csv", header = TRUE, 
                      stringsAsFactors = FALSE) # Same as DxAllBack
DxMalesDMRs <- read.csv("DMRs/Discovery/Diagnosis Males 50/DMRs_Dx_Discovery50_males.csv", header = TRUE, 
                        stringsAsFactors = FALSE)
DxMalesCand <- read.csv("DMRs/Discovery/Diagnosis Males 50/CandidateRegions_Dx_Discovery50_males.csv", header = TRUE, 
                        stringsAsFactors = FALSE)
DxMalesBack <- read.csv("DMRs/Discovery/Diagnosis Males 50/bsseq_background_Discovery50_males.csv", header = TRUE, 
                        stringsAsFactors = FALSE)
DxFemalesDMRs <- read.csv("DMRs/Discovery/Diagnosis Females 50/DMRs_Dx_Discovery50_females.csv", header = TRUE, 
                          stringsAsFactors = FALSE)
DxFemalesCand <- read.csv("DMRs/Discovery/Diagnosis Females 50/CandidateRegions_Dx_Discovery50_females.csv", header = TRUE, 
                          stringsAsFactors = FALSE)
DxFemalesBack <- read.csv("DMRs/Discovery/Diagnosis Females 50/bsseq_background_Discovery50_females.csv", header = TRUE, 
                          stringsAsFactors = FALSE)

# Region Stats
regions <- list(DxAllDMRs, DxAllCand, DxAllBack, DxSexAllDMRs, DxSexAllCand, DxSexAllBack, DxMalesDMRs, DxMalesCand, DxMalesBack, 
                DxFemalesDMRs, DxFemalesCand, DxFemalesBack)
names(regions) <- c("DxAllDMRs", "DxAllCand", "DxAllBack", "DxSexAllDMRs", "DxSexAllCand", "DxSexAllBack", 
                    "DxMalesDMRs", "DxMalesCand", "DxMalesBack", "DxFemalesDMRs", "DxFemalesCand", "DxFemalesBack")
regions <- sapply(regions, function(x){
        colnames(x)[1] <- "chr"
        return(x)
        })

getRegionStats <- function(regionData){
        regionStats <- sapply(regionData, function(x){
                c("Number" = nrow(x), "Width" = sum(x$width), "CpGs" = sum(x[,colnames(x) %in% c("n", "L")]))
                })
        regionStats <- regionStats %>% t %>% as.data.frame
        regionStats$Regions <- row.names(regionStats)
        row.names(regionStats) <- 1:nrow(regionStats)
        regionStats$TD <- rep(c(56, 56, 39, 17), each = 3)
        regionStats$ASD <- rep(c(52, 52, 37, 15), each = 3)
        regionStats$All <- regionStats$TD + regionStats$ASD
        regionStats$NumberPerBack <- round(regionStats$Number * 100 / rep(regionStats$Number[grepl("Back", regionStats$Regions)], each = 3), 3)
        regionStats$WidthPerBack <- round(regionStats$Width * 100 / rep(regionStats$Width[grepl("Back", regionStats$Regions)], each = 3), 3)
        regionStats$CpGsPerBack <- round(regionStats$CpGs * 100 / rep(regionStats$CpGs[grepl("Back", regionStats$Regions)], each = 3), 3)
        regionStats <- regionStats[,c("Regions", "TD", "ASD", "All", "Number", "NumberPerBack", "Width", "WidthPerBack", "CpGs", "CpGsPerBack")]
        return(regionStats)
}

regionStats <- getRegionStats(regions)
write.table(regionStats, "Tables/DMR Region Stats All Chr Dx Discovery 50.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Region Stats chrX
regionsX <- sapply(regions, function(x){subset(x, chr == "chrX")})
regionStatsX <- getRegionStats(regionsX)
write.table(regionStatsX, "Tables/DMR Region Stats ChrX Dx Discovery 50.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Region Overlap Venn Diagrams
library(ChIPpeakAnno)
GR_regions <- sapply(regions, function(x) {GRanges(seqnames = x$chr, ranges = IRanges(start = x$start, end = x$end))})

pdf(file="Figures/DMR Overlap Dx Discovery 50.pdf", width=10, height=8, onefile = FALSE)
venn <- suppressMessages(makeVennDiagram(Peaks = list(GR_regions$DxAllDMRs, GR_regions$DxSexAllDMRs, 
                                                      GR_regions$DxMalesDMRs, GR_regions$DxFemalesDMRs), 
                                         NameOfPeaks = c("All_DMRs", "All_DMRs_Sex", "Males_DMRs", "Females_DMRs"), 
                                         maxgap = -1, minoverlap = 1, by = "region", connectedPeaks = "min", 
                                         rotation.degree = 0, margin = 0.02, cat.cex = 2, cex = 2.5, 
                                         fill = c("lightblue", "lavender", "lightpink", "lightgreen"), 
                                         cat.pos = c(355, 10, 0, 0), cat.dist = c(0.2, 0.2, 0.1, 0.08), 
                                         fontfamily = "sans", cat.fontfamily = "sans", ext.dist = -0.4,
                                         ext.length = 0.85))
dev.off()

pdf(file="Figures/Hyper DMR Overlap Dx Discovery 50.pdf", width=10, height=8, onefile = FALSE)
venn <- suppressMessages(makeVennDiagram(Peaks = list(GR_regions$DxAllDMRs[DxAllDMRs$percentDifference > 0], 
                                                      GR_regions$DxSexAllDMRs[DxSexAllDMRs$percentDifference > 0], 
                                                      GR_regions$DxMalesDMRs[DxMalesDMRs$percentDifference > 0], 
                                                      GR_regions$DxFemalesDMRs[DxFemalesDMRs$percentDifference > 0]), 
                                         NameOfPeaks = c("All_DMRs", "All_DMRs_Sex", "Males_DMRs", "Females_DMRs"), 
                                         maxgap = -1, minoverlap = 1, by = "region", connectedPeaks = "min", 
                                         rotation.degree = 0, margin = 0.02, cat.cex = 2, cex = 2.5, 
                                         fill = c("lightblue", "lavender", "lightpink", "lightgreen"), 
                                         cat.pos = c(350, 10, 0, 0), cat.dist = c(0.2, 0.2, 0.1, 0.08), 
                                         fontfamily = "sans", cat.fontfamily = "sans", ext.dist = -0.4, 
                                         ext.length = 0.85))
dev.off()

pdf(file="Figures/Hypo DMR Overlap Dx Discovery 50.pdf", width=10, height=8, onefile = FALSE)
venn <- suppressMessages(makeVennDiagram(Peaks = list(GR_regions$DxAllDMRs[DxAllDMRs$percentDifference < 0],
                                                      GR_regions$DxSexAllDMRs[DxSexAllDMRs$percentDifference < 0],
                                                      GR_regions$DxMalesDMRs[DxMalesDMRs$percentDifference < 0], 
                                                      GR_regions$DxFemalesDMRs[DxFemalesDMRs$percentDifference < 0]), 
                                         NameOfPeaks = c("All_DMRs", "All_DMRs_Sex", "Males_DMRs", "Females_DMRs"), 
                                         maxgap = -1, minoverlap = 1, by = "region", connectedPeaks = "min", 
                                         rotation.degree = 0, margin = 0.02, cat.cex = 2, cex = 2.5, 
                                         fill = c("lightblue", "lavender", "lightpink", "lightgreen"), 
                                         cat.pos = c(350, 10, 0, 0), cat.dist = c(0.2, 0.2, 0.1, 0.08), 
                                         fontfamily = "sans", cat.fontfamily = "sans", ext.dist = -0.4, 
                                         ext.length = 0.85))
dev.off()

# Male Female Overlap
HyperMF <- intersect(x = GR_regions$DxMalesDMRs[DxMalesDMRs$percentDifference > 0],
                     y = GR_regions$DxFemalesDMRs[DxFemalesDMRs$percentDifference > 0]) %>% as.data.frame
HypoMF <- intersect(x = GR_regions$DxMalesDMRs[DxMalesDMRs$percentDifference < 0],
                    y = GR_regions$DxFemalesDMRs[DxFemalesDMRs$percentDifference < 0]) %>% as.data.frame
MF <- rbind(HyperMF, HypoMF)

HyperMalesDMRsInF <- subset(DxMalesDMRs[DxMalesDMRs$percentDifference > 0,], GR_regions$DxMalesDMRs[DxMalesDMRs$percentDifference > 0] %over%
                                    GR_regions$DxFemalesDMRs[DxFemalesDMRs$percentDifference > 0])
HypoMalesDMRsInF <- subset(DxMalesDMRs[DxMalesDMRs$percentDifference < 0,], GR_regions$DxMalesDMRs[DxMalesDMRs$percentDifference < 0] %over%
                                    GR_regions$DxFemalesDMRs[DxFemalesDMRs$percentDifference < 0])
MalesDMRsInF <- rbind(HyperMalesDMRsInF, HypoMalesDMRsInF)
write.table(MalesDMRsInF, "Tables/Dx Discovery Males DMRs Overlapping with Females.txt", sep = "\t", quote = FALSE, row.names = FALSE)

HyperFemalesDMRsInM <- subset(DxFemalesDMRs[DxFemalesDMRs$percentDifference > 0,], GR_regions$DxFemalesDMRs[DxFemalesDMRs$percentDifference > 0] %over%
                                      GR_regions$DxMalesDMRs[DxMalesDMRs$percentDifference > 0])
HypoFemalesDMRsInM <- subset(DxFemalesDMRs[DxFemalesDMRs$percentDifference < 0,], GR_regions$DxFemalesDMRs[DxFemalesDMRs$percentDifference < 0] %over%
                                      GR_regions$DxMalesDMRs[DxMalesDMRs$percentDifference < 0])
FemalesDMRsInM <- rbind(HyperFemalesDMRsInM, HypoFemalesDMRsInM)
write.table(FemalesDMRsInM, "Tables/Dx Discovery Females DMRs Overlapping with Males.txt", sep = "\t", quote = FALSE, row.names = FALSE)

rm(HyperMF, HypoMF, HyperMalesDMRsInF, HypoMalesDMRsInF, HyperFemalesDMRsInM, HypoFemalesDMRsInM, venn)









