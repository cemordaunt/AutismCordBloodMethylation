# Diagnosis and Sex DMRs Replication ####
# Autism Cord Blood Methylation
# Charles Mordaunt
# 2/5/19

# cmordaunt@epigenerate:/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Replication

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
set.seed(5)
register(MulticoreParam(1))

# DMRs with All Samples ####
# Load and Process Samples (Done)
name <- gsub( "_.*$","", list.files(path = getwd(), pattern = "*.txt.gz"))
bs.filtered <- processBismark(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                              meta = read.xlsx("sample_info.xlsx", colNames = TRUE) %>% mutate_if(is.character,as.factor),
                              groups = testCovariate, Cov = coverage, mc.cores = cores, per.Group = perGroup)
saveRDS(bs.filtered, "BSseq_05Group.rds")

# Background Regions (Done)
bs.filtered <- readRDS("BSseq_05Group.rds")
background <- getBackground(bs.filtered, minNumRegion = minCpGs, maxGap = 1000)
write.table(background, file = "bsseq_background_Replication50.csv", sep = ",", quote = FALSE, row.names = FALSE)

# Smoothed Methylation (Done)
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

saveRDS(bs.filtered.bsseq, "Filtered_Smoothed_BSseq_Replication50.rds")

# Meth ~ Diagnosis DMRs and Plots ####
# DMRs (Done)
regions <- dmrseq(bs = bs.filtered, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms,
                    testCovariate = testCovariate, adjustCovariate = NULL, matchCovariate = NULL)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "CandidateRegions_Dx_Replication50.csv")
gr2csv(sigRegions, "DMRs_Dx_Replication50.csv")

# Smoothed Methylation (Done)
annoTrack <- readRDS("hg38_annoTrack.rds")
smoothed <- getSmooth(bsseq = bs.filtered.bsseq, regions = sigRegions, 
                      out = "DMR_smoothed_methylation_Dx_Replication50.txt")

# Plots (Done)
pdf("DMRs_Dx_Replication50.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate, 
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2,
          addRegions = sigRegions, annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE,
          horizLegend = TRUE)
dev.off()

pdf("DMRs_Dx_Replication50_noPoints.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate, 
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2,
          addRegions = sigRegions, annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE,
          horizLegend = TRUE, addPoints = FALSE)
dev.off()

# Get raw DMR methylation (Done)
bs.filtered <- readRDS("BSseq_05Group.rds")
sigRegions <- read.csv("DMRs_Dx_Replication50.csv", header = TRUE, stringsAsFactors = FALSE)
raw <- as.data.frame(getMeth(BSseq = bs.filtered, regions = data.frame2GRanges(sigRegions), type = "raw", what = "perRegion"))
raw <- cbind(sigRegions, raw)
write.table(raw, "DMR_raw_methylation_Dx_Replication50.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Meth ~ Diagnosis + AdjSex DMRs and Plots ####
# DMRs (Done)
bs.filtered <- readRDS("BSseq_05Group.rds")
regions <- dmrseq(bs = bs.filtered, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms,
                  testCovariate = testCovariate, adjustCovariate = "Sex", matchCovariate = NULL)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "backgroundRegions50.csv")

# Smoothed Methylation (Done)
regions <- read.csv("backgroundRegions50.csv", header = TRUE, stringsAsFactors = FALSE)
sigRegions <- subset(regions, pval < 0.05)
write.csv(sigRegions, "DMRs_DxAdjSex_Replication50.csv", row.names = FALSE, quote = FALSE)

smoothed <- getSmooth(bsseq = bs.filtered.bsseq, regions = sigRegions, 
                      out = "DMR_smoothed_methylation_DxAdjSex_Replication50.txt")

# Plots (Done)
sigRegions <- data.frame2GRanges(sigRegions)
pdf("DMRs_DxAdjSex_Replication50.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate, 
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2,
          addRegions = sigRegions, annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE,
          horizLegend = TRUE)
dev.off()

pdf("DMRs_DxAdjSex_Replication50_noPoints.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate, 
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2,
          addRegions = sigRegions, annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE,
          horizLegend = TRUE, addPoints = FALSE)
dev.off()

# Get raw DMR methylation (Done)
bs.filtered <- readRDS("BSseq_05Group.rds")
sigRegions <- read.csv("DMRs_DxAdjSex_Replication50.csv", header = TRUE, stringsAsFactors = FALSE)
raw <- as.data.frame(getMeth(BSseq = bs.filtered, regions = data.frame2GRanges(sigRegions), type = "raw", what = "perRegion"))
raw <- cbind(sigRegions, raw)
write.table(raw, "DMR_raw_methylation_DxAdjSex_Replication50.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Meth ~ Diagnosis, exclude chrX and Y, DMRs and Plots ####
# (Done)
bs.filtered <- readRDS("BSseq_05Group.rds")
bs.filtered <- chrSelectBSseq(bs.filtered, seqnames = c(paste("chr", 1:22, sep = ""), "chrM")) # Remove chrX, chrY
regions <- dmrseq(bs = bs.filtered, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms,
                  testCovariate = testCovariate, adjustCovariate = NULL, matchCovariate = NULL)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "CandidateRegions_DxNoXY_Replication50.csv")
gr2csv(sigRegions, "DMRs_DxNoXY_Replication50.csv")

bs.filtered.bsseq <- readRDS("Filtered_Smoothed_BSseq_Replication50.rds")
annoTrack <- readRDS("hg38_annoTrack.rds")
smoothed <- getSmooth(bsseq = bs.filtered.bsseq, regions = sigRegions, 
                      out = "DMR_smoothed_methylation_DxNoXY_Replication50.txt")

pdf("DMRs_DxNoXY_Replication50.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate, 
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2,
          addRegions = sigRegions, annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE,
          horizLegend = TRUE)
dev.off()

pdf("DMRs_DxNoXY_Replication50_noPoints.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate, 
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2,
          addRegions = sigRegions, annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE,
          horizLegend = TRUE, addPoints = FALSE)
dev.off()

# Get raw DMR methylation (Done)
sigRegions <- read.csv("DMRs_DxNoXY_Replication50.csv", header = TRUE, stringsAsFactors = FALSE)
raw <- as.data.frame(getMeth(BSseq = bs.filtered, regions = data.frame2GRanges(sigRegions), type = "raw", what = "perRegion"))
raw <- cbind(sigRegions, raw)
write.table(raw, "DMR_raw_methylation_DxNoXY_Replication50.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# DMRs with Males Only ####
# Load and Process Samples (Done)
cores <- 5
name <- gsub( "_.*$","", list.files(path = getwd(), pattern = "*.txt.gz"))
bs.filtered <- processBismark(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                              meta = read.xlsx("sample_info_males.xlsx", colNames = TRUE) %>% mutate_if(is.character,as.factor),
                              groups = testCovariate, Cov = coverage, mc.cores = cores, per.Group = perGroup)
saveRDS(bs.filtered, "Filtered_BSseq_Replication50_males.rds")

# Background Regions (Done)
background <- getBackground(bs.filtered, minNumRegion = minCpGs, maxGap = 1000)
write.table(background, file = "bsseq_background_Replication50_males.csv", sep = ",", quote = FALSE, row.names = FALSE)

# Meth ~ Diagnosis DMRs and Plots (Done)
regions <- dmrseq(bs = bs.filtered, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms,
                  testCovariate = testCovariate, adjustCovariate = NULL, matchCovariate = NULL)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "CandidateRegions_Dx_Replication50_males.csv")
gr2csv(sigRegions, "DMRs_Dx_Replication50_males.csv")

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
saveRDS(bs.filtered.bsseq, "Filtered_Smoothed_BSseq_Replication50_males.rds")

smoothed <- getSmooth(bsseq = bs.filtered.bsseq, regions = sigRegions, 
                      out = "DMR_smoothed_methylation_Dx_Replication50_males.txt")

annoTrack <- readRDS("hg38_annoTrack.rds")
pdf("DMRs_Dx_Replication50_males_noPoints.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate,
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2, addRegions = sigRegions, 
          annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE, horizLegend = TRUE, addPoints = FALSE)
dev.off()

pdf("DMRs_Dx_Replication50_males.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate,
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2, addRegions = sigRegions, 
          annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE, horizLegend = TRUE)
dev.off()

# Get raw DMR methylation (Done)
bs.filtered <- readRDS("Filtered_BSseq_Replication50_males.rds")
sigRegions <- read.csv("DMRs_Dx_Replication50_males.csv", header = TRUE, stringsAsFactors = FALSE)
raw <- as.data.frame(getMeth(BSseq = bs.filtered, regions = data.frame2GRanges(sigRegions), type = "raw", what = "perRegion"))
raw <- cbind(sigRegions, raw)
write.table(raw, "DMR_raw_methylation_Dx_Replication50_males.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# DMRs with Females Only ####
# Load and Process Samples (Done)
cores <- 10
name <- gsub( "_.*$","", list.files(path = getwd(), pattern = "*.txt.gz"))
bs.filtered <- processBismark(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                              meta = read.xlsx("sample_info_females.xlsx", colNames = TRUE) %>% mutate_if(is.character,as.factor),
                              groups = testCovariate, Cov = coverage, mc.cores = cores, per.Group = perGroup)
saveRDS(bs.filtered, "Filtered_BSseq_Replication50_females.rds")

# Background Regions (Done)
background <- getBackground(bs.filtered, minNumRegion = minCpGs, maxGap = 1000)
write.table(background, file = "bsseq_background_Replication50_females.csv", sep = ",", quote = FALSE, row.names = FALSE)

# Meth ~ Diagnosis DMRs (Done)
cores <- 10
bs.filtered <- readRDS("Filtered_BSseq_Replication50_females.rds")
regions <- dmrseq(bs = bs.filtered, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms,
                  testCovariate = testCovariate, adjustCovariate = NULL, matchCovariate = NULL)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "CandidateRegions_Dx_Replication50_females.csv")
gr2csv(sigRegions, "DMRs_Dx_Replication50_females.csv")

# Smoothed Methylation (Done)
bs.filtered.bsseq <- BSmooth(bs.filtered, BPPARAM = MulticoreParam(workers = cores, progressbar = TRUE))
pData <- pData(bs.filtered.bsseq)
pData$col <- NULL
pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "#3366CC"
pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <-  "#FF3366"
pData$label <- NULL
pData$label[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "TD"
pData$label[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <-  "ASD"
pData(bs.filtered.bsseq) <- pData
saveRDS(bs.filtered.bsseq, "Filtered_Smoothed_BSseq_Replication50_females.rds")

smoothed <- getSmooth(bsseq = bs.filtered.bsseq, regions = sigRegions, 
                      out = "DMR_smoothed_methylation_Dx_Replication50_females.txt")

# Plots (Done)
annoTrack <- readRDS("hg38_annoTrack.rds")
pdf("DMRs_Dx_Replication50_females_noPoints.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate, 
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2, addRegions = sigRegions, 
          annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE, horizLegend = TRUE, 
          addPoints = FALSE)
dev.off()

pdf("DMRs_Dx_Replication50_females.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate, 
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2, addRegions = sigRegions, 
          annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE, horizLegend = TRUE, 
          addPoints = TRUE)
dev.off()

# Get raw DMR methylation (Done)
bs.filtered <- readRDS("Filtered_BSseq_Replication50_females.rds")
sigRegions <- read.csv("DMRs_Dx_Replication50_females.csv", header = TRUE, stringsAsFactors = FALSE)
raw <- as.data.frame(getMeth(BSseq = bs.filtered, regions = data.frame2GRanges(sigRegions), type = "raw", what = "perRegion"))
raw <- cbind(sigRegions, raw)
write.table(raw, "DMR_raw_methylation_Dx_Replication50_females.txt", sep = "\t", quote = FALSE, row.names = FALSE)

