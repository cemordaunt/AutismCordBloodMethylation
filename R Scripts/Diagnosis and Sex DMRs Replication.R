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

# DMRs with Females Only, 100 per group ####
# Load Data (Done)
perGroup <- (as.numeric(100)/100)
cores <- 10
name <- gsub( "_.*$","", list.files(path = getwd(), pattern = "*.txt.gz"))
bs.filtered <- processBismark(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                              meta = read.xlsx("sample_info_females.xlsx", colNames = TRUE) %>% mutate_if(is.character,as.factor),
                              groups = testCovariate, Cov = coverage, mc.cores = cores, per.Group = perGroup)
saveRDS(bs.filtered, "Filtered_BSseq_Replication100_females.rds")

# Background Regions (Done)
background <- getBackground(bs.filtered, minNumRegion = minCpGs, maxGap = 1000)
write.table(background, file = "bsseq_background_Replication100_females.csv", sep = ",", quote = FALSE, row.names = FALSE)

# Meth ~ Diagnosis DMRs (Done)
cores <- 10
regions <- dmrseq(bs = bs.filtered, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms,
                  testCovariate = testCovariate, adjustCovariate = NULL, matchCovariate = NULL)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "CandidateRegions_Dx_Replication100_females.csv")
gr2csv(sigRegions, "DMRs_Dx_Replication100_females.csv")

# Get raw DMR methylation (Done)
raw <- as.data.frame(getMeth(BSseq = bs.filtered, regions = sigRegions, type = "raw", what = "perRegion"))
raw <- cbind(sigRegions, raw)
write.table(raw, "DMR_raw_methylation_Dx_Replication100_females.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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
saveRDS(bs.filtered.bsseq, "Filtered_Smoothed_BSseq_Replication100_females.rds")

smoothed <- getSmooth(bsseq = bs.filtered.bsseq, regions = sigRegions, 
                      out = "DMR_smoothed_methylation_Dx_Replication100_females.txt")

# Plots (Done)
annoTrack <- readRDS("hg38_annoTrack.rds")
pdf("DMRs_Dx_Replication100_females_noPoints.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate, 
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2, addRegions = sigRegions, 
          annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE, horizLegend = TRUE, 
          addPoints = FALSE)
dev.off()

pdf("DMRs_Dx_Replication100_females.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate, 
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2, addRegions = sigRegions, 
          annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE, horizLegend = TRUE, 
          addPoints = TRUE)
dev.off()

# DMR Comparison ####
samples <- read.xlsx("DMRs/Replication/Diagnosis 50/sample_info.xlsx", colNames = TRUE)
table(samples$Diagnosis, samples$Sex)
#          F  M
# Ctrl_TD  3 17
# Exp_ASD  5 21

# Load DMRs ####
DxAllDMRs <- read.csv("DMRs/Replication/Diagnosis 50/DMRs_DxNoXY_Replication50.csv", header = TRUE, 
                      stringsAsFactors = FALSE) %>% subset(!seqnames %in% c("chrX", "chrY"))
DxAllCand <- read.csv("DMRs/Replication/Diagnosis 50/CandidateRegions_DxNoXY_Replication50.csv", header = TRUE, 
                      stringsAsFactors = FALSE) %>% subset(!seqnames %in% c("chrX", "chrY"))
DxAllBack <- read.csv("DMRs/Replication/Diagnosis 50/bsseq_background_Replication50.csv", header = TRUE, 
                      stringsAsFactors = FALSE) %>% subset(!chr %in% c("chrX", "chrY"))

DxSexAllDMRs <- read.csv("DMRs/Replication/Diagnosis and Sex 50/DMRs_DxAdjSex_Replication50.csv", header = TRUE, 
                         stringsAsFactors = FALSE) %>% subset(!seqnames %in% c("chrY"))
DxSexAllCand <- read.csv("DMRs/Replication/Diagnosis and Sex 50/backgroundRegions50.csv", header = TRUE, 
                         stringsAsFactors = FALSE) %>% subset(!seqnames %in% c("chrY"))
DxSexAllBack <- read.csv("DMRs/Replication/Diagnosis 50/bsseq_background_Replication50.csv", header = TRUE, 
                         stringsAsFactors = FALSE) %>% subset(!chr %in% c("chrY"))

DxMalesDMRs <- read.csv("DMRs/Replication/Diagnosis Males 50/DMRs_Dx_Replication50_males.csv", header = TRUE, 
                        stringsAsFactors = FALSE)
DxMalesCand <- read.csv("DMRs/Replication/Diagnosis Males 50/CandidateRegions_Dx_Replication50_males.csv", header = TRUE, 
                        stringsAsFactors = FALSE)
DxMalesBack <- read.csv("DMRs/Replication/Diagnosis Males 50/bsseq_background_Replication50_males.csv", header = TRUE, 
                        stringsAsFactors = FALSE)

DxFemalesDMRs <- read.csv("DMRs/Replication/Diagnosis Females 100/DMRs_Dx_Replication100_females.csv", header = TRUE, 
                          stringsAsFactors = FALSE) %>% subset(!seqnames %in% c("chrY"))
DxFemalesCand <- read.csv("DMRs/Replication/Diagnosis Females 100/CandidateRegions_Dx_Replication100_females.csv", header = TRUE, 
                          stringsAsFactors = FALSE) %>% subset(!seqnames %in% c("chrY"))
DxFemalesBack <- read.csv("DMRs/Replication/Diagnosis Females 100/bsseq_background_Replication100_females.csv", header = TRUE, 
                          stringsAsFactors = FALSE) %>% subset(!chr %in% c("chrY"))
DxFemales50DMRs <- read.csv("DMRs/Replication/Diagnosis Females 50/DMRs_Dx_Replication50_females.csv", header = TRUE, 
                          stringsAsFactors = FALSE) %>% subset(!seqnames %in% c("chrY"))

# Region Stats ####
regions <- list(DxAllDMRs, DxAllCand, DxAllBack, DxSexAllDMRs, DxSexAllCand, DxSexAllBack, DxMalesDMRs, DxMalesCand, DxMalesBack, 
                DxFemalesDMRs, DxFemalesCand, DxFemalesBack)
names(regions) <- c("DxAllDMRs", "DxAllCand", "DxAllBack", "DxSexAllDMRs", "DxSexAllCand", "DxSexAllBack", 
                    "DxMalesDMRs", "DxMalesCand", "DxMalesBack", "DxFemalesDMRs", "DxFemalesCand", "DxFemalesBack")
regions <- sapply(regions, function(x){colnames(x)[1] <- "chr"; return(x)})
regionStats <- getRegionStats(regions, TD = c(20, 20, 17, 3), ASD = c(26, 26, 21, 5))
write.table(regionStats, "Tables/DMR Region Stats All Chr Dx Replication 50 Females 100.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Region Stats chrX
regionsX <- regions[!names(regions) %in% c("DxAllDMRs", "DxAllCand", "DxAllBack")] # no chrX
regionsX <- sapply(regions, function(x){subset(x, chr == "chrX")})
regionStatsX <- getRegionStats(regionsX, TD = c(20, 20, 17, 3), ASD = c(26, 26, 21, 5))
write.table(regionStatsX, "Tables/DMR Region Stats ChrX Dx Replication 50 Females 100.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Region Overlap Venn Diagrams ####
library(ChIPpeakAnno)
GR_regions <- sapply(regions, function(x) {GRanges(seqnames = x$chr, ranges = IRanges(start = x$start, end = x$end))})

pdf(file="Figures/DMR Overlap Dx Replication 50 Females 100.pdf", width=10, height=8, onefile = FALSE)
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

pdf(file="Figures/Hyper DMR Overlap Dx Replication 50 Females 100.pdf", width=10, height=8, onefile = FALSE)
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

pdf(file="Figures/Hypo DMR Overlap Dx Replication 50 Females 100.pdf", width=10, height=8, onefile = FALSE)
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

# Diagnosis DMRs (noXY), Diagnosis + Sex DMRs
pdf(file="Figures/Hyper DMR Overlap All Dx or Sex Replication 50.pdf", width=10, height=8, onefile = FALSE)
venn <- suppressMessages(makeVennDiagram(Peaks = list(GR_regions$DxAllDMRs[DxAllDMRs$percentDifference > 0], 
                                                      GR_regions$DxSexAllDMRs[DxSexAllDMRs$percentDifference > 0]),
                                         NameOfPeaks = c("All_DMRs", "All_DMRs_Sex"), 
                                         maxgap = -1, minoverlap = 1, by = "region", connectedPeaks = "min", 
                                         rotation.degree = 0, margin = 0.02, cat.cex = 2, cex = 2.5, 
                                         fill = c("lightblue", "lightpink"), 
                                         cat.pos = c(0, 0), cat.dist = c(-0.02, -0.02), 
                                         fontfamily = "sans", cat.fontfamily = "sans", ext.text = FALSE))
dev.off()

pdf(file="Figures/Hypo DMR Overlap All Dx or Sex Replication 50.pdf", width=10, height=8, onefile = FALSE)
venn <- suppressMessages(makeVennDiagram(Peaks = list(GR_regions$DxAllDMRs[DxAllDMRs$percentDifference < 0], 
                                                      GR_regions$DxSexAllDMRs[DxSexAllDMRs$percentDifference < 0]),
                                         NameOfPeaks = c("All_DMRs", "All_DMRs_Sex"), 
                                         maxgap = -1, minoverlap = 1, by = "region", connectedPeaks = "min", 
                                         rotation.degree = 0, margin = 0.02, cat.cex = 2, cex = 2.5, 
                                         fill = c("lightblue", "lightpink"), 
                                         cat.pos = c(0, 355), cat.dist = c(-0.02, -0.025), 
                                         fontfamily = "sans", cat.fontfamily = "sans", ext.text = FALSE))
dev.off()

# Males Diagnosis DMRs, Females Diagnosis DMRs
pdf(file="Figures/Hyper DMR Overlap Males and Females Dx Replication 50 Females 100.pdf", width=10, height=8, onefile = FALSE)
venn <- suppressMessages(makeVennDiagram(Peaks = list(GR_regions$DxMalesDMRs[DxMalesDMRs$percentDifference > 0], 
                                                      GR_regions$DxFemalesDMRs[DxFemalesDMRs$percentDifference > 0]),
                                         NameOfPeaks = c("Males_DMRs", "Females_DMRs"), 
                                         maxgap = -1, minoverlap = 1, by = "region", connectedPeaks = "min", 
                                         rotation.degree = 0, margin = 0.02, cat.cex = 2, cex = 2.5, 
                                         fill = c("lightblue", "lightpink"), 
                                         cat.pos = c(0, 0), cat.dist = c(0.02, 0.02), 
                                         fontfamily = "sans", cat.fontfamily = "sans", ext.dist = -0.15, 
                                         ext.length = 0.85))
dev.off()

pdf(file="Figures/Hypo DMR Overlap Males and Females Dx Replication 50 Females 100.pdf", width=10, height=8, onefile = FALSE)
venn <- suppressMessages(makeVennDiagram(Peaks = list(GR_regions$DxMalesDMRs[DxMalesDMRs$percentDifference < 0], 
                                                      GR_regions$DxFemalesDMRs[DxFemalesDMRs$percentDifference < 0]),
                                         NameOfPeaks = c("Males_DMRs", "Females_DMRs"), 
                                         maxgap = -1, minoverlap = 1, by = "region", connectedPeaks = "min", 
                                         rotation.degree = 0, margin = 0.02, cat.cex = 2, cex = 2.5, 
                                         fill = c("lightblue", "lightpink"), 
                                         cat.pos = c(0, 0), cat.dist = c(0.02, 0.02), 
                                         fontfamily = "sans", cat.fontfamily = "sans", ext.dist = -0.2, 
                                         ext.length = 0.85))
dev.off()

# Females Diagnosis DMRs 50 vs 100
regions_females <- list(DxFemales50DMRs[DxFemales50DMRs$percentDifference > 0,],
                        DxFemales50DMRs[DxFemales50DMRs$percentDifference < 0,],
                        DxFemalesDMRs[DxFemalesDMRs$percentDifference > 0,],
                        DxFemalesDMRs[DxFemalesDMRs$percentDifference < 0,])
names(regions_females) <- c("DxFemales50DMRsHyper", "DxFemales50DMRsHypo", "DxFemales100DMRsHyper", 
                            "DxFemales100DMRsHypo")
regions_females <- lapply(regions_females, function(x){colnames(x)[1] <- "chr"; return(x)})
GR_regions_females <- sapply(regions_females, function(x) {GRanges(seqnames = x$chr, ranges = IRanges(start = x$start, end = x$end))})

pdf(file="Figures/Hyper DMR Overlap Females Dx Replication 50 vs 100.pdf", width=10, height=8, onefile = FALSE)
venn <- suppressMessages(makeVennDiagram(Peaks = list(GR_regions_females$DxFemales50DMRsHyper,
                                                      GR_regions_females$DxFemales100DMRsHyper),
                                         NameOfPeaks = c("Females50_DMRs", "Females100_DMRs"), 
                                         maxgap = -1, minoverlap = 1, by = "region", connectedPeaks = "min", 
                                         rotation.degree = 0, margin = 0.02, cat.cex = 2, cex = 2.5, 
                                         fill = c("lightblue", "lightpink"), 
                                         cat.pos = c(0, 15), cat.dist = c(0.02, 0.02), 
                                         fontfamily = "sans", cat.fontfamily = "sans", ext.dist = -0.15, 
                                         ext.length = 0.85))
dev.off()

pdf(file="Figures/Hypo DMR Overlap Females Dx Replication 50 vs 100.pdf", width=10, height=8, onefile = FALSE)
venn <- suppressMessages(makeVennDiagram(Peaks = list(GR_regions_females$DxFemales50DMRsHypo,
                                                      GR_regions_females$DxFemales100DMRsHypo),
                                         NameOfPeaks = c("Females50_DMRs", "Females100_DMRs"), 
                                         maxgap = -1, minoverlap = 1, by = "region", connectedPeaks = "min", 
                                         rotation.degree = 0, margin = 0.02, cat.cex = 2, cex = 2.5, 
                                         fill = c("lightblue", "lightpink"), 
                                         cat.pos = c(0, 15), cat.dist = c(0.02, 0.02), 
                                         fontfamily = "sans", cat.fontfamily = "sans", ext.dist = -0.15, 
                                         ext.length = 0.85))
dev.off()

# Male Female Overlap ####
HyperMF <- intersect(x = GR_regions$DxMalesDMRs[DxMalesDMRs$percentDifference > 0],
                     y = GR_regions$DxFemalesDMRs[DxFemalesDMRs$percentDifference > 0]) %>% as.data.frame
HypoMF <- intersect(x = GR_regions$DxMalesDMRs[DxMalesDMRs$percentDifference < 0],
                    y = GR_regions$DxFemalesDMRs[DxFemalesDMRs$percentDifference < 0]) %>% as.data.frame
MF <- rbind(HyperMF, HypoMF)
colnames(MF)[1] <- "chr"
MF$DMRid <- paste("DMR", 1:nrow(MF), sep = "_")
MF$percentDifference <- c(rep(1, nrow(HyperMF)), rep(-1, nrow(HypoMF)))
regDomains <- read.delim("Tables/Regulatory domains hg38.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
HyperGenes <- getDMRgeneList(MF, regDomains = regDomains, direction = "hyper", type = "gene_name")
HypoGenes <- getDMRgeneList(MF, regDomains = regDomains, direction = "hypo", type = "gene_name")

HyperMalesDMRsInF <- subset(DxMalesDMRs[DxMalesDMRs$percentDifference > 0,], GR_regions$DxMalesDMRs[DxMalesDMRs$percentDifference > 0] %over%
                                    GR_regions$DxFemalesDMRs[DxFemalesDMRs$percentDifference > 0])
HypoMalesDMRsInF <- subset(DxMalesDMRs[DxMalesDMRs$percentDifference < 0,], GR_regions$DxMalesDMRs[DxMalesDMRs$percentDifference < 0] %over%
                                   GR_regions$DxFemalesDMRs[DxFemalesDMRs$percentDifference < 0])
MalesDMRsInF <- rbind(HyperMalesDMRsInF, HypoMalesDMRsInF)
write.table(MalesDMRsInF, "Tables/Dx Replication Males DMRs Overlapping with Females 100.txt", sep = "\t", quote = FALSE, row.names = FALSE)

HyperFemalesDMRsInM <- subset(DxFemalesDMRs[DxFemalesDMRs$percentDifference > 0,], GR_regions$DxFemalesDMRs[DxFemalesDMRs$percentDifference > 0] %over%
                                      GR_regions$DxMalesDMRs[DxMalesDMRs$percentDifference > 0])
HypoFemalesDMRsInM <- subset(DxFemalesDMRs[DxFemalesDMRs$percentDifference < 0,], GR_regions$DxFemalesDMRs[DxFemalesDMRs$percentDifference < 0] %over%
                                     GR_regions$DxMalesDMRs[DxMalesDMRs$percentDifference < 0])
FemalesDMRsInM <- rbind(HyperFemalesDMRsInM, HypoFemalesDMRsInM)
write.table(FemalesDMRsInM, "Tables/Dx Replication Females 100 DMRs Overlapping with Males.txt", sep = "\t", quote = FALSE, row.names = FALSE)

rm(HyperMF, HypoMF, HyperMalesDMRsInF, HypoMalesDMRsInF, HyperFemalesDMRsInM, HypoFemalesDMRsInM, venn)

