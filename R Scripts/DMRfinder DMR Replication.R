# DMRfinder DMR Replication -----------------------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 5/11/19

# Packages ####
.libPaths("/share/lasallelab/programs/DMRfinder/R_3.5")
sapply(c("tidyverse", "BiocGenerics", "Biobase", "S4Vectors", "matrixStats", "DelayedArray", "bsseq", "DSS", "permute",
         "GenomicRanges", "scales", "optparse"), require, character.only = TRUE)

# Global Variables ####
mc.cores <- 20
set.seed(5)

# Discovery DMR Comparisons -----------------------------------------------
# Discovery Diagnosis All DMRs (Complete) ####
# Setup Comparison
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Discovery/Dx_All/")
chroms <- c(paste("chr", 1:22, sep = ""), "chrM") # Excluded chrX, chrY
numCtrl <- 56
numExp <- 52
CTRLgroup <- paste("C", 1:numCtrl, sep = "")
EXPgroup <- paste("E", 1:numExp, sep = "")
cutoff <- qt(1 - 0.05 / 2, numCtrl + numExp - 2)
outprefix <- "DMRfinder_Dx_All"

# Load Methylation Data
BSobj <- readRDS("Filtered_Smoothed_BSseq_Discovery50.rds")
BSobj <- chrSelectBSseq(BSobj, seqnames = chroms)
pData <- pData(BSobj) # check that all control samples are first
sampleNames(BSobj) <- c(CTRLgroup, EXPgroup)

# Identify DMRs
tstat <- BSmooth.tstat(BSobj, group1 = EXPgroup, group2 = CTRLgroup, estimate.var = "same", local.correct = TRUE, qSd = 0.75, 
                       k = 101, maxGap = 10^8, mc.cores = mc.cores, verbose = TRUE)  
DMRs <- dmrFinder(tstat, cutoff = c(-cutoff, cutoff), maxGap = 300, verbose = FALSE)
DMRs$end <- DMRs$end + 1 # Add one base to end to include last CpG
DMRs <- subset(DMRs, n >= 3 & abs(meanDiff) > 0.05 & invdensity <= 300)

# Write DMR Table
DMRs <- DMRs[,c("chr", "start", "end", "n", "width", "invdensity", "areaStat", "maxStat", "tstat.sd", "group2.mean", "group1.mean", 
                "meanDiff", "direction")]
colnames(DMRs) <- c("chr", "start", "end", "CpGs", "width", "invdensity", "areaStat", "maxStat", "tstat_sd", "Ctrl_mean", "Exp_mean",
                    "meanDiff", "direction")
write.csv(DMRs, paste(outprefix, "DMRs.csv", sep = "_"), quote = FALSE, row.names = FALSE)

rm(chroms, numCtrl, numExp, CTRLgroup, EXPgroup, cutoff, outprefix, BSobj, pData, tstat, DMRs)

# Discovery Diagnosis Males DMRs (Complete) ####
# Setup Comparison
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Discovery/Dx_Males/")
chroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")
numCtrl <- 39
numExp <- 37
CTRLgroup <- paste("C", 1:numCtrl, sep = "")
EXPgroup <- paste("E", 1:numExp, sep = "")
cutoff <- qt(1 - 0.05 / 2, numCtrl + numExp - 2)
outprefix <- "DMRfinder_Dx_Males"

# Load Methylation Data
BSobj <- readRDS("Filtered_Smoothed_BSseq_Discovery50_males.rds")
BSobj <- chrSelectBSseq(BSobj, seqnames = chroms)
pData <- pData(BSobj) # check that all control samples are first
sampleNames(BSobj) <- c(CTRLgroup, EXPgroup)

# Identify DMRs
tstat <- BSmooth.tstat(BSobj, group1 = EXPgroup, group2 = CTRLgroup, estimate.var = "same", local.correct = TRUE, qSd = 0.75, 
                       k = 101, maxGap = 10^8, mc.cores = mc.cores, verbose = TRUE)  
DMRs <- dmrFinder(tstat, cutoff = c(-cutoff, cutoff), maxGap = 300, verbose = FALSE)
DMRs$end <- DMRs$end + 1 # Add one base to end to include last CpG
DMRs <- subset(DMRs, n >= 3 & abs(meanDiff) > 0.05 & invdensity <= 300)

# Write DMR Table
DMRs <- DMRs[,c("chr", "start", "end", "n", "width", "invdensity", "areaStat", "maxStat", "tstat.sd", "group2.mean", "group1.mean", 
                "meanDiff", "direction")]
colnames(DMRs) <- c("chr", "start", "end", "CpGs", "width", "invdensity", "areaStat", "maxStat", "tstat_sd", "Ctrl_mean", "Exp_mean",
                    "meanDiff", "direction")
write.csv(DMRs, paste(outprefix, "DMRs.csv", sep = "_"), quote = FALSE, row.names = FALSE)
rm(chroms, numCtrl, numExp, CTRLgroup, EXPgroup, cutoff, outprefix, BSobj, pData, tstat, DMRs)

# Discovery Diagnosis Females DMRs (Complete) ####
# Setup Comparison
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Discovery/Dx_Females/")
chroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrM") # Exclude chrY
numCtrl <- 17
numExp <- 15
CTRLgroup <- paste("C", 1:numCtrl, sep = "")
EXPgroup <- paste("E", 1:numExp, sep = "")
cutoff <- qt(1 - 0.05 / 2, numCtrl + numExp - 2)
outprefix <- "DMRfinder_Dx_Females"

# Load Methylation Data
BSobj <- readRDS("Filtered_Smoothed_BSseq_Discovery50_females.rds")
BSobj <- chrSelectBSseq(BSobj, seqnames = chroms)
pData <- pData(BSobj) # check that all control samples are first
sampleNames(BSobj) <- c(CTRLgroup, EXPgroup)

# Identify DMRs
tstat <- BSmooth.tstat(BSobj, group1 = EXPgroup, group2 = CTRLgroup, estimate.var = "same", local.correct = TRUE, qSd = 0.75, 
                       k = 101, maxGap = 10^8, mc.cores = mc.cores, verbose = TRUE)  
DMRs <- dmrFinder(tstat, cutoff = c(-cutoff, cutoff), maxGap = 300, verbose = FALSE)
DMRs$end <- DMRs$end + 1 # Add one base to end to include last CpG
DMRs <- subset(DMRs, n >= 3 & abs(meanDiff) > 0.05 & invdensity <= 300)

# Write DMR Table
DMRs <- DMRs[,c("chr", "start", "end", "n", "width", "invdensity", "areaStat", "maxStat", "tstat.sd", "group2.mean", "group1.mean", 
                "meanDiff", "direction")]
colnames(DMRs) <- c("chr", "start", "end", "CpGs", "width", "invdensity", "areaStat", "maxStat", "tstat_sd", "Ctrl_mean", "Exp_mean",
                    "meanDiff", "direction")
write.csv(DMRs, paste(outprefix, "DMRs.csv", sep = "_"), quote = FALSE, row.names = FALSE)
rm(chroms, numCtrl, numExp, CTRLgroup, EXPgroup, cutoff, outprefix, BSobj, pData, tstat, DMRs)

# Replication DMR Comparisons -----------------------------------------------
# Replication Diagnosis All DMRs (Complete) ####
# Setup Comparison
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Replication/Dx_All/")
chroms <- c(paste("chr", 1:22, sep = ""), "chrM") # Excluded chrX, chrY
numCtrl <- 20
numExp <- 26
CTRLgroup <- paste("C", 1:numCtrl, sep = "")
EXPgroup <- paste("E", 1:numExp, sep = "")
cutoff <- qt(1 - 0.05 / 2, numCtrl + numExp - 2)
outprefix <- "Replication_DMRfinder_Dx_All"

# Load Methylation Data
BSobj <- readRDS("Filtered_Smoothed_BSseq_Replication50.rds")
BSobj <- chrSelectBSseq(BSobj, seqnames = chroms)
pData <- pData(BSobj) # check that all control samples are first
sampleNames(BSobj) <- c(CTRLgroup, EXPgroup)

# Identify DMRs
tstat <- BSmooth.tstat(BSobj, group1 = EXPgroup, group2 = CTRLgroup, estimate.var = "same", local.correct = TRUE, qSd = 0.75, 
                       k = 101, maxGap = 10^8, mc.cores = mc.cores, verbose = TRUE)  
DMRs <- dmrFinder(tstat, cutoff = c(-cutoff, cutoff), maxGap = 300, verbose = FALSE)
DMRs$end <- DMRs$end + 1 # Add one base to end to include last CpG
DMRs <- subset(DMRs, n >= 3 & abs(meanDiff) > 0.05 & invdensity <= 300)

# Write DMR Table
DMRs <- DMRs[,c("chr", "start", "end", "n", "width", "invdensity", "areaStat", "maxStat", "tstat.sd", "group2.mean", "group1.mean", 
                "meanDiff", "direction")]
colnames(DMRs) <- c("chr", "start", "end", "CpGs", "width", "invdensity", "areaStat", "maxStat", "tstat_sd", "Ctrl_mean", "Exp_mean",
                    "meanDiff", "direction")
write.csv(DMRs, paste(outprefix, "DMRs.csv", sep = "_"), quote = FALSE, row.names = FALSE)

rm(chroms, numCtrl, numExp, CTRLgroup, EXPgroup, cutoff, outprefix, BSobj, pData, tstat, DMRs)

# Replication Diagnosis Males DMRs (Complete) ####
# Waiting for Barbera/Epigenerate
# Setup Comparison
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Replication/Dx_Males/")
chroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")
numCtrl <- 17
numExp <- 21
CTRLgroup <- paste("C", 1:numCtrl, sep = "")
EXPgroup <- paste("E", 1:numExp, sep = "")
cutoff <- qt(1 - 0.05 / 2, numCtrl + numExp - 2)
outprefix <- "Replication_DMRfinder_Dx_Males"

# Load Methylation Data
BSobj <- readRDS("Filtered_Smoothed_BSseq_Replication50_males.rds")
BSobj <- chrSelectBSseq(BSobj, seqnames = chroms)
pData <- pData(BSobj) # check that all control samples are first
sampleNames(BSobj) <- c(CTRLgroup, EXPgroup)

# Identify DMRs
tstat <- BSmooth.tstat(BSobj, group1 = EXPgroup, group2 = CTRLgroup, estimate.var = "same", local.correct = TRUE, qSd = 0.75, 
                       k = 101, maxGap = 10^8, mc.cores = mc.cores, verbose = TRUE)  
DMRs <- dmrFinder(tstat, cutoff = c(-cutoff, cutoff), maxGap = 300, verbose = FALSE)
DMRs$end <- DMRs$end + 1 # Add one base to end to include last CpG
DMRs <- subset(DMRs, n >= 3 & abs(meanDiff) > 0.05 & invdensity <= 300)

# Write DMR Table
DMRs <- DMRs[,c("chr", "start", "end", "n", "width", "invdensity", "areaStat", "maxStat", "tstat.sd", "group2.mean", "group1.mean", 
                "meanDiff", "direction")]
colnames(DMRs) <- c("chr", "start", "end", "CpGs", "width", "invdensity", "areaStat", "maxStat", "tstat_sd", "Ctrl_mean", "Exp_mean",
                    "meanDiff", "direction")
write.csv(DMRs, paste(outprefix, "DMRs.csv", sep = "_"), quote = FALSE, row.names = FALSE)
rm(chroms, numCtrl, numExp, CTRLgroup, EXPgroup, cutoff, outprefix, BSobj, pData, tstat, DMRs)

# Replication Diagnosis Females DMRs (Complete) ####
# Setup Comparison
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Replication/Dx_Females_100/")
chroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrM") # Exclude chrY
numCtrl <- 3
numExp <- 5
CTRLgroup <- paste("C", 1:numCtrl, sep = "")
EXPgroup <- paste("E", 1:numExp, sep = "")
cutoff <- qt(1 - 0.05 / 2, numCtrl + numExp - 2)
outprefix <- "Replication_DMRfinder_Dx_Females_100"

# Load Methylation Data
BSobj <- readRDS("Filtered_Smoothed_BSseq_Replication100_females.rds")
BSobj <- chrSelectBSseq(BSobj, seqnames = chroms)
pData <- pData(BSobj) # check that all control samples are first
sampleNames(BSobj) <- c(CTRLgroup, EXPgroup)

# Identify DMRs
tstat <- BSmooth.tstat(BSobj, group1 = EXPgroup, group2 = CTRLgroup, estimate.var = "same", local.correct = TRUE, qSd = 0.75, 
                       k = 101, maxGap = 10^8, mc.cores = mc.cores, verbose = TRUE)  
DMRs <- dmrFinder(tstat, cutoff = c(-cutoff, cutoff), maxGap = 300, verbose = FALSE)
DMRs$end <- DMRs$end + 1 # Add one base to end to include last CpG
DMRs <- subset(DMRs, n >= 3 & abs(meanDiff) > 0.05 & invdensity <= 300)

# Write DMR Table
DMRs <- DMRs[,c("chr", "start", "end", "n", "width", "invdensity", "areaStat", "maxStat", "tstat.sd", "group2.mean", "group1.mean", 
                "meanDiff", "direction")]
colnames(DMRs) <- c("chr", "start", "end", "CpGs", "width", "invdensity", "areaStat", "maxStat", "tstat_sd", "Ctrl_mean", "Exp_mean",
                    "meanDiff", "direction")
write.csv(DMRs, paste(outprefix, "DMRs.csv", sep = "_"), quote = FALSE, row.names = FALSE)
rm(chroms, numCtrl, numExp, CTRLgroup, EXPgroup, cutoff, outprefix, BSobj, pData, tstat, DMRs)

# DMRfinder vs DMRichR Comparison -----------------------------------------
source("R Scripts/DMR Analysis Functions.R")
# Data ####
# Load Regions
DMRichR_DMRs <- list(Disc_All = read.csv("DMRs/Discovery/Diagnosis 50/DMRs_DxNoXY_Discovery50.csv", 
                                         header = TRUE, stringsAsFactors = FALSE),
                     Disc_Males = read.csv("DMRs/Discovery/Diagnosis Males 50/DMRs_Dx_Discovery50_males.csv", 
                                           header = TRUE, stringsAsFactors = FALSE),
                     Disc_Females = read.csv("DMRs/Discovery/Diagnosis Females 50/DMRs_Dx_Discovery50_females.csv", 
                                             header = TRUE, stringsAsFactors = FALSE),
                     Rep_All = read.csv("DMRs/Replication/Diagnosis 50/DMRs_DxNoXY_Replication50.csv",
                                        header = TRUE, stringsAsFactors = FALSE),
                     Rep_Males = read.csv("DMRs/Replication/Diagnosis Males 50/DMRs_Dx_Replication50_males.csv",
                                          header = TRUE, stringsAsFactors = FALSE),
                     Rep_Females = read.csv("DMRs/Replication/Diagnosis Females 100/DMRs_Dx_Replication100_females.csv",
                                            header = TRUE, stringsAsFactors = FALSE))
DMRfinder_DMRs <- list(Disc_All = read.csv("DMRs/Discovery/Diagnosis 50/DMRfinder_Dx_All_DMRs.csv", 
                                           header = TRUE, stringsAsFactors = FALSE),
                       Disc_Males = read.csv("DMRs/Discovery/Diagnosis Males 50/DMRfinder_Dx_Males_DMRs.csv", 
                                             header = TRUE, stringsAsFactors = FALSE),
                       Disc_Females = read.csv("DMRs/Discovery/Diagnosis Females 50/DMRfinder_Dx_Females_DMRs.csv", 
                                               header = TRUE, stringsAsFactors = FALSE),
                       Rep_All = read.csv("DMRs/Replication/Diagnosis 50/Replication_DMRfinder_Dx_All_DMRs.csv",
                                          header = TRUE, stringsAsFactors = FALSE),
                       Rep_Males = read.csv("DMRs/Replication/Diagnosis Males 50/Replication_DMRfinder_Dx_Males_DMRs.csv",
                                            header = TRUE, stringsAsFactors = FALSE),
                       Rep_Females = read.csv("DMRs/Replication/Diagnosis Females 100/Replication_DMRfinder_Dx_Females_100_DMRs.csv",
                                              header = TRUE, stringsAsFactors = FALSE))
background <- list(Disc_All = read.csv("DMRs/Discovery/Diagnosis 50/bsseq_background_Discovery50.csv", 
                                       header = TRUE, stringsAsFactors = FALSE),
                   Disc_Males = read.csv("DMRs/Discovery/Diagnosis Males 50/bsseq_background_Discovery50_males.csv", 
                                         header = TRUE, stringsAsFactors = FALSE),
                   Disc_Females = read.csv("DMRs/Discovery/Diagnosis Females 50/bsseq_background_Discovery50_females.csv", 
                                           header = TRUE, stringsAsFactors = FALSE),
                   Rep_All = read.csv("DMRs/Replication/Diagnosis 50/bsseq_background_Replication50.csv",
                                      header = TRUE, stringsAsFactors = FALSE),
                   Rep_Males = read.csv("DMRs/Replication/Diagnosis Males 50/bsseq_background_Replication50_males.csv",
                                        header = TRUE, stringsAsFactors = FALSE),
                   Rep_Females = read.csv("DMRs/Replication/Diagnosis Females 100/bsseq_background_Replication100_females.csv",
                                          header = TRUE, stringsAsFactors = FALSE))

# Convert to GRanges
GR_DMRichR_DMRs <- sapply(DMRichR_DMRs, function(x) {GRanges(seqnames = x$seqnames, ranges = IRanges(start = x$start, end = x$end))})
GR_DMRfinder_DMRs <- sapply(DMRfinder_DMRs, function(x) {GRanges(seqnames = x$chr, ranges = IRanges(start = x$start, end = x$end))})
GR_background <- sapply(background, function(x) {GRanges(seqnames = x$chr, ranges = IRanges(start = x$start, end = x$end))})

# Get Region Stats ####
DMRichR_DMRstats <- getRegionStats(DMRs = DMRichR_DMRs, background = background, n = c(108, 76, 32, 46, 38, 8))
write.csv(DMRichR_DMRstats, "Tables/DMRichR DMR Stats DMRfinder DMRichR Comparison.csv", row.names = FALSE, quote = FALSE)

DMRfinder_DMRstats <- getRegionStats(DMRs = DMRfinder_DMRs, background = background, n = c(108, 76, 32, 46, 38, 8))
write.csv(DMRfinder_DMRstats, "Tables/DMRfinder DMR Stats DMRfinder DMRichR Comparison.csv", row.names = FALSE, quote = FALSE)

cor(DMRichR_DMRstats$DMR_Number, DMRfinder_DMRstats$DMR_Number) # 0.326
cor(DMRichR_DMRstats$DMR_Width_KB, DMRfinder_DMRstats$DMR_Width_KB) # 0.167
cor(DMRichR_DMRstats$DMR_CpGs, DMRfinder_DMRstats$DMR_CpGs) # 0.226

# Intersect Regions ####
GR_DMRichR_DMRs_Hyper <- mapply(function(x, y) x[y$percentDifference > 0], x = GR_DMRichR_DMRs, y = DMRichR_DMRs)
GR_DMRichR_DMRs_Hypo <- mapply(function(x, y) x[y$percentDifference < 0], x = GR_DMRichR_DMRs, y = DMRichR_DMRs)
GR_DMRfinder_DMRs_Hyper <- mapply(function(x, y) x[y$meanDiff > 0], x = GR_DMRfinder_DMRs, y = DMRfinder_DMRs)
GR_DMRfinder_DMRs_Hypo <- mapply(function(x, y) x[y$meanDiff < 0], x = GR_DMRfinder_DMRs, y = DMRfinder_DMRs)
DMR_Overlap_Hyper <- mapply(function(x, y) sum(x %over% y), x = GR_DMRichR_DMRs_Hyper, y = GR_DMRfinder_DMRs_Hyper)
# Disc_All   Disc_Males Disc_Females      Rep_All    Rep_Males  Rep_Females 
#        9           36          260          143          218          439 
DMR_Overlap_Hypo <- mapply(function(x, y) sum(x %over% y), x = GR_DMRichR_DMRs_Hypo, y = GR_DMRfinder_DMRs_Hypo)
# Disc_All   Disc_Males Disc_Females      Rep_All    Rep_Males  Rep_Females 
#       24           68          371          373          462          454 

DMR_Overlap <- data.frame(Comparison = names(DMRichR_DMRs), DMRichR_DMRs = sapply(GR_DMRichR_DMRs, length),
                          DMRichR_Hyper_DMRs = sapply(GR_DMRichR_DMRs_Hyper, length), 
                          DMRichR_Hypo_DMRs = sapply(GR_DMRichR_DMRs_Hypo, length),
                          DMRfinder_DMRs = sapply(GR_DMRfinder_DMRs, length), 
                          DMRfinder_Hyper_DMRs = sapply(GR_DMRfinder_DMRs_Hyper, length),
                          DMRfinder_Hypo_DMRs = sapply(GR_DMRfinder_DMRs_Hypo, length),
                          All_DMR_Overlap = DMR_Overlap_Hyper + DMR_Overlap_Hypo,
                          Hyper_DMR_Overlap = DMR_Overlap_Hyper, Hypo_DMR_Overlap = DMR_Overlap_Hypo,
                          All_DMR_PerOverlap = (DMR_Overlap_Hyper + DMR_Overlap_Hypo) * 100 / sapply(GR_DMRichR_DMRs, length),
                          Hyper_DMR_PerOverlap = DMR_Overlap_Hyper * 100 / sapply(GR_DMRichR_DMRs_Hyper, length),
                          Hypo_DMR_PerOverlap = DMR_Overlap_Hypo * 100 / sapply(GR_DMRichR_DMRs_Hypo, length))
write.csv(DMR_Overlap, "Tables/Overlap Table DMRfinder DMRichR Comparison.csv", row.names = FALSE, quote = FALSE)

# Venn Diagrams ####
# Discovery All
DMRoverlapVenn(list(GR_DMRichR_DMRs_Hyper$Disc_All, GR_DMRfinder_DMRs_Hyper$Disc_All),
               NameOfPeaks = c("DMRichR_Disc_All_Hyper", "DMRfinder_Disc_All_Hyper"), 
               totalTest = length(GR_background$Disc_All),
               file = "Figures/Discovery All Hyper DMR Overlap DMRfinder DMRichR Comparison Venn.pdf",
               rotation.degree = 180, cat.pos = c(180, 180), cat.dist = c(0.03, 0.03), ext.text = FALSE)
# Hypergeometric Test p = 2.415242e-14 

DMRoverlapVenn(list(GR_DMRichR_DMRs_Hypo$Disc_All, GR_DMRfinder_DMRs_Hypo$Disc_All),
               NameOfPeaks = c("DMRichR_Disc_All_Hypo", "DMRfinder_Disc_All_Hypo"), 
               totalTest = length(GR_background$Disc_All),
               file = "Figures/Discovery All Hypo DMR Overlap DMRfinder DMRichR Comparison Venn.pdf",
               rotation.degree = 180, cat.pos = c(180, 180), cat.dist = c(0.03, 0.03), ext.text = FALSE)
# Hypergeometric Test p = 1.334224e-30 

# Discovery Males
DMRoverlapVenn(list(GR_DMRichR_DMRs_Hyper$Disc_Males, GR_DMRfinder_DMRs_Hyper$Disc_Males),
               NameOfPeaks = c("DMRichR_Disc_Males_Hyper", "DMRfinder_Disc_Males_Hyper"), 
               totalTest = length(GR_background$Disc_Males),
               file = "Figures/Discovery Males Hyper DMR Overlap DMRfinder DMRichR Comparison Venn.pdf",
               rotation.degree = 180, cat.pos = c(190, 180), cat.dist = c(0.03, 0.03), ext.text = FALSE)
#Hypergeometric Test p = 2.314635e-48 

DMRoverlapVenn(list(GR_DMRichR_DMRs_Hypo$Disc_Males, GR_DMRfinder_DMRs_Hypo$Disc_Males),
               NameOfPeaks = c("DMRichR_Disc_Males_Hypo", "DMRfinder_Disc_Males_Hypo"), 
               totalTest = length(GR_background$Disc_Males),
               file = "Figures/Discovery Males Hypo DMR Overlap DMRfinder DMRichR Comparison Venn.pdf",
               rotation.degree = 180, cat.pos = c(180, 180), cat.dist = c(0.03, 0.03), ext.text = FALSE)
# Hypergeometric Test p = 2.735767e-53 

# Discovery Females
DMRoverlapVenn(list(GR_DMRichR_DMRs_Hyper$Disc_Females, GR_DMRfinder_DMRs_Hyper$Disc_Females),
               NameOfPeaks = c("DMRichR_Disc_Females_Hyper", "DMRfinder_Disc_Females_Hyper"), 
               totalTest = length(GR_background$Disc_Females),
               file = "Figures/Discovery Females Hyper DMR Overlap DMRfinder DMRichR Comparison Venn.pdf",
               rotation.degree = 180, cat.pos = c(205, 180), cat.dist = c(0.03, 0.03), ext.text = FALSE)
# Hypergeometric Test p = 2.157491e-132 

DMRoverlapVenn(list(GR_DMRichR_DMRs_Hypo$Disc_Females, GR_DMRfinder_DMRs_Hypo$Disc_Females),
               NameOfPeaks = c("DMRichR_Disc_Females_Hypo", "DMRfinder_Disc_Females_Hypo"), 
               totalTest = length(GR_background$Disc_Females),
               file = "Figures/Discovery Females Hypo DMR Overlap DMRfinder DMRichR Comparison Venn.pdf",
               rotation.degree = 180, cat.pos = c(200, 180), cat.dist = c(0.03, 0.03), ext.text = FALSE)
# Hypergeometric Test p = 1.725105e-203 

# Replication All
DMRoverlapVenn(list(GR_DMRichR_DMRs_Hyper$Rep_All, GR_DMRfinder_DMRs_Hyper$Rep_All),
               NameOfPeaks = c("DMRichR_Rep_All_Hyper", "DMRfinder_Rep_All_Hyper"), 
               totalTest = length(GR_background$Rep_All),
               file = "Figures/Replication All Hyper DMR Overlap DMRfinder DMRichR Comparison Venn.pdf",
               rotation.degree = 180, cat.pos = c(180, 180), cat.dist = c(0.03, 0.03), ext.text = FALSE)
# Hypergeometric Test p = 9.841473e-95 

DMRoverlapVenn(list(GR_DMRichR_DMRs_Hypo$Rep_All, GR_DMRfinder_DMRs_Hypo$Rep_All),
               NameOfPeaks = c("DMRichR_Rep_All_Hypo", "DMRfinder_Rep_All_Hypo"), 
               totalTest = length(GR_background$Rep_All),
               file = "Figures/Replication All Hypo DMR Overlap DMRfinder DMRichR Comparison Venn.pdf",
               rotation.degree = 180, cat.pos = c(180, 180), cat.dist = c(0.03, 0.03), ext.text = FALSE)
# Hypergeometric Test p = 9.735604e-74 

# Replication Males
DMRoverlapVenn(list(GR_DMRichR_DMRs_Hyper$Rep_Males, GR_DMRfinder_DMRs_Hyper$Rep_Males),
               NameOfPeaks = c("DMRichR_Rep_Males_Hyper", "DMRfinder_Rep_Males_Hyper"), 
               totalTest = length(GR_background$Rep_Males),
               file = "Figures/Replication Males Hyper DMR Overlap DMRfinder DMRichR Comparison Venn.pdf",
               rotation.degree = 180, cat.pos = c(190, 180), cat.dist = c(0.03, 0.03), ext.text = FALSE)
# Hypergeometric Test p = 3.690528e-125 

DMRoverlapVenn(list(GR_DMRichR_DMRs_Hypo$Rep_Males, GR_DMRfinder_DMRs_Hypo$Rep_Males),
               NameOfPeaks = c("DMRichR_Rep_Males_Hypo", "DMRfinder_Rep_Males_Hypo"), 
               totalTest = length(GR_background$Rep_Males),
               file = "Figures/Replication Males Hypo DMR Overlap DMRfinder DMRichR Comparison Venn.pdf",
               rotation.degree = 180, cat.pos = c(180, 180), cat.dist = c(0.03, 0.03), ext.text = FALSE)
# Hypergeometric Test p = 1.032171e-100 

# Replication Females
DMRoverlapVenn(list(GR_DMRichR_DMRs_Hyper$Rep_Females, GR_DMRfinder_DMRs_Hyper$Rep_Females),
               NameOfPeaks = c("DMRichR_Rep_Females_Hyper", "DMRfinder_Rep_Females_Hyper"), 
               totalTest = length(GR_background$Rep_Females),
               file = "Figures/Replication Females Hyper DMR Overlap DMRfinder DMRichR Comparison Venn.pdf",
               rotation.degree = 180, cat.pos = c(205, 180), cat.dist = c(0.03, 0.03), ext.text = FALSE)
# Hypergeometric Test p = 1.3363e-267 

DMRoverlapVenn(list(GR_DMRichR_DMRs_Hypo$Rep_Females, GR_DMRfinder_DMRs_Hypo$Rep_Females),
               NameOfPeaks = c("DMRichR_Rep_Females_Hypo", "DMRfinder_Rep_Females_Hypo"), 
               totalTest = length(GR_background$Rep_Females),
               file = "Figures/Replication Females Hypo DMR Overlap DMRfinder DMRichR Comparison Venn.pdf",
               rotation.degree = 180, cat.pos = c(200, 180), cat.dist = c(0.03, 0.03), ext.text = FALSE)
# Hypergeometric Test p = 0 


