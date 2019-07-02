# DMRfinder DMR Replication -----------------------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 7/1/19
# Excluded JLCM032B and JLCM050B
# Removed all comparison and updated overlap stats

# Packages ####
.libPaths("/share/lasallelab/programs/DMRfinder/R_3.5")
sapply(c("tidyverse", "BiocGenerics", "Biobase", "S4Vectors", "matrixStats", "DelayedArray", "bsseq", "DSS", "permute",
         "GenomicRanges", "scales", "optparse", "regioneR", "LOLA"), require, character.only = TRUE)

# Global Variables ####
mc.cores <- 20
set.seed(5)

# Discovery DMR Comparisons -----------------------------------------------
# Discovery Diagnosis All DMRs (Rerun without JLCM032B and JLCM050B, Complete) ####
# Setup Comparison
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Discovery/Dx_All/")
chroms <- c(paste("chr", 1:22, sep = ""), "chrM") # Excluded chrX, chrY
numCtrl <- 56
numExp <- 50
CTRLgroup <- paste("C", 1:numCtrl, sep = "")
EXPgroup <- paste("E", 1:numExp, sep = "")
cutoff <- qt(1 - 0.05 / 2, numCtrl + numExp - 2)
outprefix <- "DMRfinder_Dx_All"

# Load Methylation Data
BSobj <- readRDS("Filtered_Smoothed_BSseq_Discovery50.rds")
BSobj <- chrSelectBSseq(BSobj, seqnames = chroms)
c("JLCM032B", "JLCM050B") %in% sampleNames(BSobj)
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

# Discovery Diagnosis Males DMRs (Rerun without JLCM032B and JLCM050B, Complete) ####
# Setup Comparison
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Discovery/Dx_Males/")
chroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")
numCtrl <- 39
numExp <- 35
CTRLgroup <- paste("C", 1:numCtrl, sep = "")
EXPgroup <- paste("E", 1:numExp, sep = "")
cutoff <- qt(1 - 0.05 / 2, numCtrl + numExp - 2)
outprefix <- "DMRfinder_Dx_Males"

# Load Methylation Data
BSobj <- readRDS("Filtered_Smoothed_BSseq_Discovery50_males.rds")
BSobj <- chrSelectBSseq(BSobj, seqnames = chroms)
c("JLCM032B", "JLCM050B") %in% sampleNames(BSobj)
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

# DMRfinder vs DMRichR Comparison (Rerun without JLCM032B and JLCM050B, Complete) -----------------------------------------
source("R Scripts/DMR Analysis Functions.R")
# Data ####
# Load Regions
DMRichR_DMRs <- list(Disc_Males = read.csv("DMRs/Discovery/Diagnosis Males 50/DMRs_Dx_Discovery50_males.csv", 
                                           header = TRUE, stringsAsFactors = FALSE),
                     Disc_Females = read.csv("DMRs/Discovery/Diagnosis Females 50/DMRs_Dx_Discovery50_females.csv", 
                                             header = TRUE, stringsAsFactors = FALSE),
                     Rep_Males = read.csv("DMRs/Replication/Diagnosis Males 50/DMRs_Dx_Replication50_males.csv",
                                          header = TRUE, stringsAsFactors = FALSE),
                     Rep_Females = read.csv("DMRs/Replication/Diagnosis Females 100/DMRs_Dx_Replication100_females.csv",
                                            header = TRUE, stringsAsFactors = FALSE))
DMRfinder_DMRs <- list(Disc_Males = read.csv("DMRs/Discovery/Diagnosis Males 50/DMRfinder_Dx_Males_DMRs.csv", 
                                             header = TRUE, stringsAsFactors = FALSE),
                       Disc_Females = read.csv("DMRs/Discovery/Diagnosis Females 50/DMRfinder_Dx_Females_DMRs.csv", 
                                               header = TRUE, stringsAsFactors = FALSE),
                       Rep_Males = read.csv("DMRs/Replication/Diagnosis Males 50/Replication_DMRfinder_Dx_Males_DMRs.csv",
                                            header = TRUE, stringsAsFactors = FALSE),
                       Rep_Females = read.csv("DMRs/Replication/Diagnosis Females 100/Replication_DMRfinder_Dx_Females_100_DMRs.csv",
                                              header = TRUE, stringsAsFactors = FALSE))
background <- list(Disc_Males = read.csv("DMRs/Discovery/Diagnosis Males 50/bsseq_background_Discovery50_males.csv", 
                                         header = TRUE, stringsAsFactors = FALSE),
                   Disc_Females = read.csv("DMRs/Discovery/Diagnosis Females 50/bsseq_background_Discovery50_females.csv", 
                                           header = TRUE, stringsAsFactors = FALSE),
                   Rep_Males = read.csv("DMRs/Replication/Diagnosis Males 50/bsseq_background_Replication50_males.csv",
                                        header = TRUE, stringsAsFactors = FALSE),
                   Rep_Females = read.csv("DMRs/Replication/Diagnosis Females 100/bsseq_background_Replication100_females.csv",
                                          header = TRUE, stringsAsFactors = FALSE))

# Convert to GRanges
GR_DMRichR_DMRs <- sapply(DMRichR_DMRs, makeGRange, direction = "all")
GR_DMRfinder_DMRs <- sapply(DMRfinder_DMRs, makeGRange, direction = "all")
GR_background <- sapply(background, makeGRange, direction = "all")

# Get Region Stats ####
DMRichR_DMRstats <- getRegionStats(DMRs = DMRichR_DMRs, background = background, n = c(76, 32, 38, 8))
write.csv(DMRichR_DMRstats, "Tables/DMRichR DMR Stats DMRfinder DMRichR Comparison.csv", row.names = FALSE, quote = FALSE)

DMRfinder_DMRstats <- getRegionStats(DMRs = DMRfinder_DMRs, background = background, n = c(76, 32, 38, 8))
write.csv(DMRfinder_DMRstats, "Tables/DMRfinder DMR Stats DMRfinder DMRichR Comparison.csv", row.names = FALSE, quote = FALSE)

# Intersect Regions ####
GR_DMRichR_DMRs_Hyper <- mapply(function(x, y) x[y$percentDifference > 0], x = GR_DMRichR_DMRs, y = DMRichR_DMRs)
GR_DMRichR_DMRs_Hypo <- mapply(function(x, y) x[y$percentDifference < 0], x = GR_DMRichR_DMRs, y = DMRichR_DMRs)
GR_DMRfinder_DMRs_Hyper <- mapply(function(x, y) x[y$meanDiff > 0], x = GR_DMRfinder_DMRs, y = DMRfinder_DMRs)
GR_DMRfinder_DMRs_Hypo <- mapply(function(x, y) x[y$meanDiff < 0], x = GR_DMRfinder_DMRs, y = DMRfinder_DMRs)
DMR_Overlap_Hyper <- mapply(function(x, y) sum(x %over% y), x = GR_DMRichR_DMRs_Hyper, y = GR_DMRfinder_DMRs_Hyper)
# Disc_Males Disc_Females    Rep_Males  Rep_Females 
#         29          260          218          439
DMR_Overlap_Hypo <- mapply(function(x, y) sum(x %over% y), x = GR_DMRichR_DMRs_Hypo, y = GR_DMRfinder_DMRs_Hypo)
# Disc_Males Disc_Females    Rep_Males  Rep_Females 
#         68          371          462          454 

DMR_Overlap <- data.frame(Comparison = names(DMRichR_DMRs), DMRichR_DMRs = sapply(GR_DMRichR_DMRs, length),
                          DMRichR_Hyper_DMRs = sapply(GR_DMRichR_DMRs_Hyper, length), 
                          DMRichR_Hypo_DMRs = sapply(GR_DMRichR_DMRs_Hypo, length),
                          DMRfinder_DMRs = sapply(GR_DMRfinder_DMRs, length), 
                          DMRfinder_Hyper_DMRs = sapply(GR_DMRfinder_DMRs_Hyper, length),
                          DMRfinder_Hypo_DMRs = sapply(GR_DMRfinder_DMRs_Hypo, length),
                          Hyper_DMR_Overlap = DMR_Overlap_Hyper, Hypo_DMR_Overlap = DMR_Overlap_Hypo,
                          Hyper_DMR_PerOverlap = DMR_Overlap_Hyper * 100 / sapply(GR_DMRichR_DMRs_Hyper, length),
                          Hypo_DMR_PerOverlap = DMR_Overlap_Hypo * 100 / sapply(GR_DMRichR_DMRs_Hypo, length))
write.csv(DMR_Overlap, "Tables/Overlap Table DMRfinder DMRichR Comparison.csv", row.names = FALSE, quote = FALSE)

# Venn Diagrams ####
# Discovery Males
DMRoverlapVenn(list(GR_DMRichR_DMRs_Hyper$Disc_Males, GR_DMRfinder_DMRs_Hyper$Disc_Males),
               NameOfPeaks = c("DMRichR", "DMRfinder"), 
               file = "Figures/Discovery Males Hyper DMR Overlap DMRfinder DMRichR Comparison Venn.pdf",
               rotation.degree = 180, cat.pos = c(170, 180), cat.dist = c(0.03, 0.03), ext.text = FALSE)
DMRoverlapVenn(list(GR_DMRichR_DMRs_Hypo$Disc_Males, GR_DMRfinder_DMRs_Hypo$Disc_Males),
               NameOfPeaks = c("DMRichR", "DMRfinder"), 
               file = "Figures/Discovery Males Hypo DMR Overlap DMRfinder DMRichR Comparison Venn.pdf",
               rotation.degree = 180, cat.pos = c(180, 180), cat.dist = c(0.03, 0.03), ext.text = FALSE)

# Discovery Females
DMRoverlapVenn(list(GR_DMRichR_DMRs_Hyper$Disc_Females, GR_DMRfinder_DMRs_Hyper$Disc_Females),
               NameOfPeaks = c("DMRichR", "DMRfinder"), 
               file = "Figures/Discovery Females Hyper DMR Overlap DMRfinder DMRichR Comparison Venn.pdf",
               rotation.degree = 180, cat.pos = c(155, 180), cat.dist = c(0.03, 0.03), ext.text = FALSE)
DMRoverlapVenn(list(GR_DMRichR_DMRs_Hypo$Disc_Females, GR_DMRfinder_DMRs_Hypo$Disc_Females),
               NameOfPeaks = c("DMRichR", "DMRfinder"), 
               file = "Figures/Discovery Females Hypo DMR Overlap DMRfinder DMRichR Comparison Venn.pdf",
               rotation.degree = 180, cat.pos = c(155, 180), cat.dist = c(0.03, 0.03), ext.text = FALSE)

# Replication Males
DMRoverlapVenn(list(GR_DMRichR_DMRs_Hyper$Rep_Males, GR_DMRfinder_DMRs_Hyper$Rep_Males),
               NameOfPeaks = c("DMRichR", "DMRfinder"), 
               file = "Figures/Replication Males Hyper DMR Overlap DMRfinder DMRichR Comparison Venn.pdf",
               rotation.degree = 180, cat.pos = c(170, 180), cat.dist = c(0.03, 0.03), ext.text = FALSE)
DMRoverlapVenn(list(GR_DMRichR_DMRs_Hypo$Rep_Males, GR_DMRfinder_DMRs_Hypo$Rep_Males),
               NameOfPeaks = c("DMRichR", "DMRfinder"), 
               file = "Figures/Replication Males Hypo DMR Overlap DMRfinder DMRichR Comparison Venn.pdf",
               rotation.degree = 180, cat.pos = c(180, 180), cat.dist = c(0.03, 0.03), ext.text = FALSE)

# Replication Females
DMRoverlapVenn(list(GR_DMRichR_DMRs_Hyper$Rep_Females, GR_DMRfinder_DMRs_Hyper$Rep_Females),
               NameOfPeaks = c("DMRichR", "DMRfinder"), 
               file = "Figures/Replication Females Hyper DMR Overlap DMRfinder DMRichR Comparison Venn.pdf",
               rotation.degree = 180, cat.pos = c(180, 180), cat.dist = c(0.03, 0.03), ext.text = FALSE)
DMRoverlapVenn(list(GR_DMRichR_DMRs_Hypo$Rep_Females, GR_DMRfinder_DMRs_Hypo$Rep_Females),
               NameOfPeaks = c("DMRichR", "DMRfinder"), 
               file = "Figures/Replication Females Hypo DMR Overlap DMRfinder DMRichR Comparison Venn.pdf",
               rotation.degree = 0, cat.pos = c(0, 0), cat.dist = c(0.03, 0.03), ext.text = FALSE)

# RegioneR Stats ####
options("mc.cores" = 12) # sets 12 cores for parallelization
set.seed(5)
hg38_XY <- getGenomeAndMask("hg38", mask = NULL)$genome %>%
        filterChromosomes(keep.chr = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"))
hg38_X <- getGenomeAndMask("hg38", mask = NULL)$genome %>%
        filterChromosomes(keep.chr = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"))

# Discovery Males
stats_Disc_Males_Hyper <- DMRpermTest(A = redefineUserSets(GRangesList(GR_DMRichR_DMRs_Hyper$Disc_Males), 
                                                           GR_background$Disc_Males)[[1]], 
                                      B = redefineUserSets(GRangesList(GR_DMRfinder_DMRs_Hyper$Disc_Males), 
                                                           GR_background$Disc_Males)[[1]],
                                      genome = hg38_XY, universe = GR_background$Disc_Males, Comparison = "Disc_Males_Hyper",
                                      file = "Figures/Hyper DMR Overlap Discovery Males DMRichR vs DMRfinder RegioneR Plots.pdf")
stats_Disc_Males_Hypo <- DMRpermTest(A = redefineUserSets(GRangesList(GR_DMRichR_DMRs_Hypo$Disc_Males), 
                                                          GR_background$Disc_Males)[[1]], 
                                     B = redefineUserSets(GRangesList(GR_DMRfinder_DMRs_Hypo$Disc_Males), 
                                                          GR_background$Disc_Males)[[1]],
                                     genome = hg38_XY, universe = GR_background$Disc_Males, Comparison = "Disc_Males_Hypo",
                                     file = "Figures/Hypo DMR Overlap Discovery Males DMRichR vs DMRfinder RegioneR Plots.pdf")
# Discovery Females
stats_Disc_Females_Hyper <- DMRpermTest(A = redefineUserSets(GRangesList(GR_DMRichR_DMRs_Hyper$Disc_Females), 
                                                             GR_background$Disc_Females)[[1]], 
                                        B = redefineUserSets(GRangesList(GR_DMRfinder_DMRs_Hyper$Disc_Females), 
                                                             GR_background$Disc_Females)[[1]],
                                        genome = hg38_X, universe = GR_background$Disc_Females, Comparison = "Disc_Females_Hyper",
                                        file = "Figures/Hyper DMR Overlap Discovery Females DMRichR vs DMRfinder RegioneR Plots.pdf")
stats_Disc_Females_Hypo <- DMRpermTest(A = redefineUserSets(GRangesList(GR_DMRichR_DMRs_Hypo$Disc_Females), 
                                                            GR_background$Disc_Females)[[1]], 
                                       B = redefineUserSets(GRangesList(GR_DMRfinder_DMRs_Hypo$Disc_Females), 
                                                            GR_background$Disc_Females)[[1]],
                                       genome = hg38_X, universe = GR_background$Disc_Females, Comparison = "Disc_Females_Hypo",
                                       file = "Figures/Hypo DMR Overlap Discovery Females DMRichR vs DMRfinder RegioneR Plots.pdf")
# Replication Males
stats_Rep_Males_Hyper <- DMRpermTest(A = redefineUserSets(GRangesList(GR_DMRichR_DMRs_Hyper$Rep_Males), 
                                                          GR_background$Rep_Males)[[1]], 
                                     B = redefineUserSets(GRangesList(GR_DMRfinder_DMRs_Hyper$Rep_Males), 
                                                          GR_background$Rep_Males)[[1]],
                                     genome = hg38_XY, universe = GR_background$Rep_Males, Comparison = "Rep_Males_Hyper",
                                     file = "Figures/Hyper DMR Overlap Replication Males DMRichR vs DMRfinder RegioneR Plots.pdf")
stats_Rep_Males_Hypo <- DMRpermTest(A = redefineUserSets(GRangesList(GR_DMRichR_DMRs_Hypo$Rep_Males), 
                                                         GR_background$Rep_Males)[[1]], 
                                    B = redefineUserSets(GRangesList(GR_DMRfinder_DMRs_Hypo$Rep_Males), 
                                                         GR_background$Rep_Males)[[1]],
                                    genome = hg38_XY, universe = GR_background$Rep_Males, Comparison = "Rep_Males_Hypo",
                                    file = "Figures/Hypo DMR Overlap Replication Males DMRichR vs DMRfinder RegioneR Plots.pdf")
# Replication Females
stats_Rep_Females_Hyper <- DMRpermTest(A = redefineUserSets(GRangesList(GR_DMRichR_DMRs_Hyper$Rep_Females), 
                                                            GR_background$Rep_Females)[[1]], 
                                       B = redefineUserSets(GRangesList(GR_DMRfinder_DMRs_Hyper$Rep_Females), 
                                                            GR_background$Rep_Females)[[1]],
                                       genome = hg38_X, universe = GR_background$Rep_Females, Comparison = "Rep_Females_Hyper",
                                       file = "Figures/Hyper DMR Overlap Replication Females DMRichR vs DMRfinder RegioneR Plots.pdf")
stats_Rep_Females_Hypo <- DMRpermTest(A = redefineUserSets(GRangesList(GR_DMRichR_DMRs_Hypo$Rep_Females), 
                                                           GR_background$Rep_Females)[[1]], 
                                      B = redefineUserSets(GRangesList(GR_DMRfinder_DMRs_Hypo$Rep_Females), 
                                                           GR_background$Rep_Females)[[1]],
                                      genome = hg38_X, universe = GR_background$Rep_Females, Comparison = "Rep_Females_Hypo",
                                      file = "Figures/Hypo DMR Overlap Replication Females DMRichR vs DMRfinder RegioneR Plots.pdf")
# Combine
stats_All <- rbind(stats_Disc_Males_Hyper, stats_Disc_Males_Hypo, stats_Disc_Females_Hyper, stats_Disc_Females_Hypo,
                   stats_Rep_Males_Hyper, stats_Rep_Males_Hypo, stats_Rep_Females_Hyper, stats_Rep_Females_Hypo)
write.csv(stats_All, file = "Tables/DMR Overlap DMRichR vs DMRfinder RegioneR Stats.csv", quote = FALSE, row.names = FALSE)


