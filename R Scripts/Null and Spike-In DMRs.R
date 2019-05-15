# Null and Spike-In DMRs --------------------------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 4/4/19

# Load Packages ####
.libPaths("/share/lasallelab/programs/DMRichR/R_3.5")
sapply(c("tidyverse", "openxlsx", "bsseq", "dmrseq", "DMRichR", "ChIPpeakAnno"), require, character.only = TRUE)

# Global Variables ####
genome <- as.character("hg38")
coverage <- as.numeric(1)
perGroup <- (as.numeric(50)/100)
minCpGs <- as.numeric(3)
maxPerms <- as.numeric(10)
testCovariate <- as.character("Diagnosis")
cores <- 10
set.seed(5)
register(MulticoreParam(1))

# Sample Meta Data (laptop) ####
# Get male TD samples
meta <- read.xlsx("DMRs/Discovery/Diagnosis Males 50/sample_info_males.xlsx", colNames = TRUE) 
meta <- subset(meta, Diagnosis == "Ctrl_TD" & Sex == "M", select = c("Name", "Diagnosis", "Sex"))

# Assign null diagnosis groups
meta$Diagnosis <- "Ctrl_Null"
Exp_Null <- sample(1:nrow(meta), size = round(nrow(meta)/2))
meta$Diagnosis[Exp_Null] <- "Exp_Null"
table(meta$Diagnosis)
# Ctrl_Null  Exp_Null 
#        19        20 
write.xlsx(meta, file = "DMRs/Discovery/Null and Spike In/sample_info_null_spikein.xlsx")

# Stats for Replication Males Diagnosis DMRs (True Comparison) ####
DMRs <- read.csv("DMRs/Replication/Diagnosis Males 50/DMRs_Dx_Replication50_males.csv", header = TRUE,
                 stringsAsFactors = FALSE)
nrow(DMRs) #4650
abs(DMRs$percentDifference / 100) %>% mean # 0.098
# For simulation, use round numbers: nDMRs = 5000, methDiff = 0.1
nDMRs <- 5000
methDiff <- 0.1

# Null DMRs ---------------------------------------------------------------
# Load and Process Samples (Complete)
bs.filtered <- processBismark(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                              meta = read.xlsx("sample_info_null_spikein.xlsx", colNames = TRUE) %>% 
                                      mutate_if(is.character, as.factor), groups = testCovariate, Cov = coverage, 
                              mc.cores = cores, per.Group = perGroup)
saveRDS(bs.filtered, "Filtered_BSseq_Discovery50_null_spikein.rds")

# Background Regions (Complete)
background <- getBackground(bs.filtered, minNumRegion = minCpGs)
write.csv(background, file = "bsseq_background_Discovery50_null.csv", quote = FALSE, row.names = FALSE)

# Null Meth ~ Diagnosis DMRs (Complete)
regions <- dmrseq(bs = bs.filtered, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms, testCovariate = testCovariate)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "CandidateRegions_Dx_Discovery50_null.csv")
gr2csv(sigRegions, "DMRs_Dx_Discovery50_null.csv")
rm(background, regions, sigRegions)

# Spike-In DMRs ---------------------------------------------------------------
source("R Scripts/simDMRs2.R")
bs.filtered <- readRDS("Filtered_BSseq_Discovery50_null_spikein.rds")
bs.filtered <- bs.filtered[,c(which(pData(bs.filtered)$Diagnosis == "Ctrl_Null"),
                              which(pData(bs.filtered)$Diagnosis == "Exp_Null"))] # Reorder bs.filtered so groups are correct
# simDMRs2 assigns first n = sampleSize samples as condition 1 and remaining samples as condition 2

# Simulate DMRs, Spiked-In DMRs are expected values for a true comparison based on replication males (Epigenerate 4/5)
sim <- simDMRs2(bs.filtered, num.dmrs = 5000, delta.max0 = 0.1, sampleSize = sum(pData(bs.filtered)$Diagnosis == "Ctrl_Null")) 
saveRDS(sim, "Filtered_BSseq_Discovery50_with_spikein_DMRs.rds")

sim <- readRDS("Filtered_BSseq_Discovery50_with_spikein_DMRs.rds")
bs.sim <- sim$bs
pData(bs.sim) <- pData(bs.filtered)
rm(sim, bs.filtered)

# Background Regions (Complete)
background <- getBackground(bs.sim, minNumRegion = minCpGs)
write.csv(background, file = "bsseq_background_Discovery50_spikein.csv", quote = FALSE, row.names = FALSE)

# Spike-In Meth ~ Diagnosis DMRs (Complete)
regions <- dmrseq(bs = bs.sim, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms, testCovariate = testCovariate)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "CandidateRegions_Dx_Discovery50_spikein.csv")
gr2csv(sigRegions, "DMRs_Dx_Discovery50_spikein.csv")

# Permuted Null DMRs ---------------------------------------------------------------
# Sample Meta Data (laptop) ####
set.seed(5)

# Get male samples
meta <- read.xlsx("DMRs/Discovery/Diagnosis Males 50/sample_info_males.xlsx", colNames = TRUE) 
meta <- meta[ ,c("Name", "Diagnosis", "Sex")]

# Assign permuted null diagnosis groups
meta$Diagnosis <- "Ctrl_Null"
Exp_Null <- sample(1:nrow(meta), size = round(nrow(meta)/2))
meta$Diagnosis[Exp_Null] <- "Exp_Null"
table(meta$Diagnosis)
# Ctrl_Null  Exp_Null 
#        38        38 
write.xlsx(meta, file = "DMRs/Discovery/Null and Spike In/sample_info_permuted_null_spikein.xlsx")

# Load and Process Samples (Done) ####
bs.filtered <- processBismark(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                              meta = read.xlsx("sample_info_permuted_null_spikein.xlsx", colNames = TRUE) %>% 
                                      mutate_if(is.character, as.factor), groups = testCovariate, Cov = coverage, 
                              mc.cores = cores, per.Group = perGroup)
saveRDS(bs.filtered, "Filtered_BSseq_Discovery50_permuted_null_spikein.rds")

# Background Regions (Done)
background <- getBackground(bs.filtered, minNumRegion = minCpGs)
write.csv(background, file = "bsseq_background_Discovery50_permuted_null.csv", quote = FALSE, row.names = FALSE)

# Null Meth ~ Diagnosis DMRs (Done)
regions <- dmrseq(bs = bs.filtered, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms, testCovariate = testCovariate)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "CandidateRegions_Dx_Discovery50_permuted_null.csv")
gr2csv(sigRegions, "DMRs_Dx_Discovery50_permuted_null.csv")
rm(background, regions, sigRegions)

# Permuted Spike-In DMRs ---------------------------------------------------------------
# Simulate DMRs (Done) ####
source("R Scripts/simDMRs2.R")
#bs.filtered <- readRDS("Filtered_BSseq_Discovery50_permuted_null_spikein.rds")
bs.filtered <- bs.filtered[,c(which(pData(bs.filtered)$Diagnosis == "Ctrl_Null"),
                              which(pData(bs.filtered)$Diagnosis == "Exp_Null"))] # Reorder bs.filtered so groups are correct
# simDMRs2 assigns first n = sampleSize samples as condition 1 and remaining samples as condition 2

# Simulate DMRs, Spiked-In DMRs are expected values for a true comparison based on replication males
sim <- simDMRs2(bs.filtered, num.dmrs = 5000, delta.max0 = 0.1, sampleSize = sum(pData(bs.filtered)$Diagnosis == "Ctrl_Null")) 
saveRDS(sim, "Filtered_BSseq_Discovery50_permuted_with_spikein_DMRs.rds")

#sim <- readRDS("Filtered_BSseq_Discovery50_permuted_with_spikein_DMRs.rds")
bs.sim <- sim$bs
pData(bs.sim) <- pData(bs.filtered)
rm(sim, bs.filtered)

# Background Regions (Epigenerate 4/16) ####
background <- getBackground(bs.sim, minNumRegion = minCpGs)
write.csv(background, file = "bsseq_background_Discovery50_permuted_spikein.csv", quote = FALSE, row.names = FALSE)

# Spike-In Meth ~ Diagnosis DMRs (Epigenerate 4/16) ####
regions <- dmrseq(bs = bs.sim, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms, testCovariate = testCovariate)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "CandidateRegions_Dx_Discovery50_permuted_spikein.csv")
gr2csv(sigRegions, "DMRs_Dx_Discovery50_permuted_spikein.csv")

# Sensitivity and Specificity ---------------------------------------------
# See Korthauer et al. 2018 Figure 3 and Table 1
source("R Scripts/DMR Analysis Functions.R")

# Regions ####
# All DMR Comparisons used male samples

# Get Spiked In DMRs from cluster
Null_SpikeIn <- readRDS("Filtered_BSseq_Discovery50_with_spikein_DMRs.rds")
Null_SpikeIn_DMRs <- as.data.frame(Null_SpikeIn$gr.dmrs)
Null_SpikeIn_DMRs$mncov <- Null_SpikeIn$dmr.mncov
Null_SpikeIn_DMRs$L <- Null_SpikeIn$dmr.L
Null_SpikeIn_DMRs$delta <- Null_SpikeIn$delta
write.csv(Null_SpikeIn_DMRs, "SpikedInDMRsOnly_DxDiscovery50_spikein.csv", quote = FALSE, row.names = FALSE)

Permute_SpikeIn <- readRDS("Filtered_BSseq_Discovery50_permuted_with_spikein_DMRs.rds")
Permute_SpikeIn_DMRs <- as.data.frame(Permute_SpikeIn$gr.dmrs)
Permute_SpikeIn_DMRs$mncov <- Permute_SpikeIn$dmr.mncov
Permute_SpikeIn_DMRs$L <- Permute_SpikeIn$dmr.L
Permute_SpikeIn_DMRs$delta <- Permute_SpikeIn$delta
write.csv(Permute_SpikeIn_DMRs, "SpikedInDMRsOnly_DxDiscovery50_permuted_spikein.csv", quote = FALSE, row.names = FALSE)

# DMRs
DMRs <- list(DMRs_Dx = read.csv("DMRs/Discovery/Diagnosis Males 50/DMRs_Dx_Discovery50_males.csv", header = TRUE, 
                                stringsAsFactors = FALSE), # ASD vs TD, All Samples, Discovery, 39 TD, 37 ASD
             DMRs_DxRed = read.csv("DMRs/Discovery/Diagnosis Males 50/DMRs_Dx_Discovery50_reduced_males.csv", header = TRUE, 
                                   stringsAsFactors = FALSE), # ASD vs TD, Reduced Samples, Discovery, 17 TD, 21 ASD
             DMRs_DxRep = read.csv("DMRs/Replication/Diagnosis Males 50/DMRs_Dx_Replication50_males.csv", header = TRUE,
                                   stringsAsFactors = FALSE), # ASD vs TD, All Samples, Replication, 17 TD, 21 ASD
             DMRs_Null = read.csv("DMRs/Discovery/Null and Spike In/DMRs_Dx_Discovery50_null.csv", header = TRUE,
                                  stringsAsFactors = FALSE), # Null Exp vs Null Ctrl, All TD Samples, Discovery, 19 Ctrl, 20 Exp
             DMRs_NullSpikeOnly = read.csv("DMRs/Discovery/Null and Spike In/SpikedInDMRsOnly_DxDiscovery50_spikein.csv", header = TRUE,
                                       stringsAsFactors = FALSE), # 5000 DMRs spiked in only
             DMRs_NullSpike = read.csv("DMRs/Discovery/Null and Spike In/DMRs_Dx_Discovery50_spikein.csv", header = TRUE,
                                       stringsAsFactors = FALSE), # Null Exp vs Null Ctrl, 5000 DMRs spiked in, All TD Samples, Discovery, 19 Ctrl, 20 Exp
             DMRs_Permute = read.csv("DMRs/Discovery/Null and Spike In/DMRs_Dx_Discovery50_permuted_null.csv", header = TRUE,
                                     stringsAsFactors = FALSE), # Null Exp vs Null Ctrl, All Samples, Discovery, 38 Ctrl, 38 Exp
             DMRs_PermuteSpikeOnly = read.csv("DMRs/Discovery/Null and Spike In/SpikedInDMRsOnly_DxDiscovery50_permuted_spikein.csv", header = TRUE,
                                          stringsAsFactors = FALSE), # 5000 DMRs spiked in only
             DMRs_PermuteSpike = read.csv("DMRs/Discovery/Null and Spike In/DMRs_Dx_Discovery50_permuted_spikein.csv", header = TRUE,
                                          stringsAsFactors = FALSE)) # Null Exp vs Null Ctrl, 5000 DMRs spiked in, All Samples, Discovery, 38 Ctrl, 38 Exp

# Background
background <- list(background_Dx = read.csv("DMRs/Discovery/Diagnosis Males 50/bsseq_background_Discovery50_males.csv", header = TRUE, 
                                stringsAsFactors = FALSE), # ASD vs TD, All Samples, Discovery, 39 TD, 37 ASD
             background_DxRed = read.csv("DMRs/Discovery/Diagnosis Males 50/bsseq_background_Discovery50_reduced_males.csv", header = TRUE, 
                                   stringsAsFactors = FALSE), # ASD vs TD, Reduced Samples, Discovery, 17 TD, 21 ASD
             background_DxRep = read.csv("DMRs/Replication/Diagnosis Males 50/bsseq_background_Replication50_males.csv", header = TRUE,
                                   stringsAsFactors = FALSE), # ASD vs TD, All Samples, Replication, 17 TD, 21 ASD
             background_Null = read.csv("DMRs/Discovery/Null and Spike In/bsseq_background_Discovery50_null.csv", header = TRUE,
                                  stringsAsFactors = FALSE), # Null Exp vs Null Ctrl, All TD Samples, Discovery, 19 Ctrl, 20 Exp
             background_NullSpikeOnly = read.csv("DMRs/Discovery/Null and Spike In/bsseq_background_Discovery50_spikein.csv", header = TRUE,
                                             stringsAsFactors = FALSE), # Null Exp vs Null Ctrl, 5000 background spiked in, All TD Samples, Discovery, 19 Ctrl, 20 Exp
             background_NullSpike = read.csv("DMRs/Discovery/Null and Spike In/bsseq_background_Discovery50_spikein.csv", header = TRUE,
                                       stringsAsFactors = FALSE), # Null Exp vs Null Ctrl, 5000 background spiked in, All TD Samples, Discovery, 19 Ctrl, 20 Exp
             background_Permute = read.csv("DMRs/Discovery/Null and Spike In/bsseq_background_Discovery50_permuted_null.csv", header = TRUE,
                                     stringsAsFactors = FALSE), # Null Exp vs Null Ctrl, All Samples, Discovery, 38 Ctrl, 38 Exp
             background_PermuteSpikeOnly = read.csv("DMRs/Discovery/Null and Spike In/bsseq_background_Discovery50_permuted_spikein.csv", header = TRUE,
                                                stringsAsFactors = FALSE), # Null Exp vs Null Ctrl, 5000 background spiked in, All Samples, Discovery, 38 Ctrl, 38 Exp
             background_PermuteSpike = read.csv("DMRs/Discovery/Null and Spike In/bsseq_background_Discovery50_permuted_spikein.csv", header = TRUE,
                                          stringsAsFactors = FALSE)) # Null Exp vs Null Ctrl, 5000 background spiked in, All Samples, Discovery, 38 Ctrl, 38 Exp

# DMR Stats ####
# Filter out Spike In DMRs with delta < 0.05
DMRs$DMRs_NullSpikeOnly <- subset(DMRs$DMRs_NullSpikeOnly, delta < -0.05 | delta > 0.05)
DMRs$DMRs_PermuteSpikeOnly <- subset(DMRs$DMRs_PermuteSpikeOnly, delta < -0.05 | delta > 0.05)

# Run Stats
DMRstats <- getRegionStats(DMRs = DMRs, background = background, n = c(76, 38, 38, 39, 39, 39, 76, 76, 76))
write.csv(DMRstats, "Tables/Null and Spike-In DMR and Background Region Stats.csv", row.names = FALSE, quote = FALSE)

# Overlap DMRs ####
GR_DMRs <- sapply(DMRs, function(x) {GRanges(seqnames = x$seqnames, ranges = IRanges(start = x$start, end = x$end))})
GR_background <- sapply(background, function(x) {GRanges(seqnames = x$chr, ranges = IRanges(start = x$start, end = x$end))})

# Overlap Diagnosis All Males with Diagnosis Reduced Males DMRs
DMRoverlapVenn(list(GR_DMRs$DMRs_Dx[DMRs$DMRs_Dx$percentDifference > 0], 
                    GR_DMRs$DMRs_DxRed[DMRs$DMRs_DxRed$percentDifference > 0]), 
               NameOfPeaks = c("Dx_All_Hyper", "Dx_Reduced_Hyper"), 
               file = "Figures/DMR Overlap Diagnosis All Males with Diagnosis Reduced Males Hyper DMRs Venn.pdf",
               rotation.degree = 180, cat.pos = c(170, 180), cat.dist = c(0.03, 0.03))
DMRoverlapVenn(list(GR_DMRs$DMRs_Dx[DMRs$DMRs_Dx$percentDifference < 0], 
                    GR_DMRs$DMRs_DxRed[DMRs$DMRs_DxRed$percentDifference < 0]), 
               NameOfPeaks = c("Dx_All_Hypo", "Dx_Reduced_Hypo"), 
               file = "Figures/DMR Overlap Diagnosis All Males with Diagnosis Reduced Males Hypo DMRs Venn.pdf",
               rotation.degree = 180, cat.pos = c(170, 180), cat.dist = c(0.03, 0.03))
DMRoverlapVenn(list(GR_background$background_Dx, GR_background$background_DxRed), 
               NameOfPeaks = c("Dx_All_Background", "Dx_Reduced_Background"), 
               file = "Figures/DMR Overlap Diagnosis All Males with Diagnosis Reduced Males Background Venn.pdf",
               rotation.degree = 180, cat.pos = c(155, 205), cat.dist = c(0.05, 0.05))
GR_background$background_Dx[!GR_background$background_Dx %over% GR_background$background_DxRed] #2666 only in background_Dx

# Overlap Diagnosis Reduced Males with Null TD Males DMRs
DMRoverlapVenn(list(GR_DMRs$DMRs_DxRed[DMRs$DMRs_DxRed$percentDifference > 0], 
                    GR_DMRs$DMRs_Null[DMRs$DMRs_Null$percentDifference > 0]), 
               NameOfPeaks = c("Dx_Reduced_Hyper", "Null_Hyper"), 
               file = "Figures/DMR Overlap Diagnosis Reduced Males with Null TD Males Hyper DMRs Venn.pdf",
               rotation.degree = 180, cat.pos = c(170, 180), cat.dist = c(0.03, 0.03))
DMRoverlapVenn(list(GR_DMRs$DMRs_DxRed[DMRs$DMRs_DxRed$percentDifference < 0], 
                    GR_DMRs$DMRs_Null[DMRs$DMRs_Null$percentDifference < 0]), 
               NameOfPeaks = c("Dx_Reduced_Hypo", "Null_Hypo"), 
               file = "Figures/DMR Overlap Diagnosis Reduced Males with Null TD Males Hypo DMRs Venn.pdf",
               rotation.degree = 180, cat.pos = c(170, 180), cat.dist = c(0.03, 0.03))
DMRoverlapVenn(list(GR_background$background_DxRed, GR_background$background_Null), 
               NameOfPeaks = c("Dx_Reduced_Background", "Dx_Null_Background"), 
               file = "Figures/DMR Overlap Diagnosis Reduced Males with Null TD Males Background Venn.pdf",
               rotation.degree = 180, cat.pos = c(155, 205), cat.dist = c(0.05, 0.05))
GR_background$background_DxRed[!GR_background$background_DxRed %over% GR_background$background_Null] #1720 only in background_DxRed

# Overlap Diagnosis All Males with Permute TD ASD DMRs
DMRoverlapVenn(list(GR_DMRs$DMRs_Dx[DMRs$DMRs_Dx$percentDifference > 0], 
                    GR_DMRs$DMRs_Permute[DMRs$DMRs_Permute$percentDifference > 0]), 
               NameOfPeaks = c("Dx_All_Hyper", "Permute_Hyper"), 
               file = "Figures/DMR Overlap Diagnosis All Males with Permute TD ASD Males Hyper DMRs Venn.pdf",
               rotation.degree = 180, cat.pos = c(170, 180), cat.dist = c(0.03, 0.03))
DMRoverlapVenn(list(GR_DMRs$DMRs_Dx[DMRs$DMRs_Dx$percentDifference < 0], 
                    GR_DMRs$DMRs_Permute[DMRs$DMRs_Permute$percentDifference < 0]), 
               NameOfPeaks = c("Dx_All_Hypo", "Permute_Hypo"), 
               file = "Figures/DMR Overlap Diagnosis All Males with Permute TD ASD Males Hypo DMRs Venn.pdf",
               rotation.degree = 180, cat.pos = c(170, 180), cat.dist = c(0.03, 0.03))
DMRoverlapVenn(list(GR_background$background_Dx, GR_background$background_Permute), 
               NameOfPeaks = c("Dx_All_Background", "Permute_Background"), 
               file = "Figures/DMR Overlap Diagnosis All Males with Permute TD ASD Males Background Venn.pdf",
               rotation.degree = 180, cat.pos = c(155, 205), cat.dist = c(0.05, 0.05))
GR_background$background_Dx[!GR_background$background_Dx %over% GR_background$background_Permute] #690 only in background_Dx

# Overlap Null TD Males with Null TD Spike In Males DMRs
DMRoverlapVenn(list(GR_DMRs$DMRs_Null[DMRs$DMRs_Null$percentDifference > 0], 
                    GR_DMRs$DMRs_NullSpike[DMRs$DMRs_NullSpike$percentDifference > 0]), 
               NameOfPeaks = c("Null_Hyper", "Null_SpikeIn_Hyper"), 
               file = "Figures/DMR Overlap Null TD Males with Null Spike In TD Males Hyper DMRs Venn.pdf",
               rotation.degree = 180, cat.pos = c(190, 180), cat.dist = c(0.03, 0.03))
DMRoverlapVenn(list(GR_DMRs$DMRs_Null[DMRs$DMRs_Null$percentDifference < 0], 
                    GR_DMRs$DMRs_NullSpike[DMRs$DMRs_NullSpike$percentDifference < 0]), 
               NameOfPeaks = c("Null_Hypo", "Null_SpikeIn_Hypo"), 
               file = "Figures/DMR Overlap Null TD Males with Null Spike In TD Males Hypo DMRs Venn.pdf",
               rotation.degree = 180, cat.pos = c(195, 180), cat.dist = c(0.03, 0.03))
# Background is the same

# Overlap Null TD Spike In Males with Null TD Spike In Only Males DMRs
DMRoverlapVenn(list(GR_DMRs$DMRs_NullSpike[DMRs$DMRs_NullSpike$percentDifference > 0],
                    GR_DMRs$DMRs_NullSpikeOnly[DMRs$DMRs_NullSpikeOnly$delta > 0.05]), 
               NameOfPeaks = c("Null_SpikeIn_Hyper", "Null_SpikeIn_Only_Hyper"), 
               file = "Figures/DMR Overlap Null Spike In TD Males with Null Spike In Only TD Males Hyper DMRs Venn.pdf",
               rotation.degree = 180, cat.pos = c(190, 180), cat.dist = c(0.03, 0.03))
DMRoverlapVenn(list(GR_DMRs$DMRs_NullSpike[DMRs$DMRs_NullSpike$percentDifference < 0],
                    GR_DMRs$DMRs_NullSpikeOnly[DMRs$DMRs_NullSpikeOnly$delta < -0.05]), 
               NameOfPeaks = c("Null_SpikeIn_Hypo", "Null_SpikeIn_Only_Hypo"), 
               file = "Figures/DMR Overlap Null Spike In TD Males with Null Spike In Only TD Males Hypo DMRs Venn.pdf",
               rotation.degree = 180, cat.pos = c(190, 180), cat.dist = c(0.03, 0.03))

pdf("Figures/Null Spike In Only DMR Delta Histogram.pdf", height = 5, width = 5)
hist(DMRs$DMRs_NullSpikeOnly$delta, breaks = 100)
abline(v = 0.05)
abline(v = -0.05)
dev.off()

table(DMRs$DMRs_NullSpikeOnly$delta > 0.05) # 3633 DMRs
table(DMRs$DMRs_NullSpikeOnly$delta < -0.05) # 31 DMRs
table(DMRs$DMRs_NullSpikeOnly$delta <= 0.05 & DMRs$DMRs_NullSpikeOnly$delta >= -0.05) # 1336 DMRs

# Null Sensitivity and Specificity
NullSensitivity <- (sum(GR_DMRs$DMRs_NullSpikeOnly[DMRs$DMRs_NullSpikeOnly$delta > 0.05] %over%
                                GR_DMRs$DMRs_NullSpike[DMRs$DMRs_NullSpike$percentDifference > 0]) +
                            sum(GR_DMRs$DMRs_NullSpikeOnly[DMRs$DMRs_NullSpikeOnly$delta < -0.05] %over% 
                                        GR_DMRs$DMRs_NullSpike[DMRs$DMRs_NullSpike$percentDifference < 0])) * 100 / 
        sum(DMRs$DMRs_NullSpikeOnly$delta > 0.05 | DMRs$DMRs_NullSpikeOnly$delta < -0.05) # 969 / 3664 = 26.4%
NullTruePositives <- sum(GR_DMRs$DMRs_NullSpike[DMRs$DMRs_NullSpike$percentDifference > 0] %over%
                                  GR_DMRs$DMRs_NullSpikeOnly[DMRs$DMRs_NullSpikeOnly$delta > 0.05]) +
        sum(GR_DMRs$DMRs_NullSpike[DMRs$DMRs_NullSpike$percentDifference < 0] %over%
                    GR_DMRs$DMRs_NullSpikeOnly[DMRs$DMRs_NullSpikeOnly$delta < -0.05]) # 1055
NullFalsePositives <- sum(!GR_DMRs$DMRs_NullSpike[DMRs$DMRs_NullSpike$percentDifference > 0] %over%
                                  GR_DMRs$DMRs_NullSpikeOnly[DMRs$DMRs_NullSpikeOnly$delta > 0.05]) +
        sum(!GR_DMRs$DMRs_NullSpike[DMRs$DMRs_NullSpike$percentDifference < 0] %over%
                    GR_DMRs$DMRs_NullSpikeOnly[DMRs$DMRs_NullSpikeOnly$delta < -0.05]) # 3410
NullTrueNegatives <- sum(!GR_background$background_Null %over% GR_DMRs$DMRs_Null) # 198584
NullSpikeTrueNegatives <- sum((!GR_background$background_NullSpike %over% GR_DMRs$DMRs_NullSpike) &
                                      (!GR_background$background_NullSpike %over% GR_DMRs$DMRs_NullSpikeOnly)) # 195029
NullSpecificity <- NullTrueNegatives / length(GR_background$background_Null) # 0.9895
NullSpikeSpecificity <- NullSpikeTrueNegatives / sum(!GR_background$background_NullSpike %over% GR_DMRs$DMRs_NullSpikeOnly) # 0.9894

# Background is the same

# Overlap Permute TD ASD Males with Permute TD ASD Spike In Males DMRs
DMRoverlapVenn(list(GR_DMRs$DMRs_Permute[DMRs$DMRs_Permute$percentDifference > 0], GR_DMRs$DMRs_PermuteSpike[DMRs$DMRs_PermuteSpike$percentDifference > 0]), 
               NameOfPeaks = c("Permute_Hyper", "Permute_SpikeIn_Hyper"), 
               file = "Figures/DMR Overlap Permute TD ASD Males with Permute Spike In TD ASD Males Hyper DMRs Venn.pdf",
               rotation.degree = 180, cat.pos = c(200, 180), cat.dist = c(0.03, 0.03))
DMRoverlapVenn(list(GR_DMRs$DMRs_Permute[DMRs$DMRs_Permute$percentDifference < 0], GR_DMRs$DMRs_PermuteSpike[DMRs$DMRs_PermuteSpike$percentDifference < 0]), 
               NameOfPeaks = c("Permute_Hypo", "Permute_SpikeIn_Hypo"), 
               file = "Figures/DMR Overlap Permute TD ASD Males with Permute Spike In TD ASD Males Hypo DMRs Venn.pdf",
               rotation.degree = 180, cat.pos = c(200, 180), cat.dist = c(0.03, 0.03))
# Background is the same

# Overlap Permute TD ASD Spike In Males with Permute TD ASD Spike In Only Males DMRs
DMRoverlapVenn(list(GR_DMRs$DMRs_PermuteSpike[DMRs$DMRs_PermuteSpike$percentDifference > 0],
                    GR_DMRs$DMRs_PermuteSpikeOnly[DMRs$DMRs_PermuteSpikeOnly$delta > 0.05]), 
               NameOfPeaks = c("Permute_SpikeIn_Hyper", "Permute_SpikeIn_Only_Hyper"), 
               file = "Figures/DMR Overlap Permute Spike In TD ASD Males with Permute Spike In Only TD ASD Males Hyper DMRs Venn.pdf",
               rotation.degree = 180, cat.pos = c(190, 180), cat.dist = c(0.03, 0.03))
DMRoverlapVenn(list(GR_DMRs$DMRs_PermuteSpike[DMRs$DMRs_PermuteSpike$percentDifference < 0],
                    GR_DMRs$DMRs_PermuteSpikeOnly[DMRs$DMRs_PermuteSpikeOnly$delta < -0.05]), 
               NameOfPeaks = c("Permute_SpikeIn_Hypo", "Permute_SpikeIn_Only_Hypo"), 
               file = "Figures/DMR Overlap Permute Spike In TD ASD Males with Permute Spike In Only TD ASD Males Hypo DMRs Venn.pdf",
               rotation.degree = 180, cat.pos = c(180, 130), cat.dist = c(0.03, 0.05))
pdf("Figures/Permute Spike In Only DMR Delta Histogram.pdf", height = 5, width = 5)
hist(DMRs$DMRs_PermuteSpikeOnly$delta, breaks = 100)
abline(v = 0.05)
abline(v = -0.05)
dev.off()

table(DMRs$DMRs_PermuteSpikeOnly$delta > 0.05) # 3522 DMRs
table(DMRs$DMRs_PermuteSpikeOnly$delta < -0.05) # 36 DMRs
table(DMRs$DMRs_PermuteSpikeOnly$delta <= 0.05 & DMRs$DMRs_PermuteSpikeOnly$delta >= -0.05) # 1442 DMRs

# Background is the same

# Permute Sensitivity and Specificity
PermuteSensitivity <- (sum(GR_DMRs$DMRs_PermuteSpikeOnly[DMRs$DMRs_PermuteSpikeOnly$delta > 0.05] %over%
                                   GR_DMRs$DMRs_PermuteSpike[DMRs$DMRs_PermuteSpike$percentDifference > 0]) +
                            sum(GR_DMRs$DMRs_PermuteSpikeOnly[DMRs$DMRs_PermuteSpikeOnly$delta < -0.05] %over%
                                        GR_DMRs$DMRs_PermuteSpike[DMRs$DMRs_PermuteSpike$percentDifference < 0])) * 100 / 
        sum(DMRs$DMRs_PermuteSpikeOnly$delta > 0.05 | DMRs$DMRs_PermuteSpikeOnly$delta < -0.05) # 1136 / 3558 = 31.9%
PermuteTruePositives <- sum(GR_DMRs$DMRs_PermuteSpike[DMRs$DMRs_PermuteSpike$percentDifference > 0] %over%
                                 GR_DMRs$DMRs_PermuteSpikeOnly[DMRs$DMRs_PermuteSpikeOnly$delta > 0.05]) +
        sum(GR_DMRs$DMRs_PermuteSpike[DMRs$DMRs_PermuteSpike$percentDifference < 0] %over%
                    GR_DMRs$DMRs_PermuteSpikeOnly[DMRs$DMRs_PermuteSpikeOnly$delta < -0.05]) # 1214
PermuteFalsePositives <- sum(!GR_DMRs$DMRs_PermuteSpike[DMRs$DMRs_PermuteSpike$percentDifference > 0] %over%
                                  GR_DMRs$DMRs_PermuteSpikeOnly[DMRs$DMRs_PermuteSpikeOnly$delta > 0.05]) +
        sum(!GR_DMRs$DMRs_PermuteSpike[DMRs$DMRs_PermuteSpike$percentDifference < 0] %over%
                    GR_DMRs$DMRs_PermuteSpikeOnly[DMRs$DMRs_PermuteSpikeOnly$delta < -0.05]) # 1774
PermuteTrueNegatives <- sum(!GR_background$background_Permute %over% GR_DMRs$DMRs_Permute) # 191693
PermuteSpikeTrueNegatives <- sum((!GR_background$background_PermuteSpike %over% GR_DMRs$DMRs_PermuteSpike) &
                                         (!GR_background$background_PermuteSpike %over% GR_DMRs$DMRs_PermuteSpikeOnly)) #188192
PermuteSpecificity <- PermuteTrueNegatives / length(GR_background$background_Permute) # 0.9972
PermuteSpikeSpecificity <- PermuteSpikeTrueNegatives / sum(!GR_background$background_PermuteSpike %over% GR_DMRs$DMRs_PermuteSpikeOnly) # 0.9970


