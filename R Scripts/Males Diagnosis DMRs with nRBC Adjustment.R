# Males Diagnosis DMRs with nRBC Adjustment --------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 3/4/20

# Load Packages and Setup R Session ####
.libPaths("/share/lasallelab/programs/DMRichR/R_3.6")
sapply(c("tidyverse", "openxlsx", "bsseq", "dmrseq", "DMRichR", "regioneR", "LOLA"), require, character.only = TRUE)
set.seed(5)
register(MulticoreParam(1))
source("R Scripts/DMR Analysis Functions.R")

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

# Males Replication DMRs with nRBC Adjustment ------------------------
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

# Overlap with Diagnosis DMRs ----------------------------------------
# Load DMRs ####
chroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")
Dx_DMRs <- list(DiscMales = loadRegions("DMRs/Discovery/Diagnosis Males 50/DMRs_Dx_Discovery50_males.csv", chroms = chroms,
                                        DMRid = TRUE),
                RepMales = loadRegions("DMRs/Replication/Diagnosis Males 50/DMRs_Dx_Replication50_males.csv", chroms = chroms,
                                       DMRid = TRUE))
Dx_nRBC_DMRs <- list(DiscMales = loadRegions("DMRs/Discovery/Diagnosis and nRBCs Males 50/DMRs_Dx_nRBC_Discovery50_males.csv", 
                                             chroms = chroms, DMRid = TRUE),
                     RepMales = loadRegions("DMRs/Replication/Diagnosis and nRBCs Males 50/DMRs_Dx_nRBC_Replication50_males.csv",
                                            chroms = chroms, DMRid = TRUE))
background <- list(DiscMales = loadRegions("DMRs/Discovery/Diagnosis Males 50/bsseq_background_Discovery50_males.csv", 
                                           chroms = chroms, DMRid = FALSE),
                   RepMales = loadRegions("DMRs/Replication/Diagnosis Males 50/bsseq_background_Replication50_males.csv", 
                                          chroms = chroms, DMRid = FALSE))

# Get Region Stats ####
DMRs <- c(Dx_DMRs, lapply(Dx_DMRs, subset, percentDifference > 0), lapply(Dx_DMRs, subset, percentDifference < 0), Dx_nRBC_DMRs,
          lapply(Dx_nRBC_DMRs, subset, percentDifference > 0), lapply(Dx_nRBC_DMRs, subset, percentDifference < 0))
names(DMRs) <- paste(rep(c("Dx", "Dx_nRBC"), each = 6), rep(c("All", "Hyper", "Hypo"), each = 2), names(DMRs), sep = "_")
regionStats <- getRegionStats(DMRs = DMRs, background = background, n = c(76,38))
write.csv(regionStats, "Tables/Region Stats Diagnosis nRBC DMR Comparison.csv", row.names = FALSE)

# Convert to GRanges ####
GR_Dx_Hyper_DMRs <- sapply(Dx_DMRs, makeGRange, direction = "hyper")
GR_Dx_Hypo_DMRs <- sapply(Dx_DMRs, makeGRange, direction = "hypo")
GR_Dx_nRBC_Hyper_DMRs <- sapply(Dx_nRBC_DMRs, makeGRange, direction = "hyper")
GR_Dx_nRBC_Hypo_DMRs <- sapply(Dx_nRBC_DMRs, makeGRange, direction = "hypo")
GR_background <- sapply(background, makeGRange)

# Overlap ####
hyperDMRoverlaps <- mapply(function(x, y) sum(x %over% y), x = GR_Dx_Hyper_DMRs, y = GR_Dx_nRBC_Hyper_DMRs)
hypoDMRoverlaps <- mapply(function(x, y) sum(x %over% y), x = GR_Dx_Hypo_DMRs, y = GR_Dx_nRBC_Hypo_DMRs)

# Direction DiscMales RepMales 
#     Hyper       122      696 
#      Hypo       293     1619 

DMRoverlap <- data.frame(Comparison = paste(rep(c("Hyper", "Hypo"), each = 2), names(Dx_DMRs), sep = "_"), 
                         Dx_DMRs = sapply(c(GR_Dx_Hyper_DMRs, GR_Dx_Hypo_DMRs), length),
                         Dx_nRBC_DMRs = sapply(c(GR_Dx_nRBC_Hyper_DMRs, GR_Dx_nRBC_Hypo_DMRs), length),
                         DMR_Overlap = c(hyperDMRoverlaps, hypoDMRoverlaps),
                         DMR_PerOverlap = c(hyperDMRoverlaps, hypoDMRoverlaps) * 100 / 
                                 sapply(c(GR_Dx_Hyper_DMRs, GR_Dx_Hypo_DMRs), length))
write.csv(DMRoverlap, "Tables/DMR Overlap Diagnosis nRBC DMR Comparison.csv", row.names = FALSE)


# Overlap Stats with regioneR ----------------------------------------
# Setup ####
options("mc.cores" = 12) # sets 12 cores for parallelization
set.seed(5)
hg38_XY <- getGenomeAndMask("hg38", mask = NULL)$genome %>% filterChromosomes(keep.chr = chroms)

# Discovery Males Overlaps ####
stats_Disc_Males_Hyper <- DMRpermTest(A = redefineUserSets(GRangesList(GR_Dx_Hyper_DMRs$DiscMales), 
                                                           GR_background$DiscMales)[[1]], 
                                      B = redefineUserSets(GRangesList(GR_Dx_nRBC_Hyper_DMRs$DiscMales), 
                                                           GR_background$DiscMales)[[1]],
                                      genome = hg38_XY, universe = GR_background$DiscMales, Comparison = "Disc_Males_Hyper",
                                      file = "Figures/Hyper DMR Overlap Discovery Males Diagnosis vs nRBC RegioneR Plots.pdf")
stats_Disc_Males_Hypo <- DMRpermTest(A = redefineUserSets(GRangesList(GR_Dx_Hypo_DMRs$DiscMales), 
                                                           GR_background$DiscMales)[[1]], 
                                      B = redefineUserSets(GRangesList(GR_Dx_nRBC_Hypo_DMRs$DiscMales), 
                                                           GR_background$DiscMales)[[1]],
                                      genome = hg38_XY, universe = GR_background$DiscMales, Comparison = "Disc_Males_Hypo",
                                      file = "Figures/Hypo DMR Overlap Discovery Males Diagnosis vs nRBC RegioneR Plots.pdf")

# Replication Males Overlaps ####
stats_Rep_Males_Hyper <- DMRpermTest(A = redefineUserSets(GRangesList(GR_Dx_Hyper_DMRs$RepMales), 
                                                           GR_background$RepMales)[[1]], 
                                      B = redefineUserSets(GRangesList(GR_Dx_nRBC_Hyper_DMRs$RepMales), 
                                                           GR_background$RepMales)[[1]],
                                      genome = hg38_XY, universe = GR_background$RepMales, Comparison = "Rep_Males_Hyper",
                                      file = "Figures/Hyper DMR Overlap Replication Males Diagnosis vs nRBC RegioneR Plots.pdf")
stats_Rep_Males_Hypo <- DMRpermTest(A = redefineUserSets(GRangesList(GR_Dx_Hypo_DMRs$RepMales), 
                                                          GR_background$RepMales)[[1]], 
                                     B = redefineUserSets(GRangesList(GR_Dx_nRBC_Hypo_DMRs$RepMales), 
                                                          GR_background$RepMales)[[1]],
                                     genome = hg38_XY, universe = GR_background$RepMales, Comparison = "Rep_Males_Hypo",
                                     file = "Figures/Hypo DMR Overlap Replication Males Diagnosis vs nRBC RegioneR Plots.pdf")

# Combine ####
stats <- rbind(stats_Disc_Males_Hyper, stats_Disc_Males_Hypo, stats_Rep_Males_Hyper, stats_Rep_Males_Hypo)
write.csv(stats, "Tables/DMR Overlap Diagnosis vs nRBC RegioneR Stats.csv", row.names = FALSE)
