# LOLA Enrichment Diagnosis and Sex ChrX DMRs -----------------------------
# ASD Cord Blood Methylation
# Charles Mordaunt
# 7/28/19

# Packages ####
sapply(c("tidyverse", "LOLA", "simpleCache", "GenomicRanges", "qvalue", "RColorBrewer", "scales", "reshape2"), 
       require, character.only = TRUE)

# Functions ####
# Epigenerate
loadRegions <- function(file, chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), sort = TRUE){
        if(grepl("txt", file, fixed = TRUE)){
                regions <- read.delim(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
        }
        else{
                regions <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
        }
        if("seqnames" %in% colnames(regions)){
                colnames(regions)[colnames(regions) == "seqnames"] <- "chr"
        }
        regions <- subset(regions, chr %in% chroms)
        regions$chr <- factor(regions$chr, levels = chroms)
        if(sort){
                regions <- regions[order(regions$chr, regions$start),]
        }
        return(regions)
}
makeGRange <- function(DMRs, direction = c("all", "hyper", "hypo")){
        if(direction == "hyper"){DMRs <- subset(DMRs, percentDifference > 0)}
        if(direction == "hypo"){DMRs <- subset(DMRs, percentDifference < 0)}
        GR <- GRanges(seqnames = DMRs$chr, ranges = IRanges(start = DMRs$start, end = DMRs$end))
}

# Laptop
source("R Scripts/DMR Analysis Functions.R")

# Get ChrX LOLA Enrichments ----------------------------------------------------
# Data ####
regionDB <- loadRegionDB(dbLocation = "/share/lasallelab/programs/LOLA/hg38", useCache = TRUE, limit = NULL, 
                         collections = c("Roadmap_ChromHMM", "roadmap_epigenomics"))

# Discovery Diagnosis Males ChrX DMRs ####
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Discovery/")
DMRs <- loadRegions(file = "Dx_Males/DMRs_Dx_Discovery50_males.csv", chroms = "chrX", sort = TRUE)
DMRlist <- list("AllDMRs" = makeGRange(DMRs = DMRs, direction = "all"),
                "HyperDMRs" = makeGRange(DMRs = DMRs, direction = "hyper"),
                "HypoDMRs" = makeGRange(DMRs = DMRs, direction = "hypo"))
Background <- loadRegions(file = "Dx_Males/bsseq_background_Discovery50_males.csv", chroms = "chrX", sort = TRUE) %>% 
        makeGRange(direction = "all")
Results <- runLOLA(userSets = DMRlist, userUniverse = Background, regionDB = regionDB, cores = 2, redefineUserSets = TRUE)
writeCombinedEnrichment(combinedResults = Results, outFolder = "LOLA", includeSplits = FALSE)
file.copy(from = "LOLA/allEnrichments.tsv", to = "LOLA/LOLA_Dx_Discovery50_males_chrX_DMRs.tsv", overwrite = TRUE)
rm(DMRs, DMRlist, Background, Results)

# Discovery Diagnosis Females ChrX DMRs ####
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Discovery/")
DMRs <- loadRegions(file = "Dx_Females/DMRs_Dx_Discovery50_females.csv", chroms = "chrX", sort = TRUE)
DMRlist <- list("AllDMRs" = makeGRange(DMRs = DMRs, direction = "all"),
                "HyperDMRs" = makeGRange(DMRs = DMRs, direction = "hyper"),
                "HypoDMRs" = makeGRange(DMRs = DMRs, direction = "hypo"))
Background <- loadRegions(file = "Dx_Females/bsseq_background_Discovery50_females.csv", chroms = "chrX", sort = TRUE) %>% 
        makeGRange(direction = "all")
Results <- runLOLA(userSets = DMRlist, userUniverse = Background, regionDB = regionDB, cores = 2, redefineUserSets = TRUE)
writeCombinedEnrichment(combinedResults = Results, outFolder = "LOLA", includeSplits = FALSE)
file.copy(from = "LOLA/allEnrichments.tsv", to = "LOLA/LOLA_Dx_Discovery50_females_chrX_DMRs.tsv", overwrite = TRUE)
rm(DMRs, DMRlist, Background, Results)

# Replication Diagnosis Males ChrX DMRs ####
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Replication/")
DMRs <- loadRegions(file = "Dx_Males/DMRs_Dx_Replication50_males.csv", chroms = "chrX", sort = TRUE)
DMRlist <- list("AllDMRs" = makeGRange(DMRs = DMRs, direction = "all"),
                "HyperDMRs" = makeGRange(DMRs = DMRs, direction = "hyper"),
                "HypoDMRs" = makeGRange(DMRs = DMRs, direction = "hypo"))
Background <- loadRegions(file = "Dx_Males/bsseq_background_Replication50_males.csv", chroms = "chrX", sort = TRUE) %>% 
        makeGRange(direction = "all")
Results <- runLOLA(userSets = DMRlist, userUniverse = Background, regionDB = regionDB, cores = 2, redefineUserSets = TRUE)
writeCombinedEnrichment(combinedResults = Results, outFolder = "LOLA", includeSplits = FALSE)
file.copy(from = "LOLA/allEnrichments.tsv", to = "LOLA/LOLA_Dx_Replication50_males_chrX_DMRs.tsv", overwrite = TRUE)
rm(DMRs, DMRlist, Background, Results)

# Replication Diagnosis Females ChrX DMRs ####
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Replication/")
DMRs <- loadRegions(file = "Dx_Females_100/DMRs_Dx_Replication100_females.csv", chroms = "chrX", sort = TRUE)
DMRlist <- list("AllDMRs" = makeGRange(DMRs = DMRs, direction = "all"),
                "HyperDMRs" = makeGRange(DMRs = DMRs, direction = "hyper"),
                "HypoDMRs" = makeGRange(DMRs = DMRs, direction = "hypo"))
Background <- loadRegions(file = "Dx_Females_100/bsseq_background_Replication100_females.csv", chroms = "chrX", sort = TRUE) %>% 
        makeGRange(direction = "all")
Results <- runLOLA(userSets = DMRlist, userUniverse = Background, regionDB = regionDB, cores = 2, redefineUserSets = TRUE)
writeCombinedEnrichment(combinedResults = Results, outFolder = "LOLA", includeSplits = FALSE)
file.copy(from = "LOLA/allEnrichments.tsv", to = "LOLA/LOLA_Dx_Replication100_females_chrX_DMRs.tsv", overwrite = TRUE)
rm(DMRs, DMRlist, Background, Results)

# Analyze Discovery Diagnosis Males ChrX DMRs LOLA ------------------------------------
# Load Data ####
lola <- read.delim("Tables/LOLA_Dx_Discovery50_males_chrX_DMRs.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
index <- read.delim("Tables/LOLA Roadmap ChromHMM index.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Histone Mark Enrichment ####
# Prep Data
histone <- subset(lola, collection == "roadmap_epigenomics")
histone_hyper <- prepLOLAhistone(histone, index = index, regions = "HyperDMRs", 
                                 file = "Tables/LOLA Histone Dx Discovery 50 Males chrX HyperDMRs.csv")
histone_hypo <- prepLOLAhistone(histone, index = index, regions = "HypoDMRs", 
                                file = "Tables/LOLA Histone Dx Discovery 50 Males chrX HypoDMRs.csv")
histone_combined <- rbind(histone_hyper, histone_hypo)
histone_combined$userSet <- ifelse(histone_combined$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))

# Plot Odds Ratio
hm.max <- quantile(histone_combined$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_combined, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                width = 9, height = 7, legend.position = c(1.11, 0.855),
                file = "Figures/LOLA Histone Dx Discovery 50 Males chrX OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(histone_combined$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_combined, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                width = 9, height = 7, legend.position = c(1.12, 0.855),
                file = "Figures/LOLA Histone Dx Discovery 50 Males chrX log qvalue Heatmap.png")

# Plot Legend
plotLOLAhistone(histone = histone_hyper, type = "legend", file = "Figures/LOLA chrX DMRs Heatmap Legend.png")

# ChromHMM Enrichment ####
# Prep Data
chromHMM <- subset(lola, collection == "Roadmap_ChromHMM")
chromHMM_hyper <- prepLOLAchromHMM(chromHMM, index = index, regions = "HyperDMRs", 
                                   file = "Tables/LOLA chromHMM Dx Discovery 50 Males chrX HyperDMRs.csv")
chromHMM_hypo <- prepLOLAchromHMM(chromHMM, index = index, regions = "HypoDMRs", 
                                  file = "Tables/LOLA chromHMM Dx Discovery 50 Males chrX HypoDMRs.csv")
chromHMM_combined <- rbind(chromHMM_hyper, chromHMM_hypo)
chromHMM_combined$userSet <- ifelse(chromHMM_combined$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))

# Plot Odds Ratio
hm.max <- quantile(chromHMM_combined$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_combined, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                 width = 9, height = 7, legend.position = c(1.11, 0.855),
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Males chrX OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(chromHMM_combined$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_combined, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                 width = 9, height = 7, legend.position = c(1.12, 0.855),
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Males chrX log qvalue Heatmap.png")

# Analyze Discovery Diagnosis Females ChrX DMRs LOLA ------------------------------------
# Load Data ####
lola <- read.delim("Tables/LOLA_Dx_Discovery50_females_chrX_DMRs.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Histone Mark Enrichment ####
# Prep Data
histone <- subset(lola, collection == "roadmap_epigenomics")
histone_hyper <- prepLOLAhistone(histone, index = index, regions = "HyperDMRs", 
                                 file = "Tables/LOLA Histone Dx Discovery 50 Females chrX HyperDMRs.csv")
histone_hypo <- prepLOLAhistone(histone, index = index, regions = "HypoDMRs", 
                                file = "Tables/LOLA Histone Dx Discovery 50 Females chrX HypoDMRs.csv")
histone_combined <- rbind(histone_hyper, histone_hypo)
histone_combined$userSet <- ifelse(histone_combined$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))

# Plot Odds Ratio
hm.max <- quantile(histone_combined$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_combined, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                width = 9, height = 7, legend.position = c(1.11, 0.855),
                file = "Figures/LOLA Histone Dx Discovery 50 Females chrX OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(histone_combined$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_combined, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                width = 9, height = 7, legend.position = c(1.12, 0.855),
                file = "Figures/LOLA Histone Dx Discovery 50 Females chrX log qvalue Heatmap.png")

# ChromHMM Enrichment ####
# Prep Data
chromHMM <- subset(lola, collection == "Roadmap_ChromHMM")
chromHMM_hyper <- prepLOLAchromHMM(chromHMM, index = index, regions = "HyperDMRs", 
                                   file = "Tables/LOLA chromHMM Dx Discovery 50 Females chrX HyperDMRs.csv")
chromHMM_hypo <- prepLOLAchromHMM(chromHMM, index = index, regions = "HypoDMRs", 
                                  file = "Tables/LOLA chromHMM Dx Discovery 50 Females chrX HypoDMRs.csv")
chromHMM_combined <- rbind(chromHMM_hyper, chromHMM_hypo)
chromHMM_combined$userSet <- ifelse(chromHMM_combined$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))

# Plot Odds Ratio
hm.max <- quantile(chromHMM_combined$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_combined, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                 width = 9, height = 7, legend.position = c(1.11, 0.855),
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Females chrX OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(chromHMM_combined$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_combined, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                 width = 9, height = 7, legend.position = c(1.12, 0.855),
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Females chrX log qvalue Heatmap.png")

# Analyze Replication Diagnosis Males ChrX DMRs LOLA ------------------------------------
# Load Data ####
lola <- read.delim("Tables/LOLA_Dx_Replication50_males_chrX_DMRs.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Histone Mark Enrichment ####
# Prep Data
histone <- subset(lola, collection == "roadmap_epigenomics")
histone_hyper <- prepLOLAhistone(histone, index = index, regions = "HyperDMRs", 
                                 file = "Tables/LOLA Histone Dx Replication 50 Males chrX HyperDMRs.csv")
histone_hypo <- prepLOLAhistone(histone, index = index, regions = "HypoDMRs", 
                                file = "Tables/LOLA Histone Dx Replication 50 Males chrX HypoDMRs.csv")
histone_combined <- rbind(histone_hyper, histone_hypo)
histone_combined$userSet <- ifelse(histone_combined$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))

# Plot Odds Ratio
hm.max <- quantile(histone_combined$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_combined, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                width = 9, height = 7, legend.position = c(1.11, 0.855),
                file = "Figures/LOLA Histone Dx Replication 50 Males chrX OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(histone_combined$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_combined, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                width = 9, height = 7, legend.position = c(1.12, 0.855),
                file = "Figures/LOLA Histone Dx Replication 50 Males chrX log qvalue Heatmap.png")

# ChromHMM Enrichment ####
# Prep Data
chromHMM <- subset(lola, collection == "Roadmap_ChromHMM")
chromHMM_hyper <- prepLOLAchromHMM(chromHMM, index = index, regions = "HyperDMRs", 
                                   file = "Tables/LOLA chromHMM Dx Replication 50 Males chrX HyperDMRs.csv")
chromHMM_hypo <- prepLOLAchromHMM(chromHMM, index = index, regions = "HypoDMRs", 
                                  file = "Tables/LOLA chromHMM Dx Replication 50 Males chrX HypoDMRs.csv")
chromHMM_combined <- rbind(chromHMM_hyper, chromHMM_hypo)
chromHMM_combined$userSet <- ifelse(chromHMM_combined$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))

# Plot Odds Ratio
hm.max <- quantile(chromHMM_combined$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_combined, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                 width = 9, height = 7, legend.position = c(1.11, 0.855),
                 file = "Figures/LOLA chromHMM Dx Replication 50 Males chrX OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(chromHMM_combined$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_combined, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                 width = 9, height = 7, legend.position = c(1.12, 0.855),
                 file = "Figures/LOLA chromHMM Dx Replication 50 Males chrX log qvalue Heatmap.png")

# Analyze Replication Diagnosis Females ChrX DMRs LOLA ------------------------------------
# Load Data ####
lola <- read.delim("Tables/LOLA_Dx_Replication100_females_chrX_DMRs.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Histone Mark Enrichment ####
# Prep Data
histone <- subset(lola, collection == "roadmap_epigenomics")
histone_hyper <- prepLOLAhistone(histone, index = index, regions = "HyperDMRs", 
                                 file = "Tables/LOLA Histone Dx Replication 100 Females chrX HyperDMRs.csv")
histone_hypo <- prepLOLAhistone(histone, index = index, regions = "HypoDMRs", 
                                file = "Tables/LOLA Histone Dx Replication 100 Females chrX HypoDMRs.csv")
histone_combined <- rbind(histone_hyper, histone_hypo)
histone_combined$userSet <- ifelse(histone_combined$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))

# Plot Odds Ratio
hm.max <- quantile(histone_combined$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_combined, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                width = 9, height = 7, legend.position = c(1.11, 0.855),
                file = "Figures/LOLA Histone Dx Replication 100 Females chrX OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(histone_combined$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_combined, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                width = 9, height = 7, legend.position = c(1.12, 0.855),
                file = "Figures/LOLA Histone Dx Replication 100 Females chrX log qvalue Heatmap.png")

# ChromHMM Enrichment ####
# Prep Data
chromHMM <- subset(lola, collection == "Roadmap_ChromHMM")
chromHMM_hyper <- prepLOLAchromHMM(chromHMM, index = index, regions = "HyperDMRs", 
                                   file = "Tables/LOLA chromHMM Dx Replication 100 Females chrX HyperDMRs.csv")
chromHMM_hypo <- prepLOLAchromHMM(chromHMM, index = index, regions = "HypoDMRs", 
                                  file = "Tables/LOLA chromHMM Dx Replication 100 Females chrX HypoDMRs.csv")
chromHMM_combined <- rbind(chromHMM_hyper, chromHMM_hypo)
chromHMM_combined$userSet <- ifelse(chromHMM_combined$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))

# Plot Odds Ratio
hm.max <- quantile(chromHMM_combined$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_combined, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                 width = 9, height = 7, legend.position = c(1.11, 0.855),
                 file = "Figures/LOLA chromHMM Dx Replication 100 Females chrX OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(chromHMM_combined$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_combined, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                 width = 9, height = 7, legend.position = c(1.12, 0.855),
                 file = "Figures/LOLA chromHMM Dx Replication 100 Females chrX log qvalue Heatmap.png")

# Get Autosome LOLA Enrichments ----------------------------------------------------
# Data ####
regionDB <- loadRegionDB(dbLocation = "/share/lasallelab/programs/LOLA/hg38", useCache = TRUE, limit = NULL, 
                         collections = c("Roadmap_ChromHMM", "roadmap_epigenomics"))
chroms <- paste("chr", 1:22, sep = "")

# Discovery Diagnosis Males Autosome DMRs ####
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Discovery/")
DMRs <- loadRegions(file = "Dx_Males/DMRs_Dx_Discovery50_males.csv", chroms = chroms, sort = TRUE)
DMRlist <- list("AllDMRs" = makeGRange(DMRs = DMRs, direction = "all"),
                "HyperDMRs" = makeGRange(DMRs = DMRs, direction = "hyper"),
                "HypoDMRs" = makeGRange(DMRs = DMRs, direction = "hypo"))
Background <- loadRegions(file = "Dx_Males/bsseq_background_Discovery50_males.csv", chroms = chroms, sort = TRUE) %>% 
        makeGRange(direction = "all")
Results <- runLOLA(userSets = DMRlist, userUniverse = Background, regionDB = regionDB, cores = 2, redefineUserSets = TRUE)
writeCombinedEnrichment(combinedResults = Results, outFolder = "LOLA", includeSplits = FALSE)
file.copy(from = "LOLA/allEnrichments.tsv", to = "LOLA/LOLA_Dx_Discovery50_males_auto_DMRs.tsv", overwrite = TRUE)
rm(DMRs, DMRlist, Background, Results)

# Discovery Diagnosis Females Autosome DMRs ####
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Discovery/")
DMRs <- loadRegions(file = "Dx_Females/DMRs_Dx_Discovery50_females.csv", chroms = chroms, sort = TRUE)
DMRlist <- list("AllDMRs" = makeGRange(DMRs = DMRs, direction = "all"),
                "HyperDMRs" = makeGRange(DMRs = DMRs, direction = "hyper"),
                "HypoDMRs" = makeGRange(DMRs = DMRs, direction = "hypo"))
Background <- loadRegions(file = "Dx_Females/bsseq_background_Discovery50_females.csv", chroms = chroms, sort = TRUE) %>% 
        makeGRange(direction = "all")
Results <- runLOLA(userSets = DMRlist, userUniverse = Background, regionDB = regionDB, cores = 2, redefineUserSets = TRUE)
writeCombinedEnrichment(combinedResults = Results, outFolder = "LOLA", includeSplits = FALSE)
file.copy(from = "LOLA/allEnrichments.tsv", to = "LOLA/LOLA_Dx_Discovery50_females_auto_DMRs.tsv", overwrite = TRUE)
rm(DMRs, DMRlist, Background, Results)

# Replication Diagnosis Males Autosome DMRs ####
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Replication/")
DMRs <- loadRegions(file = "Dx_Males/DMRs_Dx_Replication50_males.csv", chroms = chroms, sort = TRUE)
DMRlist <- list("AllDMRs" = makeGRange(DMRs = DMRs, direction = "all"),
                "HyperDMRs" = makeGRange(DMRs = DMRs, direction = "hyper"),
                "HypoDMRs" = makeGRange(DMRs = DMRs, direction = "hypo"))
Background <- loadRegions(file = "Dx_Males/bsseq_background_Replication50_males.csv", chroms = chroms, sort = TRUE) %>% 
        makeGRange(direction = "all")
Results <- runLOLA(userSets = DMRlist, userUniverse = Background, regionDB = regionDB, cores = 2, redefineUserSets = TRUE)
writeCombinedEnrichment(combinedResults = Results, outFolder = "LOLA", includeSplits = FALSE)
file.copy(from = "LOLA/allEnrichments.tsv", to = "LOLA/LOLA_Dx_Replication50_males_auto_DMRs.tsv", overwrite = TRUE)
rm(DMRs, DMRlist, Background, Results)

# Replication Diagnosis Females Autosome DMRs ####
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Replication/")
DMRs <- loadRegions(file = "Dx_Females_100/DMRs_Dx_Replication100_females.csv", chroms = chroms, sort = TRUE)
DMRlist <- list("AllDMRs" = makeGRange(DMRs = DMRs, direction = "all"),
                "HyperDMRs" = makeGRange(DMRs = DMRs, direction = "hyper"),
                "HypoDMRs" = makeGRange(DMRs = DMRs, direction = "hypo"))
Background <- loadRegions(file = "Dx_Females_100/bsseq_background_Replication100_females.csv", chroms = chroms, sort = TRUE) %>% 
        makeGRange(direction = "all")
Results <- runLOLA(userSets = DMRlist, userUniverse = Background, regionDB = regionDB, cores = 2, redefineUserSets = TRUE)
writeCombinedEnrichment(combinedResults = Results, outFolder = "LOLA", includeSplits = FALSE)
file.copy(from = "LOLA/allEnrichments.tsv", to = "LOLA/LOLA_Dx_Replication100_females_auto_DMRs.tsv", overwrite = TRUE)
rm(DMRs, DMRlist, Background, Results)

# Analyze Discovery Diagnosis Males Autosome DMRs LOLA ------------------------------------
# Load Data ####
lola <- read.delim("Tables/LOLA_Dx_Discovery50_males_auto_DMRs.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
index <- read.delim("Tables/LOLA Roadmap ChromHMM index.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Histone Mark Enrichment ####
# Prep Data
histone <- subset(lola, collection == "roadmap_epigenomics")
histone_hyper <- prepLOLAhistone(histone, index = index, regions = "HyperDMRs", 
                                 file = "Tables/LOLA Histone Dx Discovery 50 Males autosome HyperDMRs.csv")
histone_hypo <- prepLOLAhistone(histone, index = index, regions = "HypoDMRs", 
                                file = "Tables/LOLA Histone Dx Discovery 50 Males autosome HypoDMRs.csv")
histone_combined <- rbind(histone_hyper, histone_hypo)
histone_combined$userSet <- ifelse(histone_combined$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))

# Plot Odds Ratio
hm.max <- quantile(histone_combined$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_combined, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                width = 9, height = 7, legend.position = c(1.11, 0.855),
                file = "Figures/LOLA Histone Dx Discovery 50 Males autosome OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(histone_combined$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_combined, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                width = 9, height = 7, legend.position = c(1.12, 0.855),
                file = "Figures/LOLA Histone Dx Discovery 50 Males autosome log qvalue Heatmap.png")

# ChromHMM Enrichment ####
# Prep Data
chromHMM <- subset(lola, collection == "Roadmap_ChromHMM")
chromHMM_hyper <- prepLOLAchromHMM(chromHMM, index = index, regions = "HyperDMRs", 
                                   file = "Tables/LOLA chromHMM Dx Discovery 50 Males autosome HyperDMRs.csv")
chromHMM_hypo <- prepLOLAchromHMM(chromHMM, index = index, regions = "HypoDMRs", 
                                  file = "Tables/LOLA chromHMM Dx Discovery 50 Males autosome HypoDMRs.csv")
chromHMM_combined <- rbind(chromHMM_hyper, chromHMM_hypo)
chromHMM_combined$userSet <- ifelse(chromHMM_combined$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))

# Plot Odds Ratio
hm.max <- quantile(chromHMM_combined$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_combined, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                 width = 9, height = 7, legend.position = c(1.11, 0.855),
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Males autosome OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(chromHMM_combined$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_combined, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                 width = 9, height = 7, legend.position = c(1.12, 0.855),
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Males autosome log qvalue Heatmap.png")

# Analyze Discovery Diagnosis Females Autosome DMRs LOLA ------------------------------------
# Load Data ####
lola <- read.delim("Tables/LOLA_Dx_Discovery50_females_auto_DMRs.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Histone Mark Enrichment ####
# Prep Data
histone <- subset(lola, collection == "roadmap_epigenomics")
histone_hyper <- prepLOLAhistone(histone, index = index, regions = "HyperDMRs", 
                                 file = "Tables/LOLA Histone Dx Discovery 50 Females autosome HyperDMRs.csv")
histone_hypo <- prepLOLAhistone(histone, index = index, regions = "HypoDMRs", 
                                file = "Tables/LOLA Histone Dx Discovery 50 Females autosome HypoDMRs.csv")
histone_combined <- rbind(histone_hyper, histone_hypo)
histone_combined$userSet <- ifelse(histone_combined$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))

# Plot Odds Ratio
hm.max <- quantile(histone_combined$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_combined, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                width = 9, height = 7, legend.position = c(1.11, 0.855),
                file = "Figures/LOLA Histone Dx Discovery 50 Females autosome OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(histone_combined$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_combined, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                width = 9, height = 7, legend.position = c(1.12, 0.855),
                file = "Figures/LOLA Histone Dx Discovery 50 Females autosome log qvalue Heatmap.png")

# ChromHMM Enrichment ####
# Prep Data
chromHMM <- subset(lola, collection == "Roadmap_ChromHMM")
chromHMM_hyper <- prepLOLAchromHMM(chromHMM, index = index, regions = "HyperDMRs", 
                                   file = "Tables/LOLA chromHMM Dx Discovery 50 Females autosome HyperDMRs.csv")
chromHMM_hypo <- prepLOLAchromHMM(chromHMM, index = index, regions = "HypoDMRs", 
                                  file = "Tables/LOLA chromHMM Dx Discovery 50 Females autosome HypoDMRs.csv")
chromHMM_combined <- rbind(chromHMM_hyper, chromHMM_hypo)
chromHMM_combined$userSet <- ifelse(chromHMM_combined$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))

# Plot Odds Ratio
hm.max <- quantile(chromHMM_combined$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_combined, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                 width = 9, height = 7, legend.position = c(1.11, 0.855),
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Females autosome OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(chromHMM_combined$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_combined, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                 width = 9, height = 7, legend.position = c(1.12, 0.855),
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Females autosome log qvalue Heatmap.png")

# Analyze Replication Diagnosis Males Autosome DMRs LOLA ------------------------------------
# Load Data ####
lola <- read.delim("Tables/LOLA_Dx_Replication50_males_auto_DMRs.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Histone Mark Enrichment ####
# Prep Data
histone <- subset(lola, collection == "roadmap_epigenomics")
histone_hyper <- prepLOLAhistone(histone, index = index, regions = "HyperDMRs", 
                                 file = "Tables/LOLA Histone Dx Replication 50 Males autosome HyperDMRs.csv")
histone_hypo <- prepLOLAhistone(histone, index = index, regions = "HypoDMRs", 
                                file = "Tables/LOLA Histone Dx Replication 50 Males autosome HypoDMRs.csv")
histone_combined <- rbind(histone_hyper, histone_hypo)
histone_combined$userSet <- ifelse(histone_combined$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))

# Plot Odds Ratio
hm.max <- quantile(histone_combined$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_combined, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                width = 9, height = 7, legend.position = c(1.11, 0.855),
                file = "Figures/LOLA Histone Dx Replication 50 Males autosome OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(histone_combined$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_combined, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                width = 9, height = 7, legend.position = c(1.12, 0.855),
                file = "Figures/LOLA Histone Dx Replication 50 Males autosome log qvalue Heatmap.png")

# ChromHMM Enrichment ####
# Prep Data
chromHMM <- subset(lola, collection == "Roadmap_ChromHMM")
chromHMM_hyper <- prepLOLAchromHMM(chromHMM, index = index, regions = "HyperDMRs", 
                                   file = "Tables/LOLA chromHMM Dx Replication 50 Males autosome HyperDMRs.csv")
chromHMM_hypo <- prepLOLAchromHMM(chromHMM, index = index, regions = "HypoDMRs", 
                                  file = "Tables/LOLA chromHMM Dx Replication 50 Males autosome HypoDMRs.csv")
chromHMM_combined <- rbind(chromHMM_hyper, chromHMM_hypo)
chromHMM_combined$userSet <- ifelse(chromHMM_combined$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))

# Plot Odds Ratio
hm.max <- quantile(chromHMM_combined$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_combined, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                 width = 9, height = 7, legend.position = c(1.11, 0.855),
                 file = "Figures/LOLA chromHMM Dx Replication 50 Males autosome OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(chromHMM_combined$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_combined, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                 width = 9, height = 7, legend.position = c(1.12, 0.855),
                 file = "Figures/LOLA chromHMM Dx Replication 50 Males autosome log qvalue Heatmap.png")

# Analyze Replication Diagnosis Females Autosome DMRs LOLA ------------------------------------
# Load Data ####
lola <- read.delim("Tables/LOLA_Dx_Replication100_females_auto_DMRs.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Histone Mark Enrichment ####
# Prep Data
histone <- subset(lola, collection == "roadmap_epigenomics")
histone_hyper <- prepLOLAhistone(histone, index = index, regions = "HyperDMRs", 
                                 file = "Tables/LOLA Histone Dx Replication 100 Females autosome HyperDMRs.csv")
histone_hypo <- prepLOLAhistone(histone, index = index, regions = "HypoDMRs", 
                                file = "Tables/LOLA Histone Dx Replication 100 Females autosome HypoDMRs.csv")
histone_combined <- rbind(histone_hyper, histone_hypo)
histone_combined$userSet <- ifelse(histone_combined$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))

# Plot Odds Ratio
hm.max <- quantile(histone_combined$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_combined, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                width = 9, height = 7, legend.position = c(1.11, 0.855),
                file = "Figures/LOLA Histone Dx Replication 100 Females autosome OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(histone_combined$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_combined, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                width = 9, height = 7, legend.position = c(1.12, 0.855),
                file = "Figures/LOLA Histone Dx Replication 100 Females autosome log qvalue Heatmap.png")

# ChromHMM Enrichment ####
# Prep Data
chromHMM <- subset(lola, collection == "Roadmap_ChromHMM")
chromHMM_hyper <- prepLOLAchromHMM(chromHMM, index = index, regions = "HyperDMRs", 
                                   file = "Tables/LOLA chromHMM Dx Replication 100 Females autosome HyperDMRs.csv")
chromHMM_hypo <- prepLOLAchromHMM(chromHMM, index = index, regions = "HypoDMRs", 
                                  file = "Tables/LOLA chromHMM Dx Replication 100 Females autosome HypoDMRs.csv")
chromHMM_combined <- rbind(chromHMM_hyper, chromHMM_hypo)
chromHMM_combined$userSet <- ifelse(chromHMM_combined$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))

# Plot Odds Ratio
hm.max <- quantile(chromHMM_combined$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_combined, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                 width = 9, height = 7, legend.position = c(1.11, 0.855),
                 file = "Figures/LOLA chromHMM Dx Replication 100 Females autosome OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(chromHMM_combined$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_combined, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                 width = 9, height = 7, legend.position = c(1.12, 0.855),
                 file = "Figures/LOLA chromHMM Dx Replication 100 Females autosome log qvalue Heatmap.png")
rm(chromHMM, chromHMM_hyper, chromHMM_hypo, chromHMM_combined, histone, histone_hyper, histone_hypo, histone_combined, lola, hm.max)

# Compare Discovery Males Autosome and ChrX LOLA ------------------------------
# Histone Mark Enrichment ####
histone_auto_hyper <- read.csv("Tables/LOLA Histone Dx Discovery 50 Males autosome HyperDMRs.csv", header = TRUE, 
                               stringsAsFactors = FALSE)
histone_auto_hypo <- read.csv("Tables/LOLA Histone Dx Discovery 50 Males autosome HypoDMRs.csv", header = TRUE, 
                              stringsAsFactors = FALSE)
histone_chrX_hyper <- read.csv("Tables/LOLA Histone Dx Discovery 50 Males chrX HyperDMRs.csv", header = TRUE, 
                               stringsAsFactors = FALSE)
histone_chrX_hypo <- read.csv("Tables/LOLA Histone Dx Discovery 50 Males chrX HypoDMRs.csv", header = TRUE, 
                              stringsAsFactors = FALSE)
histone <- rbind(histone_auto_hyper, histone_auto_hypo, histone_chrX_hyper, histone_chrX_hypo)
histone$chroms <- c(rep("Autosomes", nrow(histone_auto_hyper) + nrow(histone_auto_hypo)),
                    rep("ChrX", nrow(histone_chrX_hyper) + nrow(histone_chrX_hypo))) %>%
        factor(levels = c("Autosomes", "ChrX"))
histone$userSet <- ifelse(histone$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))
histone$antibody <- factor(histone$antibody, levels = c("H3K4me1", "H3K4me3", "H3K9me3", "H3K27me3", "H3K36me3"))

plotLOLAhistoneBox(histone, file = "Figures/LOLA Histone Discovery Males Autosome ChrX OR Mean CL.png", facet = vars(userSet))

# ChromHMM Enrichment ####
chromHMM_auto_hyper <- read.csv("Tables/LOLA ChromHMM Dx Discovery 50 Males autosome HyperDMRs.csv", header = TRUE, 
                               stringsAsFactors = FALSE)
chromHMM_auto_hypo <- read.csv("Tables/LOLA ChromHMM Dx Discovery 50 Males autosome HypoDMRs.csv", header = TRUE, 
                              stringsAsFactors = FALSE)
chromHMM_chrX_hyper <- read.csv("Tables/LOLA ChromHMM Dx Discovery 50 Males chrX HyperDMRs.csv", header = TRUE, 
                               stringsAsFactors = FALSE)
chromHMM_chrX_hypo <- read.csv("Tables/LOLA ChromHMM Dx Discovery 50 Males chrX HypoDMRs.csv", header = TRUE, 
                              stringsAsFactors = FALSE)
chromHMM <- rbind(chromHMM_auto_hyper, chromHMM_auto_hypo, chromHMM_chrX_hyper, chromHMM_chrX_hypo)
chromHMM$chroms <- c(rep("Autosomes", nrow(chromHMM_auto_hyper) + nrow(chromHMM_auto_hypo)),
                    rep("ChrX", nrow(chromHMM_chrX_hyper) + nrow(chromHMM_chrX_hypo))) %>%
        factor(levels = c("Autosomes", "ChrX"))
chromHMM$userSet <- ifelse(chromHMM$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))
chromHMM$chromState <- factor(chromHMM$chromState, levels = unique(chromHMM$chromState))

plotLOLAchromHMMBox(chromHMM, file = "Figures/LOLA ChromHMM Discovery Males Autosome ChrX OR Mean CL.png", facet = vars(userSet))

# Differential Enrichment Stats ####
histone$chroms <- relevel(histone$chroms, ref = "ChrX")
histone_stats <- split(histone, list(histone$userSet, histone$antibody)) %>% 
        sapply(function(x){
                t.test(oddsRatio ~ chroms, data = x, paired = TRUE)[c("estimate", "conf.int", "statistic", "p.value")] %>% 
                        unlist()
        }) %>% t() %>% as.data.frame()
colnames(histone_stats) <- c("Estimate", "ConfIntLow", "ConfIntHigh", "tstat", "pvalue")
histone_stats$Direction <- strsplit(rownames(histone_stats), split = ".", fixed = TRUE) %>% sapply(function(x) x[1])
histone_stats$Antibody <- strsplit(rownames(histone_stats), split = ".", fixed = TRUE) %>% sapply(function(x) x[2])
rownames(histone_stats) <- 1:nrow(histone_stats)
histone_stats$qvalue <- p.adjust(histone_stats$pvalue, method = "fdr")

chromHMM$chroms <- relevel(chromHMM$chroms, ref = "ChrX")
chromHMM_stats <- split(chromHMM, list(chromHMM$userSet, chromHMM$chromState)) %>% 
        sapply(function(x){
                t.test(oddsRatio ~ chroms, data = x, paired = TRUE)[c("estimate", "conf.int", "statistic", "p.value")] %>% 
                        unlist()
        }) %>% t() %>% as.data.frame()
colnames(chromHMM_stats) <- c("Estimate", "ConfIntLow", "ConfIntHigh", "tstat", "pvalue")
chromHMM_stats$Direction <- strsplit(rownames(chromHMM_stats), split = ".", fixed = TRUE) %>% sapply(function(x) x[1])
chromHMM_stats$Antibody <- strsplit(rownames(chromHMM_stats), split = ".", fixed = TRUE) %>% sapply(function(x) x[2])
rownames(chromHMM_stats) <- 1:nrow(chromHMM_stats)
chromHMM_stats$qvalue <- p.adjust(chromHMM_stats$pvalue, method = "fdr")

stats <- rbind(histone_stats, chromHMM_stats)
stats <- stats[,c("Direction", "Antibody", "Estimate", "ConfIntLow", "ConfIntHigh", "tstat", "pvalue", "qvalue")]
write.csv(stats, "Tables/LOLA Discovery Males Autosome vs ChrX OR Stats.csv", quote = FALSE, row.names = FALSE)

# Compare Discovery Females Autosome and ChrX LOLA ------------------------------
# Histone Mark Enrichment ####
histone_auto_hyper <- read.csv("Tables/LOLA Histone Dx Discovery 50 Females autosome HyperDMRs.csv", header = TRUE, 
                               stringsAsFactors = FALSE)
histone_auto_hypo <- read.csv("Tables/LOLA Histone Dx Discovery 50 Females autosome HypoDMRs.csv", header = TRUE, 
                              stringsAsFactors = FALSE)
histone_chrX_hyper <- read.csv("Tables/LOLA Histone Dx Discovery 50 Females chrX HyperDMRs.csv", header = TRUE, 
                               stringsAsFactors = FALSE)
histone_chrX_hypo <- read.csv("Tables/LOLA Histone Dx Discovery 50 Females chrX HypoDMRs.csv", header = TRUE, 
                              stringsAsFactors = FALSE)
histone <- rbind(histone_auto_hyper, histone_auto_hypo, histone_chrX_hyper, histone_chrX_hypo)
histone$chroms <- c(rep("Autosomes", nrow(histone_auto_hyper) + nrow(histone_auto_hypo)),
                    rep("ChrX", nrow(histone_chrX_hyper) + nrow(histone_chrX_hypo))) %>%
        factor(levels = c("Autosomes", "ChrX"))
histone$userSet <- ifelse(histone$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))
histone$antibody <- factor(histone$antibody, levels = c("H3K4me1", "H3K4me3", "H3K9me3", "H3K27me3", "H3K36me3"))

plotLOLAhistoneBox(histone, file = "Figures/LOLA Histone Discovery Females Autosome ChrX OR Mean CL.png", facet = vars(userSet),
                   plot.margin = c(1.5, 0.5, 0.5, 0.5), ylim = c(0,13))

# ChromHMM Enrichment ####
chromHMM_auto_hyper <- read.csv("Tables/LOLA ChromHMM Dx Discovery 50 Females autosome HyperDMRs.csv", header = TRUE, 
                                stringsAsFactors = FALSE)
chromHMM_auto_hypo <- read.csv("Tables/LOLA ChromHMM Dx Discovery 50 Females autosome HypoDMRs.csv", header = TRUE, 
                               stringsAsFactors = FALSE)
chromHMM_chrX_hyper <- read.csv("Tables/LOLA ChromHMM Dx Discovery 50 Females chrX HyperDMRs.csv", header = TRUE, 
                                stringsAsFactors = FALSE)
chromHMM_chrX_hypo <- read.csv("Tables/LOLA ChromHMM Dx Discovery 50 Females chrX HypoDMRs.csv", header = TRUE, 
                               stringsAsFactors = FALSE)
chromHMM <- rbind(chromHMM_auto_hyper, chromHMM_auto_hypo, chromHMM_chrX_hyper, chromHMM_chrX_hypo)
chromHMM$chroms <- c(rep("Autosomes", nrow(chromHMM_auto_hyper) + nrow(chromHMM_auto_hypo)),
                     rep("ChrX", nrow(chromHMM_chrX_hyper) + nrow(chromHMM_chrX_hypo))) %>%
        factor(levels = c("Autosomes", "ChrX"))
chromHMM$userSet <- ifelse(chromHMM$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))
chromHMM$chromState <- factor(chromHMM$chromState, levels = unique(chromHMM$chromState))

plotLOLAchromHMMBox(chromHMM, file = "Figures/LOLA ChromHMM Discovery Females Autosome ChrX OR Mean CL.png", 
                    facet = vars(userSet))

# Differential Enrichment Stats ####
histone$chroms <- relevel(histone$chroms, ref = "ChrX")
histone_stats <- split(histone, list(histone$userSet, histone$antibody)) %>% 
        sapply(function(x){
                t.test(oddsRatio ~ chroms, data = x, paired = TRUE)[c("estimate", "conf.int", "statistic", "p.value")] %>% 
                        unlist()
        }) %>% t() %>% as.data.frame()
colnames(histone_stats) <- c("Estimate", "ConfIntLow", "ConfIntHigh", "tstat", "pvalue")
histone_stats$Direction <- strsplit(rownames(histone_stats), split = ".", fixed = TRUE) %>% sapply(function(x) x[1])
histone_stats$Antibody <- strsplit(rownames(histone_stats), split = ".", fixed = TRUE) %>% sapply(function(x) x[2])
rownames(histone_stats) <- 1:nrow(histone_stats)
histone_stats$qvalue <- p.adjust(histone_stats$pvalue, method = "fdr")

chromHMM$chroms <- relevel(chromHMM$chroms, ref = "ChrX")
chromHMM_stats <- split(chromHMM, list(chromHMM$userSet, chromHMM$chromState)) %>% 
        sapply(function(x){
                t.test(oddsRatio ~ chroms, data = x, paired = TRUE)[c("estimate", "conf.int", "statistic", "p.value")] %>% 
                        unlist()
        }) %>% t() %>% as.data.frame()
colnames(chromHMM_stats) <- c("Estimate", "ConfIntLow", "ConfIntHigh", "tstat", "pvalue")
chromHMM_stats$Direction <- strsplit(rownames(chromHMM_stats), split = ".", fixed = TRUE) %>% sapply(function(x) x[1])
chromHMM_stats$Antibody <- strsplit(rownames(chromHMM_stats), split = ".", fixed = TRUE) %>% sapply(function(x) x[2])
rownames(chromHMM_stats) <- 1:nrow(chromHMM_stats)
chromHMM_stats$qvalue <- p.adjust(chromHMM_stats$pvalue, method = "fdr")

stats <- rbind(histone_stats, chromHMM_stats)
stats <- stats[,c("Direction", "Antibody", "Estimate", "ConfIntLow", "ConfIntHigh", "tstat", "pvalue", "qvalue")]
write.csv(stats, "Tables/LOLA Discovery Females Autosome vs ChrX OR Stats.csv", quote = FALSE, row.names = FALSE)

# Compare Replication Males Autosome and ChrX LOLA ------------------------------
# Histone Mark Enrichment ####
histone_auto_hyper <- read.csv("Tables/LOLA Histone Dx Replication 50 Males autosome HyperDMRs.csv", header = TRUE, 
                               stringsAsFactors = FALSE)
histone_auto_hypo <- read.csv("Tables/LOLA Histone Dx Replication 50 Males autosome HypoDMRs.csv", header = TRUE, 
                              stringsAsFactors = FALSE)
histone_chrX_hyper <- read.csv("Tables/LOLA Histone Dx Replication 50 Males chrX HyperDMRs.csv", header = TRUE, 
                               stringsAsFactors = FALSE)
histone_chrX_hypo <- read.csv("Tables/LOLA Histone Dx Replication 50 Males chrX HypoDMRs.csv", header = TRUE, 
                              stringsAsFactors = FALSE)
histone <- rbind(histone_auto_hyper, histone_auto_hypo, histone_chrX_hyper, histone_chrX_hypo)
histone$chroms <- c(rep("Autosomes", nrow(histone_auto_hyper) + nrow(histone_auto_hypo)),
                    rep("ChrX", nrow(histone_chrX_hyper) + nrow(histone_chrX_hypo))) %>%
        factor(levels = c("Autosomes", "ChrX"))
histone$userSet <- ifelse(histone$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))
histone$antibody <- factor(histone$antibody, levels = c("H3K4me1", "H3K4me3", "H3K9me3", "H3K27me3", "H3K36me3"))

plotLOLAhistoneBox(histone, file = "Figures/LOLA Histone Replication Males Autosome ChrX OR Mean CL.png", facet = vars(userSet),
                   ylim = c(0, 7.5))

# ChromHMM Enrichment ####
chromHMM_auto_hyper <- read.csv("Tables/LOLA ChromHMM Dx Replication 50 Males autosome HyperDMRs.csv", header = TRUE, 
                                stringsAsFactors = FALSE)
chromHMM_auto_hypo <- read.csv("Tables/LOLA ChromHMM Dx Replication 50 Males autosome HypoDMRs.csv", header = TRUE, 
                               stringsAsFactors = FALSE)
chromHMM_chrX_hyper <- read.csv("Tables/LOLA ChromHMM Dx Replication 50 Males chrX HyperDMRs.csv", header = TRUE, 
                                stringsAsFactors = FALSE)
chromHMM_chrX_hypo <- read.csv("Tables/LOLA ChromHMM Dx Replication 50 Males chrX HypoDMRs.csv", header = TRUE, 
                               stringsAsFactors = FALSE)
chromHMM <- rbind(chromHMM_auto_hyper, chromHMM_auto_hypo, chromHMM_chrX_hyper, chromHMM_chrX_hypo)
chromHMM$chroms <- c(rep("Autosomes", nrow(chromHMM_auto_hyper) + nrow(chromHMM_auto_hypo)),
                     rep("ChrX", nrow(chromHMM_chrX_hyper) + nrow(chromHMM_chrX_hypo))) %>%
        factor(levels = c("Autosomes", "ChrX"))
chromHMM$userSet <- ifelse(chromHMM$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))
chromHMM$chromState <- factor(chromHMM$chromState, levels = unique(chromHMM$chromState))

plotLOLAchromHMMBox(chromHMM, file = "Figures/LOLA ChromHMM Replication Males Autosome ChrX OR Mean CL.png", 
                    facet = vars(userSet))

# Differential Enrichment Stats ####
histone$chroms <- relevel(histone$chroms, ref = "ChrX")
histone_stats <- split(histone, list(histone$userSet, histone$antibody)) %>% 
        sapply(function(x){
                t.test(oddsRatio ~ chroms, data = x, paired = TRUE)[c("estimate", "conf.int", "statistic", "p.value")] %>% 
                        unlist()
        }) %>% t() %>% as.data.frame()
colnames(histone_stats) <- c("Estimate", "ConfIntLow", "ConfIntHigh", "tstat", "pvalue")
histone_stats$Direction <- strsplit(rownames(histone_stats), split = ".", fixed = TRUE) %>% sapply(function(x) x[1])
histone_stats$Antibody <- strsplit(rownames(histone_stats), split = ".", fixed = TRUE) %>% sapply(function(x) x[2])
rownames(histone_stats) <- 1:nrow(histone_stats)
histone_stats$qvalue <- p.adjust(histone_stats$pvalue, method = "fdr")

chromHMM$chroms <- relevel(chromHMM$chroms, ref = "ChrX")
chromHMM_stats <- split(chromHMM, list(chromHMM$userSet, chromHMM$chromState)) %>% 
        sapply(function(x){
                t.test(oddsRatio ~ chroms, data = x, paired = TRUE)[c("estimate", "conf.int", "statistic", "p.value")] %>% 
                        unlist()
        }) %>% t() %>% as.data.frame()
colnames(chromHMM_stats) <- c("Estimate", "ConfIntLow", "ConfIntHigh", "tstat", "pvalue")
chromHMM_stats$Direction <- strsplit(rownames(chromHMM_stats), split = ".", fixed = TRUE) %>% sapply(function(x) x[1])
chromHMM_stats$Antibody <- strsplit(rownames(chromHMM_stats), split = ".", fixed = TRUE) %>% sapply(function(x) x[2])
rownames(chromHMM_stats) <- 1:nrow(chromHMM_stats)
chromHMM_stats$qvalue <- p.adjust(chromHMM_stats$pvalue, method = "fdr")

stats <- rbind(histone_stats, chromHMM_stats)
stats <- stats[,c("Direction", "Antibody", "Estimate", "ConfIntLow", "ConfIntHigh", "tstat", "pvalue", "qvalue")]
write.csv(stats, "Tables/LOLA Replication Males Autosome vs ChrX OR Stats.csv", quote = FALSE, row.names = FALSE)

# Compare Replication Females Autosome and ChrX LOLA ------------------------------
# Histone Mark Enrichment ####
histone_auto_hyper <- read.csv("Tables/LOLA Histone Dx Replication 100 Females autosome HyperDMRs.csv", header = TRUE, 
                               stringsAsFactors = FALSE)
histone_auto_hypo <- read.csv("Tables/LOLA Histone Dx Replication 100 Females autosome HypoDMRs.csv", header = TRUE, 
                              stringsAsFactors = FALSE)
histone_chrX_hyper <- read.csv("Tables/LOLA Histone Dx Replication 100 Females chrX HyperDMRs.csv", header = TRUE, 
                               stringsAsFactors = FALSE)
histone_chrX_hypo <- read.csv("Tables/LOLA Histone Dx Replication 100 Females chrX HypoDMRs.csv", header = TRUE, 
                              stringsAsFactors = FALSE)
histone <- rbind(histone_auto_hyper, histone_auto_hypo, histone_chrX_hyper, histone_chrX_hypo)
histone$chroms <- c(rep("Autosomes", nrow(histone_auto_hyper) + nrow(histone_auto_hypo)),
                    rep("ChrX", nrow(histone_chrX_hyper) + nrow(histone_chrX_hypo))) %>%
        factor(levels = c("Autosomes", "ChrX"))
histone$userSet <- ifelse(histone$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))
histone$antibody <- factor(histone$antibody, levels = c("H3K4me1", "H3K4me3", "H3K9me3", "H3K27me3", "H3K36me3"))

plotLOLAhistoneBox(histone, file = "Figures/LOLA Histone Replication Females Autosome ChrX OR Mean CL.png", facet = vars(userSet),
                   plot.margin = c(1.5, 0.5, 0.5, 0.5), ylim = c(0,12))

# ChromHMM Enrichment ####
chromHMM_auto_hyper <- read.csv("Tables/LOLA ChromHMM Dx Replication 100 Females autosome HyperDMRs.csv", header = TRUE, 
                                stringsAsFactors = FALSE)
chromHMM_auto_hypo <- read.csv("Tables/LOLA ChromHMM Dx Replication 100 Females autosome HypoDMRs.csv", header = TRUE, 
                               stringsAsFactors = FALSE)
chromHMM_chrX_hyper <- read.csv("Tables/LOLA ChromHMM Dx Replication 100 Females chrX HyperDMRs.csv", header = TRUE, 
                                stringsAsFactors = FALSE)
chromHMM_chrX_hypo <- read.csv("Tables/LOLA ChromHMM Dx Replication 100 Females chrX HypoDMRs.csv", header = TRUE, 
                               stringsAsFactors = FALSE)
chromHMM <- rbind(chromHMM_auto_hyper, chromHMM_auto_hypo, chromHMM_chrX_hyper, chromHMM_chrX_hypo)
chromHMM$chroms <- c(rep("Autosomes", nrow(chromHMM_auto_hyper) + nrow(chromHMM_auto_hypo)),
                     rep("ChrX", nrow(chromHMM_chrX_hyper) + nrow(chromHMM_chrX_hypo))) %>%
        factor(levels = c("Autosomes", "ChrX"))
chromHMM$userSet <- ifelse(chromHMM$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))
chromHMM$chromState <- factor(chromHMM$chromState, levels = unique(chromHMM$chromState))

plotLOLAchromHMMBox(chromHMM, file = "Figures/LOLA ChromHMM Replication Females Autosome ChrX OR Mean CL.png", 
                    facet = vars(userSet))

# Differential Enrichment Stats ####
histone$chroms <- relevel(histone$chroms, ref = "ChrX")
histone_stats <- split(histone, list(histone$userSet, histone$antibody)) %>% 
        sapply(function(x){
                t.test(oddsRatio ~ chroms, data = x, paired = TRUE)[c("estimate", "conf.int", "statistic", "p.value")] %>% 
                        unlist()
        }) %>% t() %>% as.data.frame()
colnames(histone_stats) <- c("Estimate", "ConfIntLow", "ConfIntHigh", "tstat", "pvalue")
histone_stats$Direction <- strsplit(rownames(histone_stats), split = ".", fixed = TRUE) %>% sapply(function(x) x[1])
histone_stats$Antibody <- strsplit(rownames(histone_stats), split = ".", fixed = TRUE) %>% sapply(function(x) x[2])
rownames(histone_stats) <- 1:nrow(histone_stats)
histone_stats$qvalue <- p.adjust(histone_stats$pvalue, method = "fdr")

chromHMM$chroms <- relevel(chromHMM$chroms, ref = "ChrX")
chromHMM_stats <- split(chromHMM, list(chromHMM$userSet, chromHMM$chromState)) %>% 
        sapply(function(x){
                t.test(oddsRatio ~ chroms, data = x, paired = TRUE)[c("estimate", "conf.int", "statistic", "p.value")] %>% 
                        unlist()
        }) %>% t() %>% as.data.frame()
colnames(chromHMM_stats) <- c("Estimate", "ConfIntLow", "ConfIntHigh", "tstat", "pvalue")
chromHMM_stats$Direction <- strsplit(rownames(chromHMM_stats), split = ".", fixed = TRUE) %>% sapply(function(x) x[1])
chromHMM_stats$Antibody <- strsplit(rownames(chromHMM_stats), split = ".", fixed = TRUE) %>% sapply(function(x) x[2])
rownames(chromHMM_stats) <- 1:nrow(chromHMM_stats)
chromHMM_stats$qvalue <- p.adjust(chromHMM_stats$pvalue, method = "fdr")

stats <- rbind(histone_stats, chromHMM_stats)
stats <- stats[,c("Direction", "Antibody", "Estimate", "ConfIntLow", "ConfIntHigh", "tstat", "pvalue", "qvalue")]
write.csv(stats, "Tables/LOLA Replication Females Autosome vs ChrX OR Stats.csv", quote = FALSE, row.names = FALSE)
rm(chromHMM, chromHMM_auto_hyper, chromHMM_auto_hypo, chromHMM_chrX_hyper, chromHMM_chrX_hypo, gg, histone, 
   histone_auto_hyper, histone_auto_hypo, histone_chrX_hyper, histone_chrX_hypo, stats, chromHMM_stats, histone_stats)

# Differential Enrichment Summary -----------------------------------------
# Load Data ####
stats_malesDisc <- read.csv("Tables/LOLA Discovery Males Autosome vs ChrX OR Stats.csv", header = TRUE, stringsAsFactors = FALSE)
stats_malesRep <- read.csv("Tables/LOLA Replication Males Autosome vs ChrX OR Stats.csv", header = TRUE, stringsAsFactors = FALSE)
stats_femalesDisc <- read.csv("Tables/LOLA Discovery Females Autosome vs ChrX OR Stats.csv", header = TRUE, stringsAsFactors = FALSE)
stats_femalesRep <- read.csv("Tables/LOLA Replication Females Autosome vs ChrX OR Stats.csv", header = TRUE, stringsAsFactors = FALSE)
stats <- rbind(stats_malesDisc, stats_malesRep, stats_femalesDisc, stats_femalesRep)
stats$DMRs <- rep(c("Discovery\nMales", "Replication\nMales", "Discovery\nFemales", "Replication\nFemales"), each = 40) %>%
        factor(levels = c("Discovery\nMales", "Replication\nMales", "Discovery\nFemales", "Replication\nFemales"))
stats$Direction <- factor(stats$Direction, levels = c("Hypermethylated", "Hypomethylated"))
stats$Significant <- (stats$qvalue < 0.05) %>% factor(levels = c("TRUE", "FALSE"))
stats$Antibody <- factor(stats$Antibody, levels = rev(unique(stats$Antibody)))
hm.lim <- max(abs(stats$Estimate))

# Plot Heatmap ####
gg <- ggplot(data = stats)
gg +
        geom_tile(aes(x = DMRs, y = Antibody, fill = Estimate)) +
        geom_text(aes(x = DMRs, y = Antibody, alpha = Significant), label = "*", color = "white", size = 12, nudge_y = -0.35) +
        facet_grid(cols = vars(Direction)) +
        scale_fill_gradientn("Mean\nOdds Ratio\nDifference", colors = c("#0000FF", "black", "#FF0000"), values = c(0, 0.5, 1), 
                             na.value = "black", limits = c(-hm.lim, hm.lim), 
                             breaks = pretty_breaks(n = 3)) +
        scale_alpha_manual(breaks = c("TRUE", "FALSE"), values = c(1, 0), guide = FALSE) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              panel.background = element_rect(fill = "black"), axis.ticks.x = element_line(size = 1.25), 
              axis.ticks.y = element_line(size = 1.25), legend.key = element_blank(),  legend.position = c(1.14, 0.84), 
              legend.background = element_blank(), legend.title = element_text(size = 20), legend.text = element_text(size = 18),
              plot.margin = unit(c(0, 8, 0.5, 0.5), "lines"), axis.text.y = element_text(size = 18, color = "black"), 
              axis.text.x = element_text(size = 18, color = "black", angle = 90, hjust = 1, vjust = 0.5), 
              axis.title = element_blank(), plot.title = element_text(size = 18, hjust = 0.5, vjust = 0),
              strip.background = element_blank(), strip.text = element_text(size = 20), legend.key.size = unit(1, "lines")) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0))
ggsave("Figures/LOLA Histone and ChromHMM Autosome vs ChrX Enrichment Heatmap.png", dpi = 600, width = 9, height = 8, 
       units = "in")


