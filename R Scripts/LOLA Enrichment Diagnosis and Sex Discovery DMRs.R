# LOLA Enrichment Diagnosis and Sex Discovery DMRs ------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 2/26/19

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

# Get LOLA Enrichments ----------------------------------------------------
# Discovery Diagnosis DMRs ####
DMRs <- loadRegions(file = "DMRs_DxNoXY_Discovery50.csv", chroms = c(paste("chr", 1:22, sep = ""), "chrM"), sort = TRUE)
DMRlist <- list("AllDMRs" = makeGRange(DMRs = DMRs, direction = "all"),
                "HyperDMRs" = makeGRange(DMRs = DMRs, direction = "hyper"),
                "HypoDMRs" = makeGRange(DMRs = DMRs, direction = "hypo"))
Background <- loadRegions(file = "bsseq_background_Discovery50.csv", chroms = c(paste("chr", 1:22, sep = ""), "chrM"), 
                          sort = TRUE) %>% makeGRange(direction = "all")
regionDB <- loadRegionDB(dbLocation = "/share/lasallelab/programs/LOLA/hg38", useCache = TRUE, limit = NULL, 
                         collections = c("Roadmap_ChromHMM", "roadmap_epigenomics"))
Results <- runLOLA(userSets = DMRList, userUniverse = Background, regionDB = regionDB, cores = 2, redefineUserSets = TRUE)
writeCombinedEnrichment(combinedResults = Results, outFolder = "LOLA", includeSplits = FALSE)
file.copy(from = "LOLA/allEnrichments.tsv", to = "LOLA/LOLA_DxNoXY_Discovery50_DMRs.tsv")
rm(DMRs, DMRlist, Background, Results)

# Discovery and Sex Diagnosis DMRs ####
DMRs <- loadRegions(file = "DMRs_DxAdjSex_Discovery50.csv", chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)
DMRlist <- list("AllDMRs" = makeGRange(DMRs = DMRs, direction = "all"),
                "HyperDMRs" = makeGRange(DMRs = DMRs, direction = "hyper"),
                "HypoDMRs" = makeGRange(DMRs = DMRs, direction = "hypo"))
Background <- loadRegions(file = "bsseq_background_Discovery50.csv", chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), 
                          sort = TRUE) %>% makeGRange(direction = "all")
Results <- runLOLA(userSets = DMRlist, userUniverse = Background, regionDB = regionDB, cores = 2, redefineUserSets = TRUE)
writeCombinedEnrichment(combinedResults = Results, outFolder = "LOLA", includeSplits = FALSE)
file.copy(from = "LOLA/allEnrichments.tsv", to = "LOLA/LOLA_DxAdjSex_Discovery50_DMRs.tsv")
rm(DMRs, DMRlist, Background, Results)

# Discovery Diagnosis Males DMRs ####
DMRs <- loadRegions(file = "DMRs_Dx_Discovery50_males.csv", chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), sort = TRUE)
DMRlist <- list("AllDMRs" = makeGRange(DMRs = DMRs, direction = "all"),
                "HyperDMRs" = makeGRange(DMRs = DMRs, direction = "hyper"),
                "HypoDMRs" = makeGRange(DMRs = DMRs, direction = "hypo"))
Background <- loadRegions(file = "bsseq_background_Discovery50_males.csv", chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), 
                          sort = TRUE) %>% makeGRange(direction = "all")
Results <- runLOLA(userSets = DMRlist, userUniverse = Background, regionDB = regionDB, cores = 2, redefineUserSets = TRUE)
writeCombinedEnrichment(combinedResults = Results, outFolder = "LOLA", includeSplits = FALSE)
file.copy(from = "LOLA/allEnrichments.tsv", to = "LOLA/LOLA_Dx_Discovery50_males_DMRs.tsv")
rm(DMRs, DMRlist, Background, Results)

# Discovery Diagnosis Females DMRs ####
DMRs <- loadRegions(file = "DMRs_Dx_Discovery50_females.csv", chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)
DMRlist <- list("AllDMRs" = makeGRange(DMRs = DMRs, direction = "all"),
                "HyperDMRs" = makeGRange(DMRs = DMRs, direction = "hyper"),
                "HypoDMRs" = makeGRange(DMRs = DMRs, direction = "hypo"))
Background <- loadRegions(file = "bsseq_background_Discovery50_females.csv", chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), 
                          sort = TRUE) %>% makeGRange(direction = "all")
Results <- runLOLA(userSets = DMRlist, userUniverse = Background, regionDB = regionDB, cores = 2, redefineUserSets = TRUE)
writeCombinedEnrichment(combinedResults = Results, outFolder = "LOLA", includeSplits = FALSE)
file.copy(from = "LOLA/allEnrichments.tsv", to = "LOLA/LOLA_Dx_Discovery50_females_DMRs.tsv")
rm(DMRs, DMRlist, Background, Results)

# Analyze Discovery Diagnosis DMRs LOLA ------------------------------------
# Load Data ####
lola <- read.delim("Tables/LOLA_DxNoXY_Discovery50_DMRs.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
index <- read.delim("Tables/LOLA Roadmap ChromHMM index.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Histone Mark Enrichment ####
# Prep Data
histone <- subset(lola, collection == "roadmap_epigenomics")
histone_all <- prepLOLAhistone(histone, index = index, regions = "AllDMRs", 
                               file = "Tables/LOLA Histone DxNoXY Discovery 50 AllDMRs.csv")
histone_hyper <- prepLOLAhistone(histone, index = index, regions = "HyperDMRs", 
                               file = "Tables/LOLA Histone DxNoXY Discovery 50 HyperDMRs.csv")
histone_hypo <- prepLOLAhistone(histone, index = index, regions = "HypoDMRs", 
                               file = "Tables/LOLA Histone DxNoXY Discovery 50 HypoDMRs.csv")

# Plot Odds Ratio
hm.max <- quantile(c(histone_all$oddsRatio, histone_hyper$oddsRatio, histone_hypo$oddsRatio), probs = 0.975, 
                   names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_all, title = "All DMRs", type = "oddsRatio", hm.max = hm.max, 
                file = "Figures/LOLA Histone DxNoXY Discovery 50 All DMRs OR Heatmap.png")
plotLOLAhistone(histone = histone_hyper, title = "Hyper DMRs", type = "oddsRatio", hm.max = hm.max, 
                file = "Figures/LOLA Histone DxNoXY Discovery 50 Hyper DMRs OR Heatmap.png")
plotLOLAhistone(histone = histone_hypo, title = "Hypo DMRs", type = "oddsRatio", hm.max = hm.max, 
                file = "Figures/LOLA Histone DxNoXY Discovery 50 Hypo DMRs OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(c(histone_all$qValueLog, histone_hyper$qValueLog, histone_hypo$qValueLog), probs = 0.975, 
                   names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_all, title = "All DMRs", type = "qValueLog", hm.max = hm.max, 
                file = "Figures/LOLA Histone DxNoXY Discovery 50 All DMRs log qvalue Heatmap.png")
plotLOLAhistone(histone = histone_hyper, title = "Hyper DMRs", type = "qValueLog", hm.max = hm.max, 
                file = "Figures/LOLA Histone DxNoXY Discovery 50 Hyper DMRs log qvalue Heatmap.png")
plotLOLAhistone(histone = histone_hypo, title = "Hypo DMRs", type = "qValueLog", hm.max = hm.max, 
                file = "Figures/LOLA Histone DxNoXY Discovery 50 Hypo DMRs log qvalue Heatmap.png")

# Plot Legend
plotLOLAhistone(histone = histone_all, type = "legend", file = "Figures/LOLA Discovery 50 DMRs Heatmap Legend.png")

# ChromHMM Enrichment ####
# Prep Data
chromHMM <- subset(lola, collection == "Roadmap_ChromHMM")
chromHMM_all <- prepLOLAchromHMM(chromHMM, index = index, regions = "AllDMRs", 
                                 file = "Tables/LOLA chromHMM DxNoXY Discovery 50 AllDMRs.csv")
chromHMM_hyper <- prepLOLAchromHMM(chromHMM, index = index, regions = "HyperDMRs", 
                                   file = "Tables/LOLA chromHMM DxNoXY Discovery 50 HyperDMRs.csv")
chromHMM_hypo <- prepLOLAchromHMM(chromHMM, index = index, regions = "HypoDMRs", 
                                  file = "Tables/LOLA chromHMM DxNoXY Discovery 50 HypoDMRs.csv")

# Plot Odds Ratio
hm.max <- quantile(c(chromHMM_all$oddsRatio, chromHMM_hyper$oddsRatio, chromHMM_hypo$oddsRatio), probs = 0.975, 
                   names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_all, title = "All DMRs", type = "oddsRatio", hm.max = hm.max, 
                file = "Figures/LOLA chromHMM DxNoXY Discovery 50 All DMRs OR Heatmap.png")
plotLOLAchromHMM(chromHMM = chromHMM_hyper, title = "Hyper DMRs", type = "oddsRatio", hm.max = hm.max, 
                file = "Figures/LOLA chromHMM DxNoXY Discovery 50 Hyper DMRs OR Heatmap.png")
plotLOLAchromHMM(chromHMM = chromHMM_hypo, title = "Hypo DMRs", type = "oddsRatio", hm.max = hm.max, 
                file = "Figures/LOLA chromHMM DxNoXY Discovery 50 Hypo DMRs OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(c(chromHMM_all$qValueLog, chromHMM_hyper$qValueLog, chromHMM_hypo$qValueLog), probs = 0.975, 
                   names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_all, title = "All DMRs", type = "qValueLog", hm.max = hm.max, 
                file = "Figures/LOLA chromHMM DxNoXY Discovery 50 All DMRs log qvalue Heatmap.png")
plotLOLAchromHMM(chromHMM = chromHMM_hyper, title = "Hyper DMRs", type = "qValueLog", hm.max = hm.max, 
                file = "Figures/LOLA chromHMM DxNoXY Discovery 50 Hyper DMRs log qvalue Heatmap.png")
plotLOLAchromHMM(chromHMM = chromHMM_hypo, title = "Hypo DMRs", type = "qValueLog", hm.max = hm.max, 
                file = "Figures/LOLA chromHMM DxNoXY Discovery 50 Hypo DMRs log qvalue Heatmap.png")

# Analyze Discovery Diagnosis and Sex DMRs LOLA ------------------------------------
# Load Data ####
lola <- read.delim("Tables/LOLA_DxAdjSex_Discovery50_DMRs.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
index <- read.delim("Tables/LOLA Roadmap ChromHMM index.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Histone Mark Enrichment ####
# Prep Data
histone <- subset(lola, collection == "roadmap_epigenomics")
histone_all <- prepLOLAhistone(histone, index = index, regions = "AllDMRs", 
                               file = "Tables/LOLA Histone DxAdjSex Discovery 50 AllDMRs.csv")
histone_hyper <- prepLOLAhistone(histone, index = index, regions = "HyperDMRs", 
                                 file = "Tables/LOLA Histone DxAdjSex Discovery 50 HyperDMRs.csv")
histone_hypo <- prepLOLAhistone(histone, index = index, regions = "HypoDMRs", 
                                file = "Tables/LOLA Histone DxAdjSex Discovery 50 HypoDMRs.csv")

# Plot Odds Ratio
hm.max <- quantile(c(histone_all$oddsRatio, histone_hyper$oddsRatio, histone_hypo$oddsRatio), probs = 0.975, 
                   names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_all, title = "All DMRs", type = "oddsRatio", hm.max = hm.max, 
                file = "Figures/LOLA Histone DxAdjSex Discovery 50 All DMRs OR Heatmap.png")
plotLOLAhistone(histone = histone_hyper, title = "Hyper DMRs", type = "oddsRatio", hm.max = hm.max, 
                file = "Figures/LOLA Histone DxAdjSex Discovery 50 Hyper DMRs OR Heatmap.png")
plotLOLAhistone(histone = histone_hypo, title = "Hypo DMRs", type = "oddsRatio", hm.max = hm.max, 
                file = "Figures/LOLA Histone DxAdjSex Discovery 50 Hypo DMRs OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(c(histone_all$qValueLog, histone_hyper$qValueLog, histone_hypo$qValueLog), probs = 0.975, 
                   names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_all, title = "All DMRs", type = "qValueLog", hm.max = hm.max, 
                file = "Figures/LOLA Histone DxAdjSex Discovery 50 All DMRs log qvalue Heatmap.png")
plotLOLAhistone(histone = histone_hyper, title = "Hyper DMRs", type = "qValueLog", hm.max = hm.max, 
                file = "Figures/LOLA Histone DxAdjSex Discovery 50 Hyper DMRs log qvalue Heatmap.png")
plotLOLAhistone(histone = histone_hypo, title = "Hypo DMRs", type = "qValueLog", hm.max = hm.max, 
                file = "Figures/LOLA Histone DxAdjSex Discovery 50 Hypo DMRs log qvalue Heatmap.png")

# ChromHMM Enrichment ####
# Prep Data
chromHMM <- subset(lola, collection == "Roadmap_ChromHMM")
chromHMM_all <- prepLOLAchromHMM(chromHMM, index = index, regions = "AllDMRs", 
                                 file = "Tables/LOLA chromHMM DxAdjSex Discovery 50 AllDMRs.csv")
chromHMM_hyper <- prepLOLAchromHMM(chromHMM, index = index, regions = "HyperDMRs", 
                                   file = "Tables/LOLA chromHMM DxAdjSex Discovery 50 HyperDMRs.csv")
chromHMM_hypo <- prepLOLAchromHMM(chromHMM, index = index, regions = "HypoDMRs", 
                                  file = "Tables/LOLA chromHMM DxAdjSex Discovery 50 HypoDMRs.csv")

# Plot Odds Ratio
hm.max <- quantile(c(chromHMM_all$oddsRatio, chromHMM_hyper$oddsRatio, chromHMM_hypo$oddsRatio), probs = 0.975, 
                   names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_all, title = "All DMRs", type = "oddsRatio", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM DxAdjSex Discovery 50 All DMRs OR Heatmap.png")
plotLOLAchromHMM(chromHMM = chromHMM_hyper, title = "Hyper DMRs", type = "oddsRatio", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM DxAdjSex Discovery 50 Hyper DMRs OR Heatmap.png")
plotLOLAchromHMM(chromHMM = chromHMM_hypo, title = "Hypo DMRs", type = "oddsRatio", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM DxAdjSex Discovery 50 Hypo DMRs OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(c(chromHMM_all$qValueLog, chromHMM_hyper$qValueLog, chromHMM_hypo$qValueLog), probs = 0.975, 
                   names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_all, title = "All DMRs", type = "qValueLog", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM DxAdjSex Discovery 50 All DMRs log qvalue Heatmap.png")
plotLOLAchromHMM(chromHMM = chromHMM_hyper, title = "Hyper DMRs", type = "qValueLog", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM DxAdjSex Discovery 50 Hyper DMRs log qvalue Heatmap.png")
plotLOLAchromHMM(chromHMM = chromHMM_hypo, title = "Hypo DMRs", type = "qValueLog", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM DxAdjSex Discovery 50 Hypo DMRs log qvalue Heatmap.png")

# Analyze Discovery Diagnosis Males DMRs LOLA ------------------------------------
# Load Data ####
lola <- read.delim("Tables/LOLA_Dx_Discovery50_males_DMRs.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
index <- read.delim("Tables/LOLA Roadmap ChromHMM index.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Histone Mark Enrichment ####
# Prep Data
histone <- subset(lola, collection == "roadmap_epigenomics")
histone_all <- prepLOLAhistone(histone, index = index, regions = "AllDMRs", 
                               file = "Tables/LOLA Histone Dx Discovery 50 Males AllDMRs.csv")
histone_hyper <- prepLOLAhistone(histone, index = index, regions = "HyperDMRs", 
                                 file = "Tables/LOLA Histone Dx Discovery 50 Males HyperDMRs.csv")
histone_hypo <- prepLOLAhistone(histone, index = index, regions = "HypoDMRs", 
                                file = "Tables/LOLA Histone Dx Discovery 50 Males HypoDMRs.csv")

# Plot Odds Ratio
hm.max <- quantile(c(histone_all$oddsRatio, histone_hyper$oddsRatio, histone_hypo$oddsRatio), probs = 0.975, 
                   names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_all, title = "All DMRs", type = "oddsRatio", hm.max = hm.max, 
                file = "Figures/LOLA Histone Dx Discovery 50 Males All DMRs OR Heatmap.png")
plotLOLAhistone(histone = histone_hyper, title = "Hyper DMRs", type = "oddsRatio", hm.max = hm.max, 
                file = "Figures/LOLA Histone Dx Discovery 50 Males Hyper DMRs OR Heatmap.png")
plotLOLAhistone(histone = histone_hypo, title = "Hypo DMRs", type = "oddsRatio", hm.max = hm.max, 
                file = "Figures/LOLA Histone Dx Discovery 50 Males Hypo DMRs OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(c(histone_all$qValueLog, histone_hyper$qValueLog, histone_hypo$qValueLog), probs = 0.975, 
                   names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_all, title = "All DMRs", type = "qValueLog", hm.max = hm.max, 
                file = "Figures/LOLA Histone Dx Discovery 50 Males All DMRs log qvalue Heatmap.png")
plotLOLAhistone(histone = histone_hyper, title = "Hyper DMRs", type = "qValueLog", hm.max = hm.max, 
                file = "Figures/LOLA Histone Dx Discovery 50 Males Hyper DMRs log qvalue Heatmap.png")
plotLOLAhistone(histone = histone_hypo, title = "Hypo DMRs", type = "qValueLog", hm.max = hm.max, 
                file = "Figures/LOLA Histone Dx Discovery 50 Males Hypo DMRs log qvalue Heatmap.png")

# ChromHMM Enrichment ####
# Prep Data
chromHMM <- subset(lola, collection == "Roadmap_ChromHMM")
chromHMM_all <- prepLOLAchromHMM(chromHMM, index = index, regions = "AllDMRs", 
                                 file = "Tables/LOLA chromHMM Dx Discovery 50 Males AllDMRs.csv")
chromHMM_hyper <- prepLOLAchromHMM(chromHMM, index = index, regions = "HyperDMRs", 
                                   file = "Tables/LOLA chromHMM Dx Discovery 50 Males HyperDMRs.csv")
chromHMM_hypo <- prepLOLAchromHMM(chromHMM, index = index, regions = "HypoDMRs", 
                                  file = "Tables/LOLA chromHMM Dx Discovery 50 Males HypoDMRs.csv")

# Plot Odds Ratio
hm.max <- quantile(c(chromHMM_all$oddsRatio, chromHMM_hyper$oddsRatio, chromHMM_hypo$oddsRatio), probs = 0.975, 
                   names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_all, title = "All DMRs", type = "oddsRatio", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Males All DMRs OR Heatmap.png")
plotLOLAchromHMM(chromHMM = chromHMM_hyper, title = "Hyper DMRs", type = "oddsRatio", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Males Hyper DMRs OR Heatmap.png")
plotLOLAchromHMM(chromHMM = chromHMM_hypo, title = "Hypo DMRs", type = "oddsRatio", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Males Hypo DMRs OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(c(chromHMM_all$qValueLog, chromHMM_hyper$qValueLog, chromHMM_hypo$qValueLog), probs = 0.975, 
                   names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_all, title = "All DMRs", type = "qValueLog", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Males All DMRs log qvalue Heatmap.png")
plotLOLAchromHMM(chromHMM = chromHMM_hyper, title = "Hyper DMRs", type = "qValueLog", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Males Hyper DMRs log qvalue Heatmap.png")
plotLOLAchromHMM(chromHMM = chromHMM_hypo, title = "Hypo DMRs", type = "qValueLog", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Males Hypo DMRs log qvalue Heatmap.png")

# Analyze Discovery Diagnosis Females DMRs LOLA ------------------------------------
# Load Data ####
lola <- read.delim("Tables/LOLA_Dx_Discovery50_females_DMRs.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
index <- read.delim("Tables/LOLA Roadmap ChromHMM index.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Histone Mark Enrichment ####
# Prep Data
histone <- subset(lola, collection == "roadmap_epigenomics")
histone_all <- prepLOLAhistone(histone, index = index, regions = "AllDMRs", 
                               file = "Tables/LOLA Histone Dx Discovery 50 Females AllDMRs.csv")
histone_hyper <- prepLOLAhistone(histone, index = index, regions = "HyperDMRs", 
                                 file = "Tables/LOLA Histone Dx Discovery 50 Females HyperDMRs.csv")
histone_hypo <- prepLOLAhistone(histone, index = index, regions = "HypoDMRs", 
                                file = "Tables/LOLA Histone Dx Discovery 50 Females HypoDMRs.csv")

# Plot Odds Ratio
hm.max <- quantile(c(histone_all$oddsRatio, histone_hyper$oddsRatio, histone_hypo$oddsRatio), probs = 0.975, 
                   names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_all, title = "All DMRs", type = "oddsRatio", hm.max = hm.max, 
                file = "Figures/LOLA Histone Dx Discovery 50 Females All DMRs OR Heatmap.png")
plotLOLAhistone(histone = histone_hyper, title = "Hyper DMRs", type = "oddsRatio", hm.max = hm.max, 
                file = "Figures/LOLA Histone Dx Discovery 50 Females Hyper DMRs OR Heatmap.png")
plotLOLAhistone(histone = histone_hypo, title = "Hypo DMRs", type = "oddsRatio", hm.max = hm.max, 
                file = "Figures/LOLA Histone Dx Discovery 50 Females Hypo DMRs OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(c(histone_all$qValueLog, histone_hyper$qValueLog, histone_hypo$qValueLog), probs = 0.975, 
                   names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_all, title = "All DMRs", type = "qValueLog", hm.max = hm.max, 
                file = "Figures/LOLA Histone Dx Discovery 50 Females All DMRs log qvalue Heatmap.png")
plotLOLAhistone(histone = histone_hyper, title = "Hyper DMRs", type = "qValueLog", hm.max = hm.max, 
                file = "Figures/LOLA Histone Dx Discovery 50 Females Hyper DMRs log qvalue Heatmap.png")
plotLOLAhistone(histone = histone_hypo, title = "Hypo DMRs", type = "qValueLog", hm.max = hm.max, 
                file = "Figures/LOLA Histone Dx Discovery 50 Females Hypo DMRs log qvalue Heatmap.png")

# ChromHMM Enrichment ####
# Prep Data
chromHMM <- subset(lola, collection == "Roadmap_ChromHMM")
chromHMM_all <- prepLOLAchromHMM(chromHMM, index = index, regions = "AllDMRs", 
                                 file = "Tables/LOLA chromHMM Dx Discovery 50 Females AllDMRs.csv")
chromHMM_hyper <- prepLOLAchromHMM(chromHMM, index = index, regions = "HyperDMRs", 
                                   file = "Tables/LOLA chromHMM Dx Discovery 50 Females HyperDMRs.csv")
chromHMM_hypo <- prepLOLAchromHMM(chromHMM, index = index, regions = "HypoDMRs", 
                                  file = "Tables/LOLA chromHMM Dx Discovery 50 Females HypoDMRs.csv")

# Plot Odds Ratio
hm.max <- quantile(c(chromHMM_all$oddsRatio, chromHMM_hyper$oddsRatio, chromHMM_hypo$oddsRatio), probs = 0.975, 
                   names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_all, title = "All DMRs", type = "oddsRatio", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Females All DMRs OR Heatmap.png")
plotLOLAchromHMM(chromHMM = chromHMM_hyper, title = "Hyper DMRs", type = "oddsRatio", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Females Hyper DMRs OR Heatmap.png")
plotLOLAchromHMM(chromHMM = chromHMM_hypo, title = "Hypo DMRs", type = "oddsRatio", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Females Hypo DMRs OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(c(chromHMM_all$qValueLog, chromHMM_hyper$qValueLog, chromHMM_hypo$qValueLog), probs = 0.975, 
                   names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_all, title = "All DMRs", type = "qValueLog", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Females All DMRs log qvalue Heatmap.png")
plotLOLAchromHMM(chromHMM = chromHMM_hyper, title = "Hyper DMRs", type = "qValueLog", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Females Hyper DMRs log qvalue Heatmap.png")
plotLOLAchromHMM(chromHMM = chromHMM_hypo, title = "Hypo DMRs", type = "qValueLog", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Females Hypo DMRs log qvalue Heatmap.png")












