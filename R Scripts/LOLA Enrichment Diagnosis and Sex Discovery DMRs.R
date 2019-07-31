# LOLA Enrichment Diagnosis and Sex DMRs ------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 6/17/19
# Excluded JLCM032B and JLCM050B

# Packages ####
sapply(c("tidyverse", "LOLA", "simpleCache", "GenomicRanges", "qvalue", "RColorBrewer", "scales", "reshape2"), 
       require, character.only = TRUE)

# Functions ####
# Epigenerate
# loadRegions <- function(file, chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), sort = TRUE){
#         if(grepl("txt", file, fixed = TRUE)){
#                 regions <- read.delim(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#         }
#         else{
#                 regions <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
#         }
#         if("seqnames" %in% colnames(regions)){
#                 colnames(regions)[colnames(regions) == "seqnames"] <- "chr"
#         }
#         regions <- subset(regions, chr %in% chroms)
#         regions$chr <- factor(regions$chr, levels = chroms)
#         if(sort){
#                 regions <- regions[order(regions$chr, regions$start),]
#         }
#         return(regions)
# }
# makeGRange <- function(DMRs, direction = c("all", "hyper", "hypo")){
#         if(direction == "hyper"){DMRs <- subset(DMRs, percentDifference > 0)}
#         if(direction == "hypo"){DMRs <- subset(DMRs, percentDifference < 0)}
#         GR <- GRanges(seqnames = DMRs$chr, ranges = IRanges(start = DMRs$start, end = DMRs$end))
# }

# Laptop
source("R Scripts/DMR Analysis Functions.R")

# Get LOLA Enrichments ----------------------------------------------------
# Data ####
regionDB <- loadRegionDB(dbLocation = "/share/lasallelab/programs/LOLA/hg38", useCache = TRUE, limit = NULL, 
                         collections = c("Roadmap_ChromHMM", "roadmap_epigenomics"))
maleChroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")
femaleChroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrM")

# Discovery Diagnosis Males DMRs ####
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Discovery/")
DMRs <- loadRegions(file = "Dx_Males/DMRs_Dx_Discovery50_males.csv", chroms = maleChroms, sort = TRUE)
DMRlist <- list("AllDMRs" = makeGRange(DMRs = DMRs, direction = "all"),
                "HyperDMRs" = makeGRange(DMRs = DMRs, direction = "hyper"),
                "HypoDMRs" = makeGRange(DMRs = DMRs, direction = "hypo"))
Background <- loadRegions(file = "Dx_Males/bsseq_background_Discovery50_males.csv", chroms = maleChroms, sort = TRUE) %>% 
        makeGRange(direction = "all")
Results <- runLOLA(userSets = DMRlist, userUniverse = Background, regionDB = regionDB, cores = 2, redefineUserSets = TRUE)
writeCombinedEnrichment(combinedResults = Results, outFolder = "LOLA", includeSplits = FALSE)
file.copy(from = "LOLA/allEnrichments.tsv", to = "LOLA/LOLA_Dx_Discovery50_males_DMRs.tsv", overwrite = TRUE)
rm(DMRs, DMRlist, Background, Results)

# Discovery Diagnosis Females DMRs ####
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Discovery/")
DMRs <- loadRegions(file = "Dx_Females/DMRs_Dx_Discovery50_females.csv", chroms = femaleChroms, sort = TRUE)
DMRlist <- list("AllDMRs" = makeGRange(DMRs = DMRs, direction = "all"),
                "HyperDMRs" = makeGRange(DMRs = DMRs, direction = "hyper"),
                "HypoDMRs" = makeGRange(DMRs = DMRs, direction = "hypo"))
Background <- loadRegions(file = "Dx_Females/bsseq_background_Discovery50_females.csv", chroms = femaleChroms, sort = TRUE) %>% 
        makeGRange(direction = "all")
Results <- runLOLA(userSets = DMRlist, userUniverse = Background, regionDB = regionDB, cores = 2, redefineUserSets = TRUE)
writeCombinedEnrichment(combinedResults = Results, outFolder = "LOLA", includeSplits = FALSE)
file.copy(from = "LOLA/allEnrichments.tsv", to = "LOLA/LOLA_Dx_Discovery50_females_DMRs.tsv", overwrite = TRUE)
rm(DMRs, DMRlist, Background, Results)

# Replication Diagnosis Males DMRs ####
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Replication/")
DMRs <- loadRegions(file = "Dx_Males/DMRs_Dx_Replication50_males.csv", chroms = maleChroms, sort = TRUE)
DMRlist <- list("AllDMRs" = makeGRange(DMRs = DMRs, direction = "all"),
                "HyperDMRs" = makeGRange(DMRs = DMRs, direction = "hyper"),
                "HypoDMRs" = makeGRange(DMRs = DMRs, direction = "hypo"))
Background <- loadRegions(file = "Dx_Males/bsseq_background_Replication50_males.csv", chroms = maleChroms, sort = TRUE) %>% 
        makeGRange(direction = "all")
Results <- runLOLA(userSets = DMRlist, userUniverse = Background, regionDB = regionDB, cores = 2, redefineUserSets = TRUE)
writeCombinedEnrichment(combinedResults = Results, outFolder = "LOLA", includeSplits = FALSE)
file.copy(from = "LOLA/allEnrichments.tsv", to = "LOLA/LOLA_Dx_Replication50_males_DMRs.tsv", overwrite = TRUE)
rm(DMRs, DMRlist, Background, Results)

# Replication Diagnosis Females DMRs ####
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Replication/")
DMRs <- loadRegions(file = "Dx_Females_100/DMRs_Dx_Replication100_females.csv", chroms = femaleChroms, sort = TRUE)
DMRlist <- list("AllDMRs" = makeGRange(DMRs = DMRs, direction = "all"),
                "HyperDMRs" = makeGRange(DMRs = DMRs, direction = "hyper"),
                "HypoDMRs" = makeGRange(DMRs = DMRs, direction = "hypo"))
Background <- loadRegions(file = "Dx_Females_100/bsseq_background_Replication100_females.csv", chroms = femaleChroms, sort = TRUE) %>% 
        makeGRange(direction = "all")
Results <- runLOLA(userSets = DMRlist, userUniverse = Background, regionDB = regionDB, cores = 2, redefineUserSets = TRUE)
writeCombinedEnrichment(combinedResults = Results, outFolder = "LOLA", includeSplits = FALSE)
file.copy(from = "LOLA/allEnrichments.tsv", to = "LOLA/LOLA_Dx_Replication100_females_DMRs.tsv", overwrite = TRUE)
rm(DMRs, DMRlist, Background, Results)

# Analyze Discovery Diagnosis Males DMRs LOLA ------------------------------------
# Load Data ####
lola <- read.delim("Tables/LOLA_Dx_Discovery50_males_DMRs.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
index <- read.delim("Tables/LOLA Roadmap ChromHMM index.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Histone Mark Enrichment ####
# Prep Data
histone <- subset(lola, collection == "roadmap_epigenomics")
histone_hyper <- prepLOLAhistone(histone, index = index, regions = "HyperDMRs", 
                                 file = "Tables/LOLA Histone Dx Discovery 50 Males HyperDMRs.csv")
histone_hypo <- prepLOLAhistone(histone, index = index, regions = "HypoDMRs", 
                                file = "Tables/LOLA Histone Dx Discovery 50 Males HypoDMRs.csv")
histone_combined <- rbind(histone_hyper, histone_hypo)
histone_combined$userSet <- ifelse(histone_combined$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))

# Plot Odds Ratio
hm.max <- quantile(histone_combined$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_combined, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                width = 9, height = 7, legend.position = c(1.11, 0.855),
                file = "Figures/LOLA Histone Dx Discovery 50 Males OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(histone_combined$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_combined, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                width = 9, height = 7, legend.position = c(1.12, 0.855),
                file = "Figures/LOLA Histone Dx Discovery 50 Males log qvalue Heatmap.png")

# Plot Legend
plotLOLAhistone(histone = histone_hyper, type = "legend", file = "Figures/LOLA DMRs Heatmap Legend.png")

# ChromHMM Enrichment ####
# Prep Data
chromHMM <- subset(lola, collection == "Roadmap_ChromHMM")
chromHMM_hyper <- prepLOLAchromHMM(chromHMM, index = index, regions = "HyperDMRs", 
                                   file = "Tables/LOLA chromHMM Dx Discovery 50 Males HyperDMRs.csv")
chromHMM_hypo <- prepLOLAchromHMM(chromHMM, index = index, regions = "HypoDMRs", 
                                  file = "Tables/LOLA chromHMM Dx Discovery 50 Males HypoDMRs.csv")
chromHMM_combined <- rbind(chromHMM_hyper, chromHMM_hypo)
chromHMM_combined$userSet <- ifelse(chromHMM_combined$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))

# Plot Odds Ratio
hm.max <- quantile(chromHMM_combined$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_combined, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                 width = 9, height = 7, legend.position = c(1.11, 0.855),
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Males OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(chromHMM_combined$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_combined, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                 width = 9, height = 7, legend.position = c(1.12, 0.855),
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Males log qvalue Heatmap.png")

# Top Enriched ####
# Top enriched if regions in at least 50% of cell types are above the median of odds ratio and log(qvalue) for that DMR set
table(histone_hyper$oddsRatio >= median(histone_hyper$oddsRatio) &
              histone_hyper$qValueLog >= median(histone_hyper$qValueLog) &
              histone_hyper$qValue < 0.05,
      histone_hyper$antibody)["TRUE",] > length(unique(histone_hyper$cellType)) * 0.5
# H3K4me1  H3K4me3  H3K9me3 H3K27me3 H3K36me3 
#    TRUE     TRUE    FALSE     TRUE    FALSE

table(histone_hypo$oddsRatio >= median(histone_hypo$oddsRatio) &
              histone_hypo$qValueLog >= median(histone_hypo$qValueLog) &
              histone_hypo$qValue < 0.05,
      histone_hypo$antibody)["TRUE",] > length(unique(histone_hypo$cellType)) * 0.5
# H3K4me1  H3K4me3  H3K9me3 H3K27me3 H3K36me3 
#   FALSE     TRUE    FALSE     TRUE    FALSE 

table(chromHMM_hyper$oddsRatio >= median(chromHMM_hyper$oddsRatio) &
              chromHMM_hyper$qValueLog >= median(chromHMM_hyper$qValueLog) &
              chromHMM_hyper$qValue < 0.05,
      chromHMM_hyper$chromState)["TRUE",] > length(unique(chromHMM_hyper$cellType)) * 0.5
#     TssA TssAFlnk   TxFlnk       Tx     TxWk     EnhG      Enh  ZnfRpts      Het   TssBiv  BivFlnk   EnhBiv   ReprPC 
#     TRUE     TRUE     TRUE    FALSE    FALSE     TRUE    FALSE    FALSE    FALSE     TRUE     TRUE     TRUE     TRUE 
# ReprPCwk    Quies 
#    FALSE    FALSE 

table(chromHMM_hypo$oddsRatio >= median(chromHMM_hypo$oddsRatio) &
              chromHMM_hypo$qValueLog >= median(chromHMM_hypo$qValueLog) &
              chromHMM_hypo$qValue < 0.05,
      chromHMM_hypo$chromState)["TRUE",] > length(unique(chromHMM_hypo$cellType)) * 0.5
#     TssA TssAFlnk   TxFlnk       Tx     TxWk     EnhG      Enh  ZnfRpts      Het   TssBiv  BivFlnk   EnhBiv   ReprPC 
#     TRUE     TRUE    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE     TRUE     TRUE     TRUE     TRUE 
# ReprPCwk    Quies 
#    FALSE    FALSE 

# Analyze Discovery Diagnosis Females DMRs LOLA ------------------------------------
# Load Data ####
lola <- read.delim("Tables/LOLA_Dx_Discovery50_females_DMRs.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Histone Mark Enrichment ####
# Prep Data
histone <- subset(lola, collection == "roadmap_epigenomics")
histone_hyper <- prepLOLAhistone(histone, index = index, regions = "HyperDMRs", 
                                 file = "Tables/LOLA Histone Dx Discovery 50 Females HyperDMRs.csv")
histone_hypo <- prepLOLAhistone(histone, index = index, regions = "HypoDMRs", 
                                file = "Tables/LOLA Histone Dx Discovery 50 Females HypoDMRs.csv")
histone_combined <- rbind(histone_hyper, histone_hypo)
histone_combined$userSet <- ifelse(histone_combined$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))

# Plot Odds Ratio
hm.max <- quantile(histone_combined$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_combined, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                width = 9, height = 7, legend.position = c(1.11, 0.855),
                file = "Figures/LOLA Histone Dx Discovery 50 Females OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(histone_combined$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_combined, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                width = 9, height = 7, legend.position = c(1.12, 0.855),
                file = "Figures/LOLA Histone Dx Discovery 50 Females log qvalue Heatmap.png")

# ChromHMM Enrichment ####
# Prep Data
chromHMM <- subset(lola, collection == "Roadmap_ChromHMM")
chromHMM_hyper <- prepLOLAchromHMM(chromHMM, index = index, regions = "HyperDMRs", 
                                   file = "Tables/LOLA chromHMM Dx Discovery 50 Females HyperDMRs.csv")
chromHMM_hypo <- prepLOLAchromHMM(chromHMM, index = index, regions = "HypoDMRs", 
                                  file = "Tables/LOLA chromHMM Dx Discovery 50 Females HypoDMRs.csv")
chromHMM_combined <- rbind(chromHMM_hyper, chromHMM_hypo)
chromHMM_combined$userSet <- ifelse(chromHMM_combined$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))

# Plot Odds Ratio
hm.max <- quantile(chromHMM_combined$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_combined, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                 width = 9, height = 7, legend.position = c(1.11, 0.855),
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Females OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(chromHMM_combined$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_combined, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                 width = 9, height = 7, legend.position = c(1.12, 0.855),
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Females log qvalue Heatmap.png")

# Top Enriched ####
# Top enriched if regions in at least 50% of cell types are above the median of odds ratio and log(qvalue) for that DMR set
table(histone_hyper$oddsRatio >= median(histone_hyper$oddsRatio) &
              histone_hyper$qValueLog >= median(histone_hyper$qValueLog) &
              histone_hyper$qValue < 0.05,
      histone_hyper$antibody)["TRUE",] > length(unique(histone_hyper$cellType)) * 0.5
# H3K4me1  H3K4me3  H3K9me3 H3K27me3 H3K36me3 
#    TRUE     TRUE    FALSE     TRUE    FALSE

table(histone_hypo$oddsRatio >= median(histone_hypo$oddsRatio) &
              histone_hypo$qValueLog >= median(histone_hypo$qValueLog) &
              histone_hypo$qValue < 0.05,
      histone_hypo$antibody)["TRUE",] > length(unique(histone_hypo$cellType)) * 0.5
# H3K4me1  H3K4me3  H3K9me3 H3K27me3 H3K36me3 
#    TRUE     TRUE    FALSE     TRUE    FALSE 

table(chromHMM_hyper$oddsRatio >= median(chromHMM_hyper$oddsRatio) &
              chromHMM_hyper$qValueLog >= median(chromHMM_hyper$qValueLog) &
              chromHMM_hyper$qValue < 0.05,
      chromHMM_hyper$chromState)["TRUE",] > length(unique(chromHMM_hyper$cellType)) * 0.5
#     TssA TssAFlnk   TxFlnk       Tx     TxWk     EnhG      Enh  ZnfRpts      Het   TssBiv  BivFlnk   EnhBiv   ReprPC 
#     TRUE     TRUE    FALSE    FALSE    FALSE     TRUE    FALSE    FALSE    FALSE    FALSE     TRUE     TRUE     TRUE 
# ReprPCwk    Quies 
#    FALSE    FALSE 

table(chromHMM_hypo$oddsRatio >= median(chromHMM_hypo$oddsRatio) &
              chromHMM_hypo$qValueLog >= median(chromHMM_hypo$qValueLog) &
              chromHMM_hypo$qValue < 0.05,
      chromHMM_hypo$chromState)["TRUE",] > length(unique(chromHMM_hypo$cellType)) * 0.5
#     TssA TssAFlnk   TxFlnk       Tx     TxWk     EnhG      Enh  ZnfRpts      Het   TssBiv  BivFlnk   EnhBiv   ReprPC 
#     TRUE     TRUE    FALSE    FALSE    FALSE     TRUE    FALSE    FALSE    FALSE    FALSE     TRUE     TRUE     TRUE 
# ReprPCwk    Quies 
#    FALSE    FALSE 

# Analyze Replication Diagnosis Males DMRs LOLA ------------------------------------
# Load Data ####
lola <- read.delim("Tables/LOLA_Dx_Replication50_males_DMRs.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Histone Mark Enrichment ####
# Prep Data
histone <- subset(lola, collection == "roadmap_epigenomics")
histone_hyper <- prepLOLAhistone(histone, index = index, regions = "HyperDMRs", 
                                 file = "Tables/LOLA Histone Dx Replication 50 Males HyperDMRs.csv")
histone_hypo <- prepLOLAhistone(histone, index = index, regions = "HypoDMRs", 
                                file = "Tables/LOLA Histone Dx Replication 50 Males HypoDMRs.csv")
histone_combined <- rbind(histone_hyper, histone_hypo)
histone_combined$userSet <- ifelse(histone_combined$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))

# Plot Odds Ratio
hm.max <- quantile(histone_combined$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_combined, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                width = 9, height = 7, legend.position = c(1.11, 0.855),
                file = "Figures/LOLA Histone Dx Replication 50 Males OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(histone_combined$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_combined, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                width = 9, height = 7, legend.position = c(1.12, 0.855),
                file = "Figures/LOLA Histone Dx Replication 50 Males log qvalue Heatmap.png")

# ChromHMM Enrichment ####
# Prep Data
chromHMM <- subset(lola, collection == "Roadmap_ChromHMM")
chromHMM_hyper <- prepLOLAchromHMM(chromHMM, index = index, regions = "HyperDMRs", 
                                   file = "Tables/LOLA chromHMM Dx Replication 50 Males HyperDMRs.csv")
chromHMM_hypo <- prepLOLAchromHMM(chromHMM, index = index, regions = "HypoDMRs", 
                                  file = "Tables/LOLA chromHMM Dx Replication 50 Males HypoDMRs.csv")
chromHMM_combined <- rbind(chromHMM_hyper, chromHMM_hypo)
chromHMM_combined$userSet <- ifelse(chromHMM_combined$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))

# Plot Odds Ratio
hm.max <- quantile(chromHMM_combined$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_combined, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                 width = 9, height = 7, legend.position = c(1.11, 0.855),
                 file = "Figures/LOLA chromHMM Dx Replication 50 Males OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(chromHMM_combined$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_combined, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                 width = 9, height = 7, legend.position = c(1.12, 0.855),
                 file = "Figures/LOLA chromHMM Dx Replication 50 Males log qvalue Heatmap.png")

# Top Enriched ####
# Top enriched if regions in at least 50% of cell types are above the median of odds ratio and log(qvalue) for that DMR set
table(histone_hyper$oddsRatio >= median(histone_hyper$oddsRatio) &
              histone_hyper$qValueLog >= median(histone_hyper$qValueLog) &
              histone_hyper$qValue < 0.05,
      histone_hyper$antibody)["TRUE",] > length(unique(histone_hyper$cellType)) * 0.5
# H3K4me1  H3K4me3  H3K9me3 H3K27me3 H3K36me3 
#    TRUE     TRUE    FALSE    FALSE    FALSE 

table(histone_hypo$oddsRatio >= median(histone_hypo$oddsRatio) &
              histone_hypo$qValueLog >= median(histone_hypo$qValueLog) &
              histone_hypo$qValue < 0.05,
      histone_hypo$antibody)["TRUE",] > length(unique(histone_hypo$cellType)) * 0.5
# H3K4me1  H3K4me3  H3K9me3 H3K27me3 H3K36me3 
#   FALSE     TRUE    FALSE     TRUE    FALSE 

table(chromHMM_hyper$oddsRatio >= median(chromHMM_hyper$oddsRatio) &
              chromHMM_hyper$qValueLog >= median(chromHMM_hyper$qValueLog) &
              chromHMM_hyper$qValue < 0.05,
      chromHMM_hyper$chromState)["TRUE",] > length(unique(chromHMM_hyper$cellType)) * 0.5
#     TssA TssAFlnk   TxFlnk       Tx     TxWk     EnhG      Enh  ZnfRpts      Het   TssBiv  BivFlnk   EnhBiv   ReprPC 
#     TRUE     TRUE    FALSE    FALSE    FALSE     TRUE    FALSE    FALSE    FALSE    FALSE    FALSE     TRUE    FALSE 
# ReprPCwk    Quies 
#    FALSE    FALSE

table(chromHMM_hypo$oddsRatio >= median(chromHMM_hypo$oddsRatio) &
              chromHMM_hypo$qValueLog >= median(chromHMM_hypo$qValueLog) &
              chromHMM_hypo$qValue < 0.05,
      chromHMM_hypo$chromState)["TRUE",] > length(unique(chromHMM_hypo$cellType)) * 0.5
#     TssA TssAFlnk   TxFlnk       Tx     TxWk     EnhG      Enh  ZnfRpts      Het   TssBiv  BivFlnk   EnhBiv   ReprPC 
#     TRUE     TRUE    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE     TRUE     TRUE     TRUE 
# ReprPCwk    Quies 
#    FALSE    FALSE 

# Analyze Replication Diagnosis Females DMRs LOLA ------------------------------------
# Load Data ####
lola <- read.delim("Tables/LOLA_Dx_Replication100_females_DMRs.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Histone Mark Enrichment ####
# Prep Data
histone <- subset(lola, collection == "roadmap_epigenomics")
histone_hyper <- prepLOLAhistone(histone, index = index, regions = "HyperDMRs", 
                                 file = "Tables/LOLA Histone Dx Replication 100 Females HyperDMRs.csv")
histone_hypo <- prepLOLAhistone(histone, index = index, regions = "HypoDMRs", 
                                file = "Tables/LOLA Histone Dx Replication 100 Females HypoDMRs.csv")
histone_combined <- rbind(histone_hyper, histone_hypo)
histone_combined$userSet <- ifelse(histone_combined$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))

# Plot Odds Ratio
hm.max <- quantile(histone_combined$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_combined, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                width = 9, height = 7, legend.position = c(1.11, 0.855),
                file = "Figures/LOLA Histone Dx Replication 100 Females OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(histone_combined$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_combined, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                width = 9, height = 7, legend.position = c(1.12, 0.855),
                file = "Figures/LOLA Histone Dx Replication 100 Females log qvalue Heatmap.png")

# ChromHMM Enrichment ####
# Prep Data
chromHMM <- subset(lola, collection == "Roadmap_ChromHMM")
chromHMM_hyper <- prepLOLAchromHMM(chromHMM, index = index, regions = "HyperDMRs", 
                                   file = "Tables/LOLA chromHMM Dx Replication 100 Females HyperDMRs.csv")
chromHMM_hypo <- prepLOLAchromHMM(chromHMM, index = index, regions = "HypoDMRs", 
                                  file = "Tables/LOLA chromHMM Dx Replication 100 Females HypoDMRs.csv")
chromHMM_combined <- rbind(chromHMM_hyper, chromHMM_hypo)
chromHMM_combined$userSet <- ifelse(chromHMM_combined$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))

# Plot Odds Ratio
hm.max <- quantile(chromHMM_combined$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_combined, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                 width = 9, height = 7, legend.position = c(1.11, 0.855),
                 file = "Figures/LOLA chromHMM Dx Replication 100 Females OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(chromHMM_combined$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_combined, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                 width = 9, height = 7, legend.position = c(1.12, 0.855),
                 file = "Figures/LOLA chromHMM Dx Replication 100 Females log qvalue Heatmap.png")

# Top Enriched ####
# Top enriched if regions in at least 50% of cell types are above the median of odds ratio and log(qvalue) for that DMR set
table(histone_hyper$oddsRatio >= median(histone_hyper$oddsRatio) &
              histone_hyper$qValueLog >= median(histone_hyper$qValueLog) &
              histone_hyper$qValue < 0.05,
      histone_hyper$antibody)["TRUE",] > length(unique(histone_hyper$cellType)) * 0.5
# H3K4me1  H3K4me3  H3K9me3 H3K27me3 H3K36me3 
#    TRUE     TRUE    FALSE    FALSE    FALSE 

table(histone_hypo$oddsRatio >= median(histone_hypo$oddsRatio) &
              histone_hypo$qValueLog >= median(histone_hypo$qValueLog) &
              histone_hypo$qValue < 0.05,
      histone_hypo$antibody)["TRUE",] > length(unique(histone_hypo$cellType)) * 0.5
# H3K4me1  H3K4me3  H3K9me3 H3K27me3 H3K36me3 
#    TRUE     TRUE    FALSE     TRUE    FALSE 

table(chromHMM_hyper$oddsRatio >= median(chromHMM_hyper$oddsRatio) &
              chromHMM_hyper$qValueLog >= median(chromHMM_hyper$qValueLog) &
              chromHMM_hyper$qValue < 0.05,
      chromHMM_hyper$chromState)["TRUE",] > length(unique(chromHMM_hyper$cellType)) * 0.5
#     TssA TssAFlnk   TxFlnk       Tx     TxWk     EnhG      Enh  ZnfRpts      Het   TssBiv  BivFlnk   EnhBiv   ReprPC 
#     TRUE     TRUE    FALSE    FALSE    FALSE     TRUE    FALSE    FALSE    FALSE    FALSE     TRUE     TRUE     TRUE 
# ReprPCwk    Quies 
#    FALSE    FALSE 

table(chromHMM_hypo$oddsRatio >= median(chromHMM_hypo$oddsRatio) &
              chromHMM_hypo$qValueLog >= median(chromHMM_hypo$qValueLog) &
              chromHMM_hypo$qValue < 0.05,
      chromHMM_hypo$chromState)["TRUE",] > length(unique(chromHMM_hypo$cellType)) * 0.5
#     TssA TssAFlnk   TxFlnk       Tx     TxWk     EnhG      Enh  ZnfRpts      Het   TssBiv  BivFlnk   EnhBiv   ReprPC 
#     TRUE     TRUE    FALSE    FALSE    FALSE     TRUE    FALSE    FALSE    FALSE    FALSE     TRUE     TRUE     TRUE 
# ReprPCwk    Quies 
#    FALSE    FALSE 
