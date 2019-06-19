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

# Plot Odds Ratio
hm.max <- quantile(c(histone_hyper$oddsRatio, histone_hypo$oddsRatio), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_hyper, title = "Hypermethylated", type = "oddsRatio", hm.max = hm.max, 
                file = "Figures/LOLA Histone Dx Discovery 50 Males Hypermethylated OR Heatmap.png")
plotLOLAhistone(histone = histone_hypo, title = "Hypomethylated", type = "oddsRatio", hm.max = hm.max, 
                file = "Figures/LOLA Histone Dx Discovery 50 Males Hypomethylated OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(c(histone_hyper$qValueLog, histone_hypo$qValueLog), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_hyper, title = "Hypermethylated", type = "qValueLog", hm.max = hm.max, 
                file = "Figures/LOLA Histone Dx Discovery 50 Males Hypermethylated log qvalue Heatmap.png")
plotLOLAhistone(histone = histone_hypo, title = "Hypomethylated", type = "qValueLog", hm.max = hm.max, 
                file = "Figures/LOLA Histone Dx Discovery 50 Males Hypomethylated log qvalue Heatmap.png")

# Plot Legend
plotLOLAhistone(histone = histone_hyper, type = "legend", file = "Figures/LOLA DMRs Heatmap Legend.png")

# Plot Odds Ratio for Blood
hm.sub <- rbind(histone_hyper, histone_hypo) %>% subset(tissue %in% c("HSC & B-cell", "Blood & T-cell"))
hm.sub$cellType[grepl(pattern = "E038", x = hm.sub$filename, fixed = TRUE)] <- "Primary T helper naive cells from peripheral blood 1"
hm.sub$cellType[grepl(pattern = "E039", x = hm.sub$filename, fixed = TRUE)] <- "Primary T helper naive cells from peripheral blood 2"
hm.sub$cellType <- iconv(hm.sub$cellType, from = 'UTF-8', to = 'ASCII//TRANSLIT') %>% # remove special characters
        str_replace_all(c("Primary " = "")) %>% str_to_sentence %>% 
        str_replace_all(c(" from peripheral blood" = "", "cd" = "CD", "pma-i" = "PMA-I", "g-csf" = "G-CSF"))

hm.sub$userSet <- ifelse(hm.sub$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))
hm.max <- quantile(hm.sub$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = hm.sub, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                axis.text.y = element_text(size = 13, color = "Black"), axis.ticks.y = element_line(size = 1.25), 
                labels = unique(hm.sub$cellType), width = 12, height = 7, legend.position = c(1.13, 0.855),
                file = "Figures/LOLA Histone Blood Dx Discovery 50 Males OR Heatmap.png")

# Plot log q-value for Blood
hm.max <- quantile(hm.sub$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = hm.sub, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                axis.text.y = element_text(size = 13, color = "Black"), axis.ticks.y = element_line(size = 1.25), 
                labels = unique(hm.sub$cellType), width = 12, height = 7, legend.position = c(1.14, 0.855),
                file = "Figures/LOLA Histone Blood Dx Discovery 50 Males log qvalue Heatmap.png")

# ChromHMM Enrichment ####
# Prep Data
chromHMM <- subset(lola, collection == "Roadmap_ChromHMM")
chromHMM_hyper <- prepLOLAchromHMM(chromHMM, index = index, regions = "HyperDMRs", 
                                   file = "Tables/LOLA chromHMM Dx Discovery 50 Males HyperDMRs.csv")
chromHMM_hypo <- prepLOLAchromHMM(chromHMM, index = index, regions = "HypoDMRs", 
                                  file = "Tables/LOLA chromHMM Dx Discovery 50 Males HypoDMRs.csv")

# Plot Odds Ratio
hm.max <- quantile(c(chromHMM_hyper$oddsRatio, chromHMM_hypo$oddsRatio), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_hyper, title = "Hypermethylated", type = "oddsRatio", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Males Hypermethylated OR Heatmap.png")
plotLOLAchromHMM(chromHMM = chromHMM_hypo, title = "Hypomethylated", type = "oddsRatio", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Males Hypomethylated OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(c(chromHMM_hyper$qValueLog, chromHMM_hypo$qValueLog), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_hyper, title = "Hypermethylated", type = "qValueLog", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Males Hypermethylated log qvalue Heatmap.png")
plotLOLAchromHMM(chromHMM = chromHMM_hypo, title = "Hypomethylated", type = "qValueLog", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Males Hypomethylated log qvalue Heatmap.png")

# Plot Odds Ratio for Blood
hm.sub <- rbind(chromHMM_hyper, chromHMM_hypo) %>% subset(tissue %in% c("HSC & B-cell", "Blood & T-cell"))
hm.sub$cellType[grepl(pattern = "E038", x = hm.sub$filename, fixed = TRUE)] <- "Primary T helper naive cells from peripheral blood 1"
hm.sub$cellType[grepl(pattern = "E039", x = hm.sub$filename, fixed = TRUE)] <- "Primary T helper naive cells from peripheral blood 2"
hm.sub$cellType <- iconv(hm.sub$cellType, from = 'UTF-8', to = 'ASCII//TRANSLIT') %>% # remove special characters
        str_replace_all(c("Primary " = "")) %>% str_to_sentence %>% 
        str_replace_all(c(" from peripheral blood" = "", "cd" = "CD", "pma-i" = "PMA-I", "g-csf" = "G-CSF"))
hm.sub$userSet <- ifelse(hm.sub$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))
hm.max <- quantile(hm.sub$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = hm.sub, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                axis.text.y = element_text(size = 13, color = "Black"), axis.ticks.y = element_line(size = 1.25), 
                labels = unique(hm.sub$cellType), width = 12, height = 7, legend.position = c(1.13, 0.855),
                file = "Figures/LOLA chromHMM Blood Dx Discovery 50 Males OR Heatmap.png")

# Plot log q-value for Blood
hm.max <- quantile(hm.sub$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = hm.sub, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                axis.text.y = element_text(size = 13, color = "Black"), axis.ticks.y = element_line(size = 1.25), 
                labels = unique(hm.sub$cellType), width = 12, height = 7, legend.position = c(1.14, 0.855),
                file = "Figures/LOLA chromHMM Blood Dx Discovery 50 Males log qvalue Heatmap.png")

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

# Plot Odds Ratio
hm.max <- quantile(c(histone_hyper$oddsRatio, histone_hypo$oddsRatio), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_hyper, title = "Hypermethylated", type = "oddsRatio", hm.max = hm.max, 
                file = "Figures/LOLA Histone Dx Discovery 50 Females Hypermethylated OR Heatmap.png")
plotLOLAhistone(histone = histone_hypo, title = "Hypomethylated", type = "oddsRatio", hm.max = hm.max, 
                file = "Figures/LOLA Histone Dx Discovery 50 Females Hypomethylated OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(c(histone_hyper$qValueLog, histone_hypo$qValueLog), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_hyper, title = "Hypermethylated", type = "qValueLog", hm.max = hm.max, 
                file = "Figures/LOLA Histone Dx Discovery 50 Females Hypermethylated log qvalue Heatmap.png")
plotLOLAhistone(histone = histone_hypo, title = "Hypomethylated", type = "qValueLog", hm.max = hm.max, 
                file = "Figures/LOLA Histone Dx Discovery 50 Females Hypomethylated log qvalue Heatmap.png")

# Plot Odds Ratio for Blood
hm.sub <- rbind(histone_hyper, histone_hypo) %>% subset(tissue %in% c("HSC & B-cell", "Blood & T-cell"))
hm.sub$cellType[grepl(pattern = "E038", x = hm.sub$filename, fixed = TRUE)] <- "Primary T helper naive cells from peripheral blood 1"
hm.sub$cellType[grepl(pattern = "E039", x = hm.sub$filename, fixed = TRUE)] <- "Primary T helper naive cells from peripheral blood 2"
hm.sub$cellType <- iconv(hm.sub$cellType, from = 'UTF-8', to = 'ASCII//TRANSLIT') %>% # remove special characters
        str_replace_all(c("Primary " = "")) %>% str_to_sentence %>% 
        str_replace_all(c(" from peripheral blood" = "", "cd" = "CD", "pma-i" = "PMA-I", "g-csf" = "G-CSF"))

hm.sub$userSet <- ifelse(hm.sub$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))
hm.max <- quantile(hm.sub$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = hm.sub, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                axis.text.y = element_text(size = 13, color = "Black"), axis.ticks.y = element_line(size = 1.25), 
                labels = unique(hm.sub$cellType), width = 12, height = 7, legend.position = c(1.13, 0.855),
                file = "Figures/LOLA Histone Blood Dx Discovery 50 Females OR Heatmap.png")

# Plot log q-value for Blood
hm.max <- quantile(hm.sub$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = hm.sub, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                axis.text.y = element_text(size = 13, color = "Black"), axis.ticks.y = element_line(size = 1.25), 
                labels = unique(hm.sub$cellType), width = 12, height = 7, legend.position = c(1.14, 0.855),
                file = "Figures/LOLA Histone Blood Dx Discovery 50 Females log qvalue Heatmap.png")

# ChromHMM Enrichment ####
# Prep Data
chromHMM <- subset(lola, collection == "Roadmap_ChromHMM")
chromHMM_hyper <- prepLOLAchromHMM(chromHMM, index = index, regions = "HyperDMRs", 
                                   file = "Tables/LOLA chromHMM Dx Discovery 50 Females HyperDMRs.csv")
chromHMM_hypo <- prepLOLAchromHMM(chromHMM, index = index, regions = "HypoDMRs", 
                                  file = "Tables/LOLA chromHMM Dx Discovery 50 Females HypoDMRs.csv")

# Plot Odds Ratio
hm.max <- quantile(c(chromHMM_hyper$oddsRatio, chromHMM_hypo$oddsRatio), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_hyper, title = "Hypermethylated", type = "oddsRatio", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Females Hypermethylated OR Heatmap.png")
plotLOLAchromHMM(chromHMM = chromHMM_hypo, title = "Hypomethylated", type = "oddsRatio", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Females Hypomethylated OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(c(chromHMM_hyper$qValueLog, chromHMM_hypo$qValueLog), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_hyper, title = "Hypermethylated", type = "qValueLog", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Females Hypermethylated log qvalue Heatmap.png")
plotLOLAchromHMM(chromHMM = chromHMM_hypo, title = "Hypomethylated", type = "qValueLog", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM Dx Discovery 50 Females Hypomethylated log qvalue Heatmap.png")

# Plot Odds Ratio for Blood
hm.sub <- rbind(chromHMM_hyper, chromHMM_hypo) %>% subset(tissue %in% c("HSC & B-cell", "Blood & T-cell"))
hm.sub$cellType[grepl(pattern = "E038", x = hm.sub$filename, fixed = TRUE)] <- "Primary T helper naive cells from peripheral blood 1"
hm.sub$cellType[grepl(pattern = "E039", x = hm.sub$filename, fixed = TRUE)] <- "Primary T helper naive cells from peripheral blood 2"
hm.sub$cellType <- iconv(hm.sub$cellType, from = 'UTF-8', to = 'ASCII//TRANSLIT') %>% # remove special characters
        str_replace_all(c("Primary " = "")) %>% str_to_sentence %>% 
        str_replace_all(c(" from peripheral blood" = "", "cd" = "CD", "pma-i" = "PMA-I", "g-csf" = "G-CSF"))
hm.sub$userSet <- ifelse(hm.sub$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))
hm.max <- quantile(hm.sub$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = hm.sub, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                 axis.text.y = element_text(size = 13, color = "Black"), axis.ticks.y = element_line(size = 1.25), 
                 labels = unique(hm.sub$cellType), width = 12, height = 7, legend.position = c(1.13, 0.855),
                 file = "Figures/LOLA chromHMM Blood Dx Discovery 50 Females OR Heatmap.png")

# Plot log q-value for Blood
hm.max <- quantile(hm.sub$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = hm.sub, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                 axis.text.y = element_text(size = 13, color = "Black"), axis.ticks.y = element_line(size = 1.25), 
                 labels = unique(hm.sub$cellType), width = 12, height = 7, legend.position = c(1.14, 0.855),
                 file = "Figures/LOLA chromHMM Blood Dx Discovery 50 Females log qvalue Heatmap.png")

# Analyze Replication Diagnosis Males DMRs LOLA ------------------------------------
# Load Data ####
index <- read.delim("Tables/LOLA Roadmap ChromHMM index.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Histone Mark Enrichment ####
# Prep Data
histone <- subset(lola, collection == "roadmap_epigenomics")
histone_hyper <- prepLOLAhistone(histone, index = index, regions = "HyperDMRs", 
                                 file = "Tables/LOLA Histone Dx Replication 50 Males HyperDMRs.csv")
histone_hypo <- prepLOLAhistone(histone, index = index, regions = "HypoDMRs", 
                                file = "Tables/LOLA Histone Dx Replication 50 Males HypoDMRs.csv")

# Plot Odds Ratio
hm.max <- quantile(c(histone_hyper$oddsRatio, histone_hypo$oddsRatio), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_hyper, title = "Hypermethylated", type = "oddsRatio", hm.max = hm.max, 
                file = "Figures/LOLA Histone Dx Replication 50 Males Hypermethylated OR Heatmap.png")
plotLOLAhistone(histone = histone_hypo, title = "Hypomethylated", type = "oddsRatio", hm.max = hm.max, 
                file = "Figures/LOLA Histone Dx Replication 50 Males Hypomethylated OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(c(histone_hyper$qValueLog, histone_hypo$qValueLog), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_hyper, title = "Hypermethylated", type = "qValueLog", hm.max = hm.max, 
                file = "Figures/LOLA Histone Dx Replication 50 Males Hypermethylated log qvalue Heatmap.png")
plotLOLAhistone(histone = histone_hypo, title = "Hypomethylated", type = "qValueLog", hm.max = hm.max, 
                file = "Figures/LOLA Histone Dx Replication 50 Males Hypomethylated log qvalue Heatmap.png")

# Plot Odds Ratio for Blood
hm.sub <- rbind(histone_hyper, histone_hypo) %>% subset(tissue %in% c("HSC & B-cell", "Blood & T-cell"))
hm.sub$cellType[grepl(pattern = "E038", x = hm.sub$filename, fixed = TRUE)] <- "Primary T helper naive cells from peripheral blood 1"
hm.sub$cellType[grepl(pattern = "E039", x = hm.sub$filename, fixed = TRUE)] <- "Primary T helper naive cells from peripheral blood 2"
hm.sub$cellType <- iconv(hm.sub$cellType, from = 'UTF-8', to = 'ASCII//TRANSLIT') %>% # remove special characters
        str_replace_all(c("Primary " = "")) %>% str_to_sentence %>% 
        str_replace_all(c(" from peripheral blood" = "", "cd" = "CD", "pma-i" = "PMA-I", "g-csf" = "G-CSF"))

hm.sub$userSet <- ifelse(hm.sub$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))
hm.max <- quantile(hm.sub$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = hm.sub, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                axis.text.y = element_text(size = 13, color = "Black"), axis.ticks.y = element_line(size = 1.25), 
                labels = unique(hm.sub$cellType), width = 12, height = 7, legend.position = c(1.13, 0.855),
                file = "Figures/LOLA Histone Blood Dx Replication 50 Males OR Heatmap.png")

# Plot log q-value for Blood
hm.max <- quantile(hm.sub$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = hm.sub, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                axis.text.y = element_text(size = 13, color = "Black"), axis.ticks.y = element_line(size = 1.25), 
                labels = unique(hm.sub$cellType), width = 12, height = 7, legend.position = c(1.14, 0.855),
                file = "Figures/LOLA Histone Blood Dx Replication 50 Males log qvalue Heatmap.png")

# ChromHMM Enrichment ####
# Prep Data
chromHMM <- subset(lola, collection == "Roadmap_ChromHMM")
chromHMM_hyper <- prepLOLAchromHMM(chromHMM, index = index, regions = "HyperDMRs", 
                                   file = "Tables/LOLA chromHMM Dx Replication 50 Males HyperDMRs.csv")
chromHMM_hypo <- prepLOLAchromHMM(chromHMM, index = index, regions = "HypoDMRs", 
                                  file = "Tables/LOLA chromHMM Dx Replication 50 Males HypoDMRs.csv")

# Plot Odds Ratio
hm.max <- quantile(c(chromHMM_hyper$oddsRatio, chromHMM_hypo$oddsRatio), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_hyper, title = "Hypermethylated", type = "oddsRatio", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM Dx Replication 50 Males Hypermethylated OR Heatmap.png")
plotLOLAchromHMM(chromHMM = chromHMM_hypo, title = "Hypomethylated", type = "oddsRatio", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM Dx Replication 50 Males Hypomethylated OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(c(chromHMM_hyper$qValueLog, chromHMM_hypo$qValueLog), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_hyper, title = "Hypermethylated", type = "qValueLog", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM Dx Replication 50 Males Hypermethylated log qvalue Heatmap.png")
plotLOLAchromHMM(chromHMM = chromHMM_hypo, title = "Hypomethylated", type = "qValueLog", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM Dx Replication 50 Males Hypomethylated log qvalue Heatmap.png")

# Plot Odds Ratio for Blood
hm.sub <- rbind(chromHMM_hyper, chromHMM_hypo) %>% subset(tissue %in% c("HSC & B-cell", "Blood & T-cell"))
hm.sub$cellType[grepl(pattern = "E038", x = hm.sub$filename, fixed = TRUE)] <- "Primary T helper naive cells from peripheral blood 1"
hm.sub$cellType[grepl(pattern = "E039", x = hm.sub$filename, fixed = TRUE)] <- "Primary T helper naive cells from peripheral blood 2"
hm.sub$cellType <- iconv(hm.sub$cellType, from = 'UTF-8', to = 'ASCII//TRANSLIT') %>% # remove special characters
        str_replace_all(c("Primary " = "")) %>% str_to_sentence %>% 
        str_replace_all(c(" from peripheral blood" = "", "cd" = "CD", "pma-i" = "PMA-I", "g-csf" = "G-CSF"))
hm.sub$userSet <- ifelse(hm.sub$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))
hm.max <- quantile(hm.sub$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = hm.sub, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                 axis.text.y = element_text(size = 13, color = "Black"), axis.ticks.y = element_line(size = 1.25), 
                 labels = unique(hm.sub$cellType), width = 12, height = 7, legend.position = c(1.13, 0.855),
                 file = "Figures/LOLA chromHMM Blood Dx Replication 50 Males OR Heatmap.png")

# Plot log q-value for Blood
hm.max <- quantile(hm.sub$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = hm.sub, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                 axis.text.y = element_text(size = 13, color = "Black"), axis.ticks.y = element_line(size = 1.25), 
                 labels = unique(hm.sub$cellType), width = 12, height = 7, legend.position = c(1.14, 0.855),
                 file = "Figures/LOLA chromHMM Blood Dx Replication 50 Males log qvalue Heatmap.png")

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

# Plot Odds Ratio
hm.max <- quantile(c(histone_hyper$oddsRatio, histone_hypo$oddsRatio), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_hyper, title = "Hypermethylated", type = "oddsRatio", hm.max = hm.max, 
                file = "Figures/LOLA Histone Dx Replication 100 Females Hypermethylated OR Heatmap.png")
plotLOLAhistone(histone = histone_hypo, title = "Hypomethylated", type = "oddsRatio", hm.max = hm.max, 
                file = "Figures/LOLA Histone Dx Replication 100 Females Hypomethylated OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(c(histone_hyper$qValueLog, histone_hypo$qValueLog), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = histone_hyper, title = "Hypermethylated", type = "qValueLog", hm.max = hm.max, 
                file = "Figures/LOLA Histone Dx Replication 100 Females Hypermethylated log qvalue Heatmap.png")
plotLOLAhistone(histone = histone_hypo, title = "Hypomethylated", type = "qValueLog", hm.max = hm.max, 
                file = "Figures/LOLA Histone Dx Replication 100 Females Hypomethylated log qvalue Heatmap.png")

# Plot Odds Ratio for Blood
hm.sub <- rbind(histone_hyper, histone_hypo) %>% subset(tissue %in% c("HSC & B-cell", "Blood & T-cell"))
hm.sub$cellType[grepl(pattern = "E038", x = hm.sub$filename, fixed = TRUE)] <- "Primary T helper naive cells from peripheral blood 1"
hm.sub$cellType[grepl(pattern = "E039", x = hm.sub$filename, fixed = TRUE)] <- "Primary T helper naive cells from peripheral blood 2"
hm.sub$cellType <- iconv(hm.sub$cellType, from = 'UTF-8', to = 'ASCII//TRANSLIT') %>% # remove special characters
        str_replace_all(c("Primary " = "")) %>% str_to_sentence %>% 
        str_replace_all(c(" from peripheral blood" = "", "cd" = "CD", "pma-i" = "PMA-I", "g-csf" = "G-CSF"))

hm.sub$userSet <- ifelse(hm.sub$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))
hm.max <- quantile(hm.sub$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = hm.sub, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                axis.text.y = element_text(size = 13, color = "Black"), axis.ticks.y = element_line(size = 1.25), 
                labels = unique(hm.sub$cellType), width = 12, height = 7, legend.position = c(1.13, 0.855),
                file = "Figures/LOLA Histone Blood Dx Replication 100 Females OR Heatmap.png")

# Plot log q-value for Blood
hm.max <- quantile(hm.sub$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAhistone(histone = hm.sub, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                axis.text.y = element_text(size = 13, color = "Black"), axis.ticks.y = element_line(size = 1.25), 
                labels = unique(hm.sub$cellType), width = 12, height = 7, legend.position = c(1.14, 0.855),
                file = "Figures/LOLA Histone Blood Dx Replication 100 Females log qvalue Heatmap.png")

# ChromHMM Enrichment ####
# Prep Data
chromHMM <- subset(lola, collection == "Roadmap_ChromHMM")
chromHMM_hyper <- prepLOLAchromHMM(chromHMM, index = index, regions = "HyperDMRs", 
                                   file = "Tables/LOLA chromHMM Dx Replication 100 Females HyperDMRs.csv")
chromHMM_hypo <- prepLOLAchromHMM(chromHMM, index = index, regions = "HypoDMRs", 
                                  file = "Tables/LOLA chromHMM Dx Replication 100 Females HypoDMRs.csv")

# Plot Odds Ratio
hm.max <- quantile(c(chromHMM_hyper$oddsRatio, chromHMM_hypo$oddsRatio), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_hyper, title = "Hypermethylated", type = "oddsRatio", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM Dx Replication 100 Females Hypermethylated OR Heatmap.png")
plotLOLAchromHMM(chromHMM = chromHMM_hypo, title = "Hypomethylated", type = "oddsRatio", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM Dx Replication 100 Females Hypomethylated OR Heatmap.png")

# Plot log q-value
hm.max <- quantile(c(chromHMM_hyper$qValueLog, chromHMM_hypo$qValueLog), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = chromHMM_hyper, title = "Hypermethylated", type = "qValueLog", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM Dx Replication 100 Females Hypermethylated log qvalue Heatmap.png")
plotLOLAchromHMM(chromHMM = chromHMM_hypo, title = "Hypomethylated", type = "qValueLog", hm.max = hm.max, 
                 file = "Figures/LOLA chromHMM Dx Replication 100 Females Hypomethylated log qvalue Heatmap.png")

# Plot Odds Ratio for Blood
hm.sub <- rbind(chromHMM_hyper, chromHMM_hypo) %>% subset(tissue %in% c("HSC & B-cell", "Blood & T-cell"))
hm.sub$cellType[grepl(pattern = "E038", x = hm.sub$filename, fixed = TRUE)] <- "Primary T helper naive cells from peripheral blood 1"
hm.sub$cellType[grepl(pattern = "E039", x = hm.sub$filename, fixed = TRUE)] <- "Primary T helper naive cells from peripheral blood 2"
hm.sub$cellType <- iconv(hm.sub$cellType, from = 'UTF-8', to = 'ASCII//TRANSLIT') %>% # remove special characters
        str_replace_all(c("Primary " = "")) %>% str_to_sentence %>% 
        str_replace_all(c(" from peripheral blood" = "", "cd" = "CD", "pma-i" = "PMA-I", "g-csf" = "G-CSF"))
hm.sub$userSet <- ifelse(hm.sub$userSet == "HyperDMRs", yes = "Hypermethylated", no = "Hypomethylated") %>% 
        factor(levels = c("Hypermethylated", "Hypomethylated"))
hm.max <- quantile(hm.sub$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = hm.sub, facet = vars(userSet), title = NULL, type = "oddsRatio", hm.max = hm.max, 
                 axis.text.y = element_text(size = 13, color = "Black"), axis.ticks.y = element_line(size = 1.25), 
                 labels = unique(hm.sub$cellType), width = 12, height = 7, legend.position = c(1.13, 0.855),
                 file = "Figures/LOLA chromHMM Blood Dx Replication 100 Females OR Heatmap.png")

# Plot log q-value for Blood
hm.max <- quantile(hm.sub$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
plotLOLAchromHMM(chromHMM = hm.sub, facet = vars(userSet), title = NULL, type = "qValueLog", hm.max = hm.max, 
                 axis.text.y = element_text(size = 13, color = "Black"), axis.ticks.y = element_line(size = 1.25), 
                 labels = unique(hm.sub$cellType), width = 12, height = 7, legend.position = c(1.14, 0.855),
                 file = "Figures/LOLA chromHMM Blood Dx Replication 100 Females log qvalue Heatmap.png")


