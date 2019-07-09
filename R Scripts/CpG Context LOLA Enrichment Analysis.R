# CpG Context LOLA Enrichment Analysis ---------------------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 7/8/19

# Packages ####
.libPaths("/share/lasallelab/Charles/R")
sapply(c("tidyverse", "LOLA", "simpleCache", "GenomicRanges", "qvalue", "annotatr", "scales"), require, character.only = TRUE)

# Functions ####
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

# Get CpG Context Annotation ####
CpGs <- build_annotations(genome = "hg38", annotations = "hg38_cpgs") %>%
        GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") %>% as.data.frame()
CpGs <- CpGs[,c("seqnames", "start", "end", "type")]
CpGs$type <- str_replace_all(CpGs$type, pattern = c("hg38_cpg_islands" = "CpG_Island", "hg38_cpg_shores" = "CpG_Shore",
                                                    "hg38_cpg_shelves" = "CpG_Shelf", "hg38_cpg_inter" = "CpG_Open_Sea"))
write.table(CpGs, "UCSC Tracks/CpG_Context.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
splitFileIntoCollection("UCSC Tracks/CpG_Context.bed", splitCol = 4, collectionFolder = "UCSC Tracks")
index <- data.frame(filename = c("CpG_Island.bed", "CpG_Shore.bed", "CpG_Shelf.bed", "CpG_Open_Sea.bed"),
                    description = c("CpG Island", "CpG Shore", "CpG Shelf", "CpG Open Sea"))
write.table(index, "Tables/index.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
rm(CpGs, index)

# Get LOLA Enrichments ----------------------------------------------------
# Data ####
regionDB <- loadRegionDB(dbLocation = "/share/lasallelab/programs/LOLA/hg38", useCache = TRUE, limit = NULL, 
                         collections = c("CpG_context"))
maleChroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")
femaleChroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrM")

# Discovery Males DMRs ####
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Discovery/")
DMRs <- loadRegions(file = "Dx_Males/DMRs_Dx_Discovery50_males.csv", chroms = maleChroms, sort = TRUE)
DMRlist <- list("AllDMRs" = makeGRange(DMRs = DMRs, direction = "all"),
                "HyperDMRs" = makeGRange(DMRs = DMRs, direction = "hyper"),
                "HypoDMRs" = makeGRange(DMRs = DMRs, direction = "hypo"))
Background <- loadRegions(file = "Dx_Males/bsseq_background_Discovery50_males.csv", chroms = maleChroms, sort = TRUE) %>% 
        makeGRange(direction = "all")
Results <- runLOLA(userSets = DMRlist, userUniverse = Background, regionDB = regionDB, cores = 2, redefineUserSets = TRUE)
writeCombinedEnrichment(combinedResults = Results, outFolder = "LOLA", includeSplits = FALSE)
file.copy(from = "LOLA/allEnrichments.tsv", to = "LOLA/LOLA_CpG_Context_Dx_Discovery50_males_DMRs.tsv", overwrite = TRUE)
rm(DMRs, DMRlist, Background, Results)

# Discovery Females DMRs ####
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Discovery/")
DMRs <- loadRegions(file = "Dx_Females/DMRs_Dx_Discovery50_females.csv", chroms = femaleChroms, sort = TRUE)
DMRlist <- list("AllDMRs" = makeGRange(DMRs = DMRs, direction = "all"),
                "HyperDMRs" = makeGRange(DMRs = DMRs, direction = "hyper"),
                "HypoDMRs" = makeGRange(DMRs = DMRs, direction = "hypo"))
Background <- loadRegions(file = "Dx_Females/bsseq_background_Discovery50_females.csv", chroms = femaleChroms, sort = TRUE) %>% 
        makeGRange(direction = "all")
Results <- runLOLA(userSets = DMRlist, userUniverse = Background, regionDB = regionDB, cores = 2, redefineUserSets = TRUE)
writeCombinedEnrichment(combinedResults = Results, outFolder = "LOLA", includeSplits = FALSE)
file.copy(from = "LOLA/allEnrichments.tsv", to = "LOLA/LOLA_CpG_Context_Dx_Discovery50_females_DMRs.tsv", overwrite = TRUE)
rm(DMRs, DMRlist, Background, Results)

# Replication Males DMRs ####
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Replication/")
DMRs <- loadRegions(file = "Dx_Males/DMRs_Dx_Replication50_males.csv", chroms = maleChroms, sort = TRUE)
DMRlist <- list("AllDMRs" = makeGRange(DMRs = DMRs, direction = "all"),
                "HyperDMRs" = makeGRange(DMRs = DMRs, direction = "hyper"),
                "HypoDMRs" = makeGRange(DMRs = DMRs, direction = "hypo"))
Background <- loadRegions(file = "Dx_Males/bsseq_background_Replication50_males.csv", chroms = maleChroms, sort = TRUE) %>% 
        makeGRange(direction = "all")
Results <- runLOLA(userSets = DMRlist, userUniverse = Background, regionDB = regionDB, cores = 2, redefineUserSets = TRUE)
writeCombinedEnrichment(combinedResults = Results, outFolder = "LOLA", includeSplits = FALSE)
file.copy(from = "LOLA/allEnrichments.tsv", to = "LOLA/LOLA_CpG_Context_Dx_Replication50_males_DMRs.tsv", overwrite = TRUE)
rm(DMRs, DMRlist, Background, Results)

# Replication Females DMRs ####
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Replication/")
DMRs <- loadRegions(file = "Dx_Females_100/DMRs_Dx_Replication100_females.csv", chroms = femaleChroms, sort = TRUE)
DMRlist <- list("AllDMRs" = makeGRange(DMRs = DMRs, direction = "all"),
                "HyperDMRs" = makeGRange(DMRs = DMRs, direction = "hyper"),
                "HypoDMRs" = makeGRange(DMRs = DMRs, direction = "hypo"))
Background <- loadRegions(file = "Dx_Females_100/bsseq_background_Replication100_females.csv", chroms = femaleChroms, sort = TRUE) %>% 
        makeGRange(direction = "all")
Results <- runLOLA(userSets = DMRlist, userUniverse = Background, regionDB = regionDB, cores = 2, redefineUserSets = TRUE)
writeCombinedEnrichment(combinedResults = Results, outFolder = "LOLA", includeSplits = FALSE)
file.copy(from = "LOLA/allEnrichments.tsv", to = "LOLA/LOLA_CpG_Context_Dx_Replication100_females_DMRs.tsv", overwrite = TRUE)
rm(DMRs, DMRlist, Background, Results)

# Analyze LOLA Enrichments -----------------------------------------
# Load Data and Combine ####
lola <- rbind(read.delim("Tables/LOLA_CpG_Context_Dx_Discovery50_males_DMRs.tsv", sep = "\t", header = TRUE, 
                         stringsAsFactors = FALSE),
              read.delim("Tables/LOLA_CpG_Context_Dx_Discovery50_females_DMRs.tsv", sep = "\t", header = TRUE, 
                         stringsAsFactors = FALSE),
              read.delim("Tables/LOLA_CpG_Context_Dx_Replication50_males_DMRs.tsv", sep = "\t", header = TRUE, 
                         stringsAsFactors = FALSE),
              read.delim("Tables/LOLA_CpG_Context_Dx_Replication100_females_DMRs.tsv", sep = "\t", header = TRUE, 
                         stringsAsFactors = FALSE))
lola$DMRs <- c(rep(c("Discovery\nMales", "Discovery\nFemales", "Replication\nMales", "Replication\nFemales"), each = 12)) %>%
        factor(levels = c("Discovery\nMales", "Replication\nMales", "Discovery\nFemales", "Replication\nFemales"))
table(lola$support < 5) # All FALSE
lola$pct_DMRs <- lola$support * 100 / (lola$support + lola$c)
lola$pValueLog[is.infinite(lola$pValueLog)] <- NA
lola$pValueLog[is.na(lola$pValueLog)] <- 1.5 * max(lola$pValueLog, na.rm = TRUE)
lola$pValue <- 10^(-lola$pValueLog)
lola$qValueLog <- -log10(lola$qValue)
lola$qValueLog[is.infinite(lola$qValueLog)] <- NA
lola$qValueLog[is.na(lola$qValueLog)] <- 1.5 * max(lola$qValueLog, na.rm = TRUE)
lola <- lola[,c("DMRs", "userSet", "description", "pValue", "qValue", "pValueLog", "qValueLog", "oddsRatio", "support", "pct_DMRs", 
                "b", "c", "d", "size")]
lola$userSet <- factor(lola$userSet, levels = c("AllDMRs", "HyperDMRs", "HypoDMRs"))
lola$description <- gsub("CpG ", replacement = "", x = lola$description, fixed = TRUE) %>% 
        factor(levels = rev(c("Island", "Shore", "Shelf", "Open Sea")))
lola$Significant <- (lola$qValue < 0.05) %>% factor(levels = c("TRUE", "FALSE"))
lola <- lola[order(lola$DMRs, lola$userSet, lola$description),]
write.csv(lola, "Tables/LOLA CpG Context Enrichment Results.csv", row.names = FALSE)

# Plot Odds Ratio Heatmap ####
lola <- subset(lola, !userSet == "AllDMRs")
gg <- ggplot(data = lola)
gg +
        geom_tile(aes(x = userSet, y = description, fill = oddsRatio)) +
        geom_text(aes(x = userSet, y = description, alpha = Significant), label = "*", color = "white", size = 14, nudge_y = -0.2) +
        facet_grid(cols = vars(DMRs)) +
        scale_fill_gradientn("Odds Ratio", colors = c("black", "#FF0000"), values = c(0, 1), 
                             na.value = "#FF0000", limits = c(0, max(lola$oddsRatio)), breaks = pretty_breaks(n = 3)) +
        scale_alpha_manual(breaks = c("TRUE", "FALSE"), values = c(1, 0), guide = FALSE) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              panel.background = element_rect(fill = "black"), axis.ticks.x = element_line(size = 1.25), 
              axis.ticks.y = element_line(size = 1.25), legend.key = element_blank(),  legend.position = c(1.115, 0.86), 
              legend.background = element_blank(), legend.title = element_text(size = 15), legend.text = element_text(size = 14),
              plot.margin = unit(c(0, 6, 0.5, 0.5), "lines"), axis.text.y = element_text(size = 14, color = "black"), 
              axis.text.x = element_text(size = 14, color = "black", angle = 90, hjust = 1, vjust = 0.5), 
              axis.title = element_blank(), plot.title = element_text(size = 18, hjust = 0.5, vjust = 0),
              strip.background = element_blank(), strip.text = element_text(size = 15), legend.key.size = unit(1, "lines")) +
        scale_x_discrete(expand = c(0, 0), labels = c("Hyper", "Hypo")) +
        scale_y_discrete(expand = c(0, 0))
ggsave("Figures/LOLA CpG Context Enrichment Odds Ratio Heatmap.png", dpi = 600, width = 8, height = 4, units = "in")



