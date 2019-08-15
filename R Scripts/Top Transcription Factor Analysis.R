# Top Transcription Factor Analysis ---------------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 8/9/19

# Packages ####
.libPaths("/share/lasallelab/Charles/R")
sapply(c("tidyverse", "GenomicRanges", "LOLA", "reshape2", "qvalue", "simpleCache", "rlist", "scales", "data.table"), require, character.only = TRUE)

# Functions ####
source("R Scripts/DMR Analysis Functions.R")

# Get TF Expression in Fetal Brain ----------------------------------------
# Data ####
TFs <- read.csv("Tables/HOMER Top Enriched TF Motif Info.csv", header = TRUE, stringsAsFactors = FALSE)
brainExp <- read.csv("Tables/Allen_developmental_transcriptome_expression_matrix.csv", header = FALSE, stringsAsFactors = FALSE)
brainGenes <- read.csv("Tables/Allen_developmental_transcriptome_rows_metadata.csv", header = TRUE, stringsAsFactors = FALSE)
brainSamples <- read.csv("Tables/Allen_developmental_transcriptome_columns_metadata.csv", header = TRUE, stringsAsFactors = FALSE)

# Subset Samples in Brain Expression Matrix ####
brainExp <- brainExp[,2:ncol(brainExp)]
sampleCols <- brainSamples$column_num[brainSamples$structure_acronym == "DFC" & brainSamples$age == "13 pcw"] # 88  96 107
brainSamples_sub <- brainSamples[sampleCols,]
#     column_num donor_id   donor_name    age gender structure_id structure_acronym                 structure_name
# 88          88    12820 H376.IIIA.50 13 pcw      M        10173               DFC dorsolateral prefrontal cortex
# 96          96    12834 H376.IIIA.51 13 pcw      F        10173               DFC dorsolateral prefrontal cortex
# 107        107    12888 H376.IIIA.52 13 pcw      M        10173               DFC dorsolateral prefrontal cortex

brainExp <- brainExp[,sampleCols]
colnames(brainExp) <- apply(brainSamples_sub[,c("gender", "structure_acronym", "age", "donor_id")], 1, paste, collapse = "_") %>%
        gsub(pattern = " ", replacement = "", fixed = TRUE) %>% paste(., "RPKM", sep = "_")
cor(brainExp$M_DFC_13pcw_12820_RPKM, brainExp$F_DFC_13pcw_12834_RPKM) # 0.9375689 1st male is most similar to female
cor(brainExp$M_DFC_13pcw_12820_RPKM, brainExp$M_DFC_13pcw_12888_RPKM) # 0.7725906 
cor(brainExp$F_DFC_13pcw_12834_RPKM, brainExp$M_DFC_13pcw_12888_RPKM) # 0.8378264 
brainExp <- brainExp[, c("M_DFC_13pcw_12820_RPKM", "F_DFC_13pcw_12834_RPKM")] # Exclude 2nd male
brainExp$Ensembl_ID <- brainGenes$ensembl_gene_id

# Merge with TF Info ####
TFs <- TFs[!TFs$Motif == "", colnames(TFs)[!colnames(TFs) %in% c("X", "X.1")]]
TFs <- merge(x = TFs, y = brainExp, by = "Ensembl_ID", all.x = TRUE, all.y = FALSE, sort = FALSE)
TFs <- TFs[order(TFs$Motif),]
rm(brainExp, brainGenes, brainSamples, brainSamples_sub, sampleCols)
write.csv(TFs, "Tables/HOMER Top Enriched TF Motif Info with Brain Expression.csv", row.names = FALSE)

# Motif ChromHMM Enrichment -----------------------------------------------
# Get Motif Locations from HOMER (Epigenerate) ####
# See HOMER_TF_scanMotifGenomeWide.sh
# BED files in hg38
# Background?

# Get Enrichment with LOLA (Epigenerate) ####
# Functions
makeGRange <- function(DMRs, direction = c("all", "hyper", "hypo")){
         if(direction == "hyper"){DMRs <- subset(DMRs, percentDifference > 0)}
         if(direction == "hypo"){DMRs <- subset(DMRs, percentDifference < 0)}
         GR <- GRanges(seqnames = DMRs$chr, ranges = IRanges(start = DMRs$start, end = DMRs$end))
}

# Load Region Data
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/HOMER")
regionDB <- loadRegionDB(dbLocation = "/share/lasallelab/programs/LOLA/hg38", useCache = TRUE, limit = NULL, 
                         collections = "Roadmap_ChromHMM")
files <- list.files("BED") %>% paste("BED/", ., sep = "")
motifs <- lapply(files, fread, sep = "\t", header = FALSE, stringsAsFactors = FALSE, verbose = FALSE, data.table = FALSE,
                 nThread = 48)
names(motifs) <- str_replace_all(files, pattern = c("BED/" = "", "_motif_locations.bed" = ""))
motifs <- lapply(motifs, function(x){
        x <- x[,1:3]
        colnames(x) <- c("chr", "start", "end")
        x <- subset(x, chr %in% c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")) %>% makeGRange(direction = "all")
        return(x)
})
Background <- fread(file = "homer_knownMotifs_hg38.bed", sep = "\t", header = FALSE, stringsAsFactors = FALSE, 
                    verbose = FALSE, data.table = FALSE, nThread = 48) # data.table::fread() for much faster reading of 23GB file
Background <- Background[,1:3]
colnames(Background) <- c("chr", "start", "end")
Background <- subset(Background, chr %in% c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")) %>%
        makeGRange(direction = "all") %>% reduce() # Reduce regions to remove overlaps and merge adjacent sites
length(Background) # 88 854 866
summary(width(Background))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 8.00    11.00    16.00    20.34    25.00 20564.00 
rm(files, makeGRange)

# Run LOLA
Results <- runLOLA(userSets = motifs, userUniverse = Background, regionDB = regionDB, cores = 1, redefineUserSets = TRUE)
writeCombinedEnrichment(combinedResults = Results, outFolder = "LOLA", includeSplits = FALSE)
file.copy(from = "LOLA/allEnrichments.tsv", to = "LOLA/Top_TF_Motif_LOLA_ChromHMM_Enrichments.tsv", overwrite = TRUE)

# Enrichment Design ####
# support is the number of regions in the userSet that overlap a region in the testSet
# b is the number of regions in the userUniverse that overlap a region in the testSet, 
#     minus the number of regions in the userSet that overlap a region in the testSet
# c is the number of regions in the userSet that don't overlap a region in the testSet
# d is the number of regions in the userUniverse, 
#     minus the number of regions in the userUniverse that overlap a region in the testSet,
#     minus the number of regions in the userSet that don't overlap a region in the testSet
# Assumes all regions in the userSet are in the userUniverse
# Also assumes all regions in the userUniverse only overlap one region in the userSet, that they're comparable
# Multiple regions in the userSet overlapping one userUniverse region would decrease b, and decrease d
# userUniverse must be comparable to userSet

#' userSet: TF motif locations
#' testSet: chromHMM chromatin state segmentations
#' userUniverse: restricted universe (regions covered in at least one set, combines all userSets and disjoins them to fix overlapping)
#' 
#' Example:
#' userSet: ARE motif locations
#' testSet: Enhancer regions in male fetal brain
#' userUniverse: All motif locations of these 37 factors
#' Redefine ARE motif locations in terms of restricted universe
#' Question: Do ARE motifs overlap enhancers in male fetal brain more than would be expected for top motifs?
#' Would miss a common enrichment among all top motifs...
#' Use locations for all known motifs to create a restricted universe, removes bias for homer tfs








# Analysis ####
# Load Results
lola <- read.delim("Tables/Top_TF_Motif_LOLA_ChromHMM_Enrichments.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
index <- read.delim("Tables/LOLA Roadmap ChromHMM index.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
chromHMM <- split(lola, f = lola$userSets) %>%
        lapply(function(x) prepLOLAchromHMM(x, index = index, regions = unique(x$userSet)))
names(chromHMM) <- unique(lola$userSet)
chromHMM <- list.rbind(chromHMM)
write.csv(chromHMM, file = "Tables/Top TF Motif LOLA ChromHMM Enrichments.csv", quote = FALSE, row.names = FALSE)

# Subset Fetal Brain Enrichments
brain <- subset(chromHMM, cellType %in% c("Fetal Brain Male", "Fetal Brain Female"))

# Get Top Ranked
motifs <- unique(lola$userSet)
top_chromStates <- NULL
for(i in 1:length(motifs)){
        sub <- subset(lola, userSet == motifs[i] & cellType == "Fetal Brain Male")
        maleTop <- sub$chromState[sub$maxRnk == min(sub$maxRnk)]
        sub <- subset(lola, userSet == motifs[i] & cellType == "Fetal Brain Female")
        femaleTop <- sub$chromState[sub$maxRnk == min(sub$maxRnk)]
        temp <- c(motifs[i], maleTop, femaleTop)
        top_chromStates <- rbind(top_chromStates, temp)
}
colnames(top_chromStates) <- c("Motif", "Top_FetalBrainMale", "Top_FetalBrainFemale")

# Plot Odds Ratio
hm.max <- quantile(brain$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
gg <- ggplot(data = brain)
gg +
        geom_tile(aes(x = chromState, y = userSet, fill = oddsRatio)) +
        facet_grid(cols = vars(cellType)) +
        scale_fill_gradientn("Odds Ratio", colors = c("black", "#FF0000"), values = c(0, 1), 
                             na.value = "#FF0000", limits = c(0, hm.max), breaks = pretty_breaks(n = 3)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              panel.background = element_rect(fill = "black"), axis.ticks.x = element_line(size = 1.25), 
              axis.ticks.y = element_line(size = 1.25), legend.key = element_blank(), legend.position = c(1.21, 0.855), 
              legend.background = element_blank(), legend.title = element_text(size = 18), 
              plot.margin = unit(c(0.5, 8, 0.5, 0.5), "lines"), axis.text.y = element_text(size = 15, color = "black"), 
              axis.text.x = element_text(size = 15, color = "black", angle = 90, hjust = 1, vjust = 0.5), 
              axis.title = element_blank(), plot.title = element_text(size = 18, hjust = 0.5, vjust = 0),
              strip.background = element_blank(), strip.text = element_text(size = 18)) +
        scale_x_discrete(expand = c(0, 0))
ggsave("Figures/Top TF Motif LOLA Fetal Brain ChromHMM Enrichments Odds Ratio Heatmap.png", dpi = 600, width = 5.5, height = 7, units = "in")

# Plot log q-value
hm.max <- quantile(brain$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
gg <- ggplot(data = brain)
gg +
        geom_tile(aes(x = chromState, y = userSet, fill = oddsRatio)) +
        facet_grid(cols = vars(cellType)) +
        scale_fill_gradientn("Odds Ratio", colors = c("black", "#FF0000"), values = c(0, 1), 
                             na.value = "#FF0000", limits = c(0, hm.max), breaks = pretty_breaks(n = 3)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              panel.background = element_rect(fill = "black"), axis.ticks.x = element_line(size = 1.25), 
              axis.ticks.y = element_line(size = 1.25), legend.key = element_blank(), legend.position = c(1.21, 0.855), 
              legend.background = element_blank(), legend.title = element_text(size = 18), 
              plot.margin = unit(c(0.5, 8, 0.5, 0.5), "lines"), axis.text.y = element_text(size = 15, color = "black"), 
              axis.text.x = element_text(size = 15, color = "black", angle = 90, hjust = 1, vjust = 0.5), 
              axis.title = element_blank(), plot.title = element_text(size = 18, hjust = 0.5, vjust = 0),
              strip.background = element_blank(), strip.text = element_text(size = 18)) +
        scale_x_discrete(expand = c(0, 0))
ggsave("Figures/Top TF Motif LOLA Fetal Brain ChromHMM Enrichments log qvalue Heatmap.png", dpi = 600, width = 5.5, height = 7, units = "in")










