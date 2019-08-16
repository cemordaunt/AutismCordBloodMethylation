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
# Get Enrichment with LOLA (Epigenerate) ####
# Functions
makeGRange <- function(DMRs, direction = c("all", "hyper", "hypo")){
         if(direction == "hyper"){DMRs <- subset(DMRs, percentDifference > 0)}
         if(direction == "hypo"){DMRs <- subset(DMRs, percentDifference < 0)}
         GR <- GRanges(seqnames = DMRs$chr, ranges = IRanges(start = DMRs$start, end = DMRs$end))
}

# Made new RegionDB: FetalBrain_ChromHMM

# Load Region Data
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/HOMER")
regionDB <- loadRegionDB(dbLocation = "/share/lasallelab/programs/LOLA/hg38", useCache = TRUE, limit = NULL, 
                         collections = "FetalBrain_ChromHMM")
files <- list.files("BED") %>% paste("BED/", ., sep = "")
motifs <- lapply(files, fread, sep = "\t", header = FALSE, stringsAsFactors = FALSE, verbose = FALSE, data.table = FALSE,
                 nThread = 60)
names(motifs) <- str_replace_all(files, pattern = c("BED/" = "", "_motif_locations.bed" = ""))
motifs <- lapply(motifs, function(x){
        x <- x[,1:3]
        colnames(x) <- c("chr", "start", "end")
        x <- subset(x, chr %in% c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")) %>% makeGRange(direction = "all")
        return(x)
})
Background <- fread(file = "homer_knownMotifs_hg38.bed", sep = "\t", header = FALSE, stringsAsFactors = FALSE, 
                    verbose = FALSE, data.table = FALSE, nThread = 60) # data.table::fread() for much faster reading of 23GB file
Background <- Background[,1:3]
colnames(Background) <- c("chr", "start", "end")
Background <- subset(Background, chr %in% c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")) %>%
        makeGRange(direction = "all") %>% reduce() # Reduce regions to remove overlaps and merge adjacent sites
length(Background) # 88 854 866
summary(width(Background))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 8.00    11.00    16.00    20.34    25.00 20564.00 
Background_df <- as.data.frame(Background)
Background_df <- Background_df[,c("seqnames", "start", "end")]
colnames(Background_df)[colnames(Background_df) == "seqnames"] <- "chr"
fwrite(Background_df, file = "homer_knownMotifs_hg38_merged_background.bed", sep = "\t", quote = FALSE, row.names = FALSE,
            col.names = TRUE, nThread = 60, verbose = FALSE)
rm(files, makeGRange, Background_df)

# Run LOLA
Results <- runLOLA(userSets = motifs, userUniverse = Background, regionDB = regionDB, cores = 1, redefineUserSets = TRUE)
writeCombinedEnrichment(combinedResults = Results, outFolder = "LOLA", includeSplits = FALSE)
file.copy(from = "LOLA/allEnrichments.tsv", to = "LOLA/Top_TF_Motif_LOLA_ChromHMM_Enrichments.tsv", overwrite = TRUE)

# LOLA Enrichment Stats Info ####
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

# Analysis ####
# Load Results
lola <- read.delim("Tables/Top_TF_Motif_LOLA_ChromHMM_Enrichments.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
index <- read.delim("Tables/LOLA Roadmap ChromHMM index.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
chromHMM <- split(lola, f = lola$userSet) %>%
        lapply(function(x) prepLOLAchromHMM(x, index = index, regions = unique(x$userSet)))
chromHMM <- list.rbind(chromHMM)
rownames(chromHMM) <- 1:nrow(chromHMM)
write.csv(chromHMM, file = "Tables/Top TF Motif LOLA ChromHMM Enrichments.csv", quote = FALSE, row.names = FALSE)

# Get Top Ranked
motifs <- unique(chromHMM$userSet)
top_chromStates <- NULL
for(i in 1:length(motifs)){
        sub <- subset(chromHMM, userSet == motifs[i] & cellType == "Fetal Brain Male")
        maleTop <- sub$chromState[sub$maxRnk == min(sub$maxRnk) & sub$qValue < 0.05] %>% as.character() %>% 
                paste(collapse = "_")
        sub <- subset(chromHMM, userSet == motifs[i] & cellType == "Fetal Brain Female")
        femaleTop <- sub$chromState[sub$maxRnk == min(sub$maxRnk) & sub$qValue < 0.05] %>% as.character() %>% 
                paste(collapse = "_")
        temp <- c(motifs[i], maleTop, femaleTop)
        top_chromStates <- rbind(top_chromStates, temp)
}
top_chromStates <- as.data.frame(top_chromStates)
rownames(top_chromStates) <- 1:nrow(top_chromStates)
colnames(top_chromStates) <- c("Motif", "Top_FetalBrainMale", "Top_FetalBrainFemale")

top_chromStates$Motif[grepl("_", top_chromStates$Top_FetalBrainMale, fixed = TRUE)] %>% as.character() 
# None with multiple top
top_chromStates$Motif[grepl("_", top_chromStates$Top_FetalBrainFemale, fixed = TRUE)] %>% as.character()
# "ebf"      "ets-ebox" "pbx3"
# ebf: ReprPC has higher odds ratio
# ets-ebox: Enh has higher odds ratio
# pbx3: Enh has higher odds ratio

top_chromStates$Top_FetalBrainFemale[top_chromStates$Motif == "ebf"] <- "ReprPC"
top_chromStates$Top_FetalBrainFemale[top_chromStates$Motif == "ets-ebox"] <- "Enh"
top_chromStates$Top_FetalBrainFemale[top_chromStates$Motif == "pbx3"] <- "Enh"
top_chromStates$Same <- as.character(top_chromStates$Top_FetalBrainMale) == 
        as.character(top_chromStates$Top_FetalBrainFemale)
write.csv(top_chromStates, file = "Tables/Top TF Motif Chromatin States LOLA ChromHMM Enrichments.csv", quote = FALSE, 
          row.names = FALSE)

# Plot Odds Ratio
hm.max <- quantile(chromHMM$oddsRatio, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
chromHMM$userSet <- factor(chromHMM$userSet, levels = rev(unique(chromHMM$userSet)))
chromHMM$cellType <- factor(chromHMM$cellType, levels = c("Fetal Brain Male", "Fetal Brain Female"))
gg <- ggplot(data = chromHMM)
gg +
        geom_tile(aes(x = chromState, y = userSet, fill = oddsRatio)) +
        facet_grid(cols = vars(cellType)) +
        scale_fill_gradientn("Odds Ratio", colors = c("black", "#FF0000"), values = c(0, 1), 
                             na.value = "#FF0000", limits = c(0, hm.max), breaks = pretty_breaks(n = 3)) +
        scale_y_discrete(labels = function(x) str_to_upper(x)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              panel.background = element_rect(fill = "black"), axis.ticks.x = element_line(size = 1.25), 
              axis.ticks.y = element_line(size = 1.25), legend.key = element_blank(), legend.position = c(1.13, 0.89), 
              legend.background = element_blank(), legend.title = element_text(size = 18), 
              plot.margin = unit(c(0.5, 8, 0.5, 0.5), "lines"), axis.text.y = element_text(size = 15, color = "black"), 
              axis.text.x = element_text(size = 15, color = "black", angle = 90, hjust = 1, vjust = 0.5), 
              axis.title = element_blank(), plot.title = element_text(size = 18, hjust = 0.5, vjust = 0),
              strip.background = element_blank(), strip.text = element_text(size = 18)) +
        scale_x_discrete(expand = c(0, 0))
ggsave("Figures/Top TF Motif LOLA Fetal Brain ChromHMM Enrichments Odds Ratio Heatmap.png", dpi = 600, width = 9, height = 9, 
       units = "in")

# Plot log q-value
hm.max <- quantile(chromHMM$qValueLog, probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
gg <- ggplot(data = chromHMM)
gg +
        geom_tile(aes(x = chromState, y = userSet, fill = qValueLog)) +
        facet_grid(cols = vars(cellType)) +
        scale_fill_gradientn("-log(q-value)", colors = c("black", "#FF0000"), values = c(0, 1), 
                             na.value = "#FF0000", limits = c(0, hm.max), breaks = pretty_breaks(n = 3)) +
        scale_y_discrete(labels = function(x) str_to_upper(x)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              panel.background = element_rect(fill = "black"), axis.ticks.x = element_line(size = 1.25), 
              axis.ticks.y = element_line(size = 1.25), legend.key = element_blank(), legend.position = c(1.14, 0.89), 
              legend.background = element_blank(), legend.title = element_text(size = 18), 
              plot.margin = unit(c(0.5, 8, 0.5, 0.5), "lines"), axis.text.y = element_text(size = 15, color = "black"), 
              axis.text.x = element_text(size = 15, color = "black", angle = 90, hjust = 1, vjust = 0.5), 
              axis.title = element_blank(), plot.title = element_text(size = 18, hjust = 0.5, vjust = 0),
              strip.background = element_blank(), strip.text = element_text(size = 18)) +
        scale_x_discrete(expand = c(0, 0))
ggsave("Figures/Top TF Motif LOLA Fetal Brain ChromHMM Enrichments log qvalue Heatmap.png", dpi = 600, width = 9, height = 9, 
       units = "in")
rm(gg, index, lola, sub, femaleTop, hm.max, i, maleTop, motifs, temp, chromHMM)

# Combined Transcription Factor Info Heatmap ------------------------------
# Merge Tables and Setup for Plot ####
TFs$Motif_Match <- gsub(".motif", replacement = "", x = TFs$Motif_File, fixed = TRUE)
TFs <- merge(x = TFs, y = top_chromStates, by.x = "Motif_Match", by.y = "Motif", all = TRUE, sort = FALSE)
write.csv(TFs, file = "Tables/HOMER Top Enriched TF Motif Info with Brain Expression and ChromHMM.csv", quote = FALSE, 
          row.names = FALSE)

TFplot <- TFs[,c("Motif", "Gene_Name", "M_DFC_13pcw_12820_RPKM", "F_DFC_13pcw_12834_RPKM", "Unmethylated_Percent",
                 "HeterogeneouslyMethylated_Percent", "Methylated_Percent", "Top_FetalBrainMale", "Top_FetalBrainFemale")]
colnames(TFplot) <- str_replace_all(colnames(TFplot), pattern = c("M_DFC_13pcw_12820_RPKM" = "Exp_MaleFetalBrain",
                                                                  "F_DFC_13pcw_12834_RPKM" = "Exp_FemaleFetalBrain",
                                                                  "Unmethylated_Percent" = "Unmethylated",
                                                                  "HeterogeneouslyMethylated_Percent" = "Partially_Methylated",
                                                                  "Methylated_Percent" = "Methylated",
                                                                  "Top_FetalBrainMale" = "ChromState_MaleFetalBrain",
                                                                  "Top_FetalBrainFemale" = "ChromState_FemaleFetalBrain"))
TFplot$ChromState_MaleFetalBrain <- as.character(TFplot$ChromState_MaleFetalBrain)
TFplot$ChromState_FemaleFetalBrain <- as.character(TFplot$ChromState_FemaleFetalBrain)
TFplot[TFplot == ""] <- NA
TFplot$Gene_Name <- str_replace_all(TFplot$Gene_Name, pattern = c("NR3C1" = "GR", "HIF1A" = "HIF-1A", "NR1D2" = "REVERB"))
TFplot$Motif_Gene <- mapply(function(x,y){ 
        if(x == y){ return(x) 
        } else { return(paste(x, " (", y, ")", sep = ""))}}, 
        x = TFplot$Motif, y = TFplot$Gene_Name) %>%
        factor(levels = sort(unique(.), decreasing = TRUE))
TFplot <- TFplot[,c("Motif_Gene", "Exp_MaleFetalBrain", "Exp_FemaleFetalBrain", "Unmethylated", "Partially_Methylated",
                    "Methylated", "ChromState_MaleFetalBrain", "ChromState_FemaleFetalBrain")]
TFplot <- melt(TFplot, id.vars = "Motif_Gene")
TFplot$Datatype <- c(rep("Expression", 
                         length(TFplot$variable[TFplot$variable %in% c("Exp_MaleFetalBrain", "Exp_FemaleFetalBrain")])),
                     rep("Methylation", 
                         length(TFplot$variable[TFplot$variable %in% c("Unmethylated", "Partially_Methylated", "Methylated")])),
                     rep("ChromState", 
                         length(TFplot$variable[TFplot$variable %in% c("ChromState_MaleFetalBrain", "ChromState_FemaleFetalBrain")]))) %>%
        factor(levels = c("Expression", "Methylation", "ChromState"))
TFplot$Category = c(rep("Male", length(TFplot$variable[TFplot$variable == "Exp_MaleFetalBrain"])),
                    rep("Female", length(TFplot$variable[TFplot$variable == "Exp_FemaleFetalBrain"])),
                    rep("Unmethylated", length(TFplot$variable[TFplot$variable == "Unmethylated"])),
                    rep("Partially_Methylated", length(TFplot$variable[TFplot$variable == "Partially_Methylated"])),
                    rep("Methylated", length(TFplot$variable[TFplot$variable == "Methylated"])),
                    rep("Male", length(TFplot$variable[TFplot$variable == "ChromState_MaleFetalBrain"])),
                    rep("Female", length(TFplot$variable[TFplot$variable == "ChromState_FemaleFetalBrain"]))) %>%
        factor(levels = c("Male", "Female", "Unmethylated", "Partially_Methylated", "Methylated"))

# Heatmap ####
# Expression
TFplot_exp <- subset(TFplot, Datatype == "Expression")
TFplot_exp$value <- as.numeric(TFplot_exp$value)
TFplot_exp$value <- log2(TFplot_exp$value + 1)
gg <- ggplot(data = TFplot_exp)
gg +
        geom_tile(aes(x = Category, y = Motif_Gene, fill = value)) +
        scale_fill_gradientn(expression(atop("Expression", "log"[2]*"(RPKM+1)")), 
                             colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", 
                             breaks = pretty_breaks(n = 4)) +
        scale_x_discrete(expand = c(0.25, 0)) +
        scale_y_discrete(expand = c(0.01, 0)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              plot.margin = unit(c(0.5, 10.25, 4.05, 1), "lines"), axis.ticks = element_line(size = 1), 
              panel.background = element_rect(fill = "black"),
              axis.text.x = element_text(size = 18, color = "black", angle = 90, hjust = 1, vjust = 0.5),
              axis.text.y = element_text(size = 14, color = "black"), axis.title = element_blank(), 
              legend.key = element_blank(), legend.position = c(1.65, 0.89), legend.background = element_blank(), 
              legend.key.size = unit(1, "lines"), legend.title = element_text(size = 18), 
              legend.text = element_text(size = 16), legend.title.align = 0)
ggsave("Figures/HOMER Top Enriched TF Motif Info Expression Heatmap.png", dpi = 600, width = 5.75, height = 10, units = "in")

# Methylation
TFplot_meth <- subset(TFplot, Datatype == "Methylation")
TFplot_meth$value <- as.numeric(TFplot_meth$value)
gg <- ggplot(data = TFplot_meth)
gg +
        geom_tile(aes(x = Category, y = Motif_Gene, fill = value)) +
        scale_fill_gradientn("Binding\nSites (%)", colors = c("Black", "#3366CC"), values = c(0,1), na.value = "Black", 
                             breaks = pretty_breaks(n = 3), limits = c(0, 100)) +
        scale_x_discrete(expand = c(0.25, 0), labels = function(x) str_replace_all(x, pattern = c("_" = "\n"))) +
        scale_y_discrete(expand = c(0.01, 0)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              plot.margin = unit(c(0.5, 6, 0.5, 1), "lines"), axis.ticks.x = element_line(size = 1), 
              panel.background = element_rect(fill = "black"), axis.ticks.y = element_blank(),
              axis.text.x = element_text(size = 18, color = "black", angle = 90, hjust = 1, vjust = 0.5),
              axis.text.y = element_text(size = 14, color = "black"), axis.title = element_blank(), 
              legend.key = element_blank(), legend.position = c(1.25, 0.9), legend.background = element_blank(), 
              legend.key.size = unit(1, "lines"), legend.title = element_text(size = 18), 
              legend.text = element_text(size = 16), legend.title.align = 0)
ggsave("Figures/HOMER Top Enriched TF Motif Info Methylation Heatmap.png", dpi = 600, width = 5.75, height = 10, units = "in")

# Chromatin State
TFplot_chrom <- subset(TFplot, Datatype == "ChromState")
TFplot_chrom$value <- factor(TFplot_chrom$value, levels = c("TssA", "Tx", "TxWk", "Enh", "Het", "ReprPC", "ReprPCwk", "Quies"))
chromColors <- c(rgb(255,0,0, maxColorValue = 255), rgb(0,128,0, maxColorValue = 255), rgb(0,100,0, maxColorValue = 255), 
                 rgb(255,255,0, maxColorValue = 255), rgb(138,145,208, maxColorValue = 255), 
                 rgb(128,128,128, maxColorValue = 255), rgb(192,192,192, maxColorValue = 255), 
                 rgb(255,255,255, maxColorValue = 255))

gg <- ggplot(data = TFplot_chrom)
gg +
        geom_tile(aes(x = Category, y = Motif_Gene, fill = value)) +
        scale_fill_manual("Chromatin State", values = chromColors, na.value = "Black") +
        scale_x_discrete(expand = c(0.25, 0)) +
        scale_y_discrete(expand = c(0.01, 0)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              plot.margin = unit(c(0.5, 10.25, 4.05, 1), "lines"), axis.ticks.x = element_line(size = 1), 
              panel.background = element_rect(fill = "black"), axis.ticks.y = element_blank(),
              axis.text.x = element_text(size = 18, color = "black", angle = 90, hjust = 1, vjust = 0.5),
              axis.text.y = element_text(size = 14, color = "black"), axis.title = element_blank(), 
              legend.key = element_blank(), legend.position = c(1.65, 0.87), legend.background = element_blank(), 
              legend.key.size = unit(1, "lines"), legend.title = element_text(size = 18), 
              legend.text = element_text(size = 16), legend.title.align = 0)
ggsave("Figures/HOMER Top Enriched TF Motif Info Chromatin State Heatmap.png", dpi = 600, width = 5.75, height = 10, units = "in")

