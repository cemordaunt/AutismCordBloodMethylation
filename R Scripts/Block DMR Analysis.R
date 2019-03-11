# Block DMR Analysis ------------------------------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 3/4/19

# Packages ####
sapply(c("tidyverse", "ggdendro", "scales", "ggplot2", "ggbiplot", "reshape", "grid", "RColorBrewer", "CMplot", "rlist",
         "annotatr", "GenomicRanges", "LOLA", "rtracklayer", "R.utils", "rGREAT", "limma", "DMRichR",
         "ChIPpeakAnno"), require, character.only = TRUE)

# Functions ####
source("R Scripts/DMR Analysis Functions.R")

# Discovery Block DMRs ----------------------------------------------------
# Diagnosis ####
# Data
DMRs <- loadRegions("DMRs/Discovery/Diagnosis 50/DifferentialBlocks_DxNoXY_Discovery50.csv",
                    chroms = c(paste("chr", 1:22, sep = ""), "chrM"), sort = TRUE)
DMRs$DMRid <- paste("Block", 1:nrow(DMRs), sep = "_")
DMRs <- DMRs[,c("chr", "start", "end", "DMRid", "width", "L", "percentDifference", "stat", "pval", "qval")]

candidates <- loadRegions("DMRs/Discovery/Diagnosis 50/CandidateBlocks_DxNoXY_Discovery50.csv",
                          chroms = c(paste("chr", 1:22, sep = ""), "chrM"), sort = TRUE)
candidates <- candidates[,c("chr", "start", "end", "width", "L", "percentDifference", "stat", "pval", "qval")]

# Bed file
DMRs_bed <- data.frame("chr" = DMRs$chr, "start" = DMRs$start, "end" = DMRs$end, "name" = DMRs$DMRid, "score" = 0, "strand" = ".", 
                       "thickStart" = 0, "thickEnd" = 0, "RGB" = ifelse(DMRs_bed$percentDifference > 0, "255,0,0", "0,0,255"))
write.table(DMRs_bed, file = "Tables/Discovery Diagnosis Blocks.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Manhattan plot
CMplotDMR(candidates, prefix = "Figures/Discovery Diagnosis Blocks ", plot.type = "m", m.cex = 0.5, bin.max = 100)

# QQ plot
CMplotDMR(candidates, prefix = "Figures/Discovery Diagnosis Blocks ", plot.type = "q")

# DMR Annotation
regDomains <- read.delim("Tables/Regulatory domains hg38.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
DMRs_anno <- getDMRanno(DMRstats = DMRs, regDomains = regDomains, file = "Tables/Discovery Diagnosis Blocks Annotation.txt")

# Diagnosis and Sex ####
# Data
DMRs <- loadRegions("DMRs/Discovery/Diagnosis and Sex 50/DifferentialBlocks_DxAdjSex_Discovery50.csv",
                    chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)
DMRs$DMRid <- paste("Block", 1:nrow(DMRs), sep = "_")
DMRs <- DMRs[,c("chr", "start", "end", "DMRid", "width", "L", "percentDifference", "stat", "pval", "qval")]

candidates <- loadRegions("DMRs/Discovery/Diagnosis and Sex 50/CandidateBlocks_DxAdjSex_Discovery50.csv",
                          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)
candidates <- candidates[,c("chr", "start", "end", "width", "L", "percentDifference", "stat", "pval", "qval")]

# Bed file
DMRs_bed <- data.frame("chr" = DMRs$chr, "start" = DMRs$start, "end" = DMRs$end, "name" = DMRs$DMRid, "score" = 0, "strand" = ".", 
                       "thickStart" = 0, "thickEnd" = 0, "RGB" = ifelse(DMRs$percentDifference > 0, "255,0,0", "0,0,255"))
write.table(DMRs_bed, file = "Tables/Discovery Diagnosis and Sex Blocks.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Manhattan plot
CMplotDMR(candidates, prefix = "Figures/Discovery Diagnosis and Sex Blocks ", plot.type = "m", m.cex = 0.5, bin.max = 100)

# QQ plot
CMplotDMR(candidates, prefix = "Figures/Discovery Diagnosis and Sex Blocks ", plot.type = "q")

# DMR Annotation
DMRs_anno <- getDMRanno(DMRstats = DMRs, regDomains = regDomains, file = "Tables/Discovery Diagnosis and Sex Blocks Annotation.txt")

# Diagnosis Males ####
# Data
DMRs <- loadRegions("DMRs/Discovery/Diagnosis Males 50/DifferentialBlocks_Dx_Discovery50_males.csv",
                    chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)
DMRs$DMRid <- paste("Block", 1:nrow(DMRs), sep = "_")
DMRs <- DMRs[,c("chr", "start", "end", "DMRid", "width", "L", "percentDifference", "stat", "pval", "qval")]

candidates <- loadRegions("DMRs/Discovery/Diagnosis Males 50/CandidateBlocks_Dx_Discovery50_males.csv",
                          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)
candidates <- candidates[,c("chr", "start", "end", "width", "L", "percentDifference", "stat", "pval", "qval")]

# Bed file
DMRs_bed <- data.frame("chr" = DMRs$chr, "start" = DMRs$start, "end" = DMRs$end, "name" = DMRs$DMRid, "score" = 0, "strand" = ".", 
                       "thickStart" = 0, "thickEnd" = 0, "RGB" = ifelse(DMRs$percentDifference > 0, "255,0,0", "0,0,255"))
write.table(DMRs_bed, file = "Tables/Discovery Diagnosis Males Blocks.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Manhattan plot
CMplotDMR(candidates, prefix = "Figures/Discovery Diagnosis Males Blocks ", plot.type = "m", m.cex = 0.5, bin.max = 90)

# QQ plot
CMplotDMR(candidates, prefix = "Figures/Discovery Diagnosis Males Blocks ", plot.type = "q")

# DMR Annotation
DMRs_anno <- getDMRanno(DMRstats = DMRs, regDomains = regDomains, file = "Tables/Discovery Diagnosis Males Blocks Annotation.txt")

# Diagnosis Females ####
# Data
DMRs <- loadRegions("DMRs/Discovery/Diagnosis Females 50/DifferentialBlocks_Dx_Discovery50_females.csv",
                    chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)
DMRs$DMRid <- paste("Block", 1:nrow(DMRs), sep = "_")
DMRs <- DMRs[,c("chr", "start", "end", "DMRid", "width", "L", "percentDifference", "stat", "pval", "qval")]

candidates <- loadRegions("DMRs/Discovery/Diagnosis Females 50/CandidateBlocks_Dx_Discovery50_females.csv",
                          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)
candidates <- candidates[,c("chr", "start", "end", "width", "L", "percentDifference", "stat", "pval", "qval")]

# Bed file
DMRs_bed <- data.frame("chr" = DMRs$chr, "start" = DMRs$start, "end" = DMRs$end, "name" = DMRs$DMRid, "score" = 0, "strand" = ".", 
                       "thickStart" = 0, "thickEnd" = 0, "RGB" = ifelse(DMRs$percentDifference > 0, "255,0,0", "0,0,255"))
write.table(DMRs_bed, file = "Tables/Discovery Diagnosis Females Blocks.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Manhattan plot
CMplotDMR(candidates, prefix = "Figures/Discovery Diagnosis Females Blocks ", plot.type = "m", m.cex = 0.5, bin.max = 90)
CMplotDMR(subset(candidates, chr == "chrX"), prefix = "Figures/Discovery Diagnosis Females chrX Blocks ", 
          plot.type = "m", m.cex = 1.5, bin.max = 45)

# QQ plot
CMplotDMR(candidates, prefix = "Figures/Discovery Diagnosis Females Blocks ", plot.type = "q")

# DMR Annotation
DMRs_anno <- getDMRanno(DMRstats = DMRs, regDomains = regDomains, file = "Tables/Discovery Diagnosis Females Blocks Annotation.txt")

rm(candidates, DMRs, DMRs_anno, DMRs_bed, regDomains)

# Discovery Block Comparison ----------------------------------------------
DxAll <- read.delim("Tables/Discovery Diagnosis Blocks Annotation.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
DxSexAll <- read.delim("Tables/Discovery Diagnosis and Sex Blocks Annotation.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
DxMales <- read.delim("Tables/Discovery Diagnosis Males Blocks Annotation.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
DxFemales <- read.delim("Tables/Discovery Diagnosis Females Blocks Annotation.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
blocks <- list("DxAll" = DxAll, "DxSexAll" = DxSexAll, "DxMales" = DxMales, "DxFemales" = DxFemales)

# Region Stats
regionStats <- sapply(blocks, function(x){
        c("Number" = nrow(x), "Width" = sum(x$width), "CpGs" = sum(x$L), "NumberHyper" = nrow(x[x$percentDifference > 0,]), 
          "WidthHyper" = sum(x[x$percentDifference > 0, "width"]), "CpGsHyper" = sum(x[x$percentDifference > 0, "L"]),
          "NumberHypo" = nrow(x[x$percentDifference < 0,]), "WidthHypo" = sum(x[x$percentDifference < 0, "width"]),
          "CpGsHypo" = sum(x[x$percentDifference < 0, "L"]))
})
regionStats <- regionStats %>% t %>% as.data.frame
regionStats$Regions <- row.names(regionStats)
row.names(regionStats) <- 1:nrow(regionStats)
regionStats <- regionStats[,c("Regions", "Number", "Width", "CpGs", "NumberHyper", "WidthHyper", "CpGsHyper",
                              "NumberHypo", "WidthHypo", "CpGsHypo")]
write.csv(regionStats, "Tables/Discovery Block DMRs Region Stats.csv", row.names = FALSE, quote = FALSE)

# Region Stats chrX
blocks_chrX <- lapply(blocks, subset, chr == "chrX")
regionStats <- sapply(blocks_chrX, function(x){
        c("Number" = nrow(x), "Width" = sum(x$width), "CpGs" = sum(x$L), "NumberHyper" = nrow(x[x$percentDifference > 0,]), 
          "WidthHyper" = sum(x[x$percentDifference > 0, "width"]), "CpGsHyper" = sum(x[x$percentDifference > 0, "L"]),
          "NumberHypo" = nrow(x[x$percentDifference < 0,]), "WidthHypo" = sum(x[x$percentDifference < 0, "width"]),
          "CpGsHypo" = sum(x[x$percentDifference < 0, "L"]))
})
regionStats <- regionStats %>% t %>% as.data.frame
regionStats$Regions <- row.names(regionStats)
row.names(regionStats) <- 1:nrow(regionStats)
regionStats <- regionStats[,c("Regions", "Number", "Width", "CpGs", "NumberHyper", "WidthHyper", "CpGsHyper",
                              "NumberHypo", "WidthHypo", "CpGsHypo")]
write.csv(regionStats, "Tables/Discovery chrX Block DMRs Region Stats.csv", row.names = FALSE, quote = FALSE)

# Block Overlap
GR_blocks <- sapply(blocks, function(x) {GRanges(seqnames = x$chr, ranges = IRanges(start = x$start, end = x$end))})

pdf(file = "Figures/Hyper Block Overlap Dx Discovery.pdf", width = 10, height = 8, onefile = FALSE)
venn <- suppressMessages(makeVennDiagram(Peaks = list(GR_blocks$DxAll[DxAll$percentDifference > 0], 
                                                      GR_blocks$DxSexAll[DxSexAll$percentDifference > 0], 
                                                      GR_blocks$DxMales[DxMales$percentDifference > 0], 
                                                      GR_blocks$DxFemales[DxFemales$percentDifference > 0]), 
                                         NameOfPeaks = c("All_Blocks", "All_Blocks_Sex", "Males_Blocks", "Females_Blocks"), 
                                         maxgap = -1, minoverlap = 1, by = "region", connectedPeaks = "min", 
                                         rotation.degree = 0, margin = 0.02, cat.cex = 2, cex = 2.5, 
                                         fill = c("lightblue", "lavender", "lightpink", "lightgreen"), 
                                         cat.pos = c(350, 10, 0, 0), cat.dist = c(0.2, 0.2, 0.1, 0.08), 
                                         fontfamily = "sans", cat.fontfamily = "sans", ext.dist = -0.4, 
                                         ext.length = 0.85))
dev.off()

pdf(file = "Figures/Hypo Block Overlap Dx Discovery.pdf", width = 10, height = 8, onefile = FALSE)
venn <- suppressMessages(makeVennDiagram(Peaks = list(GR_blocks$DxAll[DxAll$percentDifference < 0], 
                                                      GR_blocks$DxSexAll[DxSexAll$percentDifference < 0], 
                                                      GR_blocks$DxMales[DxMales$percentDifference < 0], 
                                                      GR_blocks$DxFemales[DxFemales$percentDifference < 0]), 
                                         NameOfPeaks = c("All_Blocks", "All_Blocks_Sex", "Males_Blocks", "Females_Blocks"), 
                                         maxgap = -1, minoverlap = 1, by = "region", connectedPeaks = "min", 
                                         rotation.degree = 0, margin = 0.02, cat.cex = 2, cex = 2.5, 
                                         fill = c("lightblue", "lavender", "lightpink", "lightgreen"), 
                                         cat.pos = c(350, 10, 0, 0), cat.dist = c(0.2, 0.2, 0.1, 0.08), 
                                         fontfamily = "sans", cat.fontfamily = "sans", ext.dist = -0.4, 
                                         ext.length = 0.85))
dev.off()

# Male Female Overlap
intersect(GR_blocks$DxMales[DxMales$percentDifference > 0],
          GR_blocks$DxFemales[DxFemales$percentDifference > 0]) # No hyper regions overlap
intersect(GR_blocks$DxMales[DxMales$percentDifference < 0],
          GR_blocks$DxFemales[DxFemales$percentDifference < 0])
# Hypo in males and females chrX 58026216-58094942 chrXp11.21 pericentromeric





