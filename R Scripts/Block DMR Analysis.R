# Block DMR Analysis ------------------------------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 3/4/19

# Packages ####
sapply(c("tidyverse", "ggdendro", "scales", "ggplot2", "ggbiplot", "reshape", "grid", "RColorBrewer", "CMplot", "rlist",
         "annotatr", "GenomicRanges", "LOLA", "rtracklayer", "R.utils", "rGREAT", "limma", "DMRichR"), require, character.only = TRUE)

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
