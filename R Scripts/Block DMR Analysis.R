# Differentially Methylated Block Analysis ---------------------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 7/4/19
# Updated with misclassified females excluded and new covariate filtering

# Packages ####
sapply(c("tidyverse", "ggdendro", "scales", "ggplot2", "ggbiplot", "reshape", "grid", "RColorBrewer", "CMplot", "rlist",
         "annotatr", "GenomicRanges", "LOLA", "rtracklayer", "R.utils", "rGREAT", "limma", "DMRichR", "AnnotationDbi",
         "ChIPpeakAnno", "GeneOverlap", "rJava", "RDAVIDWebService"), require, character.only = TRUE)

# Functions ####
source("R Scripts/DMR Analysis Functions.R")
chroms <- c(paste("chr", 1:22, sep = ""), "chrM")
chromsX <- c(paste("chr", 1:22, sep = ""), "chrX", "chrM")
chromsXY <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")

# Discovery DMBs ----------------------------------------------------
# Get Regions and Genes ####
# Load Regions
DiscDMBs <- list(All = loadRegions("DMRs/Discovery/Diagnosis 50/DifferentialBlocks_DxNoXY_Discovery50.csv",
                                   chroms = chroms, sort = TRUE, DMBid = TRUE),
                 #SexAll = loadRegions("DMRs/Discovery/Diagnosis and Sex 50/DifferentialBlocks_DxAdjSex_Discovery50.csv",
                                      #chroms = chromsX, sort = TRUE, DMBid = TRUE), # No DMBs found
                 Males = loadRegions("DMRs/Discovery/Diagnosis Males 50/DifferentialBlocks_Dx_Discovery50_males.csv",
                                     chroms = chromsXY, sort = TRUE, DMBid = TRUE),
                 Females = loadRegions("DMRs/Discovery/Diagnosis Females 50/DifferentialBlocks_Dx_Discovery50_females.csv",
                                       chroms = chromsX, sort = TRUE, DMBid = TRUE))
DiscDMBs <- lapply(DiscDMBs, function(x){
        x$percentDifference <- x$beta * 100 / pi
        return(x)
})
DiscDMBs <- lapply(DiscDMBs, function(x) x[, c("chr", "start", "end", "DMBid", "width", "L", "percentDifference", "stat", "pval", 
                                               "qval")])
DiscBackground <- list(All = loadRegions("DMRs/Discovery/Diagnosis 50/bsseq_block_background_Discovery50.csv",
                                         chroms = chroms, sort = TRUE),
                       Males = loadRegions("DMRs/Discovery/Diagnosis Males 50/bsseq_block_background_Discovery50_males.csv",
                                           chroms = chromsXY, sort = TRUE),
                       Females = loadRegions("DMRs/Discovery/Diagnosis Females 50/bsseq_block_background_Discovery50_females.csv",
                                             chroms = chromsX, sort = TRUE))

DiscSexAllBackground <- loadRegions("DMRs/Discovery/Diagnosis and Sex 50/bsseq_block_background_Discovery50_DxAdjSex.csv",
                                    chroms = chromsX, sort = TRUE)
nrow(DiscSexAllBackground) # 2928 background regions
sum(DiscSexAllBackground$width) / 10^9 # 2.804715 Gb
sum(DiscSexAllBackground$n) / 10^6 # 25.22035 M CpGs
rm(DiscSexAllBackground)

# Write DMB Bed Files
mapply(writeBED, regions = DiscDMBs, file = c("UCSC Tracks/Discovery Diagnosis DMBs.bed", 
                                              "UCSC Tracks/Discovery Diagnosis Males DMBs.bed", 
                                              "UCSC Tracks/Discovery Diagnosis Females DMBs.bed"))

# Annotate DMBs
regDomains <- read.delim("Tables/Regulatory domains hg38.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
DiscDMBsAnno <- mapply(getDMRanno, DMRstats = DiscDMBs, file = c("Tables/Discovery Diagnosis DMBs Annotation.txt", 
                                                                 "Tables/Discovery Diagnosis Males DMBs Annotation.txt", 
                                                                 "Tables/Discovery Diagnosis Females DMBs Annotation.txt"),
                       MoreArgs = list(regDomains = regDomains), SIMPLIFY = FALSE)

# Get DMB Genes
DiscDMBsGenes <- lapply(DiscDMBs, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_name")
DiscBackgroundGenes <- lapply(DiscBackground, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_name")

# Replication DMBs ----------------------------------------------------
# Get Regions and Genes ####
# Load Regions
RepDMBs <- list(All = loadRegions("DMRs/Replication/Diagnosis 50/DifferentialBlocks_DxNoXY_Replication50.csv",
                                   chroms = chroms, sort = TRUE, DMBid = TRUE),
                SexAll = loadRegions("DMRs/Replication/Diagnosis and Sex 50/DifferentialBlocks_DxAdjSex_Replication50.csv",
                                     chroms = chromsX, sort = TRUE, DMBid = TRUE),
                Males = loadRegions("DMRs/Replication/Diagnosis Males 50/DifferentialBlocks_Dx_Replication50_males.csv",
                                    chroms = chromsXY, sort = TRUE, DMBid = TRUE),
                Females = loadRegions("DMRs/Replication/Diagnosis Females 100/DifferentialBlocks_Dx_Replication100_females.csv",
                                      chroms = chromsX, sort = TRUE, DMBid = TRUE))
RepDMBs <- lapply(RepDMBs, function(x){
        x$percentDifference <- x$beta * 100 / pi
        return(x)
})
RepDMBs <- lapply(RepDMBs, function(x) x[, c("chr", "start", "end", "DMBid", "width", "L", "percentDifference", "stat", "pval", 
                                               "qval")])
RepBackground <- list(All = loadRegions("DMRs/Replication/Diagnosis 50/bsseq_block_background_Replication50.csv",
                                        chroms = chroms, sort = TRUE),
                      SexAll = loadRegions("DMRs/Replication/Diagnosis and Sex 50/bsseq_block_background_Replication50_DxAdjSex.csv",
                                           chroms = chromsX, sort = TRUE),
                      Males = loadRegions("DMRs/Replication/Diagnosis Males 50/bsseq_block_background_Replication50_males.csv",
                                           chroms = chromsXY, sort = TRUE),
                      Females = loadRegions("DMRs/Replication/Diagnosis Females 100/bsseq_block_background_Replication100_females.csv",
                                            chroms = chromsX, sort = TRUE))
# Write DMB Bed Files
mapply(writeBED, regions = RepDMBs, file = c("UCSC Tracks/Replication Diagnosis DMBs.bed", 
                                             "UCSC Tracks/Replication Diagnosis and Sex DMBs.bed",
                                             "UCSC Tracks/Replication Diagnosis Males DMBs.bed",
                                             "UCSC Tracks/Replication Diagnosis Females DMBs.bed"))
# Annotate DMBs
RepDMBsAnno <- mapply(getDMRanno, DMRstats = RepDMBs, file = c("Tables/Replication Diagnosis DMBs Annotation.txt", 
                                                               "Tables/Replication Diagnosis Males DMBs Annotation.txt", 
                                                               "Tables/Replication Diagnosis Females DMBs Annotation.txt"), 
                      MoreArgs = list(regDomains = regDomains), SIMPLIFY = FALSE)
# Get DMB Genes
RepDMBsGenes <- lapply(RepDMBs, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_name")
RepBackgroundGenes <- lapply(RepBackground, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_name")

# Discovery vs Replication Block Comparison ----------------------------------------------
# Region Stats ####
samples <- read.csv("Tables/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo and Cell Type.csv", header = TRUE,
                    stringsAsFactors = FALSE)
table(samples$Sex, samples$Platform)
regionStats <- getRegionStats(DMRs = append(DiscDMBs, RepDMBs), background = append(DiscBackground, RepBackground),
                              n = c(106, 74, 32, 46, 46, 38, 8))
regionStats$Comparison <- c("Disc_All", "Disc_Males", "Disc_Females", "Rep_All", "Rep_SexAll", "Rep_Males", "Rep_Females")
write.csv(regionStats, "Tables/Differentially Methylated Blocks Region Stats.csv", row.names = FALSE, quote = FALSE)

# Region Stats chrX
DMBs_chrX <- lapply(append(DiscDMBs, RepDMBs), subset, chr == "chrX")
names(DMBs_chrX) <- c("Disc_All", "Disc_Males", "Disc_Females", "Rep_All", "Rep_SexAll", "Rep_Males", "Rep_Females")
Background_chrX <- lapply(append(DiscBackground, RepBackground), subset, chr == "chrX")
names(Background_chrX) <- c("Disc_All", "Disc_Males", "Disc_Females", "Rep_All", "Rep_SexAll", "Rep_Males", "Rep_Females")
regionStats_chrX <- getRegionStats(DMRs = DMBs_chrX, background = Background_chrX, n = c(106, 74, 32, 46, 46, 38, 8))
write.csv(regionStats_chrX, "Tables/Differentially Methylated Blocks chrX Only Region Stats.csv", row.names = FALSE, quote = FALSE)

# Percent of Background by Direction Stacked Barplots ####
DiscDMBs <- DiscDMBs[c("Males", "Females")]
RepDMBs <- RepDMBs[c("Males", "Females")]
DiscBackground <- DiscBackground[c("Males", "Females")]
RepBackground <- RepBackground[c("Males", "Females")]

# Percent of Background Hyper and Hypomethylated All Chroms
backPercent <- mapply(function(x, y){c("Hyper  " = sum(x$width[x$percentDifference > 0]) * 100 / sum(y$width), 
                                       "Hypo" = sum(x$width[x$percentDifference < 0]) * 100 / sum(y$width))},
                      x = append(DiscDMBs, RepDMBs), y = append(DiscBackground, RepBackground))
colnames(backPercent) <- c("Discovery Males", "Discovery Females", "Replication Males", "Replication Females")
backPercent <- melt(backPercent)
colnames(backPercent) <- c("Direction", "Comparison", "Percent")
backPercent$Comparison <- factor(backPercent$Comparison, levels = c("Discovery Males", "Discovery Females", 
                                                                    "Replication Males", "Replication Females"))
gg <- ggplot(backPercent, aes(x = Comparison, y = Percent, fill = Direction, color = Direction))
gg + 
        geom_bar(stat="identity") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), panel.grid.minor = element_blank(),
              legend.position = c(0.7, 1.05), legend.background = element_blank(), 
              legend.key.size = unit(0.7, "cm"), legend.spacing = unit(4, "lines"),
              axis.ticks = element_line(size = 1.25, color = "black"), 
              legend.text = element_text(size = 20, margin = unit(c(0, 0, 0, 0.5), "lines")),
              strip.background = element_blank(), legend.direction = "horizontal", 
              panel.spacing.y = unit(0, "lines"), 
              axis.text.x = element_text(size = 18, color = "black", angle = 45, hjust = 1),
              plot.margin = unit(c(2, 1, 1, 1), "lines"), axis.title.x = element_blank(), 
              axis.text.y = element_text(size = 20, color = "black"), legend.title = element_blank()) +
        ylab("DMB Width (% of Background)") +
        scale_fill_manual(breaks = c("Hyper  ", "Hypo"), values = c("#FF3366", "#3366CC")) +
        scale_color_manual(breaks = c("Hyper  ", "Hypo"), values = c("#FF3366", "#3366CC")) +
        coord_cartesian(ylim = c(0, 6)) +
        scale_y_continuous(expand = c(0.004, 0), breaks = pretty_breaks(n = 4))
ggsave("Figures/DMB Percent of Background by Direction Stacked Barplot.png", dpi = 600, width = 5, height = 7, units = "in")

# Percent of Background Hyper and Hypomethylated Autosomes
DMBs_auto <- lapply(append(DiscDMBs, RepDMBs), subset, chr %in% paste("chr", 1:22, sep = ""))
names(DMBs_auto) <- c("Disc_Males", "Disc_Females", "Rep_Males", "Rep_Females")
Background_auto <- lapply(append(DiscBackground, RepBackground), subset, chr %in% paste("chr", 1:22, sep = ""))
names(Background_auto) <- c("Disc_Males", "Disc_Females", "Rep_Males", "Rep_Females")

backPercent <- mapply(function(x, y){c("Hyper  " = sum(x$width[x$percentDifference > 0]) * 100 / sum(y$width), 
                                       "Hypo" = sum(x$width[x$percentDifference < 0]) * 100 / sum(y$width))},
                      x = DMBs_auto, y = Background_auto)
colnames(backPercent) <- c("Discovery Males", "Discovery Females", "Replication Males", "Replication Females")
backPercent <- melt(backPercent)
colnames(backPercent) <- c("Direction", "Comparison", "Percent")
backPercent$Comparison <- factor(backPercent$Comparison, levels = c("Discovery Males", "Discovery Females", 
                                                                    "Replication Males", "Replication Females"))
gg <- ggplot(backPercent, aes(x = Comparison, y = Percent, fill = Direction, color = Direction))
gg + 
        geom_bar(stat="identity") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), panel.grid.minor = element_blank(),
              legend.position = c(0.7, 1.05), legend.background = element_blank(), 
              legend.key.size = unit(0.7, "cm"), legend.spacing = unit(4, "lines"),
              axis.ticks = element_line(size = 1.25, color = "black"), 
              legend.text = element_text(size = 20, margin = unit(c(0, 0, 0, 0.5), "lines")),
              strip.background = element_blank(), legend.direction = "horizontal", 
              panel.spacing.y = unit(0, "lines"), 
              axis.text.x = element_text(size = 18, color = "black", angle = 45, hjust = 1),
              plot.margin = unit(c(2, 1, 1, 1), "lines"), axis.title.x = element_blank(), 
              axis.text.y = element_text(size = 20, color = "black"), legend.title = element_blank()) +
        ylab("DMB Width (% of Background)") +
        scale_fill_manual(breaks = c("Hyper  ", "Hypo"), values = c("#FF3366", "#3366CC")) +
        scale_color_manual(breaks = c("Hyper  ", "Hypo"), values = c("#FF3366", "#3366CC")) +
        coord_cartesian(ylim = c(0, 6)) +
        scale_y_continuous(expand = c(0.004, 0), breaks = pretty_breaks(n = 4))
ggsave("Figures/DMB Percent of Background by Direction Autosomes Only Stacked Barplot.png", dpi = 600, width = 5, height = 7, 
       units = "in")

# Percent of Background Hyper and Hypomethylated chrX
DMBs_chrX <- DMBs_chrX[c("Disc_Males", "Disc_Females", "Rep_Males", "Rep_Females")]
Background_chrX <- Background_chrX[c("Disc_Males", "Disc_Females", "Rep_Males", "Rep_Females")]
backPercent <- mapply(function(x, y){c("Hyper  " = sum(x$width[x$percentDifference > 0]) * 100 / sum(y$width), 
                                       "Hypo" = sum(x$width[x$percentDifference < 0]) * 100 / sum(y$width))},
                      x = DMBs_chrX, y = Background_chrX)
colnames(backPercent) <- c("Discovery Males", "Discovery Females", "Replication Males", "Replication Females")
backPercent <- melt(backPercent)
colnames(backPercent) <- c("Direction", "Comparison", "Percent")
backPercent$Comparison <- factor(backPercent$Comparison, levels = c("Discovery Males", "Discovery Females", 
                                                                    "Replication Males", "Replication Females"))
gg <- ggplot(backPercent, aes(x = Comparison, y = Percent, fill = Direction, color = Direction))
gg + 
        geom_bar(stat="identity") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), panel.grid.minor = element_blank(),
              legend.position = c(0.7, 1.05), legend.background = element_blank(), 
              legend.key.size = unit(0.7, "cm"), legend.spacing = unit(4, "lines"),
              axis.ticks = element_line(size = 1.25, color = "black"), 
              legend.text = element_text(size = 20, margin = unit(c(0, 0, 0, 0.5), "lines")),
              strip.background = element_blank(), legend.direction = "horizontal", 
              panel.spacing.y = unit(0, "lines"), 
              axis.text.x = element_text(size = 18, color = "black", angle = 45, hjust = 1),
              plot.margin = unit(c(2, 1, 1, 1), "lines"), axis.title.x = element_blank(), 
              axis.text.y = element_text(size = 20, color = "black"), legend.title = element_blank()) +
        ylab("DMB Width (% of Background)") +
        scale_fill_manual(breaks = c("Hyper  ", "Hypo"), values = c("#FF3366", "#3366CC")) +
        scale_color_manual(breaks = c("Hyper  ", "Hypo"), values = c("#FF3366", "#3366CC")) +
        coord_cartesian(ylim = c(0, 6)) +
        scale_y_continuous(expand = c(0.004, 0), breaks = pretty_breaks(n = 4))
ggsave("Figures/DMB Percent of Background by Direction chrX Only Stacked Barplot.png", dpi = 600, width = 5, height = 7, 
       units = "in")

# Block Overlap ####
# Make GRanges
GR_DiscDMBs_Hyper <- lapply(DiscDMBs, makeGRange, direction = "hyper")
GR_DiscDMBs_Hypo <- lapply(DiscDMBs, makeGRange, direction = "hypo")
GR_RepDMBs_Hyper <- lapply(RepDMBs, makeGRange, direction = "hyper")
GR_RepDMBs_Hypo <- lapply(RepDMBs, makeGRange, direction = "hypo")
GR_DiscBackground <- lapply(DiscBackground, makeGRange, direction = "all")
GR_RepBackground <- lapply(RepBackground, makeGRange, direction = "all")
GR_IntBackground <- mapply(intersect, x = GR_DiscBackground, y = GR_RepBackground)

# Overlap within Discovery
DMRoverlapVenn(GR_DiscDMBs_Hyper, NameOfPeaks = c("Males", "Females"), cex = 3,
               file = "Figures/Hyper DMB Overlap Discovery.pdf", ext.dist = -0.1,
               rotation.degree = 180, cat.pos = c(0, 0), cat.dist = c(0.03, 0.03))
DMRoverlapVenn(GR_DiscDMBs_Hypo, NameOfPeaks = c("Males", "Females"), cex = 3,
               file = "Figures/Hypo DMB Overlap Discovery.pdf", ext.dist = -0.1,
               rotation.degree = 0, cat.pos = c(0, 0), cat.dist = c(0.03, 0.03))

# Overlap within Replication
DMRoverlapVenn(GR_RepDMBs_Hyper, NameOfPeaks = c("Males", "Females"), cex = 3,
               file = "Figures/Hyper DMB Overlap Replication.pdf", ext.dist = -0.1,
               rotation.degree = 180, cat.pos = c(0, 0), cat.dist = c(0.03, 0.03))
DMRoverlapVenn(GR_RepDMBs_Hypo, NameOfPeaks = c("Males", "Females"), cex = 3,
               file = "Figures/Hypo DMB Overlap Replication.pdf", ext.dist = -0.1,
               rotation.degree = 0, cat.pos = c(0, 0), cat.dist = c(0.03, 0.03))

# Discovery vs Replication Males Overlap
DMRoverlapVenn(list(GR_DiscDMBs_Hyper$Males, GR_RepDMBs_Hyper$Males), cex = 3,
               NameOfPeaks = c("Discovery", "Replication"),
               file = "Figures/Hyper DMB Overlap Discovery vs Replication Males.pdf", ext.dist = -0.1,
               rotation.degree = 180, cat.pos = c(0, 0), cat.dist = c(0.04, 0.04))
DMRoverlapVenn(list(GR_DiscDMBs_Hypo$Males, GR_RepDMBs_Hypo$Males), cex = 3, margin = 0.07,
               NameOfPeaks = c("Discovery", "Replication"), ext.percent = c(0.01, 0.01, 0.01),
               file = "Figures/Hypo DMB Overlap Discovery vs Replication Males.pdf", ext.dist = -0.12,
               rotation.degree = 180, cat.pos = c(150, 180), cat.dist = c(0.04, 0.03))

intDMBs_Males_Hypo <- intersect(GR_DiscDMBs_Hypo$Males, GR_RepDMBs_Hypo$Males) %>% as.data.frame
colnames(intDMBs_Males_Hypo)[colnames(intDMBs_Males_Hypo) == "seqnames"] <- "chr"
intDMBs_Males_Hypo$DMBid <- paste("DMB", 1:nrow(intDMBs_Males_Hypo), sep = "_")
intDMBs_Males_Hypo <- getDMRanno(DMRstats = intDMBs_Males_Hypo, regDomains = regDomains,
                                 file = "Tables/Hypo DMB Intersect Discovery vs Replication Males.txt")

# Discovery vs Replication Females Overlap
DMRoverlapVenn(list(GR_DiscDMBs_Hyper$Females, GR_RepDMBs_Hyper$Females), cex = 3,
               NameOfPeaks = c("Discovery", "Replication"),
               file = "Figures/Hyper DMB Overlap Discovery vs Replication Females.pdf", ext.dist = -0.1,
               rotation.degree = 180, cat.pos = c(180, 180), cat.dist = c(0.04, 0.04), ext.text = FALSE)

DMRoverlapVenn(list(GR_DiscDMBs_Hypo$Females, GR_RepDMBs_Hypo$Females), cex = 3,
               NameOfPeaks = c("Discovery", "Replication"), margin = 0.07,
               file = "Figures/Hypo DMB Overlap Discovery vs Replication Females.pdf", ext.dist = -0.12,
               rotation.degree = 180, cat.pos = c(145, 180), cat.dist = c(0.04, 0.03))

intDMBs_Females_Hyper <- intersect(GR_DiscDMBs_Hyper$Females, GR_RepDMBs_Hyper$Females) %>% as.data.frame
colnames(intDMBs_Females_Hyper)[colnames(intDMBs_Females_Hyper) == "seqnames"] <- "chr"
intDMBs_Females_Hyper$DMBid <- paste("DMB", 1:nrow(intDMBs_Females_Hyper), sep = "_")
intDMBs_Females_Hyper <- getDMRanno(DMRstats = intDMBs_Females_Hyper, regDomains = regDomains,
                                    file = "Tables/Hyper DMB Intersect Discovery vs Replication Females.txt")

intDMBs_Females_Hypo <- intersect(GR_DiscDMBs_Hypo$Females, GR_RepDMBs_Hypo$Females) %>% as.data.frame
colnames(intDMBs_Females_Hypo)[colnames(intDMBs_Females_Hypo) == "seqnames"] <- "chr"
intDMBs_Females_Hypo$DMBid <- paste("DMB", 1:nrow(intDMBs_Females_Hypo), sep = "_")
intDMBs_Females_Hypo <- getDMRanno(DMRstats = intDMBs_Females_Hypo, regDomains = regDomains,
                                    file = "Tables/Hypo DMB Intersect Discovery vs Replication Females.txt")

rm(Background_auto, Background_chrX, backPercent, DiscDMBsAnno, DMBs_auto,
   DMBs_chrX, gg, GR_DiscBackground, GR_DiscDMBs_Hyper, GR_DiscDMBs_Hypo, GR_IntBackground,
   GR_RepBackground, GR_RepDMBs_Hyper, GR_RepDMBs_Hypo, intDMBs_Females_Hyper,
   intDMBs_Females_Hypo, intDMBs_Males_Hypo, regionStats, regionStats_chrX, RepDMBsAnno)

# Gene Overlap ####
DiscDMBsGenes <- DiscDMBsGenes[c("Males", "Females")]
RepDMBsGenes <- RepDMBsGenes[c("Males", "Females")]

# Males
geneOverlapVenn(list("Discovery" = DiscDMBsGenes$Males, "Replication" = RepDMBsGenes$Males),
                file = "Figures/DMB Gene Overlap Discovery vs Replication Males Overlap Venn.png",
                cat.pos = c(145, 180), cat.dist = c(0.05, 0.03), cex = 2.75, cat.cex = 3, margin = 0.07)

# Females
geneOverlapVenn(list("Discovery" = DiscDMBsGenes$Females, "Replication" = RepDMBsGenes$Females),
                file = "Figures/DMB Gene Overlap Discovery vs Replication Females Overlap Venn.png",
                cat.pos = c(145, 180), cat.dist = c(0.05, 0.03), cex = 2.75, cat.cex = 3, margin = 0.07)

# Intersecting Genes
intDMBsGenes <- mapply(intersect, x = DiscDMBsGenes, y = RepDMBsGenes)
mapply(write.table, x = intDMBsGenes, file = c("Tables/Discovery and Replication Males DMB Genes.txt", 
                                               "Tables/Discovery and Replication Females DMB Genes.txt"),
       MoreArgs = list(sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE))
intBackgroundGenes <- intersect(DiscBackgroundGenes$Males, DiscBackgroundGenes$Females) %>% 
        intersect(RepBackgroundGenes$Males) %>% intersect(RepBackgroundGenes$Females) # 25995 genes
intersect(intDMBsGenes$Males, intDMBsGenes$Females) # ALG10

# Stats
gom <- newGOM(gsetA = DiscDMBsGenes, gsetB = RepDMBsGenes, genome.size = length(intBackgroundGenes))
gomResults <- cbind(getMatrix(gom, "intersection") %>% melt, getMatrix(gom, "odds.ratio") %>% melt %>% .[,"value"],
                    getMatrix(gom, "pval") %>% melt %>% .[,"value"])
colnames(gomResults) <- c("Discovery", "Replication", "Intersection", "OddsRatio", "pValue")
gomResults <- subset(gomResults, Discovery == Replication)
gomResults$DiscoveryLength <- sapply(DiscDMBsGenes, length)
gomResults$ReplicationLength <- sapply(RepDMBsGenes, length)
gomResults$PerDiscovery <- gomResults$Intersection * 100 / gomResults$DiscoveryLength
gomResults$PerReplication <- gomResults$Intersection * 100 / gomResults$ReplicationLength
gomResults$qValue <- p.adjust(gomResults$pValue, method = "fdr")
gomResults <- gomResults[,c("Discovery", "Replication", "DiscoveryLength", "ReplicationLength", "Intersection",
                            "PerDiscovery", "PerReplication", "OddsRatio", "pValue", "qValue")]
write.csv(gomResults, file = "Tables/DMB Gene Overlap Discovery vs Replication Stats.csv", quote = FALSE, row.names = FALSE)

# DAVID Enrichment ####
# Get entrez IDs
DiscDMBsIDs <- lapply(DiscDMBs, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_entrezID")
RepDMBsIDs <- lapply(RepDMBs, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_entrezID")
intDMBsIDs <- mapply(intersect, x = DiscDMBsIDs, y = RepDMBsIDs)

DiscBackgroundIDs <- lapply(DiscBackground, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_entrezID")
RepBackgroundIDs <- lapply(RepBackground, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_entrezID")
intBackgroundIDs <- mapply(intersect, x = DiscBackgroundIDs, y = RepBackgroundIDs)

# Run DAVID
david <- DAVIDWebService$new(url = "https://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/",
                             email = "cemordaunt@ucdavis.edu")
categories <- getAllAnnotationCategoryNames(david)
categories <- c("BBID", "BIOCARTA", "BIOGRID_INTERACTION", "CGAP_EST_QUARTILE", "CGAP_SAGE_QUARTILE", "CHROMOSOME",                
                "COG_ONTOLOGY", "CYTOBAND", "DIP", "EC_NUMBER", "GAD_DISEASE", "GENE3D", "GNF_U133A_QUARTILE", 
                "GENERIF_SUMMARY", "GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY", "INTERPRO", 
                "INTACT", "OMIM_DISEASE", "MINT", "PIR_SEQ_FEATURE", "PIR_SUMMARY", "PIR_SUPERFAMILY", "PFAM", 
                "REACTOME_PATHWAY", "SMART", "PRINTS", "PRODOM", "PROSITE", "UP_KEYWORDS", "UNIGENE_EST_QUARTILE", 
                "TIGRFAMS", "SUPFAM", "UP_TISSUE")  

DiscDAVID <- mapply(getDAVID, genes = DiscDMBsIDs, background = DiscBackgroundIDs, 
                   file = c("Tables/Discovery Males DMB Genes DAVID Results.txt",
                            "Tables/Discovery Females DMB Genes DAVID Results.txt"),
                   MoreArgs = list(categories = categories), SIMPLIFY = FALSE)
RepDAVID <- mapply(getDAVID, genes = RepDMBsIDs, background = RepBackgroundIDs, 
                    file = c("Tables/Replication Males DMB Genes DAVID Results.txt",
                             "Tables/Replication Females DMB Genes DAVID Results.txt"),
                    MoreArgs = list(categories = categories), SIMPLIFY = FALSE)
intDAVID <- mapply(getDAVID, genes = intDMBsIDs, background = intBackgroundIDs, 
                    file = c("Tables/Discovery and Replication Males DMB Genes DAVID Results.txt",
                             "Tables/Discovery and Replication Females DMB Genes DAVID Results.txt"),
                    MoreArgs = list(categories = categories), SIMPLIFY = FALSE)

# Add Matching Term
DiscDAVID <- lapply(DiscDAVID, function(x){
        x$Category_Term <- paste(x$Category, x$Term, sep = " ")
        return(x)})
RepDAVID <- lapply(RepDAVID, function(x){
        x$Category_Term <- paste(x$Category, x$Term, sep = " ")
        return(x)})

# Merge and Write Files
DiscRepDAVID <- mapply(merge, x = DiscDAVID, y = RepDAVID, MoreArgs = list(by = "Category_Term", all = FALSE), SIMPLIFY = FALSE)
mapply(write.table, x = DiscRepDAVID, file = c("Tables/Overlapping DAVID Terms in Males DMB Genes Discovery and Replication.txt",
                                           "Tables/Overlapping DAVID Terms in Females DMB Genes Discovery and Replication.txt"),
       MoreArgs = list(sep = "\t", quote = FALSE, row.names = FALSE))

# Overlapping DAVID Terms Plots ####
plotDAVID <- read.delim("Tables/Discovery and Replication Males Females DMB Genes DAVID Results for Plot.txt", sep = "\t",
                        header = TRUE, stringsAsFactors = FALSE)
plotDAVID$Comparison <- str_replace_all(plotDAVID$Comparison, pattern = c(" " = "\n"))
plotDAVID$Comparison <- factor(plotDAVID$Comparison, levels = c("Males\nDiscovery", "Males\nReplication",
                                                                "Females\nDiscovery", "Females\nReplication"))
plotDAVIDall <- plotDAVID
plotDAVIDexp <- subset(plotDAVID, grepl("expression", plotDAVID$Term, fixed = TRUE))
plotDAVID <- subset(plotDAVID, !grepl("expression", plotDAVID$Term, fixed = TRUE))

# All Terms
plotDAVIDall <- plotDAVIDall[order(plotDAVIDall$Term),]
pvals <- reshape::cast(plotDAVIDall[,c("Term", "Comparison", "log10Benjamini")], formula = Term ~ Comparison, 
                       fun.aggregate = mean, value = "log10Benjamini", add.missing = TRUE, fill = 0)
TermOrder <- hclust(dist(pvals[,2:ncol(pvals)], method = "euclidean"), method = "ward.D")$order
plotDAVIDall$Term <- factor(plotDAVIDall$Term, levels = unique(plotDAVIDall$Term)[rev(TermOrder)], ordered = TRUE)
gg <- ggplot()
gg <- gg +
        geom_tile(data = plotDAVIDall, aes(y = Term, x = Comparison, fill = log10Benjamini)) +
        scale_fill_gradientn("-log(q-value)", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000",
                             breaks = pretty_breaks(n = 4), limits = c(0, 9)) +
        scale_x_discrete(expand = c(0,0), drop = FALSE) +
        scale_y_discrete(expand = c(0,0), drop = FALSE) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              plot.margin = unit(c(0.5, 6.5, 2, 0.5), "lines"), axis.ticks = element_line(size = 0.7), 
              axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1.1, vjust = 1.15), 
              axis.text.y = element_text(size = 14, color = "black"),
              axis.title = element_blank(), legend.key = element_blank(),  
              legend.position = c(1.24, 0.87), legend.background = element_blank(), 
              legend.key.size = unit(1, "lines"), legend.title = element_text(size = 15), 
              legend.text = element_text(size = 14), panel.background = element_rect(fill = "black"))
ggsave("Figures/Discovery and Replication Males Females DMB Genes DAVID All logp Heatmap.png", plot = gg, dpi = 600, width = 8, 
       height = 6.5, units = "in")

# All Terms, Manual Order
plotDAVIDall$Term <- factor(plotDAVIDall$Term, levels = rev(c("Chromosome 5", "Glycoprotein",
                                                              "Cell membrane", "Plasma membrane",                         
                                                              "Embryo development expression", "Cadherin-like",
                                                              "Cell adhesion", "Membrane", "Transmembrane",                           
                                                              "Transmembrane helix", "Calcium", "Brain expression", 
                                                              "Bone marrow endothelial cells expression", 
                                                              "Integral component of membrane",
                                                              "Cadherin, N-terminal", "Cadherins", "Homophilic cell adhesion",
                                                              "Cadherin signature", "Cadherin domain", "Cadherin",                                
                                                              "Cadherin conserved site")))
gg <- ggplot()
gg <- gg +
        geom_tile(data = plotDAVIDall, aes(y = Term, x = Comparison, fill = log10Benjamini)) +
        geom_text(data = plotDAVIDall, aes(y = Term, x = Comparison), label = "*", color = "white", size = 10, nudge_y = -0.38) +
        scale_fill_gradientn("-log(q-value)", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000",
                             breaks = pretty_breaks(n = 4), limits = c(0, 8)) +
        scale_x_discrete(expand = c(0,0), drop = FALSE) +
        scale_y_discrete(expand = c(0,0), drop = FALSE) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              plot.margin = unit(c(0.5, 6.5, 2, 0.5), "lines"), axis.ticks = element_line(size = 0.7), 
              axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1.1, vjust = 1.15), 
              axis.text.y = element_text(size = 14, color = "black"),
              axis.title = element_blank(), legend.key = element_blank(),  
              legend.position = c(1.24, 0.87), legend.background = element_blank(), 
              legend.key.size = unit(1, "lines"), legend.title = element_text(size = 15), 
              legend.text = element_text(size = 14), panel.background = element_rect(fill = "black"))
ggsave("Figures/Discovery and Replication Males Females DMB Genes DAVID All logp Heatmap Manual Order.png", plot = gg, dpi = 600, width = 8, 
       height = 6.5, units = "in")

# Terms, exclude expression
plotDAVID <- plotDAVID[order(plotDAVID$Term),]
pvals <- reshape::cast(plotDAVID[,c("Term", "Comparison", "log10Benjamini")], formula = Term ~ Comparison, 
                       fun.aggregate = mean, value = "log10Benjamini", add.missing = TRUE, fill = 0)
TermOrder <- hclust(dist(pvals[,2:ncol(pvals)], method = "euclidean"), method = "ward.D")$order
plotDAVID$Term <- factor(plotDAVID$Term, levels = unique(plotDAVID$Term)[rev(TermOrder)], ordered = TRUE)
gg <- ggplot()
gg <- gg +
        geom_tile(data = plotDAVID, aes(y = Term, x = Comparison, fill = log10Benjamini)) +
        scale_fill_gradientn("-log(q-value)", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000",
                             breaks = pretty_breaks(n = 4), limits = c(0, 9)) +
        scale_x_discrete(expand = c(0,0), drop = FALSE) +
        scale_y_discrete(expand = c(0,0), drop = FALSE) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              plot.margin = unit(c(0.5, 6.5, 4, 4), "lines"), axis.ticks = element_line(size = 0.7), 
              axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1.1, vjust = 1.15), 
              axis.text.y = element_text(size = 14, color = "black"),
              axis.title = element_blank(), legend.key = element_blank(),  
              legend.position = c(1.24, 0.855), legend.background = element_blank(), 
              legend.key.size = unit(1, "lines"), legend.title = element_text(size = 15), 
              legend.text = element_text(size = 14), panel.background = element_rect(fill = "black"))
ggsave("Figures/Discovery and Replication Males Females DMB Genes DAVID logp Heatmap.png", plot = gg, dpi = 600, width = 8, 
       height = 6.5, units = "in")

# Expression Terms Only
plotDAVIDexp$Term <- gsub(" expression", replacement = "", plotDAVIDexp$Term, fixed = TRUE)
plotDAVIDexp <- plotDAVIDexp[order(plotDAVIDexp$Term),]
pvals <- reshape::cast(plotDAVIDexp[,c("Term", "Comparison", "log10Benjamini")], formula = Term ~ Comparison, 
                       fun.aggregate = mean, value = "log10Benjamini", add.missing = TRUE, fill = 0)
TermOrder <- hclust(dist(pvals[,2:ncol(pvals)], method = "euclidean"), method = "ward.D")$order
plotDAVIDexp$Term <- factor(plotDAVIDexp$Term, levels = unique(plotDAVIDexp$Term)[TermOrder], ordered = TRUE)
gg <- ggplot()
gg <- gg +
        geom_tile(data = plotDAVIDexp, aes(y = Term, x = Comparison, fill = log10Benjamini)) +
        scale_fill_gradientn("-log(q-value)", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000",
                             breaks = pretty_breaks(n = 4), limits = c(0, 9)) +
        scale_x_discrete(expand = c(0,0), drop = FALSE) +
        scale_y_discrete(expand = c(0,0), drop = FALSE) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              plot.margin = unit(c(0.5, 6.5, 22.75, 5.25), "lines"), axis.ticks = element_line(size = 0.7), 
              axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1.1, vjust = 1.15), 
              axis.text.y = element_text(size = 14, color = "black"),
              axis.title = element_blank(), legend.key = element_blank(),  
              legend.position = c(1.24, 0.15), legend.background = element_blank(), 
              legend.key.size = unit(1, "lines"), legend.title = element_text(size = 15), 
              legend.text = element_text(size = 14), panel.background = element_rect(fill = "black"))
ggsave("Figures/Discovery and Replication Males Females DMB Genes DAVID Expression logp Heatmap.png", plot = gg, dpi = 600, 
       width = 8, height = 6.5, units = "in")

rm(DiscBackgroundIDs, DiscDMBsIDs, gom, gomResults, intBackgroundIDs, intDMBsIDs, RepBackgroundIDs, 
   RepDMBsIDs)

# Overlap Replicated DMB Genes with Replicated DMR Genes ####
# Load DMRs
DiscDMRs <- list(Males = loadRegions("DMRs/Discovery/Diagnosis Males 50/DMRs_Dx_Discovery50_males.csv",
                                     chroms = chromsXY, sort = TRUE, DMRid = TRUE),
                 Females = loadRegions("DMRs/Discovery/Diagnosis Females 50/DMRs_Dx_Discovery50_females.csv",
                                       chroms = chromsX, sort = TRUE, DMRid = TRUE))
DiscDMRbackground <- list(Males = loadRegions("DMRs/Discovery/Diagnosis Males 50/bsseq_background_Discovery50_males.csv",
                                              chroms = chromsXY, sort = TRUE),
                          Females = loadRegions("DMRs/Discovery/Diagnosis Females 50/bsseq_background_Discovery50_females.csv",
                                                chroms = chromsX, sort = TRUE))
RepDMRs <- list(Males = loadRegions("DMRs/Replication/Diagnosis Males 50/DMRs_Dx_Replication50_males.csv",
                                    chroms = chromsXY, sort = TRUE, DMRid = TRUE),
                Females = loadRegions("DMRs/Replication/Diagnosis Females 100/DMRs_Dx_Replication100_females.csv",
                                      chroms = chromsX, sort = TRUE, DMRid = TRUE))
RepDMRbackground <- list(Males = loadRegions("DMRs/Replication/Diagnosis Males 50/bsseq_background_Replication50_males.csv",
                                             chroms = chromsXY, sort = TRUE),
                         Females = loadRegions("DMRs/Replication/Diagnosis Females 100/bsseq_background_Replication100_females.csv",
                                               chroms = chromsX, sort = TRUE))
# Get DMR Genes
DiscDMRsGenes <- lapply(DiscDMRs, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_name")
DiscDMRbackgroundGenes <- lapply(DiscDMRbackground, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_name")
RepDMRsGenes <- lapply(RepDMRs, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_name")
RepDMRbackgroundGenes <- lapply(RepDMRbackground, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_name")

# Intersect DMR Genes
intDMRsGenes <- mapply(intersect, x = DiscDMRsGenes, y = RepDMRsGenes)
intDMRbackgroundGenes <- mapply(intersect, x = DiscDMRbackgroundGenes, y = RepDMRbackgroundGenes)

# Venn Diagrams
geneOverlapVenn(list("DMRs" = intDMRsGenes$Males, "DMBs" = intDMBsGenes$Males),
                file = "Figures/DMR DMB Gene Males Overlap Venn.png", rotation.degree = 0, cat.pos = c(0,10), 
                cat.dist = c(0.03, 0.03), margin = 0.02, cex = 2.75, cat.cex = 3)
geneOverlapVenn(list("DMRs" = intDMRsGenes$Females, "DMBs" = intDMBsGenes$Females),
                file = "Figures/DMR DMB Gene Females Overlap Venn.png", rotation.degree = 0, cat.pos = c(0,45), 
                cat.dist = c(0.03, 0.05), margin = 0.08, cex = 2.75, cat.cex = 3)

# Intersecting Genes
DMRandDMBgenes <- mapply(intersect, x = intDMRsGenes, y = intDMBsGenes)
mapply(write.table, x = DMRandDMBgenes, file = c("Tables/Males Replicated DMR and DMB Genes.txt", 
                                                 "Tables/Females Replicated DMR and DMB Genes.txt"),
       MoreArgs = list(sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE))
DMRandDMBbackgroundGenes <- intersect(intBackgroundGenes, intDMRbackgroundGenes$Males) %>% 
        intersect(intDMRbackgroundGenes$Females) 
# 25978 genes in all DMR and DMB background

# Stats
gom <- newGOM(gsetA = intDMRsGenes, gsetB = intDMBsGenes, genome.size = length(DMRandDMBbackgroundGenes))
gomResults <- cbind(getMatrix(gom, "intersection") %>% melt, getMatrix(gom, "odds.ratio") %>% melt %>% .[,"value"],
                    getMatrix(gom, "pval") %>% melt %>% .[,"value"])
colnames(gomResults) <- c("DMRs", "DMBs", "Intersection", "OddsRatio", "pValue")
gomResults <- subset(gomResults, DMRs == DMBs)
gomResults$DMRlength <- sapply(intDMRsGenes, length)
gomResults$DMBlength <- sapply(intDMBsGenes, length)
gomResults$PerDMRs <- gomResults$Intersection * 100 / gomResults$DMRlength 
gomResults$PerDMBs <- gomResults$Intersection * 100 / gomResults$DMBlength 
gomResults$qValue <- p.adjust(gomResults$pValue, method = "fdr")
gomResults <- gomResults[,c("DMRs", "DMBs", "DMRlength", "DMBlength", "Intersection", "PerDMRs", "PerDMBs", "OddsRatio", 
                            "pValue", "qValue")]
write.csv(gomResults, file = "Tables/DMR DMB Gene Overlap Stats.csv", quote = FALSE, row.names = FALSE)

# Table of Replicated DMR, DMB, and Intersecting Genes
genes <- list(MalesDMRs = intDMRsGenes$Males, FemalesDMRs = intDMRsGenes$Females,
              MalesDMBs = intDMBsGenes$Males, FemalesDMBs = intDMBsGenes$Females)
geneTable <- data.frame(Gene_Symbol = unlist(genes) %>% as.character() %>% unique() %>% sort())
geneTable$Gene_Entrez_ID <- regDomains$gene_entrezID[match(geneTable$Gene_Symbol, regDomains$gene_name)]
geneTable$Males_DMR_Gene <- geneTable$Gene_Symbol %in% genes$MalesDMRs
geneTable$Females_DMR_Gene <- geneTable$Gene_Symbol %in% genes$FemalesDMRs
geneTable$Males_DMB_Gene <- geneTable$Gene_Symbol %in% genes$MalesDMBs
geneTable$Females_DMB_Gene <- geneTable$Gene_Symbol %in% genes$FemalesDMBs
write.table(geneTable, "Tables/Replicated DMR, DMB, and Intersecting Genes.txt", sep = "\t", quote = FALSE, row.names = FALSE,
            col.names = TRUE)

genesDAVID <- read.delim("Tables/Replicated DMR, DMB, and Intersecting Genes DAVID Annotation.txt", sep = "\t",
                         header = TRUE, stringsAsFactors = FALSE)

geneTable <- merge(x = geneTable, y = genesDAVID, by.x = "Gene_Entrez_ID", by.y = "ID", all.x = TRUE, all.y = FALSE, sort = FALSE)
geneTable <- geneTable[,c("Gene_Symbol", "Gene.Name", "Gene_Entrez_ID", "Males_DMR_Gene", "Females_DMR_Gene", "Males_DMB_Gene", 
                          "Females_DMB_Gene", "CYTOBAND", "ENTREZ_GENE_SUMMARY", "GOTERM_BP_DIRECT", 
                          "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "INTERPRO", "KEGG_PATHWAY", "OMIM_DISEASE", 
                          "UNIGENE_EST_QUARTILE", "UP_TISSUE", "GNF_U133A_QUARTILE")]
geneTable <- geneTable[order(geneTable$Gene_Symbol),]
geneTable$Gene.Name <- str_split(string = geneTable$Gene.Name, pattern = "\\(") %>% sapply(function(x) x[1])
geneTable[is.na(geneTable)] <- ""
write.table(geneTable, "Tables/Replicated DMR, DMB, and Intersecting Genes with DAVID Annotation Supplemental Table.txt", sep = "\t", quote = FALSE, row.names = FALSE,
            col.names = TRUE)

