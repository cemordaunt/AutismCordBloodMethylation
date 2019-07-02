# Differentially Methylated Block Analysis ---------------------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 5/15/19

# Packages ####
sapply(c("tidyverse", "ggdendro", "scales", "ggplot2", "ggbiplot", "reshape", "grid", "RColorBrewer", "CMplot", "rlist",
         "annotatr", "GenomicRanges", "LOLA", "rtracklayer", "R.utils", "rGREAT", "limma", "DMRichR", "AnnotationDbi",
         "ChIPpeakAnno", "GeneOverlap"), require, character.only = TRUE)

# Functions ####
source("R Scripts/DMR Analysis Functions.R")

# Discovery DMBs ----------------------------------------------------
# Get Regions and Genes ####
# Load Regions
DiscDMBs <- list(All = loadRegions("DMRs/Discovery/Diagnosis 50/DifferentialBlocks_DxNoXY_Discovery50.csv",
                                   chroms = c(paste("chr", 1:22, sep = ""), "chrM"), sort = TRUE, DMBid = TRUE),
                 Males = loadRegions("DMRs/Discovery/Diagnosis Males 50/DifferentialBlocks_Dx_Discovery50_males.csv",
                                     chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), sort = TRUE, DMBid = TRUE),
                 Females = loadRegions("DMRs/Discovery/Diagnosis Females 50/DifferentialBlocks_Dx_Discovery50_females.csv",
                                       chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE, DMBid = TRUE))
DiscDMBs <- lapply(DiscDMBs, function(x){
        x$percentDifference <- x$beta * 100 / pi
        return(x)
})
DiscDMBs <- lapply(DiscDMBs, function(x) x[, c("chr", "start", "end", "DMBid", "width", "L", "percentDifference", "stat", "pval", 
                                               "qval")])
DiscBackground <- list(All = loadRegions("DMRs/Discovery/Diagnosis 50/bsseq_block_background_Discovery50.csv",
                                         chroms = c(paste("chr", 1:22, sep = ""), "chrM"), sort = TRUE),
                       Males = loadRegions("DMRs/Discovery/Diagnosis Males 50/bsseq_block_background_Discovery50_males.csv",
                                           chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), sort = TRUE),
                       Females = loadRegions("DMRs/Discovery/Diagnosis Females 50/bsseq_block_background_Discovery50_females.csv",
                                             chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE))

# Write DMB Bed Files
mapply(function(x, y) writeBED(regions = x, file = y), x = DiscDMBs, 
       y = c("UCSC Tracks/Discovery Diagnosis DMBs.bed", "UCSC Tracks/Discovery Diagnosis Males DMBs.bed", 
             "UCSC Tracks/Discovery Diagnosis Females DMBs.bed"))

# Annotate DMBs
regDomains <- read.delim("Tables/Regulatory domains hg38.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
DiscDMBsAnno <- mapply(function(x, y) getDMRanno(DMRstats = x, regDomains = regDomains, file = y), x = DiscDMBs,
                       y = c("Tables/Discovery Diagnosis DMBs Annotation.txt", 
                             "Tables/Discovery Diagnosis Males DMBs Annotation.txt", 
                             "Tables/Discovery Diagnosis Females DMBs Annotation.txt"), SIMPLIFY = FALSE)

# Get DMB Genes
DiscDMBsGenes <- lapply(DiscDMBs, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_name")
DiscBackgroundGenes <- lapply(DiscBackground, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_name")

# Replication DMBs ----------------------------------------------------
# Get Regions and Genes ####
# Load Regions
RepDMBs <- list(All = loadRegions("DMRs/Replication/Diagnosis 50/DifferentialBlocks_DxNoXY_Replication50.csv",
                                   chroms = c(paste("chr", 1:22, sep = ""), "chrM"), sort = TRUE, DMBid = TRUE),
                 Males = loadRegions("DMRs/Replication/Diagnosis Males 50/DifferentialBlocks_Dx_Replication50_males.csv",
                                     chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), sort = TRUE, DMBid = TRUE),
                 Females = loadRegions("DMRs/Replication/Diagnosis Females 100/DifferentialBlocks_Dx_Replication100_females.csv",
                                       chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE, DMBid = TRUE))
RepDMBs <- lapply(RepDMBs, function(x){
        x$percentDifference <- x$beta * 100 / pi
        return(x)
})
RepDMBs <- lapply(RepDMBs, function(x) x[, c("chr", "start", "end", "DMBid", "width", "L", "percentDifference", "stat", "pval", 
                                               "qval")])
RepBackground <- list(All = loadRegions("DMRs/Replication/Diagnosis 50/bsseq_block_background_Replication50.csv",
                                         chroms = c(paste("chr", 1:22, sep = ""), "chrM"), sort = TRUE),
                       Males = loadRegions("DMRs/Replication/Diagnosis Males 50/bsseq_block_background_Replication50_males.csv",
                                           chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), sort = TRUE),
                       Females = loadRegions("DMRs/Replication/Diagnosis Females 100/bsseq_block_background_Replication100_females.csv",
                                             chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE))

# Write DMB Bed Files
mapply(function(x, y) writeBED(regions = x, file = y), x = RepDMBs, 
       y = c("UCSC Tracks/Replication Diagnosis DMBs.bed", "UCSC Tracks/Replication Diagnosis Males DMBs.bed", 
             "UCSC Tracks/Replication Diagnosis Females DMBs.bed"))

# Annotate DMBs
regDomains <- read.delim("Tables/Regulatory domains hg38.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
RepDMBsAnno <- mapply(function(x, y) getDMRanno(DMRstats = x, regDomains = regDomains, file = y), x = RepDMBs,
                       y = c("Tables/Replication Diagnosis DMBs Annotation.txt", 
                             "Tables/Replication Diagnosis Males DMBs Annotation.txt", 
                             "Tables/Replication Diagnosis Females DMBs Annotation.txt"), SIMPLIFY = FALSE)

# Get DMB Genes
RepDMBsGenes <- lapply(RepDMBs, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_name")
RepBackgroundGenes <- lapply(RepBackground, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_name")

# Discovery vs Replication Block Comparison ----------------------------------------------
# Region Stats ####
# Region Stats All Chroms
regionStats <- getRegionStats(DMRs = append(DiscDMBs, RepDMBs), background = append(DiscBackground, RepBackground),
                              n = c(108, 76, 32, 46, 38, 8))
regionStats$Comparison <- c("Disc_All", "Disc_Males", "Disc_Females", "Rep_All", "Rep_Males", "Rep_Females")
write.csv(regionStats, "Tables/Differentially Methylated Blocks Region Stats.csv", row.names = FALSE, quote = FALSE)

# Region Stats chrX
DMBs_chrX <- lapply(append(DiscDMBs, RepDMBs), subset, chr == "chrX")
names(DMBs_chrX) <- c("Disc_All", "Disc_Males", "Disc_Females", "Rep_All", "Rep_Males", "Rep_Females")
Background_chrX <- lapply(append(DiscBackground, RepBackground), subset, chr == "chrX")
names(Background_chrX) <- c("Disc_All", "Disc_Males", "Disc_Females", "Rep_All", "Rep_Males", "Rep_Females")
regionStats_chrX <- getRegionStats(DMRs = DMBs_chrX, background = Background_chrX, n = c(108, 76, 32, 46, 38, 8))
write.csv(regionStats_chrX, "Tables/Differentially Methylated Blocks chrX Only Region Stats.csv", row.names = FALSE, quote = FALSE)

# Percent of Background by Direction Stacked Barplots ####
# Percent of Background Hyper and Hypomethylated All Chroms
backPercent <- mapply(function(x, y){c("Hyper  " = sum(x$width[x$percentDifference > 0]) * 100 / sum(y$width), 
                                       "Hypo" = sum(x$width[x$percentDifference < 0]) * 100 / sum(y$width))},
                      x = append(DiscDMBs, RepDMBs), y = append(DiscBackground, RepBackground))
colnames(backPercent) <- c("Discovery All", "Discovery Males", "Discovery Females", 
                           "Replication All", "Replication Males", "Replication Females")
backPercent <- melt(backPercent)
colnames(backPercent) <- c("Direction", "Comparison", "Percent")
backPercent$Comparison <- factor(backPercent$Comparison, levels = c("Discovery All", "Discovery Males", "Discovery Females", 
                                                                    "Replication All", "Replication Males", "Replication Females"))
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
ggsave("Figures/DMB Percent of Background by Direction Stacked Barplot.png", dpi = 600, width = 5, height = 7, 
       units = "in")

# Percent of Background Hyper and Hypomethylated Autosomes
DMBs_auto <- lapply(append(DiscDMBs, RepDMBs), subset, chr %in% paste("chr", 1:22, sep = ""))
names(DMBs_auto) <- c("Disc_All", "Disc_Males", "Disc_Females", "Rep_All", "Rep_Males", "Rep_Females")
Background_auto <- lapply(append(DiscBackground, RepBackground), subset, chr %in% paste("chr", 1:22, sep = ""))
names(Background_auto) <- c("Disc_All", "Disc_Males", "Disc_Females", "Rep_All", "Rep_Males", "Rep_Females")

backPercent <- mapply(function(x, y){c("Hyper  " = sum(x$width[x$percentDifference > 0]) * 100 / sum(y$width), 
                                       "Hypo" = sum(x$width[x$percentDifference < 0]) * 100 / sum(y$width))},
                      x = DMBs_auto, y = Background_auto)
colnames(backPercent) <- c("Discovery All", "Discovery Males", "Discovery Females", 
                           "Replication All", "Replication Males", "Replication Females")
backPercent <- melt(backPercent)
colnames(backPercent) <- c("Direction", "Comparison", "Percent")
backPercent$Comparison <- factor(backPercent$Comparison, levels = c("Discovery All", "Discovery Males", "Discovery Females", 
                                                                    "Replication All", "Replication Males", "Replication Females"))
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
backPercent <- mapply(function(x, y){c("Hyper  " = sum(x$width[x$percentDifference > 0]) * 100 / sum(y$width), 
                                       "Hypo" = sum(x$width[x$percentDifference < 0]) * 100 / sum(y$width))},
                      x = DMBs_chrX, y = Background_chrX)
colnames(backPercent) <- c("Discovery All", "Discovery Males", "Discovery Females", 
                           "Replication All", "Replication Males", "Replication Females")
backPercent <- melt(backPercent)
colnames(backPercent) <- c("Direction", "Comparison", "Percent")
backPercent$Comparison <- factor(backPercent$Comparison, levels = c("Discovery All", "Discovery Males", "Discovery Females", 
                                                                    "Replication All", "Replication Males", "Replication Females"))
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
        coord_cartesian(ylim = c(0, 18)) +
        scale_y_continuous(expand = c(0.004, 0), breaks = pretty_breaks(n = 4))
ggsave("Figures/DMB Percent of Background by Direction chrX Only Stacked Barplot.png", dpi = 600, width = 5, height = 7, 
       units = "in")

# Block Overlap ####
GR_DiscDMBs_Hyper <- lapply(DiscDMBs, function(x){x <- subset(x, percentDifference > 0)
        return(GRanges(seqnames = x$chr, ranges = IRanges(start = x$start, end = x$end)))})
GR_DiscDMBs_Hypo <- lapply(DiscDMBs, function(x){x <- subset(x, percentDifference < 0)
        return(GRanges(seqnames = x$chr, ranges = IRanges(start = x$start, end = x$end)))})
GR_RepDMBs_Hyper <- lapply(RepDMBs, function(x){x <- subset(x, percentDifference > 0)
        return(GRanges(seqnames = x$chr, ranges = IRanges(start = x$start, end = x$end)))})
GR_RepDMBs_Hypo <- lapply(RepDMBs, function(x){x <- subset(x, percentDifference < 0)
        return(GRanges(seqnames = x$chr, ranges = IRanges(start = x$start, end = x$end)))})
GR_RepDMBs_Hyper <- GR_RepDMBs_Hyper[c("Males", "Females")] # No Hyper DMBs in All Comparison

GR_DiscBackground <- lapply(DiscBackground, function(x) GRanges(seqnames = x$chr, ranges = IRanges(start = x$start, end = x$end)))
GR_RepBackground <- lapply(RepBackground, function(x) GRanges(seqnames = x$chr, ranges = IRanges(start = x$start, end = x$end)))
GR_IntBackground <- mapply(intersect, x = GR_DiscBackground, y = GR_RepBackground)

# Overlap within Discovery
pdf(file = "Figures/Hyper DMB Overlap Discovery.pdf", width = 10, height = 8, onefile = FALSE)
venn <- suppressMessages(makeVennDiagram(Peaks = GR_DiscDMBs_Hyper, NameOfPeaks = c("All", "Males", "Females"),  
                                         maxgap = -1, minoverlap = 1, by = "region", connectedPeaks = "min", 
                                         rotation.degree = 0, margin = 0.02, cat.cex = 2, cex = 2.5, 
                                         fill = c("lightblue", "lightpink", "lightgreen"), 
                                         cat.pos = c(0, 0, 0), cat.dist = c(0.05, 0.05, 0.02), 
                                         fontfamily = "sans", cat.fontfamily = "sans", ext.dist = -0.4, 
                                         ext.length = 0.85))
dev.off()

pdf(file = "Figures/Hypo DMB Overlap Discovery.pdf", width = 10, height = 8, onefile = FALSE)
venn <- suppressMessages(makeVennDiagram(Peaks = GR_DiscDMBs_Hypo, NameOfPeaks = c("All", "Males", "Females"),  
                                         maxgap = -1, minoverlap = 1, by = "region", connectedPeaks = "min", 
                                         rotation.degree = 0, margin = 0.02, cat.cex = 2, cex = 2.5, 
                                         fill = c("lightblue", "lightpink", "lightgreen"), 
                                         cat.pos = c(0, 0, 180), cat.dist = c(0.05, 0.05, 0.05), 
                                         fontfamily = "sans", cat.fontfamily = "sans", ext.dist = -0.4, 
                                         ext.length = 0.85))
dev.off()

# Overlap within Replication
pdf(file = "Figures/Hyper DMB Overlap Replication.pdf", width = 10, height = 8, onefile = FALSE)
venn <- suppressMessages(makeVennDiagram(Peaks = GR_RepDMBs_Hyper, NameOfPeaks = c("Males", "Females"),  
                                         maxgap = -1, minoverlap = 1, by = "region", connectedPeaks = "min", 
                                         rotation.degree = 0, margin = 0.02, cat.cex = 2, cex = 2.5, 
                                         fill = c("lightpink", "lightgreen"), 
                                         cat.pos = c(0, 0), cat.dist = c(0.05, 0.05), 
                                         fontfamily = "sans", cat.fontfamily = "sans", ext.dist = -0.4, 
                                         ext.length = 0.85))
dev.off()

pdf(file = "Figures/Hypo DMB Overlap Replication.pdf", width = 10, height = 8, onefile = FALSE)
venn <- suppressMessages(makeVennDiagram(Peaks = GR_RepDMBs_Hypo, NameOfPeaks = c("All", "Males", "Females"),  
                                         maxgap = -1, minoverlap = 1, by = "region", connectedPeaks = "min", 
                                         rotation.degree = 0, margin = 0.02, cat.cex = 2, cex = 2.5, 
                                         fill = c("lightblue", "lightpink", "lightgreen"), 
                                         cat.pos = c(0, 0, 180), cat.dist = c(0.05, 0.05, 0.05), 
                                         fontfamily = "sans", cat.fontfamily = "sans", ext.dist = -0.4, 
                                         ext.length = 0.85))
dev.off()

# Discovery vs Replication All Overlap
# No Hyper DMRs in All Comparison for Replication
# (p-values aren't useful because background length isn't taken into account)
DMRoverlapVenn(list(GR_DiscDMBs_Hypo$All, GR_RepDMBs_Hypo$All),
               NameOfPeaks = c("Discovery", "Replication"), totalTest = length(GR_IntBackground$All),
               file = "Figures/Hypo DMB Overlap Discovery vs Replication All.pdf", ext.dist = -0.1,
               rotation.degree = 180, cat.pos = c(180, 180), cat.dist = c(0.1, 0.03))
intDMBs_All_Hypo <- intersect(GR_DiscDMBs_Hypo$All, GR_RepDMBs_Hypo$All) %>% as.data.frame

# Discovery vs Replication Males Overlap
DMRoverlapVenn(list(GR_DiscDMBs_Hyper$Males, GR_RepDMBs_Hyper$Males),
               NameOfPeaks = c("Discovery", "Replication"), totalTest = length(GR_IntBackground$Males),
               file = "Figures/Hyper DMB Overlap Discovery vs Replication Males.pdf", ext.dist = -0.1,
               rotation.degree = 0, cat.pos = c(0, 0), cat.dist = c(0.03, 0.03))

DMRoverlapVenn(list(GR_DiscDMBs_Hypo$Males, GR_RepDMBs_Hypo$Males),
               NameOfPeaks = c("Discovery", "Replication"), totalTest = length(GR_IntBackground$Males),
               file = "Figures/Hypo DMB Overlap Discovery vs Replication Males.pdf", ext.dist = -0.12,
               rotation.degree = 180, cat.pos = c(180, 180), cat.dist = c(0.03, 0.03))

intDMBs_Males_Hypo <- intersect(GR_DiscDMBs_Hypo$Males, GR_RepDMBs_Hypo$Males) %>% as.data.frame
write.csv(intDMBs_Males_Hypo, file = "Tables/Hypo DMB Intersect Discovery vs Replication Males.csv", quote = FALSE,
          row.names = FALSE)

# Discovery vs Replication Females Overlap
DMRoverlapVenn(list(GR_DiscDMBs_Hyper$Females, GR_RepDMBs_Hyper$Females),
               NameOfPeaks = c("Discovery", "Replication"), totalTest = length(GR_IntBackground$Females),
               file = "Figures/Hyper DMB Overlap Discovery vs Replication Females.pdf", ext.dist = -0.1,
               rotation.degree = 180, cat.pos = c(180, 180), cat.dist = c(0.03, 0.03), ext.text = FALSE)

DMRoverlapVenn(list(GR_DiscDMBs_Hypo$Females, GR_RepDMBs_Hypo$Females),
               NameOfPeaks = c("Discovery", "Replication"), totalTest = length(GR_IntBackground$Females),
               file = "Figures/Hypo DMB Overlap Discovery vs Replication Females.pdf", ext.dist = -0.12,
               rotation.degree = 180, cat.pos = c(180, 180), cat.dist = c(0.1, 0.03))

intDMBs_Females_Hyper <- intersect(GR_DiscDMBs_Hyper$Females, GR_RepDMBs_Hyper$Females) %>% as.data.frame
intDMBs_Females_Hypo <- intersect(GR_DiscDMBs_Hypo$Females, GR_RepDMBs_Hypo$Females) %>% as.data.frame
write.csv(intDMBs_Females_Hypo, file = "Tables/Hypo DMB Intersect Discovery vs Replication Females.csv", quote = FALSE,
          row.names = FALSE)

rm(Background_auto, Background_chrX, backPercent, DiscDMBsAnno, DiscRegionStats, DMBs_auto,
   DMBs_chrX, gg, GR_DiscBackground, GR_DiscDMBs_Hyper, GR_DiscDMBs_Hypo, GR_IntBackground,
   GR_RepBackground, GR_RepDMBs_Hyper, GR_RepDMBs_Hypo, intDMBs_All_Hypo, intDMBs_Females_Hyper,
   intDMBs_Females_Hypo, intDMBs_Males_Hypo, regionStats, regionStats_chrX, venn, RepDMBsAnno)

# Gene Overlap ####

# All
geneOverlapVenn(list("Discovery" = DiscDMBsGenes$All, "Replication" = RepDMBsGenes$All),
                file = "Figures/DMB Gene Overlap Discovery vs Replication All Overlap Venn.png")

# Males
geneOverlapVenn(list("Discovery" = DiscDMBsGenes$Males, "Replication" = RepDMBsGenes$Males),
                file = "Figures/DMB Gene Overlap Discovery vs Replication Males Overlap Venn.png",
                cat.pos = c(160, 180))

# Females
geneOverlapVenn(list("Discovery" = DiscDMBsGenes$Females, "Replication" = RepDMBsGenes$Females),
                file = "Figures/DMB Gene Overlap Discovery vs Replication Females Overlap Venn.png",
                cat.pos = c(155, 180))

# Intersecting Genes
intDMBsGenes <- mapply(intersect, x = DiscDMBsGenes, y = RepDMBsGenes)
intBackgroundGenes <- intersect(DiscBackgroundGenes$All, DiscBackgroundGenes$Males) %>% 
        intersect(DiscBackgroundGenes$Females) %>% intersect(RepBackgroundGenes$All) %>% 
        intersect(RepBackgroundGenes$Males) %>% intersect(RepBackgroundGenes$Females) # 24941 genes
intersect(intDMBsGenes$Males, intDMBsGenes$Females) # None in common between males and females

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

# GeneOverlap DAVID ####
# Get entrez IDs
DiscDMBsIDs <- lapply(DiscDMBs, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_entrezID")
RepDMBsIDs <- lapply(RepDMBs, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_entrezID")
intDMBsIDs <- mapply(intersect, x = DiscDMBsIDs, y = RepDMBsIDs)
mapply(function(x, y) write.table(x, file = y, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE),
       x = intDMBsIDs, y = c("Tables/Overlapping DMB Gene entrezIDs All Discovery vs Replication.txt",
                             "Tables/Overlapping DMB Gene entrezIDs Males Discovery vs Replication.txt",
                             "Tables/Overlapping DMB Gene entrezIDs Females Discovery vs Replication.txt"))

DiscBackgroundIDs <- lapply(DiscBackground, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_entrezID")
RepBackgroundIDs <- lapply(RepBackground, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_entrezID")
intBackgroundIDs <- mapply(intersect, x = DiscBackgroundIDs, y = RepBackgroundIDs)
mapply(function(x, y) write.table(x, file = y, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE),
       x = intBackgroundIDs, y = c("Tables/Overlapping DMB Background Gene entrezIDs All Discovery vs Replication.txt",
                                   "Tables/Overlapping DMB Background Gene entrezIDs Males Discovery vs Replication.txt",
                                   "Tables/Overlapping DMB Background Gene entrezIDs Females Discovery vs Replication.txt"))

# No Enriched Terms for Overlapping Genes from All or Females
# 5 Enriched Terms for Overlappinig Genes from Males (chrX and mental retardation)

rm(DiscBackgroundIDs, DiscDMBsIDs, gom, gomResults, intBackgroundIDs, intDMBsIDs, RepBackgroundIDs, 
   RepDMBsIDs)

# Overlap Replicated DMB Genes with Replicated DMR Genes ####
# Load DMRs
DiscDMRs <- list(All = loadRegions("DMRs/Discovery/Diagnosis 50/DMRs_DxNoXY_Discovery50.csv",
                                   chroms = c(paste("chr", 1:22, sep = ""), "chrM"), sort = TRUE, DMRid = TRUE),
                 Males = loadRegions("DMRs/Discovery/Diagnosis Males 50/DMRs_Dx_Discovery50_males.csv",
                                     chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), sort = TRUE, DMRid = TRUE),
                 Females = loadRegions("DMRs/Discovery/Diagnosis Females 50/DMRs_Dx_Discovery50_females.csv",
                                       chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE, DMRid = TRUE))
DiscDMRbackground <- list(All = loadRegions("DMRs/Discovery/Diagnosis 50/bsseq_background_Discovery50.csv",
                                            chroms = c(paste("chr", 1:22, sep = ""), "chrM"), sort = TRUE),
                          Males = loadRegions("DMRs/Discovery/Diagnosis Males 50/bsseq_background_Discovery50_males.csv",
                                              chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), sort = TRUE),
                          Females = loadRegions("DMRs/Discovery/Diagnosis Females 50/bsseq_background_Discovery50_females.csv",
                                                chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE))
RepDMRs <- list(All = loadRegions("DMRs/Replication/Diagnosis 50/DMRs_DxNoXY_Replication50.csv",
                                  chroms = c(paste("chr", 1:22, sep = ""), "chrM"), sort = TRUE, DMRid = TRUE),
                Males = loadRegions("DMRs/Replication/Diagnosis Males 50/DMRs_Dx_Replication50_males.csv",
                                    chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), sort = TRUE, DMRid = TRUE),
                Females = loadRegions("DMRs/Replication/Diagnosis Females 100/DMRs_Dx_Replication100_females.csv",
                                      chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE, DMRid = TRUE))
RepDMRbackground <- list(All = loadRegions("DMRs/Replication/Diagnosis 50/bsseq_background_Replication50.csv",
                                           chroms = c(paste("chr", 1:22, sep = ""), "chrM"), sort = TRUE),
                         Males = loadRegions("DMRs/Replication/Diagnosis Males 50/bsseq_background_Replication50_males.csv",
                                             chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), sort = TRUE),
                         Females = loadRegions("DMRs/Replication/Diagnosis Females 100/bsseq_background_Replication100_females.csv",
                                               chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE))

# Write DMR Bed Files
mapply(function(x, y) writeBED(regions = x, file = y), x = append(DiscDMRs, RepDMRs), 
       y = c("UCSC Tracks/Discovery Diagnosis DMRs.bed", "UCSC Tracks/Discovery Diagnosis Males DMRs.bed", 
             "UCSC Tracks/Discovery Diagnosis Females DMRs.bed", "UCSC Tracks/Replication Diagnosis DMRs.bed", 
             "UCSC Tracks/Replication Diagnosis Males DMRs.bed", "UCSC Tracks/Replication Diagnosis Females DMRs.bed"))

# Get DMR Genes
DiscDMRsGenes <- lapply(DiscDMRs, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_name")
DiscDMRbackgroundGenes <- lapply(DiscDMRbackground, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_name")
RepDMRsGenes <- lapply(RepDMRs, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_name")
RepDMRbackgroundGenes <- lapply(RepDMRbackground, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_name")

# Intersect DMR Genes
intDMRsGenes <- mapply(intersect, x = DiscDMRsGenes, y = RepDMRsGenes)
intDMRbackgroundGenes <- mapply(intersect, x = DiscDMRbackgroundGenes, y = RepDMRbackgroundGenes)

# Venn Diagrams
geneOverlapVenn(list("DMR_Genes" = intDMRsGenes$All, "DMB_Genes" = intDMBsGenes$All),
                file = "Figures/DMR DMB Gene All Overlap Venn.png", rotation.degree = 0, cat.pos = c(0,45), 
                cat.dist = c(0.04, 0.06), margin = 0.08)
geneOverlapVenn(list("DMR_Genes" = intDMRsGenes$Males, "DMB_Genes" = intDMBsGenes$Males),
                file = "Figures/DMR DMB Gene Males Overlap Venn.png", rotation.degree = 0, cat.pos = c(0,10), 
                cat.dist = c(0.04, 0.04), margin = 0.02)
geneOverlapVenn(list("DMR_Genes" = intDMRsGenes$Females, "DMB_Genes" = intDMBsGenes$Females),
                file = "Figures/DMR DMB Gene Females Overlap Venn.png", rotation.degree = 0, cat.pos = c(0,60), 
                cat.dist = c(0.04, 0.07), margin = 0.08)

# Intersecting Genes
DMRandDMBgenes <- mapply(intersect, x = intDMRsGenes, y = intDMBsGenes)
DMRandDMBbackgroundGenes <- intersect(intBackgroundGenes, intDMRbackgroundGenes$All) %>% 
        intersect(intDMRbackgroundGenes$Males) %>% intersect(intDMRbackgroundGenes$Females) 
# 24926 genes in all DMR and DMB background

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

# DMR DMB Gene DAVID ####
# Get DMB IDs
intMalesDMBsIDs <- intersect(getDMRgeneList(DMRstats = DiscDMBs$Males, regDomains = regDomains, direction = "all", 
                                            type = "gene_entrezID"),
                             getDMRgeneList(DMRstats = RepDMBs$Males, regDomains = regDomains, direction = "all",
                                            type = "gene_entrezID"))
intMalesDMBsBackIDs <- intersect(getDMRgeneList(DMRstats = DiscBackground$Males, regDomains = regDomains, direction = "all", 
                                                type = "gene_entrezID"),
                                 getDMRgeneList(DMRstats = RepBackground$Males, regDomains = regDomains, direction = "all",
                                                type = "gene_entrezID"))

# Get DMR IDs
intMalesDMRsIDs <- intersect(getDMRgeneList(DMRstats = DiscDMRs$Males, regDomains = regDomains, direction = "all", 
                                            type = "gene_entrezID"),
                             getDMRgeneList(DMRstats = RepDMRs$Males, regDomains = regDomains, direction = "all",
                                            type = "gene_entrezID"))
intMalesDMRsBackIDs <- intersect(getDMRgeneList(DMRstats = DiscDMRbackground$Males, regDomains = regDomains, direction = "all", 
                                                type = "gene_entrezID"),
                                 getDMRgeneList(DMRstats = RepDMRbackground$Males, regDomains = regDomains, direction = "all",
                                                type = "gene_entrezID"))

# Intersect IDs and write files
DMRandDMBIDsMales <- intersect(intMalesDMBsIDs, intMalesDMRsIDs)
write.table(DMRandDMBIDsMales, file = "Tables/Overlapping DMR and DMB Gene entrezIDs Males.txt", sep = "\t", quote = FALSE,
            col.names = FALSE, row.names = FALSE)
DMRandDMBbackgroundIDsMales <- intersect(intMalesDMBsBackIDs, intMalesDMRsBackIDs)
write.table(DMRandDMBbackgroundIDsMales, file = "Tables/Overlapping DMR and DMB Background Gene entrezIDs Males.txt", sep = "\t", quote = FALSE,
            col.names = FALSE, row.names = FALSE)



