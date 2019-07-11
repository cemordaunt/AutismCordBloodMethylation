# Array Probe Overlap -----------------------------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 7/11/19

# Packages ####
sapply(c("tidyverse", "GenomicRanges", "rtracklayer", "reshape2"), require, character.only = TRUE)

# Functions ####
source("R Scripts/DMR Analysis Functions.R")

# Load Regions ------------------------------------------------------------
# DMRs and Background ####
chromsXY <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")
chromsX <- c(paste("chr", 1:22, sep = ""), "chrX", "chrM")
DMRs <- list(MalesDisc = loadRegions("DMRs/Discovery/Diagnosis Males 50/DMRs_Dx_Discovery50_males.csv",
                                     chroms = chromsXY, sort = TRUE, DMRid = TRUE),
             FemalesDisc = loadRegions("DMRs/Discovery/Diagnosis Females 50/DMRs_Dx_Discovery50_females.csv",
                                       chroms = chromsX, sort = TRUE, DMRid = TRUE),
             MalesRep = loadRegions("DMRs/Replication/Diagnosis Males 50/DMRs_Dx_Replication50_males.csv",
                                    chroms = chromsXY, sort = TRUE, DMRid = TRUE),
             FemalesRep = loadRegions("DMRs/Replication/Diagnosis Females 100/DMRs_Dx_Replication100_females.csv",
                                      chroms = chromsX, sort = TRUE, DMRid = TRUE))
Background <- list(MalesDisc = loadRegions("DMRs/Discovery/Diagnosis Males 50/bsseq_background_Discovery50_males.csv",
                                           chroms = chromsXY, sort = TRUE),
                   FemalesDisc = loadRegions("DMRs/Discovery/Diagnosis Females 50/bsseq_background_Discovery50_females.csv",
                                             chroms = chromsX, sort = TRUE),
                   MalesRep = loadRegions("DMRs/Replication/Diagnosis Males 50/bsseq_background_Replication50_males.csv",
                                          chroms = chromsXY, sort = TRUE),
                   FemalesRep = loadRegions("DMRs/Replication/Diagnosis Females 100/bsseq_background_Replication100_females.csv",
                                            chroms = chromsX, sort = TRUE))

# Arrays ####
# 450K
array_450K <- read.delim("UCSC Tracks/HM450_hg38.bed", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(array_450K) <- c("chr", "start", "end", "probeID")

# EPIC (liftOver from hg19)
array_EPIC_hg19 <- read.csv("Tables/MethylationEPIC_v-1-0_B4.csv", header = TRUE, stringsAsFactors = FALSE, skip = 7) %>%
        .[, c("CHR", "MAPINFO", "IlmnID")]
colnames(array_EPIC_hg19) <- c("chr", "start", "probeID")
array_EPIC_hg19$chr <- paste("chr", array_EPIC_hg19$chr, sep = "")
array_EPIC_hg19$end <- array_EPIC_hg19$start + 1
array_EPIC_hg19 <- array_EPIC_hg19[,c("chr", "start", "end", "probeID")] %>%
        subset(!is.na(array_EPIC_hg19$start) & grepl("cg", array_EPIC_hg19$probeID, fixed = TRUE))
GR_array_EPIC_hg19 <- GRanges(seqnames = array_EPIC_hg19$chr, 
                              ranges = IRanges(start = array_EPIC_hg19$start, end = array_EPIC_hg19$end),
                              probeID = array_EPIC_hg19$probeID)
chain <- import.chain("Tables/hg19ToHg38.over.chain")
seqlevelsStyle(GR_array_EPIC_hg19) <- "UCSC"
GR_array_EPIC_hg38 <- liftOver(GR_array_EPIC_hg19, chain) %>% unlist() %>% resize(width = 2)

# Combine
GR_array <- list(array_450K = GRanges(seqnames = array_450K$chr, 
                                      ranges = IRanges(start = array_450K$start, end = array_450K$end),
                                      probeID = array_450K$probeID),
                 array_EPIC = GR_array_EPIC_hg38)
overlapsAny(GR_array$array_EPIC, GR_array$array_450K) %>% sum() # 437050

rm(array_450K, array_EPIC_hg19, GR_array_EPIC_hg19, GR_array_EPIC_hg38, chain, chromsX, chromsXY)

# Discovery DMR Overlap ---------------------------------------------------
GR_Regions <- list(Males_Background = makeGRange(Background$MalesDisc, direction = "all"),
                   Males_Hyper = makeGRange(DMRs$MalesDisc, direction = "hyper"),
                   Males_Hypo = makeGRange(DMRs$MalesDisc, direction = "hypo"),
                   Females_Background = makeGRange(Background$FemalesDisc, direction = "all"),
                   Females_Hyper = makeGRange(DMRs$FemalesDisc, direction = "hyper"),
                   Females_Hypo = makeGRange(DMRs$FemalesDisc, direction = "hypo"))

overlaps_450K <- sapply(GR_Regions, function(x) suppressWarnings(x %over% GR_array$array_450K))
overlaps_EPIC <- sapply(GR_Regions, function(x) suppressWarnings(x %over% GR_array$array_EPIC))
overlaps <- mapply(table, overlaps_450K, overlaps_EPIC) %>% t() %>% as.data.frame()
colnames(overlaps) <- c("None", "Only450K", "OnlyEPIC", "EPICand450K")
overlaps$Regions <- c("Males Background", "Males Hyper DMRs", "Males Hypo DMRs", "Females Background", "Females Hyper DMRs",
                      "Females Hypo DMRs")
overlaps$Regions <- factor(overlaps$Regions, levels = overlaps$Regions)
overlaps <- overlaps[,c("Regions", "None", "Only450K", "OnlyEPIC", "EPICand450K")]
overlaps$Total <- sapply(GR_Regions, length)
overlaps <- melt(overlaps, id.vars = c("Regions", "Total"))
colnames(overlaps) <- c("Regions", "Total", "Overlap", "Count")
overlaps$Percent <- overlaps$Count * 100 / overlaps$Total
overlapColors <- c("None" = "blue3", "Only450K" = "dodgerblue", "OnlyEPIC" = "goldenrod2", "EPICand450K" = "forestgreen")

gg <- ggplot(overlaps, aes(x = Regions, y = Percent, fill = Overlap, color = Overlap))
gg + 
        geom_bar(stat = "identity") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.position = c(1.23, 0.855), 
              legend.background = element_blank(), legend.key.size = unit(0.8, "cm"), axis.title.y = element_text(size = 22),
              axis.ticks = element_line(size = 1.25, color = "black"), plot.title = element_text(size = 24, vjust = 0),
              legend.text = element_text(size = 16, margin = unit(c(0, 0, 0, 0.5), "lines")),
              strip.background = element_blank(), legend.direction = "vertical", panel.spacing.y = unit(0, "lines"), 
              plot.margin = unit(c(0.5, 11, 1, 1), "lines"), axis.title.x = element_blank(), 
              axis.text.x = element_text(size = 17, color = "black", angle = 45, hjust = 1),
              axis.text.y = element_text(size = 17, color = "black"), legend.title = element_blank()) +
        ylab("Regions with at least 1 probe (%)") +
        ggtitle("Discovery") +
        scale_fill_manual(values = overlapColors) +
        scale_color_manual(values = overlapColors) +
        coord_cartesian(ylim = c(0, 100)) +
        scale_y_continuous(expand = c(0.004, 0))
ggsave("Figures/Discovery DMR Array Probe Overlap Stacked Barplot.png", dpi = 600, width = 8, height = 7, units = "in")

# Replication DMR Overlap ---------------------------------------------------
GR_Regions <- list(Males_Background = makeGRange(Background$MalesRep, direction = "all"),
                   Males_Hyper = makeGRange(DMRs$MalesRep, direction = "hyper"),
                   Males_Hypo = makeGRange(DMRs$MalesRep, direction = "hypo"),
                   Females_Background = makeGRange(Background$FemalesRep, direction = "all"),
                   Females_Hyper = makeGRange(DMRs$FemalesRep, direction = "hyper"),
                   Females_Hypo = makeGRange(DMRs$FemalesRep, direction = "hypo"))

overlaps_450K <- sapply(GR_Regions, function(x) suppressWarnings(x %over% GR_array$array_450K))
overlaps_EPIC <- sapply(GR_Regions, function(x) suppressWarnings(x %over% GR_array$array_EPIC))
overlaps <- mapply(table, overlaps_450K, overlaps_EPIC) %>% t() %>% as.data.frame()
colnames(overlaps) <- c("None", "Only450K", "OnlyEPIC", "EPICand450K")
overlaps$Regions <- c("Males Background", "Males Hyper DMRs", "Males Hypo DMRs", "Females Background", "Females Hyper DMRs",
                      "Females Hypo DMRs")
overlaps$Regions <- factor(overlaps$Regions, levels = overlaps$Regions)
overlaps <- overlaps[,c("Regions", "None", "Only450K", "OnlyEPIC", "EPICand450K")]
overlaps$Total <- sapply(GR_Regions, length)
overlaps <- melt(overlaps, id.vars = c("Regions", "Total"))
colnames(overlaps) <- c("Regions", "Total", "Overlap", "Count")
overlaps$Percent <- overlaps$Count * 100 / overlaps$Total
overlapColors <- c("None" = "blue3", "Only450K" = "dodgerblue", "OnlyEPIC" = "goldenrod2", "EPICand450K" = "forestgreen")

gg <- ggplot(overlaps, aes(x = Regions, y = Percent, fill = Overlap, color = Overlap))
gg + 
        geom_bar(stat = "identity") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.position = c(1.23, 0.855), 
              legend.background = element_blank(), legend.key.size = unit(0.8, "cm"), axis.title.y = element_text(size = 22),
              axis.ticks = element_line(size = 1.25, color = "black"), plot.title = element_text(size = 24, vjust = 0),
              legend.text = element_text(size = 16, margin = unit(c(0, 0, 0, 0.5), "lines")),
              strip.background = element_blank(), legend.direction = "vertical", panel.spacing.y = unit(0, "lines"), 
              plot.margin = unit(c(0.5, 11, 1, 1), "lines"), axis.title.x = element_blank(), 
              axis.text.x = element_text(size = 17, color = "black", angle = 45, hjust = 1),
              axis.text.y = element_text(size = 17, color = "black"), legend.title = element_blank()) +
        ylab("Regions with at least 1 probe (%)") +
        ggtitle("Replication") +
        scale_fill_manual(values = overlapColors) +
        scale_color_manual(values = overlapColors) +
        coord_cartesian(ylim = c(0, 100)) +
        scale_y_continuous(expand = c(0.004, 0))
ggsave("Figures/Replication DMR Array Probe Overlap Stacked Barplot.png", dpi = 600, width = 8, height = 7, units = "in")

