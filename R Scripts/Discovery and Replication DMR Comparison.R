# Discovery and Replication DMR Comparison --------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 2/19/19

# Packages ####
sapply(c("reshape2", "tidyverse", "ChIPpeakAnno", "annotatr", "GeneOverlap", "UpSetR", "parallel", "regioneR", "LOLA"), 
       require, character.only = TRUE)

# Functions ####
source("R Scripts/DMR Analysis Functions.R")

# Data ####
regDomains <- read.delim("Tables/Regulatory domains hg38.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Discovery DMRs
DisDxAll <- loadRegions("DMRs/Discovery/Diagnosis 50/DMRs_DxNoXY_Discovery50.csv", 
                        chroms = c(paste("chr", 1:22, sep = ""), "chrM"))
DisDxSexAll <- loadRegions("DMRs/Discovery/Diagnosis and Sex 50/DMRs_DxAdjSex_Discovery50.csv", 
                           chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"))
DisDxMales <- loadRegions("DMRs/Discovery/Diagnosis Males 50/DMRs_Dx_Discovery50_males.csv", 
                          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"))
DisDxFemales <- loadRegions("DMRs/Discovery/Diagnosis Females 50/DMRs_Dx_Discovery50_females.csv", 
                            chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"))
DisRegions <- list(DisDxAll, subset(DisDxAll, percentDifference > 0), subset(DisDxAll, percentDifference < 0),
                   DisDxSexAll, subset(DisDxSexAll, percentDifference > 0), subset(DisDxSexAll, percentDifference < 0),
                   DisDxMales, subset(DisDxMales, percentDifference > 0), subset(DisDxMales, percentDifference < 0),
                   DisDxFemales, subset(DisDxFemales, percentDifference > 0), subset(DisDxFemales, percentDifference < 0))
names(DisRegions) <- c(paste("DisDxAll", c("", "Hyper", "Hypo"), sep = ""), paste("DisDxSexAll", c("", "Hyper", "Hypo"), sep = ""),
                       paste("DisDxMales", c("", "Hyper", "Hypo"), sep = ""), paste("DisDxFemales", c("", "Hyper", "Hypo"), sep = ""))
GR_DisRegions <- sapply(DisRegions, function(x) {GRanges(seqnames = x$chr, ranges = IRanges(start = x$start, end = x$end))})

# Replication DMRs
RepDxAll <- loadRegions("DMRs/Replication/Diagnosis 50/DMRs_DxNoXY_Replication50.csv", 
                        chroms = c(paste("chr", 1:22, sep = ""), "chrM"))
RepDxSexAll <- loadRegions("DMRs/Replication/Diagnosis and Sex 50/DMRs_DxAdjSex_Replication50.csv", 
                           chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"))
RepDxMales <- loadRegions("DMRs/Replication/Diagnosis Males 50/DMRs_Dx_Replication50_males.csv", 
                          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"))
RepDxFemales <- loadRegions("DMRs/Replication/Diagnosis Females 100/DMRs_Dx_Replication100_females.csv", 
                            chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"))
RepRegions <- list(RepDxAll, subset(RepDxAll, percentDifference > 0), subset(RepDxAll, percentDifference < 0),
                   RepDxSexAll, subset(RepDxSexAll, percentDifference > 0), subset(RepDxSexAll, percentDifference < 0),
                   RepDxMales, subset(RepDxMales, percentDifference > 0), subset(RepDxMales, percentDifference < 0),
                   RepDxFemales, subset(RepDxFemales, percentDifference > 0), subset(RepDxFemales, percentDifference < 0))
names(RepRegions) <- c(paste("RepDxAll", c("", "Hyper", "Hypo"), sep = ""), paste("RepDxSexAll", c("", "Hyper", "Hypo"), sep = ""),
                       paste("RepDxMales", c("", "Hyper", "Hypo"), sep = ""), paste("RepDxFemales", c("", "Hyper", "Hypo"), sep = ""))
GR_RepRegions <- sapply(RepRegions, function(x) {GRanges(seqnames = x$chr, ranges = IRanges(start = x$start, end = x$end))})

# Discovery Background
DisDxBackRegions <- list("DisDxAllBack" = loadRegions("DMRs/Discovery/Diagnosis 50/bsseq_background_Discovery50.csv",
                                                      chroms = c(paste("chr", 1:22, sep = ""), "chrM")),
                         "DisDxSexAllBack" = loadRegions("DMRs/Discovery/Diagnosis and Sex 50/bsseq_background_Discovery50.csv",
                                                         chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM")),
                         "DisDxMalesBack" = loadRegions("DMRs/Discovery/Diagnosis Males 50/bsseq_background_Discovery50_males.csv",
                                                        chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")),
                         "DisDxFemalesBack" = loadRegions("DMRs/Discovery/Diagnosis Females 50/bsseq_background_Discovery50_females.csv",
                                                          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM")))
GR_DisDxBackRegions <- sapply(DisDxBackRegions, function(x) {GRanges(seqnames = x$chr, ranges = IRanges(start = x$start, end = x$end))})


# Replication Background
RepDxBackRegions <- list("RepDxAllBack" = loadRegions("DMRs/Replication/Diagnosis 50/bsseq_background_Replication50.csv",
                                                      chroms = c(paste("chr", 1:22, sep = ""), "chrM")),
                         "RepDxSexAllBack" = loadRegions("DMRs/Replication/Diagnosis and Sex 50/bsseq_background_Replication50.csv",
                                                         chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM")),
                         "RepDxMalesBack" = loadRegions("DMRs/Replication/Diagnosis Males 50/bsseq_background_Replication50_males.csv",
                                                        chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")),
                         "RepDxFemalesBack" = loadRegions("DMRs/Replication/Diagnosis Females 100/bsseq_background_Replication100_females.csv",
                                                          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM")))
GR_RepDxBackRegions <- sapply(RepDxBackRegions, function(x) {GRanges(seqnames = x$chr, ranges = IRanges(start = x$start, end = x$end))})

rm(DisDxAll, DisDxSexAll, DisDxMales, DisDxFemales, RepDxAll, RepDxSexAll, RepDxMales, RepDxFemales)

# Overlap by Location -----------------------------------------------------
# Venn Diagrams ####
# Diagnosis All DMRs
DMRoverlapVenn(list(GR_DisRegions$DisDxAllHyper, GR_RepRegions$RepDxAllHyper), 
               NameOfPeaks = c("Dis_All_DMRs", "Rep_All_DMRs"), file = "Figures/Hyper DMR Overlap Dx All Dis vs Rep Venn.pdf",
               rotation.degree = 180, cat.pos = c(170, 180), cat.dist = c(0.03, 0.03))
DMRoverlapVenn(list(GR_DisRegions$DisDxAllHypo, GR_RepRegions$RepDxAllHypo), 
               NameOfPeaks = c("Dis_All_DMRs", "Rep_All_DMRs"), file = "Figures/Hypo DMR Overlap Dx All Dis vs Rep Venn.pdf",
               rotation.degree = 180, cat.pos = c(0, 0), cat.dist = c(0.03, 0.03))

# Diagnosis Sex All DMRs
DMRoverlapVenn(list(GR_DisRegions$DisDxSexAllHyper, GR_RepRegions$RepDxSexAllHyper), 
               NameOfPeaks = c("Dis_All_DMRs_Sex", "Rep_All_DMRs_Sex"), file = "Figures/Hyper DMR Overlap Dx Sex All Dis vs Rep Venn.pdf",
               rotation.degree = 180, cat.pos = c(10, 0), cat.dist = c(0.03, 0.03))
DMRoverlapVenn(list(GR_DisRegions$DisDxSexAllHypo, GR_RepRegions$RepDxSexAllHypo), 
               NameOfPeaks = c("Dis_All_DMRs_Sex", "Rep_All_DMRs_Sex"), file = "Figures/Hypo DMR Overlap Dx Sex All Dis vs Rep Venn.pdf",
               rotation.degree = 180, cat.pos = c(8, 0), cat.dist = c(0.03, 0.03))

# Diagnosis Males DMRs
DMRoverlapVenn(list(GR_DisRegions$DisDxMalesHyper, GR_RepRegions$RepDxMalesHyper), 
               NameOfPeaks = c("Dis_Males_DMRs", "Rep_Males_DMRs"), file = "Figures/Hyper DMR Overlap Dx Males Dis vs Rep Venn.pdf",
               rotation.degree = 180, cat.pos = c(170, 180), cat.dist = c(0.03, 0.03))
DMRoverlapVenn(list(GR_DisRegions$DisDxMalesHypo, GR_RepRegions$RepDxMalesHypo), 
               NameOfPeaks = c("Dis_Males_DMRs", "Rep_Males_DMRs"), file = "Figures/Hypo DMR Overlap Dx Males Dis vs Rep Venn.pdf",
               rotation.degree = 180, cat.pos = c(170, 180), cat.dist = c(0.03, 0.03))

# Diagnosis Females DMRs
DMRoverlapVenn(list(GR_DisRegions$DisDxFemalesHyper, GR_RepRegions$RepDxFemalesHyper), 
               NameOfPeaks = c("Dis_Females_DMRs", "Rep_Females_DMRs"), file = "Figures/Hyper DMR Overlap Dx Females Dis vs Rep Venn.pdf",
               rotation.degree = 180, cat.pos = c(165, 180), cat.dist = c(0.03, 0.03))
DMRoverlapVenn(list(GR_DisRegions$DisDxFemalesHypo, GR_RepRegions$RepDxFemalesHypo), 
               NameOfPeaks = c("Dis_Females_DMRs", "Rep_Females_DMRs"), file = "Figures/Hypo DMR Overlap Dx Females Dis vs Rep Venn.pdf",
               rotation.degree = 180, cat.pos = c(165, 180), cat.dist = c(0.03, 0.03))

# Intersecting Regions ####
DisRepIntersect <- rbind(intersect(GR_DisRegions$DisDxAllHyper, GR_RepRegions$RepDxAllHyper) %>% as.data.frame,
                         intersect(GR_DisRegions$DisDxAllHypo, GR_RepRegions$RepDxAllHypo) %>% as.data.frame,
                         intersect(GR_DisRegions$DisDxSexAllHyper, GR_RepRegions$RepDxSexAllHyper) %>% as.data.frame,
                         intersect(GR_DisRegions$DisDxSexAllHypo, GR_RepRegions$RepDxSexAllHypo) %>% as.data.frame,
                         intersect(GR_DisRegions$DisDxMalesHyper, GR_RepRegions$RepDxMalesHyper) %>% as.data.frame,
                         intersect(GR_DisRegions$DisDxMalesHypo, GR_RepRegions$RepDxMalesHypo) %>% as.data.frame,
                         intersect(GR_DisRegions$DisDxFemalesHyper, GR_RepRegions$RepDxFemalesHyper) %>% as.data.frame,
                         intersect(GR_DisRegions$DisDxFemalesHypo, GR_RepRegions$RepDxFemalesHypo) %>% as.data.frame)
colnames(DisRepIntersect)[1] <- "chr"
DisRepIntersect$DMRid <- paste("DMR", 1:nrow(DisRepIntersect), sep = "_")
DisRepIntersect$Intersect <- c(rep("DxAllHyper", intersect(GR_DisRegions$DisDxAllHyper, GR_RepRegions$RepDxAllHyper) %>% length),
                               rep("DxAllHypo", intersect(GR_DisRegions$DisDxAllHypo, GR_RepRegions$RepDxAllHypo) %>% length),
                               rep("DxSexAllHyper", intersect(GR_DisRegions$DisDxSexAllHyper, GR_RepRegions$RepDxSexAllHyper) %>% length),
                               rep("DxSexAllHypo", intersect(GR_DisRegions$DisDxSexAllHypo, GR_RepRegions$RepDxSexAllHypo) %>% length),
                               rep("DxMalesHyper", intersect(GR_DisRegions$DisDxMalesHyper, GR_RepRegions$RepDxMalesHyper) %>% length),
                               rep("DxMalesHypo", intersect(GR_DisRegions$DisDxMalesHypo, GR_RepRegions$RepDxMalesHypo) %>% length),
                               rep("DxFemalesHyper", intersect(GR_DisRegions$DisDxFemalesHyper, GR_RepRegions$RepDxFemalesHyper) %>% length),
                               rep("DxFemalesHypo", intersect(GR_DisRegions$DisDxFemalesHypo, GR_RepRegions$RepDxFemalesHypo) %>% length))

DisRepIntersectAnno <- getDMRanno(DisRepIntersect, regDomains = regDomains, file = "Tables/Discovery Replication DMR Intersect with Anno.txt")

# Venn Diagrams with Extended DMRs ####
GR_DisRegions_extend <- lapply(GR_DisRegions, GRangeExtend, extend = 5000)
GR_RepRegions_extend <- lapply(GR_RepRegions, GRangeExtend, extend = 5000)

# Diagnosis All DMRs
DMRoverlapVenn(list(GR_DisRegions_extend$DisDxAllHyper, GR_RepRegions_extend$RepDxAllHyper), 
               NameOfPeaks = c("Dis_All_DMRs", "Rep_All_DMRs"), 
               file = "Figures/Hyper Extended DMR Overlap Dx All Dis vs Rep Venn.pdf",
               rotation.degree = 180, cat.pos = c(170, 180), cat.dist = c(0.03, 0.03))
DMRoverlapVenn(list(GR_DisRegions_extend$DisDxAllHypo, GR_RepRegions_extend$RepDxAllHypo), 
               NameOfPeaks = c("Dis_All_DMRs", "Rep_All_DMRs"), 
               file = "Figures/Hypo Extended DMR Overlap Dx All Dis vs Rep Venn.pdf",
               rotation.degree = 180, cat.pos = c(165, 180), cat.dist = c(0.03, 0.03))

# Diagnosis Sex All DMRs
DMRoverlapVenn(list(GR_DisRegions_extend$DisDxSexAllHyper, GR_RepRegions_extend$RepDxSexAllHyper), 
               NameOfPeaks = c("Dis_All_DMRs_Sex", "Rep_All_DMRs_Sex"), 
               file = "Figures/Hyper Extended DMR Overlap Dx Sex All Dis vs Rep Venn.pdf",
               rotation.degree = 180, cat.pos = c(10, 0), cat.dist = c(0.03, 0.03))
DMRoverlapVenn(list(GR_DisRegions_extend$DisDxSexAllHypo, GR_RepRegions_extend$RepDxSexAllHypo), 
               NameOfPeaks = c("Dis_All_DMRs_Sex", "Rep_All_DMRs_Sex"), 
               file = "Figures/Hypo Extended DMR Overlap Dx Sex All Dis vs Rep Venn.pdf",
               rotation.degree = 180, cat.pos = c(8, 0), cat.dist = c(0.03, 0.03))

# Diagnosis Males DMRs
DMRoverlapVenn(list(GR_DisRegions_extend$DisDxMalesHyper, GR_RepRegions_extend$RepDxMalesHyper), 
               NameOfPeaks = c("Dis_Males_DMRs", "Rep_Males_DMRs"), 
               file = "Figures/Hyper Extended DMR Overlap Dx Males Dis vs Rep Venn.pdf",
               rotation.degree = 180, cat.pos = c(170, 180), cat.dist = c(0.03, 0.03))
DMRoverlapVenn(list(GR_DisRegions_extend$DisDxMalesHypo, GR_RepRegions_extend$RepDxMalesHypo), 
               NameOfPeaks = c("Dis_Males_DMRs", "Rep_Males_DMRs"), 
               file = "Figures/Hypo Extended DMR Overlap Dx Males Dis vs Rep Venn.pdf",
               rotation.degree = 180, cat.pos = c(170, 180), cat.dist = c(0.03, 0.03))

# Diagnosis Females DMRs
DMRoverlapVenn(list(GR_DisRegions_extend$DisDxFemalesHyper, GR_RepRegions_extend$RepDxFemalesHyper), 
               NameOfPeaks = c("Dis_Females_DMRs", "Rep_Females_DMRs"), 
               file = "Figures/Hyper Extended DMR Overlap Dx Females Dis vs Rep Venn.pdf",
               rotation.degree = 180, cat.pos = c(165, 180), cat.dist = c(0.03, 0.03))
DMRoverlapVenn(list(GR_DisRegions_extend$DisDxFemalesHypo, GR_RepRegions_extend$RepDxFemalesHypo), 
               NameOfPeaks = c("Dis_Females_DMRs", "Rep_Females_DMRs"), 
               file = "Figures/Hypo Extended DMR Overlap Dx Females Dis vs Rep Venn.pdf",
               rotation.degree = 180, cat.pos = c(165, 180), cat.dist = c(0.03, 0.03))

# Intersecting Regions with Extended DMRs ####
DisRepIntersect <- rbind(intersect(GR_DisRegions_extend$DisDxAllHyper, GR_RepRegions_extend$RepDxAllHyper) %>% as.data.frame,
                         intersect(GR_DisRegions_extend$DisDxAllHypo, GR_RepRegions_extend$RepDxAllHypo) %>% as.data.frame,
                         intersect(GR_DisRegions_extend$DisDxSexAllHyper, GR_RepRegions_extend$RepDxSexAllHyper) %>% as.data.frame,
                         intersect(GR_DisRegions_extend$DisDxSexAllHypo, GR_RepRegions_extend$RepDxSexAllHypo) %>% as.data.frame,
                         intersect(GR_DisRegions_extend$DisDxMalesHyper, GR_RepRegions_extend$RepDxMalesHyper) %>% as.data.frame,
                         intersect(GR_DisRegions_extend$DisDxMalesHypo, GR_RepRegions_extend$RepDxMalesHypo) %>% as.data.frame,
                         intersect(GR_DisRegions_extend$DisDxFemalesHyper, GR_RepRegions_extend$RepDxFemalesHyper) %>% as.data.frame,
                         intersect(GR_DisRegions_extend$DisDxFemalesHypo, GR_RepRegions_extend$RepDxFemalesHypo) %>% as.data.frame)
colnames(DisRepIntersect)[1] <- "chr"
DisRepIntersect$DMRid <- paste("DMR", 1:nrow(DisRepIntersect), sep = "_")
DisRepIntersect$Intersect <- c(rep("DxAllHyper", intersect(GR_DisRegions_extend$DisDxAllHyper, GR_RepRegions_extend$RepDxAllHyper) %>% length),
                               rep("DxAllHypo", intersect(GR_DisRegions_extend$DisDxAllHypo, GR_RepRegions_extend$RepDxAllHypo) %>% length),
                               rep("DxSexAllHyper", intersect(GR_DisRegions_extend$DisDxSexAllHyper, GR_RepRegions_extend$RepDxSexAllHyper) %>% length),
                               rep("DxSexAllHypo", intersect(GR_DisRegions_extend$DisDxSexAllHypo, GR_RepRegions_extend$RepDxSexAllHypo) %>% length),
                               rep("DxMalesHyper", intersect(GR_DisRegions_extend$DisDxMalesHyper, GR_RepRegions_extend$RepDxMalesHyper) %>% length),
                               rep("DxMalesHypo", intersect(GR_DisRegions_extend$DisDxMalesHypo, GR_RepRegions_extend$RepDxMalesHypo) %>% length),
                               rep("DxFemalesHyper", intersect(GR_DisRegions_extend$DisDxFemalesHyper, GR_RepRegions_extend$RepDxFemalesHyper) %>% length),
                               rep("DxFemalesHypo", intersect(GR_DisRegions_extend$DisDxFemalesHypo, GR_RepRegions_extend$RepDxFemalesHypo) %>% length))

DisRepIntersectAnno <- getDMRanno(DisRepIntersect, regDomains = regDomains, 
                                  file = "Tables/Discovery Replication Extended DMR Intersect with Anno.txt")

# regioneR Stats ####
set.seed(5)

# DxAll
hg38_auto <- getGenomeAndMask("hg38", mask = NULL)$genome %>% 
        filterChromosomes(keep.chr = c(paste("chr", 1:22, sep = ""), "chrM"))
stats_DxAllHyper <- DMRpermTest(A = GR_DisRegions$DisDxAllHyper, B = GR_RepRegions$RepDxAllHyper, genome = hg38_auto, 
                                universe = intersect(GR_DisDxBackRegions$DisDxAllBack, GR_RepDxBackRegions$RepDxAllBack),
                                file = "Figures/Hyper DMR Overlap Dx All Dis vs Rep RegioneR Plots.pdf",
                                Comparison = "DxAllHyper")
stats_DxAllHypo <- DMRpermTest(A = GR_DisRegions$DisDxAllHypo, B = GR_RepRegions$RepDxAllHypo, genome = hg38_auto, 
                               universe = intersect(GR_DisDxBackRegions$DisDxAllBack, GR_RepDxBackRegions$RepDxAllBack),
                               file = "Figures/Hypo DMR Overlap Dx All Dis vs Rep RegioneR Plots.pdf",
                               Comparison = "DxAllHypo")

# DxSexAll
hg38_X <- getGenomeAndMask("hg38", mask = NULL)$genome %>% 
        filterChromosomes(keep.chr = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"))
stats_DxSexAllHyper <- DMRpermTest(A = GR_DisRegions$DisDxSexAllHyper, B = GR_RepRegions$RepDxSexAllHyper, genome = hg38_X, 
                                   universe = intersect(GR_DisDxBackRegions$DisDxSexAllBack, GR_RepDxBackRegions$RepDxSexAllBack),
                                   file = "Figures/Hyper DMR Overlap Dx Sex All Dis vs Rep RegioneR Plots.pdf",
                                   Comparison = "DxSexAllHyper")
stats_DxSexAllHypo <- DMRpermTest(A = GR_DisRegions$DisDxSexAllHypo, B = GR_RepRegions$RepDxSexAllHypo, genome = hg38_X, 
                                  universe = intersect(GR_DisDxBackRegions$DisDxSexAllBack, GR_RepDxBackRegions$RepDxSexAllBack),
                                  file = "Figures/Hypo DMR Overlap Dx Sex All Dis vs Rep RegioneR Plots.pdf",
                                  Comparison = "DxSexAllHypo")

# DxMales
hg38_XY <- getGenomeAndMask("hg38", mask = NULL)$genome %>%
        filterChromosomes(keep.chr = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"))
stats_DxMalesHyper <- DMRpermTest(A = GR_DisRegions$DisDxMalesHyper, B = GR_RepRegions$RepDxMalesHyper, genome = hg38_XY, 
                                    universe = intersect(GR_DisDxBackRegions$DisDxMalesBack, GR_RepDxBackRegions$RepDxMalesBack),
                                    file = "Figures/Hyper DMR Overlap Dx Males Dis vs Rep RegioneR Plots.pdf",
                                    Comparison = "DxMalesHyper")
stats_DxMalesHypo <- DMRpermTest(A = GR_DisRegions$DisDxMalesHypo, B = GR_RepRegions$RepDxMalesHypo, genome = hg38_XY, 
                                   universe = intersect(GR_DisDxBackRegions$DisDxMalesBack, GR_RepDxBackRegions$RepDxMalesBack),
                                   file = "Figures/Hypo DMR Overlap Dx Males Dis vs Rep RegioneR Plots.pdf",
                                   Comparison = "DxMalesHypo")

# DxFemales
stats_DxFemalesHyper <- DMRpermTest(A = GR_DisRegions$DisDxFemalesHyper, B = GR_RepRegions$RepDxFemalesHyper, genome = hg38_X, 
                                   universe = intersect(GR_DisDxBackRegions$DisDxFemalesBack, GR_RepDxBackRegions$RepDxFemalesBack),
                                   file = "Figures/Hyper DMR Overlap Dx Females Dis vs Rep RegioneR Plots.pdf",
                                   Comparison = "DxFemalesHyper")
stats_DxFemalesHypo <- DMRpermTest(A = GR_DisRegions$DisDxFemalesHypo, B = GR_RepRegions$RepDxFemalesHypo, genome = hg38_X, 
                                  universe = intersect(GR_DisDxBackRegions$DisDxFemalesBack, GR_RepDxBackRegions$RepDxFemalesBack),
                                  file = "Figures/Hypo DMR Overlap Dx Females Dis vs Rep RegioneR Plots.pdf",
                                  Comparison = "DxFemalesHypo")

stats_All <- rbind(stats_DxAllHyper, stats_DxAllHypo, stats_DxSexAllHyper, stats_DxSexAllHypo, stats_DxMalesHyper, 
                   stats_DxMalesHypo, stats_DxFemalesHyper, stats_DxFemalesHypo)
write.table(stats_All, file = "Tables/DMR Overlap Dis vs Rep RegioneR Stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)
rm(stats_DxAllHyper, stats_DxAllHypo, stats_DxSexAllHyper, stats_DxSexAllHypo, stats_DxMalesHyper, stats_DxMalesHypo, 
   stats_DxFemalesHyper, stats_DxFemalesHypo)

# regioneR Stats with DMRs Redefined as Background Subset ####
hg38_auto <- getGenomeAndMask("hg38", mask = NULL)$genome %>% 
        filterChromosomes(keep.chr = c(paste("chr", 1:22, sep = ""), "chrM"))
hg38_X <- getGenomeAndMask("hg38", mask = NULL)$genome %>% 
        filterChromosomes(keep.chr = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"))
hg38_XY <- getGenomeAndMask("hg38", mask = NULL)$genome %>%
        filterChromosomes(keep.chr = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"))

# DxAll
universe_back <- intersect(GR_DisDxBackRegions$DisDxAllBack, GR_RepDxBackRegions$RepDxAllBack)
A_HyperDMRs <- redefineUserSets(GRangesList(GR_DisRegions$DisDxAllHyper), universe_back)[[1]]
B_HyperDMRs <- redefineUserSets(GRangesList(GR_RepRegions$RepDxAllHyper), universe_back)[[1]]
A_HypoDMRs <- redefineUserSets(GRangesList(GR_DisRegions$DisDxAllHypo), universe_back)[[1]]
B_HypoDMRs <- redefineUserSets(GRangesList(GR_RepRegions$RepDxAllHypo), universe_back)[[1]]

stats_DxAllHyper <- DMRpermTest(A = A_HyperDMRs, B = B_HyperDMRs, genome = hg38_auto, universe = universe_back,
                                file = "Figures/Hyper DMR Redefined Overlap Dx All Dis vs Rep RegioneR Plots.pdf",
                                Comparison = "DxAllHyper")
stats_DxAllHypo <- DMRpermTest(A = A_HypoDMRs, B = B_HypoDMRs, genome = hg38_auto, universe = universe_back,
                               file = "Figures/Hypo DMR Redefined Overlap Dx All Dis vs Rep RegioneR Plots.pdf",
                               Comparison = "DxAllHypo")

# DxSexAll
universe_back <- intersect(GR_DisDxBackRegions$DisDxSexAllBack, GR_RepDxBackRegions$RepDxSexAllBack)
A_HyperDMRs <- redefineUserSets(GRangesList(GR_DisRegions$DisDxSexAllHyper), universe_back)[[1]]
B_HyperDMRs <- redefineUserSets(GRangesList(GR_RepRegions$RepDxSexAllHyper), universe_back)[[1]]
A_HypoDMRs <- redefineUserSets(GRangesList(GR_DisRegions$DisDxSexAllHypo), universe_back)[[1]]
B_HypoDMRs <- redefineUserSets(GRangesList(GR_RepRegions$RepDxSexAllHypo), universe_back)[[1]]

stats_DxSexAllHyper <- DMRpermTest(A = A_HyperDMRs, B = B_HyperDMRs, genome = hg38_X, universe = universe_back,
                                   file = "Figures/Hyper DMR Redefined Overlap Dx Sex All Dis vs Rep RegioneR Plots.pdf",
                                   Comparison = "DxSexAllHyper")
stats_DxSexAllHypo <- DMRpermTest(A = A_HypoDMRs, B = B_HypoDMRs, genome = hg38_X, universe = universe_back,
                                  file = "Figures/Hypo DMR Redefined Overlap Dx Sex All Dis vs Rep RegioneR Plots.pdf",
                                  Comparison = "DxSexAllHypo")

# DxMales
universe_back <- intersect(GR_DisDxBackRegions$DisDxMalesBack, GR_RepDxBackRegions$RepDxMalesBack)
A_HyperDMRs <- redefineUserSets(GRangesList(GR_DisRegions$DisDxMalesHyper), universe_back)[[1]]
B_HyperDMRs <- redefineUserSets(GRangesList(GR_RepRegions$RepDxMalesHyper), universe_back)[[1]]
A_HypoDMRs <- redefineUserSets(GRangesList(GR_DisRegions$DisDxMalesHypo), universe_back)[[1]]
B_HypoDMRs <- redefineUserSets(GRangesList(GR_RepRegions$RepDxMalesHypo), universe_back)[[1]]

stats_DxMalesHyper <- DMRpermTest(A = A_HyperDMRs, B = B_HyperDMRs, genome = hg38_XY, universe = universe_back,
                                  file = "Figures/Hyper DMR Redefined Overlap Dx Males Dis vs Rep RegioneR Plots.pdf",
                                  Comparison = "DxMalesHyper")
stats_DxMalesHypo <- DMRpermTest(A = A_HypoDMRs, B = B_HypoDMRs, genome = hg38_XY, universe = universe_back,
                                 file = "Figures/Hypo DMR Redefined Overlap Dx Males Dis vs Rep RegioneR Plots.pdf",
                                 Comparison = "DxMalesHypo")

# DxFemales
universe_back <- intersect(GR_DisDxBackRegions$DisDxFemalesBack, GR_RepDxBackRegions$RepDxFemalesBack)
A_HyperDMRs <- redefineUserSets(GRangesList(GR_DisRegions$DisDxFemalesHyper), universe_back)[[1]]
B_HyperDMRs <- redefineUserSets(GRangesList(GR_RepRegions$RepDxFemalesHyper), universe_back)[[1]]
A_HypoDMRs <- redefineUserSets(GRangesList(GR_DisRegions$DisDxFemalesHypo), universe_back)[[1]]
B_HypoDMRs <- redefineUserSets(GRangesList(GR_RepRegions$RepDxFemalesHypo), universe_back)[[1]]

stats_DxFemalesHyper <- DMRpermTest(A = A_HyperDMRs, B = B_HyperDMRs, genome = hg38_X, universe = universe_back,
                                    file = "Figures/Hyper DMR Redefined Overlap Dx Females Dis vs Rep RegioneR Plots.pdf",
                                    Comparison = "DxFemalesHyper")
stats_DxFemalesHypo <- DMRpermTest(A = A_HypoDMRs, B = B_HypoDMRs, genome = hg38_X, universe = universe_back,
                                   file = "Figures/Hypo DMR Redefined Overlap Dx Females Dis vs Rep RegioneR Plots.pdf",
                                   Comparison = "DxFemalesHypo")

stats_All_redefined <- rbind(stats_DxAllHyper, stats_DxAllHypo, stats_DxSexAllHyper, stats_DxSexAllHypo, stats_DxMalesHyper, 
                             stats_DxMalesHypo, stats_DxFemalesHyper, stats_DxFemalesHypo)
write.table(stats_All_redefined, file = "Tables/DMR Redefined Overlap Dis vs Rep RegioneR Stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)
rm(stats_DxAllHyper, stats_DxAllHypo, stats_DxSexAllHyper, stats_DxSexAllHypo, stats_DxMalesHyper, stats_DxMalesHypo, 
   stats_DxFemalesHyper, stats_DxFemalesHypo)

# Overlap by Gene ---------------------------------------------------------
# Background Genes ####
DisDxBackGenes <- lapply(DisDxBackRegions, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_name")
RepDxBackGenes <- lapply(RepDxBackRegions, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_name")
intersectBack <- intersect(DisDxBackGenes$DisDxAllBack, DisDxBackGenes$DisDxSexAllBack) %>% 
        intersect(DisDxBackGenes$DisDxMalesBack) %>% intersect(DisDxBackGenes$DisDxFemalesBack) %>%
        intersect(RepDxBackGenes$RepDxAllBack) %>% intersect(RepDxBackGenes$RepDxSexAllBack) %>% 
        intersect(RepDxBackGenes$RepDxMalesBack) %>% intersect(RepDxBackGenes$RepDxFemalesBack) %>% unique %>% sort
rm(DisDxBackGenes, RepDxBackGenes)

# Gene Overlap Stats ####
DisRegions_genes <- lapply(DisRegions, function(x) getDMRgeneList(x, regDomains = regDomains, direction = "all", type = "gene_name"))
RepRegions_genes <- lapply(RepRegions, function(x) getDMRgeneList(x, regDomains = regDomains, direction = "all", type = "gene_name"))
gom <- newGOM(gsetA = DisRegions_genes, gsetB = RepRegions_genes, genome.size = length(intersectBack)) # 25088
gomResults <- getMatrix(gom, "intersection") %>% melt
colnames(gomResults) <- c("DiscoveryList", "ReplicationList", "Intersection")
gomResults$OddsRatio <- getMatrix(gom, "odds.ratio") %>% melt %>% .[,"value"]
gomResults$pValue <- getMatrix(gom, "pval") %>% melt %>% .[,"value"]
gomResults$DiscoveryLevel <- as.numeric(gomResults$DiscoveryList)
gomResults$ReplicationLevel <- as.numeric(gomResults$ReplicationList)
gomResults <- subset(gomResults, DiscoveryLevel == ReplicationLevel)
gomResults$DiscoveryLength <- sapply(DisRegions_genes, length)
gomResults$ReplicationLength <- sapply(RepRegions_genes, length)
gomResults$PerDiscovery <- gomResults$Intersection * 100 / gomResults$DiscoveryLength
gomResults$PerReplication <- gomResults$Intersection * 100 / gomResults$ReplicationLength
gomResults$qValue <- p.adjust(gomResults$pValue, method = "fdr")
gomResults <- gomResults[,c("DiscoveryList", "ReplicationList", "DiscoveryLength", "ReplicationLength", "Intersection",
                            "PerDiscovery", "PerReplication", "OddsRatio", "pValue", "qValue")]
write.table(gomResults, file = "Tables/Dis vs Rep DMR Gene Overlap Stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Venn Diagrams ####
# DxAll Venn Diagrams
venn.diagram(list("Discovery" = DisRegions_genes$DisDxAll, "Replication" = RepRegions_genes$RepDxAll),
             file = "Figures/Dis vs Rep DxAll DMR Gene Overlap Venn.png", height = 8, width = 10, imagetype = "png", 
             units = "in", fontfamily = "sans", cat.fontfamily = "sans", fill = c("lightblue", "lightpink"), cex = 2.75, 
             lwd = 4, cat.cex = 3, cat.pos = c(140, 180), cat.dist = c(0.05, 0.03), rotation.degree = 180, margin = 0.06,
             ext.text = FALSE, ext.dist = -0.1, ext.length = 0.9, ext.line.lwd = 2)
venn.diagram(list("Discovery" = DisRegions_genes$DisDxAllHyper, "Replication" = RepRegions_genes$RepDxAllHyper),
             file = "Figures/Dis vs Rep DxAllHyper DMR Gene Overlap Venn.png", height = 8, width = 10, imagetype = "png", 
             units = "in", fontfamily = "sans", cat.fontfamily = "sans", fill = c("lightblue", "lightpink"), cex = 3, 
             lwd = 4, cat.cex = 3, cat.pos = c(165, 180), cat.dist = c(0.03, 0.03), rotation.degree = 180, margin = 0.03,
             ext.text = TRUE, ext.dist = -0.1, ext.length = 0.9, ext.line.lwd = 2)
venn.diagram(list("Discovery" = DisRegions_genes$DisDxAllHypo, "Replication" = RepRegions_genes$RepDxAllHypo),
             file = "Figures/Dis vs Rep DxAllHypo DMR Gene Overlap Venn.png", height = 8, width = 10, imagetype = "png", 
             units = "in", fontfamily = "sans", cat.fontfamily = "sans", fill = c("lightblue", "lightpink"), cex = 3, 
             lwd = 4, cat.cex = 3, cat.pos = c(140, 180), cat.dist = c(0.05, 0.03), rotation.degree = 180, margin = 0.07,
             ext.text = TRUE, ext.dist = -0.1, ext.length = 0.9, ext.line.lwd = 2, ext.percent = 0.01)

# DxSexAll Venn Diagrams
venn.diagram(list("Discovery" = DisRegions_genes$DisDxSexAll, "Replication" = RepRegions_genes$RepDxSexAll),
             file = "Figures/Dis vs Rep DxSexAll DMR Gene Overlap Venn.png", height = 8, width = 10, imagetype = "png", 
             units = "in", fontfamily = "sans", cat.fontfamily = "sans", fill = c("lightblue", "lightpink"), cex = 2.75, 
             lwd = 4, cat.cex = 3, cat.pos = c(140, 180), cat.dist = c(0.05, 0.03), rotation.degree = 180, margin = 0.06,
             ext.text = TRUE, ext.dist = -0.1, ext.length = 0.9, ext.line.lwd = 2, ext.percent = 0.01)
venn.diagram(list("Discovery" = DisRegions_genes$DisDxSexAllHyper, "Replication" = RepRegions_genes$RepDxSexAllHyper),
             file = "Figures/Dis vs Rep DxSexAllHyper DMR Gene Overlap Venn.png", height = 8, width = 10, imagetype = "png", 
             units = "in", fontfamily = "sans", cat.fontfamily = "sans", fill = c("lightblue", "lightpink"), cex = 3, 
             lwd = 4, cat.cex = 3, cat.pos = c(165, 180), cat.dist = c(0.03, 0.03), rotation.degree = 180, margin = 0.04,
             ext.text = TRUE, ext.dist = -0.1, ext.length = 0.9, ext.line.lwd = 2, ext.percent = 0.01)
venn.diagram(list("Discovery" = DisRegions_genes$DisDxSexAllHypo, "Replication" = RepRegions_genes$RepDxSexAllHypo),
             file = "Figures/Dis vs Rep DxSexAllHypo DMR Gene Overlap Venn.png", height = 8, width = 10, imagetype = "png", 
             units = "in", fontfamily = "sans", cat.fontfamily = "sans", fill = c("lightblue", "lightpink"), cex = 3, 
             lwd = 4, cat.cex = 3, cat.pos = c(145, 180), cat.dist = c(0.05, 0.03), rotation.degree = 180, margin = 0.06,
             ext.text = TRUE, ext.dist = -0.1, ext.length = 0.9, ext.line.lwd = 2, ext.percent = 0.01)

# Dx Males Venn Diagrams
venn.diagram(list("Discovery" = DisRegions_genes$DisDxMales, "Replication" = RepRegions_genes$RepDxMales),
             file = "Figures/Dis vs Rep DxMales DMR Gene Overlap Venn.png", height = 8, width = 10, imagetype = "png", 
             units = "in", fontfamily = "sans", cat.fontfamily = "sans", fill = c("lightblue", "lightpink"), cex = 3, 
             lwd = 4, cat.cex = 3, cat.pos = c(150, 180), cat.dist = c(0.04, 0.03), rotation.degree = 180, margin = 0.04,
             ext.text = TRUE, ext.dist = -0.1, ext.length = 0.9, ext.line.lwd = 2, ext.percent = 0.01)
venn.diagram(list("Discovery" = DisRegions_genes$DisDxMalesHyper, "Replication" = RepRegions_genes$RepDxMalesHyper),
             file = "Figures/Dis vs Rep DxMalesHyper DMR Gene Overlap Venn.png", height = 8, width = 10, imagetype = "png", 
             units = "in", fontfamily = "sans", cat.fontfamily = "sans", fill = c("lightblue", "lightpink"), cex = 3, 
             lwd = 4, cat.cex = 3, cat.pos = c(160, 180), cat.dist = c(0.03, 0.03), rotation.degree = 180, margin = 0.04,
             ext.text = TRUE, ext.dist = -0.1, ext.length = 0.9, ext.line.lwd = 2, ext.percent = 0.01)
venn.diagram(list("Discovery" = DisRegions_genes$DisDxMalesHypo, "Replication" = RepRegions_genes$RepDxMalesHypo),
             file = "Figures/Dis vs Rep DxMalesHypo DMR Gene Overlap Venn.png", height = 8, width = 10, imagetype = "png", 
             units = "in", fontfamily = "sans", cat.fontfamily = "sans", fill = c("lightblue", "lightpink"), cex = 3, 
             lwd = 4, cat.cex = 3, cat.pos = c(155, 180), cat.dist = c(0.05, 0.03), rotation.degree = 180, margin = 0.04,
             ext.text = TRUE, ext.dist = -0.1, ext.length = 0.9, ext.line.lwd = 2, ext.percent = 0.01)

# Dx Females Venn Diagrams
venn.diagram(list("Discovery" = DisRegions_genes$DisDxFemales, "Replication" = RepRegions_genes$RepDxFemales),
             file = "Figures/Dis vs Rep DxFemales DMR Gene Overlap Venn.png", height = 8, width = 10, imagetype = "png", 
             units = "in", fontfamily = "sans", cat.fontfamily = "sans", fill = c("lightblue", "lightpink"), cex = 3, 
             lwd = 4, cat.cex = 3, cat.pos = c(160, 180), cat.dist = c(0.04, 0.03), rotation.degree = 180, margin = 0.04,
             ext.text = TRUE, ext.dist = -0.1, ext.length = 0.9, ext.line.lwd = 2, ext.percent = 0.01)
venn.diagram(list("Discovery" = DisRegions_genes$DisDxFemalesHyper, "Replication" = RepRegions_genes$RepDxFemalesHyper),
             file = "Figures/Dis vs Rep DxFemalesHyper DMR Gene Overlap Venn.png", height = 8, width = 10, imagetype = "png", 
             units = "in", fontfamily = "sans", cat.fontfamily = "sans", fill = c("lightblue", "lightpink"), cex = 3, 
             lwd = 4, cat.cex = 3, cat.pos = c(160, 180), cat.dist = c(0.04, 0.03), rotation.degree = 180, margin = 0.04,
             ext.text = TRUE, ext.dist = -0.1, ext.length = 0.9, ext.line.lwd = 2, ext.percent = 0.01)
venn.diagram(list("Discovery" = DisRegions_genes$DisDxFemalesHypo, "Replication" = RepRegions_genes$RepDxFemalesHypo),
             file = "Figures/Dis vs Rep DxFemalesHypo DMR Gene Overlap Venn.png", height = 8, width = 10, imagetype = "png", 
             units = "in", fontfamily = "sans", cat.fontfamily = "sans", fill = c("lightblue", "lightpink"), cex = 3, 
             lwd = 4, cat.cex = 3, cat.pos = c(155, 180), cat.dist = c(0.04, 0.03), rotation.degree = 180, margin = 0.04,
             ext.text = TRUE, ext.dist = -0.1, ext.length = 0.9, ext.line.lwd = 2, ext.percent = 0.01)

# Overlap Replicated Genes ####
# Compare all
intersectGenes <- getNestedList(gom, "intersection")
intersectGenes <- list(intersectGenes[[1]][[1]], intersectGenes[[2]][[2]], intersectGenes[[3]][[3]], 
                       intersectGenes[[4]][[4]], intersectGenes[[5]][[5]], intersectGenes[[6]][[6]], 
                       intersectGenes[[7]][[7]], intersectGenes[[8]][[8]], intersectGenes[[9]][[9]], 
                       intersectGenes[[10]][[10]], intersectGenes[[11]][[11]], intersectGenes[[12]][[12]])
names(intersectGenes) <- c("DxAll", "DxAllHyper", "DxAllHypo", "DxSexAll", "DxSexAllHyper", "DxSexAllHypo",
                           "DxMales", "DxMalesHyper", "DxMalesHypo", "DxFemales", "DxFemalesHyper", "DxFemalesHypo")
intersect_gom <- newGOM(intersectGenes, genome.size = length(intersectBack))
intersect_gomResults <- getMatrix(intersect_gom, "intersection") %>% melt
colnames(intersect_gomResults) <- c("ListA", "ListB", "Intersection")
intersect_gomResults$OddsRatio <- getMatrix(intersect_gom, "odds.ratio") %>% melt %>% .[,"value"]
intersect_gomResults$pValue <- getMatrix(intersect_gom, "pval") %>% melt %>% .[,"value"]
intersect_gomResults$ListA_Length <- sapply(as.character(intersect_gomResults$ListA), function(x) length(intersectGenes[[x]]))
intersect_gomResults$ListB_Length <- sapply(as.character(intersect_gomResults$ListB), function(x) length(intersectGenes[[x]]))
intersect_gomResults$PerListA <- intersect_gomResults$Intersection * 100 / intersect_gomResults$ListA_Length
intersect_gomResults$PerListB <- intersect_gomResults$Intersection * 100 / intersect_gomResults$ListB_Length
intersect_gomResults$ListALevel <- as.numeric(intersect_gomResults$ListA)
intersect_gomResults$ListBLevel <- as.numeric(intersect_gomResults$ListB)
intersect_gomResults <- subset(intersect_gomResults, ListALevel <= ListBLevel) %>% # Remove comparisons that don't make sense
        subset(!(ListA == "DxAll" & ListB %in% c("DxAllHyper", "DxAllHypo"))) %>%
        subset(!(ListA == "DxSexAll" & ListB %in% c("DxSexAllHyper", "DxSexAllHypo"))) %>%
        subset(!(ListA == "DxMales" & ListB %in% c("DxMalesHyper", "DxMalesHypo"))) %>%
        subset(!(ListA == "DxFemales" & ListB %in% c("DxFemalesHyper", "DxFemalesHypo")))
intersect_gomResults$qValue <- p.adjust(intersect_gomResults$pValue, method = "fdr")
intersect_gomResults$log_qValue <- -log10(intersect_gomResults$qValue)
intersect_gomResults$Significant <- intersect_gomResults$log_qValue > -log10(0.05)
intersect_gomResults$Significant <- factor(intersect_gomResults$Significant, levels=c("TRUE", "FALSE"), ordered=TRUE)
intersect_gomResults <- intersect_gomResults[,c("ListA", "ListB", "ListA_Length", "ListB_Length", "Intersection",
                            "PerListA", "PerListB", "OddsRatio", "pValue", "qValue", "log_qValue", "Significant")]
intersect_gomResults$OddsRatio[is.infinite(intersect_gomResults$OddsRatio)] <- max(intersect_gomResults$OddsRatio[!is.infinite(intersect_gomResults$OddsRatio)])
intersect_gomResults$ListA <- factor(intersect_gomResults$ListA, levels = names(intersectGenes))
intersect_gomResults$ListB <- factor(intersect_gomResults$ListB, levels = rev(names(intersectGenes)))
write.table(intersect_gomResults, file = "Tables/Dis vs Rep Overlapping DMR Gene Overlap Stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# DMR genes in all comparisons
intersect(intersectGenes$DxAll, intersectGenes$DxSexAll) %>% intersect(intersectGenes$DxMales) %>% 
        intersect(intersectGenes$DxFemales) %>% sort
# [1] "CHST15"       "CPS1-IT1"     "EML6"         "GRIK2"        "LINC01435"    "LINC01515"    "LOC101928622"
# [8] "LOC644669"    "PPP2R3A"      "ST6GAL2"      "TIMP3" 

# DMR genes in DxAll, DxMales, and DxFemales
intersect(intersectGenes$DxAll, intersectGenes$DxMales) %>% intersect(intersectGenes$DxFemales) %>% sort
# [1] "C7orf62"      "CHST15"       "CPS1-IT1"     "EML6"         "EPHA6"        "GRIK2"        "LHFPL3"      
# [8] "LINC00377"    "LINC01378"    "LINC01435"    "LINC01515"    "LINC01549"    "LOC101928254" "LOC101928622"
# [15] "LOC101929153" "LOC105375972" "LOC441155"    "LOC644669"    "LRRC4C"       "MALRD1"       "MIR125B2"    
# [22] "PABPC4L"      "PCGEM1"       "PIK3R1"       "PPP2R3A"      "RALYL"        "SLC6A15"      "ST6GAL2"     
# [29] "TIMP3"  

# DMR genes in all comparisons by direction
intersect(intersectGenes$DxMalesHyper, intersectGenes$DxFemalesHyper) %>% intersect(intersectGenes$DxAllHyper) %>%
        intersect(intersectGenes$DxSexAllHyper) %>% sort # None
intersect(intersectGenes$DxMalesHypo, intersectGenes$DxFemalesHypo) %>% intersect(intersectGenes$DxAllHypo) %>%
        intersect(intersectGenes$DxSexAllHypo) %>% sort # "LINC01515" "ST6GAL2"  

# DMR genes in DxAll, DxMales, and DxFemales by direction
intersect(intersectGenes$DxMalesHyper, intersectGenes$DxFemalesHyper) %>% intersect(intersectGenes$DxAllHyper) %>% 
        sort # "PIK3R1" "RALYL" 
intersect(intersectGenes$DxMalesHypo, intersectGenes$DxFemalesHypo) %>% intersect(intersectGenes$DxAllHypo) %>%
       sort # "LINC01515" "MALRD1"    "PCGEM1"    "ST6GAL2"  

# DMR genes in males and females by direction
intersect(intersectGenes$DxMalesHyper, intersectGenes$DxFemalesHyper) %>% sort
# [1] "CRAT37"       "LOC101927050" "NPAS3"        "PIK3R1"       "RALYL"        "SV2B"        

intersect(intersectGenes$DxMalesHypo, intersectGenes$DxFemalesHypo) %>% sort
# [1] "ALG10"        "ANXA1"        "ARL5B"        "BRDTP1"       "C2orf27B"     "CDH18"        "COX7B2"      
# [8] "DCAF8L1"      "DIAPH2"       "DIAPH2-AS1"   "FMR1"         "FTHL17"       "FZD1"         "GABRA2"      
# [15] "GRXCR1"       "HS3ST3A1"     "IL1RAPL2"     "LINC00648"    "LINC00889"    "LINC00968"    "LINC01515"   
# [22] "LOC100507201" "LOC100996664" "LOC101927305" "LOC101927358" "LOC101927412" "LOC101928441" "LOC101928519"
# [29] "LOC102723362" "LOC102723427" "LOC102723968" "LOC102724152" "LOC645949"    "LRRTM1"       "MALRD1"      
# [36] "MC4R"         "MEAT6"        "MIR2054"      "MIR3915"      "MIR514A3"     "MIR6134"      "MIR663B"     
# [43] "MTERF1"       "OPCML"        "PCGEM1"       "PMAIP1"       "RBFOX1"       "SLCO5A1"      "SPANXN1"     
# [50] "ST6GAL2"      "TMEM47"       "XRCC6P5"      "ZIC3"     

# Dx Males vs Females Venn Diagrams
venn.diagram(list("Males" = intersectGenes$DxMales, "Females" = intersectGenes$DxFemales),
             file = "Figures/Replicated Male vs Female DMR Gene Overlap Venn.png", height = 8, width = 10, imagetype = "png", 
             units = "in", fontfamily = "sans", cat.fontfamily = "sans", fill = c("lightblue", "lightpink"), cex = 3, 
             lwd = 4, cat.cex = 3, cat.pos = c(165, 180), cat.dist = c(0.03, 0.03), rotation.degree = 180, margin = 0.04,
             ext.text = TRUE, ext.dist = -0.1, ext.length = 0.9, ext.line.lwd = 2, ext.percent = 0.01)
venn.diagram(list("Males" = intersectGenes$DxMalesHyper, "Females" = intersectGenes$DxFemalesHyper),
             file = "Figures/Replicated Male vs Female Hyper DMR Gene Overlap Venn.png", height = 8, width = 10, imagetype = "png", 
             units = "in", fontfamily = "sans", cat.fontfamily = "sans", fill = c("lightblue", "lightpink"), cex = 3, 
             lwd = 4, cat.cex = 3, cat.pos = c(160, 180), cat.dist = c(0.03, 0.03), rotation.degree = 180, margin = 0.04,
             ext.text = TRUE, ext.dist = -0.1, ext.length = 0.9, ext.line.lwd = 2, ext.percent = 0.01)
venn.diagram(list("Males" = intersectGenes$DxMalesHypo, "Females" = intersectGenes$DxFemalesHypo),
             file = "Figures/Replicated Male vs Female Hypo DMR Gene Overlap Venn.png", height = 8, width = 10, imagetype = "png", 
             units = "in", fontfamily = "sans", cat.fontfamily = "sans", fill = c("lightblue", "lightpink"), cex = 3, 
             lwd = 4, cat.cex = 3, cat.pos = c(160, 180), cat.dist = c(0.03, 0.03), rotation.degree = 180, margin = 0.04,
             ext.text = TRUE, ext.dist = -0.1, ext.length = 0.9, ext.line.lwd = 2, ext.percent = 0.01)

# Dx All, Males, Females Venn Diagram
venn.diagram(list("All" = intersectGenes$DxAll, "Males" = intersectGenes$DxMales, "Females" = intersectGenes$DxFemales),
             file = "Figures/Replicated All vs Male vs Female DMR Gene Overlap Venn.png", height = 8, width = 10, imagetype = "png", 
             units = "in", fontfamily = "sans", cat.fontfamily = "sans", fill = c("lightgreen", "lightblue", "lightpink"), cex = 3, 
             lwd = 4, cat.cex = 3, cat.pos = c(323, 0, 45), cat.dist = c(0.063, 0.05, 0.05), rotation.degree = 60, reverse = TRUE, margin = 0.04,
             ext.text = TRUE, ext.dist = -0.1, ext.length = 0.9, ext.line.lwd = 2, ext.percent = 0.01)

# Overlap Plots ####
# All, Hyper, Hypo DMR Gene Heatmaps
plotGOMheatmap(intersect_gomResults, type = "log_qValue", 
               file = "Figures/Dis vs Rep Overlapping DMR Gene Overlap log qvalue Heatmap.png")
plotGOMheatmap(intersect_gomResults, type = "log_OddsRatio", 
               file = "Figures/Dis vs Rep Overlapping DMR Gene Overlap log Odds Ratio Heatmap.png")

# Hyper, Hypo DMR Gene Heatmaps
intersect_gomResults_hyp <- subset(intersect_gomResults, !ListA %in% c("DxAll", "DxSexAll", "DxMales", "DxFemales") & 
                                           !ListB %in% c("DxAll", "DxSexAll", "DxMales", "DxFemales"))
intersect_gomResults_hyp$ListA <- factor(intersect_gomResults_hyp$ListA, 
                                         levels = c("DxAllHyper", "DxAllHypo", "DxSexAllHyper", "DxSexAllHypo", 
                                                    "DxMalesHyper", "DxMalesHypo", "DxFemalesHyper", "DxFemalesHypo"))
intersect_gomResults_hyp$ListB <- factor(intersect_gomResults_hyp$ListB, 
                                         levels = rev(c("DxAllHyper", "DxAllHypo", "DxSexAllHyper", "DxSexAllHypo", 
                                                        "DxMalesHyper", "DxMalesHypo", "DxFemalesHyper", "DxFemalesHypo")))
plotGOMheatmap(intersect_gomResults_hyp, type = "log_qValue", expand = c(0.06, 0), intersect.size = 7, sig.size = 10,
               file = "Figures/Dis vs Rep Overlapping Hyper Hypo DMR Gene Overlap log qvalue Heatmap.png")
plotGOMheatmap(intersect_gomResults_hyp, type = "log_OddsRatio", expand = c(0.06, 0), intersect.size = 7, sig.size = 10,
               file = "Figures/Dis vs Rep Overlapping Hyper Hypo DMR Gene Overlap log Odds Ratio Heatmap.png")

# All DMR Gene Heatmaps
intersect_gomResults_all <- subset(intersect_gomResults, ListA %in% c("DxAll", "DxSexAll", "DxMales", "DxFemales") & 
                                           ListB %in% c("DxAll", "DxSexAll", "DxMales", "DxFemales"))
intersect_gomResults_all$ListA <- factor(intersect_gomResults_all$ListA, 
                                         levels = c("DxAll", "DxSexAll", "DxMales", "DxFemales"))
intersect_gomResults_all$ListB <- factor(intersect_gomResults_all$ListB, 
                                         levels = rev(c("DxAll", "DxSexAll", "DxMales", "DxFemales")))
plotGOMheatmap(intersect_gomResults_all, type = "log_qValue", expand = c(0.15, 0), intersect.size = 9, sig.size = 12,
               axis.text.size = 16, file = "Figures/Dis vs Rep Overlapping All DMR Gene Overlap log qvalue Heatmap.png")
plotGOMheatmap(intersect_gomResults_all, type = "log_OddsRatio", expand = c(0.15, 0), intersect.size = 9, sig.size = 12,
               axis.text.size = 16, file = "Figures/Dis vs Rep Overlapping All DMR Gene Overlap log Odds Ratio Heatmap.png")

# All, Hyper, Hypo Upset Plot
pdf(file = "Figures/Dis vs Rep Overlapping DMR Gene Overlap Upset Plot.pdf", width = 12, height = 8, onefile = FALSE)
upset(fromList(intersectGenes), nsets = 12, nintersects = NA, order.by = "freq", sets.x.label = "Genes", 
      mainbar.y.label = "Genes", text.scale = c(2, 2, 1.75, 1.75, 1.5, 1.05), point.size = 1.5) 
dev.off()

pdf(file = "Figures/Dis vs Rep Overlapping DMR Gene Overlap Upset Plot by degree.pdf", width = 12, height = 8, onefile = FALSE)
upset(fromList(intersectGenes), nsets = 12, nintersects = NA, order.by = "degree", sets.x.label = "Genes", 
      mainbar.y.label = "Genes", text.scale = c(2, 2, 1.75, 1.75, 1.5, 1.05), point.size = 1.5) 
dev.off()

# Hyper, Hypo Upset Plot
intersectGenes_hyp <- intersectGenes[!names(intersectGenes) %in% c("DxAll", "DxSexAll", "DxMales", "DxFemales")]
pdf(file = "Figures/Dis vs Rep Overlapping Hyper Hypo DMR Gene Overlap Upset Plot.pdf", width = 12, height = 8, onefile = FALSE)
upset(fromList(intersectGenes_hyp), nsets = 8, nintersects = NA, order.by = "freq", sets.x.label = "Genes", 
      mainbar.y.label = "Genes", text.scale = c(2, 2, 1.75, 1.75, 1.75, 1.75), point.size = 2.5) 
dev.off()

pdf(file = "Figures/Dis vs Rep Overlapping Hyper Hypo DMR Gene Overlap Upset Plot by degree.pdf", width = 12, height = 8, onefile = FALSE)
upset(fromList(intersectGenes_hyp), nsets = 8, nintersects = NA, order.by = "degree", sets.x.label = "Genes", 
      mainbar.y.label = "Genes", text.scale = c(2, 2, 1.75, 1.75, 1.75, 1.75), point.size = 2.5) 
dev.off()

# All Upset Plot
intersectGenes_all <- intersectGenes[names(intersectGenes) %in% c("DxAll", "DxSexAll", "DxMales", "DxFemales")]
pdf(file = "Figures/Dis vs Rep Overlapping All DMR Gene Overlap Upset Plot.pdf", width = 12, height = 8, onefile = FALSE)
upset(fromList(intersectGenes_all), nsets = 4, nintersects = NA, order.by = "freq", sets.x.label = "Genes", 
      mainbar.y.label = "Genes", text.scale = c(2.5, 2.5, 2, 2, 2.25, 2.5), point.size = 3.5, line.size = 1.5) 
dev.off()

pdf(file = "Figures/Dis vs Rep Overlapping All DMR Gene Overlap Upset Plot by degree.pdf", width = 12, height = 8, onefile = FALSE)
upset(fromList(intersectGenes_all), nsets = 4, nintersects = NA, order.by = "degree", sets.x.label = "Genes", 
      mainbar.y.label = "Genes", text.scale = c(2.5, 2.5, 2, 2, 2.25, 2.5), point.size = 3.5, line.size = 1.5) 
dev.off()

# Hyper Upset Plot
intersectGenes_hyper <- intersectGenes[names(intersectGenes) %in% c("DxAllHyper", "DxSexAllHyper", "DxMalesHyper", "DxFemalesHyper")]
pdf(file = "Figures/Dis vs Rep Overlapping Hyper DMR Gene Overlap Upset Plot by degree.pdf", width = 12, height = 8, onefile = FALSE)
upset(fromList(intersectGenes_hyper), nsets = 4, nintersects = NA, order.by = "degree", sets.x.label = "Genes", 
      mainbar.y.label = "Genes", text.scale = c(2.5, 2.5, 2, 2, 2.25, 2.5), point.size = 3.5, line.size = 1.5) 
dev.off()

# Hypo Upset Plot
intersectGenes_hypo <- intersectGenes[names(intersectGenes) %in% c("DxAllHypo", "DxSexAllHypo", "DxMalesHypo", "DxFemalesHypo")]
pdf(file = "Figures/Dis vs Rep Overlapping Hypo DMR Gene Overlap Upset Plot by degree.pdf", width = 12, height = 8, onefile = FALSE)
upset(fromList(intersectGenes_hypo), nsets = 4, nintersects = NA, order.by = "degree", sets.x.label = "Genes", 
      mainbar.y.label = "Genes", text.scale = c(2.5, 2.5, 2, 2, 2.25, 2.5), point.size = 3.5, line.size = 1.5) 
dev.off()

# Overlap by GREAT --------------------------------------------------------
# Term Intersect ####
# No enriched terms for Discovery DxAll, DxSexAll, so can't overlap these
# Load Data
DisMalesGreat <- read.delim("Tables/Males Diagnosis DMRs Discovery GREAT Combined Results.txt", sep = "\t",
                            stringsAsFactors = FALSE)
DisFemalesGreat <- read.delim("Tables/Females Diagnosis DMRs Discovery GREAT Combined Results.txt", sep = "\t",
                              stringsAsFactors = FALSE)
RepMalesGreat <- read.delim("Tables/Replication Males Diagnosis DMRs GREAT Combined Results.txt", sep = "\t",
                            stringsAsFactors = FALSE)
RepFemalesGreat <- read.delim("Tables/Replication Females Diagnosis 100 DMRs GREAT Combined Results.txt", sep = "\t",
                          stringsAsFactors = FALSE)

# Get terms
DisMalesGreat_names <- list("All" = DisMalesGreat$name[DisMalesGreat$Direction == "All"],
                            "Hyper" = DisMalesGreat$name[DisMalesGreat$Direction == "Hyper"],
                            "Hypo" = DisMalesGreat$name[DisMalesGreat$Direction == "Hypo"])
DisFemalesGreat_names <- list("All" = DisFemalesGreat$name[DisFemalesGreat$Direction == "All"],
                              "Hyper" = DisFemalesGreat$name[DisFemalesGreat$Direction == "Hyper"],
                              "Hypo" = DisFemalesGreat$name[DisFemalesGreat$Direction == "Hypo"])
RepMalesGreat_names <- list("All" = RepMalesGreat$name[RepMalesGreat$Direction == "All"],
                            "Hyper" = RepMalesGreat$name[RepMalesGreat$Direction == "Hyper"],
                            "Hypo" = RepMalesGreat$name[RepMalesGreat$Direction == "Hypo"])
RepFemalesGreat_names <- list("All" = RepFemalesGreat$name[RepFemalesGreat$Direction == "All"],
                              "Hyper" = RepFemalesGreat$name[RepFemalesGreat$Direction == "Hyper"],
                              "Hypo" = RepFemalesGreat$name[RepFemalesGreat$Direction == "Hypo"])

# Intersect terms
MalesGreat_int <- list("All" = intersect(DisMalesGreat_names$All, RepMalesGreat_names$All) %>% sort,
                          "Hyper" = intersect(DisMalesGreat_names$Hyper, RepMalesGreat_names$Hyper) %>% sort,
                          "Hypo" = intersect(DisMalesGreat_names$Hypo, RepMalesGreat_names$Hypo) %>% sort)
# No overlapping terms in males

FemalesGreat_int <- list("All" = intersect(DisFemalesGreat_names$All, RepFemalesGreat_names$All) %>% sort,
                            "Hyper" = intersect(DisFemalesGreat_names$Hyper, RepFemalesGreat_names$Hyper) %>% sort,
                            "Hypo" = intersect(DisFemalesGreat_names$Hypo, RepFemalesGreat_names$Hypo) %>% sort)

# Link with Stats
FemalesGreatStats_int <- merge(x = DisFemalesGreat[DisFemalesGreat$Direction == "All",],
                               y = RepFemalesGreat[RepFemalesGreat$Direction == "All",], 
                               by = "ID", all.x = FALSE, all.y = FALSE, suffixes = c("_Dis", "_Rep"))
write.csv(FemalesGreatStats_int, file = "Tables/Discovery vs Replication Overlapping Females GREAT Stats.csv",
          row.names = FALSE, quote = FALSE)







