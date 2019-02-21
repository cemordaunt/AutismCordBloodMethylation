# Discovery and Replication DMR Comparison --------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 2/19/19

# Packages ####
sapply(c("reshape2", "tidyverse", "ChIPpeakAnno", "annotatr", "GeneOverlap", "UpSetR"), require, character.only = TRUE)

# Functions ####
source("R Scripts/DMR Analysis Functions.R")

# Data ####
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
DisDxAllBack <- loadRegions("DMRs/Discovery/Diagnosis 50/bsseq_background_Discovery50.csv",
                            chroms = c(paste("chr", 1:22, sep = ""), "chrM")) %>% 
        getDMRgeneList(regDomains = regDomains, direction = "all", type = "gene_name")
DisDxSexAllBack <- loadRegions("DMRs/Discovery/Diagnosis and Sex 50/bsseq_background_Discovery50.csv",
                            chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM")) %>% 
        getDMRgeneList(regDomains = regDomains, direction = "all", type = "gene_name")
DisDxMalesBack <- loadRegions("DMRs/Discovery/Diagnosis Males 50/bsseq_background_Discovery50_males.csv",
                            chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")) %>% 
        getDMRgeneList(regDomains = regDomains, direction = "all", type = "gene_name")
DisDxFemalesBack <- loadRegions("DMRs/Discovery/Diagnosis Females 50/bsseq_background_Discovery50_females.csv",
                            chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM")) %>% 
        getDMRgeneList(regDomains = regDomains, direction = "all", type = "gene_name")

# Replication Background
RepDxAllBack <- loadRegions("DMRs/Replication/Diagnosis 50/bsseq_background_Replication50.csv",
                            chroms = c(paste("chr", 1:22, sep = ""), "chrM")) %>% 
        getDMRgeneList(regDomains = regDomains, direction = "all", type = "gene_name")
RepDxSexAllBack <- loadRegions("DMRs/Replication/Diagnosis and Sex 50/bsseq_background_Replication50.csv",
                               chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM")) %>% 
        getDMRgeneList(regDomains = regDomains, direction = "all", type = "gene_name")
RepDxMalesBack <- loadRegions("DMRs/Replication/Diagnosis Males 50/bsseq_background_Replication50_males.csv",
                              chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")) %>% 
        getDMRgeneList(regDomains = regDomains, direction = "all", type = "gene_name")
RepDxFemalesBack <- loadRegions("DMRs/Replication/Diagnosis Females 100/bsseq_background_Replication100_females.csv",
                                chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM")) %>% 
        getDMRgeneList(regDomains = regDomains, direction = "all", type = "gene_name")

intersectBack <- intersect(DisDxAllBack, DisDxSexAllBack) %>% intersect(DisDxMalesBack) %>% intersect(DisDxFemalesBack) %>%
        intersect(RepDxAllBack) %>% intersect(RepDxSexAllBack) %>% intersect(RepDxMalesBack) %>% intersect(RepDxFemalesBack)
length(intersectBack)
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

regDomains <- read.delim("Tables/Regulatory domains hg38.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
DisRepIntersectAnno <- getDMRanno(DisRepIntersect, regDomains = regDomains, file = "Tables/Discovery Replication DMR Intersect with Anno.txt")

# Overlap by Gene ---------------------------------------------------------
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

# DMR genes in males and females by direction
intersect(intersectGenes$DxMalesHyper, intersectGenes$DxFemalesHyper) %>% sort
intersect(intersectGenes$DxMalesHypo, intersectGenes$DxFemalesHypo) %>% sort

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




        





