# Discovery and Replication DMR Comparison --------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 6/6/19
# Excluded JLCM032B and JLCM050B
# Removed DxAdjSex

# Packages ####
sapply(c("reshape2", "tidyverse", "ChIPpeakAnno", "annotatr", "GeneOverlap", "UpSetR", "parallel", "regioneR", "LOLA",
         "rJava", "RDAVIDWebService", "scales"), require, character.only = TRUE)

# Functions ####
source("R Scripts/DMR Analysis Functions.R")

# Data ####
regDomains <- read.delim("Tables/Regulatory domains hg38.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Discovery DMRs
DisDxAll <- loadRegions("DMRs/Discovery/Diagnosis 50/DMRs_DxNoXY_Discovery50.csv", 
                        chroms = c(paste("chr", 1:22, sep = ""), "chrM"))
DisDxMales <- loadRegions("DMRs/Discovery/Diagnosis Males 50/DMRs_Dx_Discovery50_males.csv", 
                          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"))
DisDxFemales <- loadRegions("DMRs/Discovery/Diagnosis Females 50/DMRs_Dx_Discovery50_females.csv", 
                            chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"))
DisRegions <- list(DisDxAll, subset(DisDxAll, percentDifference > 0), subset(DisDxAll, percentDifference < 0),
                   DisDxMales, subset(DisDxMales, percentDifference > 0), subset(DisDxMales, percentDifference < 0),
                   DisDxFemales, subset(DisDxFemales, percentDifference > 0), subset(DisDxFemales, percentDifference < 0))
names(DisRegions) <- c(paste("DisDxAll", c("", "Hyper", "Hypo"), sep = ""), 
                       paste("DisDxMales", c("", "Hyper", "Hypo"), sep = ""), 
                       paste("DisDxFemales", c("", "Hyper", "Hypo"), sep = ""))
GR_DisRegions <- sapply(DisRegions, function(x) {GRanges(seqnames = x$chr, ranges = IRanges(start = x$start, end = x$end))})

# Replication DMRs
RepDxAll <- loadRegions("DMRs/Replication/Diagnosis 50/DMRs_DxNoXY_Replication50.csv", 
                        chroms = c(paste("chr", 1:22, sep = ""), "chrM"))
RepDxMales <- loadRegions("DMRs/Replication/Diagnosis Males 50/DMRs_Dx_Replication50_males.csv", 
                          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"))
RepDxFemales <- loadRegions("DMRs/Replication/Diagnosis Females 100/DMRs_Dx_Replication100_females.csv", 
                            chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"))
RepRegions <- list(RepDxAll, subset(RepDxAll, percentDifference > 0), subset(RepDxAll, percentDifference < 0),
                   RepDxMales, subset(RepDxMales, percentDifference > 0), subset(RepDxMales, percentDifference < 0),
                   RepDxFemales, subset(RepDxFemales, percentDifference > 0), subset(RepDxFemales, percentDifference < 0))
names(RepRegions) <- c(paste("RepDxAll", c("", "Hyper", "Hypo"), sep = ""),
                       paste("RepDxMales", c("", "Hyper", "Hypo"), sep = ""), 
                       paste("RepDxFemales", c("", "Hyper", "Hypo"), sep = ""))
GR_RepRegions <- sapply(RepRegions, function(x) {GRanges(seqnames = x$chr, ranges = IRanges(start = x$start, end = x$end))})

# Discovery Background
DisDxBackRegions <- list("DisDxAllBack" = loadRegions("DMRs/Discovery/Diagnosis 50/bsseq_background_Discovery50.csv",
                                                      chroms = c(paste("chr", 1:22, sep = ""), "chrM")),
                         "DisDxMalesBack" = loadRegions("DMRs/Discovery/Diagnosis Males 50/bsseq_background_Discovery50_males.csv",
                                                        chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")),
                         "DisDxFemalesBack" = loadRegions("DMRs/Discovery/Diagnosis Females 50/bsseq_background_Discovery50_females.csv",
                                                          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM")))
GR_DisDxBackRegions <- sapply(DisDxBackRegions, function(x) {GRanges(seqnames = x$chr, ranges = IRanges(start = x$start, end = x$end))})


# Replication Background
RepDxBackRegions <- list("RepDxAllBack" = loadRegions("DMRs/Replication/Diagnosis 50/bsseq_background_Replication50.csv",
                                                      chroms = c(paste("chr", 1:22, sep = ""), "chrM")),
                         "RepDxMalesBack" = loadRegions("DMRs/Replication/Diagnosis Males 50/bsseq_background_Replication50_males.csv",
                                                        chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")),
                         "RepDxFemalesBack" = loadRegions("DMRs/Replication/Diagnosis Females 100/bsseq_background_Replication100_females.csv",
                                                          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM")))
GR_RepDxBackRegions <- sapply(RepDxBackRegions, function(x) {GRanges(seqnames = x$chr, ranges = IRanges(start = x$start, end = x$end))})

rm(DisDxAll, DisDxMales, DisDxFemales, RepDxAll, RepDxMales, RepDxFemales)

# Overlap by Location -----------------------------------------------------
# Venn Diagrams ####
# Diagnosis All DMRs
DMRoverlapVenn(list(GR_DisRegions$DisDxAllHyper, GR_RepRegions$RepDxAllHyper), 
               NameOfPeaks = c("Discovery", "Replication"), file = "Figures/Hyper DMR Overlap Dx All Dis vs Rep Venn.pdf",
               rotation.degree = 180, cat.pos = c(0, 0), cat.dist = c(0.03, 0.03))
DMRoverlapVenn(list(GR_DisRegions$DisDxAllHypo, GR_RepRegions$RepDxAllHypo), 
               NameOfPeaks = c("Discovery", "Replication"), file = "Figures/Hypo DMR Overlap Dx All Dis vs Rep Venn.pdf",
               rotation.degree = 180, cat.pos = c(0, 0), cat.dist = c(0.03, 0.03))

# Diagnosis Males DMRs
DMRoverlapVenn(list(GR_DisRegions$DisDxMalesHyper, GR_RepRegions$RepDxMalesHyper), 
               NameOfPeaks = c("Discovery", "Replication"), file = "Figures/Hyper DMR Overlap Dx Males Dis vs Rep Venn.pdf",
               rotation.degree = 180, cat.pos = c(180, 180), cat.dist = c(0.03, 0.03))
DMRoverlapVenn(list(GR_DisRegions$DisDxMalesHypo, GR_RepRegions$RepDxMalesHypo), 
               NameOfPeaks = c("Discovery", "Replication"), file = "Figures/Hypo DMR Overlap Dx Males Dis vs Rep Venn.pdf",
               rotation.degree = 180, cat.pos = c(180, 180), cat.dist = c(0.03, 0.03))

# Diagnosis Females DMRs
DMRoverlapVenn(list(GR_DisRegions$DisDxFemalesHyper, GR_RepRegions$RepDxFemalesHyper), 
               NameOfPeaks = c("Discovery", "Replication"), file = "Figures/Hyper DMR Overlap Dx Females Dis vs Rep Venn.pdf",
               rotation.degree = 180, cat.pos = c(180, 180), cat.dist = c(0.03, 0.03))
DMRoverlapVenn(list(GR_DisRegions$DisDxFemalesHypo, GR_RepRegions$RepDxFemalesHypo), 
               NameOfPeaks = c("Discovery", "Replication"), file = "Figures/Hypo DMR Overlap Dx Females Dis vs Rep Venn.pdf",
               rotation.degree = 180, cat.pos = c(180, 180), cat.dist = c(0.03, 0.03))

# Intersecting Regions ####
DisRepIntersect <- rbind(intersect(GR_DisRegions$DisDxAllHyper, GR_RepRegions$RepDxAllHyper) %>% as.data.frame,
                         intersect(GR_DisRegions$DisDxAllHypo, GR_RepRegions$RepDxAllHypo) %>% as.data.frame,
                         intersect(GR_DisRegions$DisDxMalesHyper, GR_RepRegions$RepDxMalesHyper) %>% as.data.frame,
                         intersect(GR_DisRegions$DisDxMalesHypo, GR_RepRegions$RepDxMalesHypo) %>% as.data.frame,
                         intersect(GR_DisRegions$DisDxFemalesHyper, GR_RepRegions$RepDxFemalesHyper) %>% as.data.frame,
                         intersect(GR_DisRegions$DisDxFemalesHypo, GR_RepRegions$RepDxFemalesHypo) %>% as.data.frame)
colnames(DisRepIntersect)[colnames(DisRepIntersect) == "seqnames"] <- "chr"
DisRepIntersect$DMRid <- paste("DMR", 1:nrow(DisRepIntersect), sep = "_")
DisRepIntersect$Intersect <- c(rep("DxAllHyper", intersect(GR_DisRegions$DisDxAllHyper, GR_RepRegions$RepDxAllHyper) %>% length),
                               rep("DxAllHypo", intersect(GR_DisRegions$DisDxAllHypo, GR_RepRegions$RepDxAllHypo) %>% length),
                               rep("DxMalesHyper", intersect(GR_DisRegions$DisDxMalesHyper, GR_RepRegions$RepDxMalesHyper) %>% length),
                               rep("DxMalesHypo", intersect(GR_DisRegions$DisDxMalesHypo, GR_RepRegions$RepDxMalesHypo) %>% length),
                               rep("DxFemalesHyper", intersect(GR_DisRegions$DisDxFemalesHyper, GR_RepRegions$RepDxFemalesHyper) %>% length),
                               rep("DxFemalesHypo", intersect(GR_DisRegions$DisDxFemalesHypo, GR_RepRegions$RepDxFemalesHypo) %>% length))

DisRepIntersectAnno <- getDMRanno(DisRepIntersect, regDomains = regDomains, file = "Tables/Discovery Replication DMR Intersect with Anno.txt")

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
hg38_X <- getGenomeAndMask("hg38", mask = NULL)$genome %>% 
        filterChromosomes(keep.chr = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"))
stats_DxFemalesHyper <- DMRpermTest(A = GR_DisRegions$DisDxFemalesHyper, B = GR_RepRegions$RepDxFemalesHyper, genome = hg38_X, 
                                   universe = intersect(GR_DisDxBackRegions$DisDxFemalesBack, GR_RepDxBackRegions$RepDxFemalesBack),
                                   file = "Figures/Hyper DMR Overlap Dx Females Dis vs Rep RegioneR Plots.pdf",
                                   Comparison = "DxFemalesHyper")
stats_DxFemalesHypo <- DMRpermTest(A = GR_DisRegions$DisDxFemalesHypo, B = GR_RepRegions$RepDxFemalesHypo, genome = hg38_X, 
                                  universe = intersect(GR_DisDxBackRegions$DisDxFemalesBack, GR_RepDxBackRegions$RepDxFemalesBack),
                                  file = "Figures/Hypo DMR Overlap Dx Females Dis vs Rep RegioneR Plots.pdf",
                                  Comparison = "DxFemalesHypo")

stats_All <- rbind(stats_DxAllHyper, stats_DxAllHypo, stats_DxMalesHyper, stats_DxMalesHypo, stats_DxFemalesHyper, 
                   stats_DxFemalesHypo)
write.table(stats_All, file = "Tables/DMR Overlap Dis vs Rep RegioneR Stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)
rm(stats_DxAllHyper, stats_DxAllHypo, stats_DxMalesHyper, stats_DxMalesHypo, stats_DxFemalesHyper, stats_DxFemalesHypo)

# regioneR Stats with DMRs Redefined as Background Subset ####
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
# Combine
stats_All_redefined <- rbind(stats_DxAllHyper, stats_DxAllHypo, stats_DxMalesHyper, stats_DxMalesHypo, 
                             stats_DxFemalesHyper, stats_DxFemalesHypo)
write.table(stats_All_redefined, file = "Tables/DMR Redefined Overlap Dis vs Rep RegioneR Stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)
rm(stats_DxAllHyper, stats_DxAllHypo, stats_DxMalesHyper, stats_DxMalesHypo, stats_DxFemalesHyper, 
   stats_DxFemalesHypo, A_HyperDMRs, A_HypoDMRs, B_HyperDMRs, B_HypoDMRs, GR_DisDxBackRegions, GR_DisRegions, 
   GR_RepDxBackRegions, GR_RepRegions, hg38_auto, hg38_X, hg38_XY, stats_All, stats_All_redefined, universe_back)

# Overlap by Gene ---------------------------------------------------------
# Background Genes ####
DisDxBackGenes <- lapply(DisDxBackRegions, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_name")
RepDxBackGenes <- lapply(RepDxBackRegions, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_name")
intersectBack <- intersect(DisDxBackGenes$DisDxAllBack, DisDxBackGenes$DisDxMalesBack) %>% 
        intersect(DisDxBackGenes$DisDxFemalesBack) %>% intersect(RepDxBackGenes$RepDxAllBack) %>% 
        intersect(RepDxBackGenes$RepDxMalesBack) %>% intersect(RepDxBackGenes$RepDxFemalesBack) %>% unique %>% sort
rm(DisDxBackGenes, RepDxBackGenes)

# Gene Overlap Stats ####
DisRegions_genes <- lapply(DisRegions[c("DisDxAll", "DisDxMales", "DisDxFemales")], function(x){
        getDMRgeneList(x, regDomains = regDomains, direction = "all", type = "gene_name")})
RepRegions_genes <- lapply(RepRegions[c("RepDxAll", "RepDxMales", "RepDxFemales")], function(x){
        getDMRgeneList(x, regDomains = regDomains, direction = "all", type = "gene_name")})
int_genes <- mapply(intersect, x = DisRegions_genes, y = RepRegions_genes)
names(int_genes) <- c("IntDxAll", "IntDxMales", "IntDxFemales")
write.table(int_genes$DisDxAll, file = "Tables/Discovery and Replication All DMR Genes.txt", sep = "\t", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)
write.table(int_genes$DisDxMales, file = "Tables/Discovery and Replication Males DMR Genes.txt", sep = "\t", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)
write.table(int_genes$DisDxFemales, file = "Tables/Discovery and Replication Females DMR Genes.txt", sep = "\t", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)
write.table(intersect(int_genes$DisDxMales, int_genes$DisDxFemales), file = "Tables/Discovery and Replication Males and Females DMR Genes.txt", sep = "\t", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)

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
# All
geneOverlapVenn(list("Discovery" = DisRegions_genes$DisDxAll, "Replication" = RepRegions_genes$RepDxAll),
                file = "Figures/Dis vs Rep DxAll DMR Gene Overlap Venn.png", cat.pos = c(150, 180), cat.dist = c(0.04, 0.03),
                rotation.degree = 180, ext.dist = -0.1, margin = 0.03, cat.cex = 3)
# Males
geneOverlapVenn(list("Discovery" = DisRegions_genes$DisDxMales, "Replication" = RepRegions_genes$RepDxMales),
                file = "Figures/Dis vs Rep DxMales DMR Gene Overlap Venn.png", cat.pos = c(155, 180), cat.dist = c(0.04, 0.03),
                rotation.degree = 180, ext.dist = -0.1, margin = 0.03, cat.cex = 3)
# Females
geneOverlapVenn(list("Discovery" = DisRegions_genes$DisDxFemales, "Replication" = RepRegions_genes$RepDxFemales),
                file = "Figures/Dis vs Rep DxFemales DMR Gene Overlap Venn.png", cat.pos = c(160, 180), cat.dist = c(0.04, 0.03),
                rotation.degree = 180, ext.dist = -0.1, margin = 0.03, cat.cex = 3)

# Replicated DMR Gene DAVID ####
# Get Replicated DMR Gene IDs
DisRegions_IDs <- lapply(DisRegions[c("DisDxAll", "DisDxMales", "DisDxFemales")], function(x){
        getDMRgeneList(x, regDomains = regDomains, direction = "all", type = "gene_entrezID")})
RepRegions_IDs <- lapply(RepRegions[c("RepDxAll", "RepDxMales", "RepDxFemales")], function(x){
        getDMRgeneList(x, regDomains = regDomains, direction = "all", type = "gene_entrezID")})
intRegions_IDs <- mapply(intersect, x = DisRegions_IDs, y = RepRegions_IDs)
mapply(function(x, y) write.table(x, file = y, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE),
       x = intRegions_IDs, y = c("Tables/Overlapping DMR Gene entrezIDs All Discovery vs Replication.txt",
                             "Tables/Overlapping DMR Gene entrezIDs Males Discovery vs Replication.txt",
                             "Tables/Overlapping DMR Gene entrezIDs Females Discovery vs Replication.txt"))

# Get Replicated DMR Gene IDs by Male-Female Overlap
MFintRegions_IDs <- intersect(intRegions_IDs$DisDxMales, intRegions_IDs$DisDxFemales)
write.table(MFintRegions_IDs, file = "Tables/Overlapping DMR Gene entrezIDs Males and Females Discovery vs Replication.txt", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
MonlyIntRegions_IDs <- intRegions_IDs$DisDxMales[!intRegions_IDs$DisDxMales %in% intRegions_IDs$DisDxFemales]
write.table(MonlyIntRegions_IDs, file = "Tables/Overlapping DMR Gene entrezIDs Males Only Discovery vs Replication.txt", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
FonlyIntRegions_IDs <- intRegions_IDs$DisDxFemales[!intRegions_IDs$DisDxFemales %in% intRegions_IDs$DisDxMales]
write.table(FonlyIntRegions_IDs, file = "Tables/Overlapping DMR Gene entrezIDs Females Only Discovery vs Replication.txt", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# Get Background IDs
DisBack_IDs <- lapply(DisDxBackRegions[c("DisDxAllBack", "DisDxMalesBack", "DisDxFemalesBack")], function(x){
        getDMRgeneList(x, regDomains = regDomains, direction = "all", type = "gene_entrezID")})
RepBack_IDs <- lapply(RepDxBackRegions[c("RepDxAllBack", "RepDxMalesBack", "RepDxFemalesBack")], function(x){
        getDMRgeneList(x, regDomains = regDomains, direction = "all", type = "gene_entrezID")})
intBack_IDs <- mapply(intersect, x = DisBack_IDs, y = RepBack_IDs)
mapply(function(x, y) write.table(x, file = y, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE),
       x = intBack_IDs, y = c("Tables/Overlapping DMR Background Gene entrezIDs All Discovery vs Replication.txt",
                                 "Tables/Overlapping DMR Background Gene entrezIDs Males Discovery vs Replication.txt",
                                 "Tables/Overlapping DMR Background Gene entrezIDs Females Discovery vs Replication.txt"))
MFintBack_IDs <- intersect(intBack_IDs$DisDxMalesBack, intBack_IDs$DisDxFemalesBack)
write.table(MFintBack_IDs, file = "Tables/Overlapping DMR Background Gene entrezIDs Males and Females Discovery vs Replication.txt", 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# Run DAVID
david <- DAVIDWebService$new(url="https://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/",
                             email="cemordaunt@ucdavis.edu")
categories <- getAllAnnotationCategoryNames(david)
categories <- c("BBID", "BIOCARTA", "BIOGRID_INTERACTION", "CGAP_EST_QUARTILE", "CGAP_SAGE_QUARTILE", "CHROMOSOME",                
                "COG_ONTOLOGY", "CYTOBAND", "DIP", "EC_NUMBER", "GAD_DISEASE", "GENE3D", "GNF_U133A_QUARTILE", 
                "GENERIF_SUMMARY", "GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY", "INTERPRO", 
                "INTACT", "OMIM_DISEASE", "MINT", "PIR_SEQ_FEATURE", "PIR_SUMMARY", "PIR_SUPERFAMILY", "PFAM", 
                "REACTOME_PATHWAY", "SMART", "PRINTS", "PRODOM", "PROSITE", "UP_KEYWORDS", "UNIGENE_EST_QUARTILE", 
                "SP_COMMENT_TYPE", "SP_COMMENT", "TIGRFAMS", "SUPFAM", "UP_TISSUE", "UP_SEQ_FEATURE")  

DisDAVID <- mapply(getDAVID, genes = DisRegions_IDs, background = DisBack_IDs, 
                   file = c("Tables/Discovery All DMR Genes DAVID Results.txt", 
                            "Tables/Discovery Males DMR Genes DAVID Results.txt",
                            "Tables/Discovery Females DMR Genes DAVID Results.txt"),
                   MoreArgs = list(categories = categories), SIMPLIFY = FALSE)

RepDAVID <- mapply(getDAVID, genes = RepRegions_IDs, background = RepBack_IDs, 
                   file = c("Tables/Replication All DMR Genes DAVID Results.txt", 
                            "Tables/Replication Males DMR Genes DAVID Results.txt",
                            "Tables/Replication Females DMR Genes DAVID Results.txt"),
                   MoreArgs = list(categories = categories), SIMPLIFY = FALSE)

intDAVID <- mapply(getDAVID, genes = intRegions_IDs, background = intBack_IDs, 
                   file = c("Tables/Overlapping All DMR Genes DAVID Results.txt", 
                            "Tables/Overlapping Males DMR Genes DAVID Results.txt",
                            "Tables/Overlapping Females DMR Genes DAVID Results.txt"),
                   MoreArgs = list(categories = categories), SIMPLIFY = FALSE)

MFoverlap_IDs <- list(MandF = MFintRegions_IDs, Monly = MonlyIntRegions_IDs, Fonly = FonlyIntRegions_IDs)
MFoverlapBack_IDs <- list(MandF = MFintBack_IDs, Monly = intBack_IDs$DisDxMalesBack, Fonly = intBack_IDs$DisDxFemalesBack)
MFoverlapDAVID <- mapply(getDAVID, genes = MFoverlap_IDs, background = MFoverlapBack_IDs, 
                         file = c("Tables/Overlapping Males and Females DMR Genes DAVID Results.txt", 
                                  "Tables/Overlapping Males Only DMR Genes DAVID Results.txt",
                                  "Tables/Overlapping Females Only DMR Genes DAVID Results.txt"),
                         MoreArgs = list(categories = categories), SIMPLIFY = FALSE)

# Overlap Replicated Genes ####
# Compare all
intersect_gom <- newGOM(int_genes, genome.size = length(intersectBack))
intersect_gomResults <- getMatrix(intersect_gom, "intersection") %>% melt
colnames(intersect_gomResults) <- c("ListA", "ListB", "Intersection")
intersect_gomResults$OddsRatio <- getMatrix(intersect_gom, "odds.ratio") %>% melt %>% .[,"value"]
intersect_gomResults$pValue <- getMatrix(intersect_gom, "pval") %>% melt %>% .[,"value"]
intersect_gomResults$ListA_Length <- sapply(as.character(intersect_gomResults$ListA), function(x) length(int_genes[[x]]))
intersect_gomResults$ListB_Length <- sapply(as.character(intersect_gomResults$ListB), function(x) length(int_genes[[x]]))
intersect_gomResults$PerListA <- intersect_gomResults$Intersection * 100 / intersect_gomResults$ListA_Length
intersect_gomResults$PerListB <- intersect_gomResults$Intersection * 100 / intersect_gomResults$ListB_Length
intersect_gomResults <- subset(intersect_gomResults, !as.character(intersect_gomResults$ListA) == 
                                       as.character(intersect_gomResults$ListB))
intersect_gomResults$qValue <- p.adjust(intersect_gomResults$pValue, method = "fdr")
intersect_gomResults$log_qValue <- -log10(intersect_gomResults$qValue)
intersect_gomResults$Significant <- intersect_gomResults$log_qValue > -log10(0.05)
intersect_gomResults$Significant <- factor(intersect_gomResults$Significant, levels=c("TRUE", "FALSE"), ordered=TRUE)
intersect_gomResults <- intersect_gomResults[,c("ListA", "ListB", "ListA_Length", "ListB_Length", "Intersection",
                            "PerListA", "PerListB", "OddsRatio", "pValue", "qValue", "log_qValue", "Significant")]

write.table(intersect_gomResults, file = "Tables/Dis vs Rep Overlapping DMR Gene Overlap Stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# DMR genes in all comparisons
intersect(int_genes$IntDxAll, int_genes$IntDxMales) %>% intersect(int_genes$IntDxFemales) %>% sort
#  [1] "ADGRL3-AS1"   "ALG10B"       "CACNA2D1"     "CPS1-IT1"     "EDIL3"        "EML6"         "EPHA6"        "FAM49A"      
#  [9] "FOXD4L5"      "GRIK2"        "LINC00550"    "LINC01378"    "LINC01435"    "LINC01515"    "LOC100505978" "LOC101928622"
# [17] "LOC101929153" "LOC105375972" "LUZP2"        "MALRD1"       "MIR548XHG"    "PABPC4L"      "PCGEM1"       "PCLO"        
# [25] "PIK3R1"       "PPP2R3A"      "RALYL"        "SLC6A15"      "ST6GAL2"      "TIMP3"       

# DMR genes in males and females
intersect(int_genes$IntDxMales, int_genes$IntDxFemales) %>% sort
#   [1] "ADGRL3-AS1"   "AFF2"         "AKAP6"        "ALG10B"       "ANXA1"        "ARL5B"        "ASTN2"        "ATXN7L3B"    
#   [9] "AUTS2"        "BCL11B"       "BLACE"        "BMS1P22"      "C14orf177"    "C2orf27B"     "C7orf62"      "CACNA2D1"    
#  [17] "CADM2-AS2"    "CBLN4"        "CCT8L2"       "CDH12"        "CDH18"        "CHST15"       "COX7B2"       "CPS1-IT1"    
#  [25] "CPXCR1"       "CPXM2"        "CSMD3"        "CSTF2T"       "DAAM1"        "EDIL3"        "EML6"         "EPHA6"       
#  [33] "EPHB1"        "FAM155B"      "FAM49A"       "FOXD4L5"      "FRG1BP"       "FZD1"         "GABRA2"       "GABRG1"      
#  [41] "GPR135"       "GRIA3"        "GRIK1-AS1"    "GRIK2"        "GRXCR1"       "HS3ST3A1"     "KC6"          "KCNC2"       
#  [49] "KDR"          "KLHL13"       "LINC00269"    "LINC00377"    "LINC00396"    "LINC00550"    "LINC00564"    "LINC00583"   
#  [57] "LINC00613"    "LINC00648"    "LINC00889"    "LINC01378"    "LINC01435"    "LINC01470"    "LINC01491"    "LINC01515"   
#  [65] "LINC01549"    "LINC01598"    "LOC100505978" "LOC100506585" "LOC100507201" "LOC100996664" "LOC101927050" "LOC101927141"
#  [73] "LOC101927305" "LOC101927358" "LOC101927359" "LOC101927406" "LOC101927657" "LOC101927907" "LOC101928137" "LOC101928201"
#  [81] "LOC101928437" "LOC101928441" "LOC101928519" "LOC101928622" "LOC101928820" "LOC101928887" "LOC101929153" "LOC102546299"
#  [89] "LOC102723362" "LOC102723427" "LOC102723968" "LOC105372038" "LOC105374693" "LOC105375972" "LOC339862"    "LOC389906"   
#  [97] "LOC441155"    "LOC644669"    "LOC645949"    "LRFN5"        "LRRC4C"       "LRRCC1"       "LRRTM1"       "LRRTM4"      
# [105] "LUZP2"        "LVCAT1"       "MAGEB16"      "MALRD1"       "MAT2B"        "MC4R"         "MIR125B2"     "MIR2053"     
# [113] "MIR2054"      "MIR3144"      "MIR3672"      "MIR378C"      "MIR3924"      "MIR4275"      "MIR4465"      "MIR548AD"    
# [121] "MIR548XHG"    "MIR595"       "MIR663B"      "MIR7157"      "MTERF1"       "NAV3"         "NPAS3"        "NUBPL"       
# [129] "PABPC4L"      "PCDH11X"      "PCGEM1"       "PCLO"         "PIK3R1"       "PMAIP1"       "PPP2R3A"      "PPP4R3CP"    
# [137] "PTPRD-AS2"    "RALYL"        "RBFOX1"       "RBMS3-AS1"    "RTN4"         "SLC2A5"       "SLC6A15"      "SNORA70C"    
# [145] "SPANXN1"      "SPRY2"        "ST6GAL2"      "SYNPR"        "TBX18"        "TDGF1P3"      "TGFBR2"       "THSD7A"      
# [153] "TIMP3"        "TMEM47"       "TOX"          "TYRP1"        "VTA1"         "XACT"         "XRCC6P5"      "ZC4H2"       
# [161] "ZIC3"         "ZWINT"       

# Dx Males vs Females Venn Diagram
geneOverlapVenn(list("Males" = int_genes$IntDxMales, "Females" = int_genes$IntDxFemales), 
                file = "Figures/Replicated Male vs Female DMR Gene Overlap Venn.png", cat.pos = c(175, 180), 
                cat.dist = c(0.03, 0.03), rotation.degree = 180, ext.dist = -0.1, margin = 0.03, cat.cex = 3)

# Dx All, Males, Females Venn Diagram
geneOverlapVenn(list("All" = int_genes$IntDxAll, "Males" = int_genes$IntDxMales, "Females" = int_genes$IntDxFemales),
             file = "Figures/Replicated All vs Male vs Female DMR Gene Overlap Venn.png", 
             fill = c("lightgreen", "lightblue", "lightpink"), cat.cex = 3, cat.pos = c(323, 0, 45), 
             cat.dist = c(0.063, 0.05, 0.05), rotation.degree = 60, reverse = TRUE, margin = 0.04,
             ext.dist = -0.1)

# Dx Males Females chrX Venn Diagram
chrX_regions <- list(DisDxMales = subset(DisRegions$DisDxMales, chr == "chrX"), 
                     DisDxFemales = subset(DisRegions$DisDxFemales, chr == "chrX"),
                     RepDxMales = subset(RepRegions$RepDxMales, chr == "chrX"), 
                     RepDxFemales = subset(RepRegions$RepDxFemales, chr == "chrX"))

chrX_genes <- lapply(chrX_regions, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_name")
chrX_int <- list(DxMales = intersect(chrX_genes$DisDxMales, chrX_genes$RepDxMales), # 55 genes in males
                 DxFemales = intersect(chrX_genes$DisDxFemales, chrX_genes$RepDxFemales)) # 173 genes in females
geneOverlapVenn(list("Males" = chrX_int$DxMales, "Females" = chrX_int$DxFemales),
             file = "Figures/Replicated Male vs Female chrX DMR Gene Overlap Venn.png", cat.cex = 3,
             cat.pos = c(175, 180), cat.dist = c(0.03, 0.03), rotation.degree = 180, margin = 0.04)

# Start here, stats for chrX intersect, need background for chrX
chrX_Back <- list(DisDxMalesBack = subset(DisDxBackRegions$DisDxMalesBack, chr == "chrX"),
                  DisDxFemalesBack = subset(DisDxBackRegions$DisDxFemalesBack, chr == "chrX"),
                  RepDxMalesBack = subset(RepDxBackRegions$RepDxMalesBack, chr == "chrX"),
                  RepDxFemalesBack = subset(RepDxBackRegions$RepDxFemalesBack, chr == "chrX"))
chrX_Back_genes <- lapply(chrX_Back, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_name")
chrX_intersectBack <- intersect(chrX_Back_genes$DisDxMalesBack, chrX_Back_genes$DisDxFemalesBack) %>%
        intersect(chrX_Back_genes$RepDxMalesBack) %>% intersect(chrX_Back_genes$RepDxFemalesBack) %>% unique %>% sort #1064 genes

chrX_intersect_gom <- newGOM(chrX_int, genome.size = length(chrX_intersectBack))
chrX_intersect_gomResults <- getMatrix(chrX_intersect_gom, "intersection") %>% melt
colnames(chrX_intersect_gomResults) <- "intersection"
chrX_intersect_gomResults$OddsRatio <- getMatrix(chrX_intersect_gom, "odds.ratio") %>% melt %>% .[,"value"]
chrX_intersect_gomResults$pValue <- getMatrix(chrX_intersect_gom, "pval") %>% melt %>% .[,"value"]
#                   intersection OddsRatio      pValue
# DxFemales.DxMales           21  3.476742 4.60506e-05

# Overlap Plots ####

# All DMR Gene Heatmaps
intersect_gomResults$ListA <- gsub("IntDx", replacement = "", x = intersect_gomResults$ListA) %>% factor(levels = c("All", "Males"))
intersect_gomResults$ListB <- gsub("IntDx", replacement = "", x = intersect_gomResults$ListB) %>% factor(levels = c("Females", "Males"))

plotGOMheatmap(intersect_gomResults, type = "log_qValue", expand = c(0.15, 0), intersect.size = 10, sig.size = 12, 
               limits = c(1, max(intersect_gomResults$log_qValue)), axis.text.size = 16, legend.text.size = 12, 
               legend.title.size = 14, plot.margin = c(1, 6, 1, 1), legend.position = c(1.25, 0.75),
               file = "Figures/Dis vs Rep Overlapping All DMR Gene Overlap log qvalue Heatmap.png",
               width = 5, height = 4)

plotGOMheatmap(intersect_gomResults, type = "log_OddsRatio", expand = c(0.15, 0), intersect.size = 10, sig.size = 12, 
               limits = c(0, max(log10(intersect_gomResults$OddsRatio))), axis.text.size = 16, legend.text.size = 12, 
               legend.title.size = 14, plot.margin = c(1, 6, 1, 1), legend.position = c(1.2, 0.75),
               file = "Figures/Dis vs Rep Overlapping All DMR Gene Overlap log Odds Ratio Heatmap.png",
               width = 5, height = 4)

# All Upset Plot
int_genes_all <- int_genes
names(int_genes_all) <- c("All", "Males", "Females")
pdf(file = "Figures/Dis vs Rep Overlapping All DMR Gene Overlap Upset Plot by degree no adjsex.pdf", 
    width = 12, height = 8, onefile = FALSE)
upset(fromList(int_genes_all), nsets = 3, nintersects = NA, order.by = "degree", 
      sets.x.label = "Total Genes", mainbar.y.label = "Genes", text.scale = c(3, 3, 2.5, 2.5, 3, 3), 
      point.size = 3.75, line.size = 1.75) 
dev.off()

intersect(intersectGenes_all$All, intersectGenes_all$Males) %>% intersect(intersectGenes_all$Females)

# Male Female Discovery Replication Upset Plot
DisRep_MF_genes <- list("Males Discovery" = DisRegions_genes$DisDxMales, 
                        "Males Replication" = RepRegions_genes$RepDxMales,
                        "Females Discovery" = DisRegions_genes$DisDxFemales,
                        "Females Replication" = RepRegions_genes$RepDxFemales)
pdf(file = "Figures/Dis vs Rep Males and Females DMR Gene Overlap Upset Plot by degree.pdf", 
    width = 12, height = 8, onefile = FALSE)
upset(fromList(DisRep_MF_genes), nsets = 4, nintersects = NA, order.by = "degree", 
      sets = rev(c("Males Discovery","Males Replication", "Females Discovery", "Females Replication")), keep.order = TRUE,
      sets.x.label = "Total Genes", mainbar.y.label = "Genes", text.scale = c(3, 3, 2.5, 2.5, 3, 2.5), 
      point.size = 3.75, line.size = 1.75) 
dev.off()

# Overlap by GREAT --------------------------------------------------------
# Term Intersect ####
# No enriched terms for Discovery DxAll, DxMales, so can't overlap these
# Load Data
DisFemalesGreat <- read.delim("Tables/Females Diagnosis DMRs Discovery GREAT Combined Results.txt", sep = "\t",
                              stringsAsFactors = FALSE)
RepFemalesGreat <- read.delim("Tables/Replication Females Diagnosis 100 DMRs GREAT Combined Results.txt", sep = "\t",
                          stringsAsFactors = FALSE)

# Get terms
DisFemalesGreat_names <- list("All" = DisFemalesGreat$name[DisFemalesGreat$Direction == "All"],
                              "Hyper" = DisFemalesGreat$name[DisFemalesGreat$Direction == "Hyper"],
                              "Hypo" = DisFemalesGreat$name[DisFemalesGreat$Direction == "Hypo"])
RepFemalesGreat_names <- list("All" = RepFemalesGreat$name[RepFemalesGreat$Direction == "All"],
                              "Hyper" = RepFemalesGreat$name[RepFemalesGreat$Direction == "Hyper"],
                              "Hypo" = RepFemalesGreat$name[RepFemalesGreat$Direction == "Hypo"])

# Intersect terms
FemalesGreat_int <- list("All" = intersect(DisFemalesGreat_names$All, RepFemalesGreat_names$All) %>% sort,
                            "Hyper" = intersect(DisFemalesGreat_names$Hyper, RepFemalesGreat_names$Hyper) %>% sort,
                            "Hypo" = intersect(DisFemalesGreat_names$Hypo, RepFemalesGreat_names$Hypo) %>% sort)

# Link with Stats
FemalesGreatStats_int <- merge(x = DisFemalesGreat[DisFemalesGreat$Direction == "All",],
                               y = RepFemalesGreat[RepFemalesGreat$Direction == "All",], 
                               by = "ID", all.x = FALSE, all.y = FALSE, suffixes = c("_Dis", "_Rep"))
write.table(FemalesGreatStats_int, file = "Tables/Discovery vs Replication Overlapping Females GREAT Stats.txt",
            row.names = FALSE, quote = FALSE, sep = "\t")

# Overlap by DAVID --------------------------------------------------------
# DAVID Results ####
DisDAVID <- list("All" = read.delim("Tables/Discovery All DMR Genes DAVID Results.txt", sep = "\t", stringsAsFactors = FALSE), 
                 "Males" = read.delim("Tables/Discovery Males DMR Genes DAVID Results.txt", sep = "\t", stringsAsFactors = FALSE),
                 "Females" = read.delim("Tables/Discovery Females DMR Genes DAVID Results.txt", sep = "\t", stringsAsFactors = FALSE))
RepDAVID <- list("All" = read.delim("Tables/Replication All DMR Genes DAVID Results.txt", sep = "\t", stringsAsFactors = FALSE), 
                 "Males" = read.delim("Tables/Replication Males DMR Genes DAVID Results.txt", sep = "\t", stringsAsFactors = FALSE),
                 "Females" = read.delim("Tables/Replication Females DMR Genes DAVID Results.txt", sep = "\t", stringsAsFactors = FALSE))

# Add Matching Term
DisDAVID <- lapply(DisDAVID, function(x){
        x$Category_Term <- paste(x$Category, x$Term, sep = " ")
        return(x)})
RepDAVID <- lapply(RepDAVID, function(x){
        x$Category_Term <- paste(x$Category, x$Term, sep = " ")
        return(x)})

# Merge and Write Files
IntDAVID <- mapply(merge, x = DisDAVID, y = RepDAVID, MoreArgs = list(by = "Category_Term", all = FALSE), SIMPLIFY = FALSE)
mapply(write.table, x = IntDAVID, file = c("Tables/Overlapping DAVID Terms in All Discovery and Replication.txt",
                                           "Tables/Overlapping DAVID Terms in Males Discovery and Replication.txt",
                                           "Tables/Overlapping DAVID Terms in Females Discovery and Replication.txt"),
       MoreArgs = list(sep = "\t", quote = FALSE, row.names = FALSE))
sapply(DisDAVID, nrow)
# All   Males Females 
#  10      39     120 
sapply(RepDAVID, nrow)
# All   Males Females 
# 260     216     182 
sapply(IntDAVID, nrow)
# All   Males Females 
#   9      36      78 

# Overlapping Terms Plots ####
plotDAVID <- read.delim("Tables/Discovery and Replication Males Females DAVID Results for Plot.txt", sep = "\t",
                        header = TRUE, stringsAsFactors = FALSE)
plotDAVID$Comparison[plotDAVID$Comparison == "Males Discovery"] <- "Males\nDiscovery"
plotDAVID$Comparison[plotDAVID$Comparison == "Males Replication"] <- "Males\nReplication"
plotDAVID$Comparison[plotDAVID$Comparison == "Females Discovery"] <- "Females\nDiscovery"
plotDAVID$Comparison[plotDAVID$Comparison == "Females Replication"] <- "Females\nReplication"
plotDAVID$Comparison <- factor(plotDAVID$Comparison, levels = c("Males\nDiscovery", "Males\nReplication",
                                                                "Females\nDiscovery", "Females\nReplication"))
plotDAVIDexp <- subset(plotDAVID, grepl("expression", plotDAVID$Term, fixed = TRUE))
plotDAVID <- subset(plotDAVID, !grepl("expression", plotDAVID$Term, fixed = TRUE))

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
              plot.margin = unit(c(0.5, 6.5, 0.5, 0.5), "lines"), axis.ticks = element_line(size = 0.7), 
              axis.text.x = element_text(size = 12, color = "black"), 
              axis.text.y = element_text(size = 12, color = "black"),
              axis.title = element_blank(), legend.key = element_blank(),  
              legend.position = c(1.2, 0.88), legend.background = element_blank(), 
              legend.key.size = unit(1, "lines"), legend.title = element_text(size = 14), 
              legend.text = element_text(size = 13), panel.background = element_rect(fill = "black"))
ggsave("Figures/Discovery and Replication Males Females DAVID logp Heatmap.png", plot = gg, dpi = 600, width = 8, 
       height = 6, units = "in")

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
              plot.margin = unit(c(0.5, 6.5, 0.5, 4.6), "lines"), axis.ticks = element_line(size = 0.7), 
              axis.text.x = element_text(size = 12, color = "black"), 
              axis.text.y = element_text(size = 12, color = "black"),
              axis.title = element_blank(), legend.key = element_blank(),  
              legend.position = c(1.2, 0.88), legend.background = element_blank(), 
              legend.key.size = unit(1, "lines"), legend.title = element_text(size = 14), 
              legend.text = element_text(size = 13), panel.background = element_rect(fill = "black"))
ggsave("Figures/Discovery and Replication Males Females DAVID Expression logp Heatmap.png", plot = gg, dpi = 600, width = 8, 
       height = 6, units = "in")

# DMR Percent of Background by Direction Stacked Barplots -----------------
# DMR Percent of Background Hyper and Hypomethylated All Chroms
DiscDMRs <- list(DisRegions$DisDxAll, DisRegions$DisDxMales, DisRegions$DisDxFemales)
RepDMRs <- list(RepRegions$RepDxAll, RepRegions$RepDxMales, RepRegions$RepDxFemales)

DMRbackPercent <- mapply(function(x, y){c("Hyper  " = sum(x$width[x$percentDifference > 0]) * 100 / sum(y$width), 
                                          "Hypo" = sum(x$width[x$percentDifference < 0]) * 100 / sum(y$width))},
                         x = append(DiscDMRs, RepDMRs), y = append(DisDxBackRegions, RepDxBackRegions))
colnames(DMRbackPercent) <- c("Discovery All", "Discovery Males", "Discovery Females", 
                              "Replication All", "Replication Males", "Replication Females")
DMRbackPercent <- melt(DMRbackPercent)
colnames(DMRbackPercent) <- c("Direction", "Comparison", "Percent")
DMRbackPercent$Comparison <- factor(DMRbackPercent$Comparison, levels = c("Discovery All", "Discovery Males", "Discovery Females", 
                                                                          "Replication All", "Replication Males", "Replication Females"))
gg <- ggplot(DMRbackPercent, aes(x = Comparison, y = Percent, fill = Direction, color = Direction))
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
        ylab("DMR Width (% of Background)") +
        scale_fill_manual(breaks = c("Hyper  ", "Hypo"), values = c("#FF3366", "#3366CC")) +
        scale_color_manual(breaks = c("Hyper  ", "Hypo"), values = c("#FF3366", "#3366CC")) +
        coord_cartesian(ylim = c(0, 0.85)) +
        scale_y_continuous(expand = c(0.004, 0), breaks = pretty_breaks(n = 4))
ggsave("Figures/DMR Percent of Background by Direction Stacked Barplot.png", dpi = 600, width = 5, height = 7, 
       units = "in")

# DMR Percent of Background Hyper and Hypomethylated Autosomes
DMRs_auto <- lapply(append(DiscDMRs, RepDMRs), subset, chr %in% paste("chr", 1:22, sep = ""))
names(DMRs_auto) <- c("Disc_All", "Disc_Males", "Disc_Females", "Rep_All", "Rep_Males", "Rep_Females")
DMRbackground_auto <- lapply(append(DisDxBackRegions, RepDxBackRegions), subset, chr %in% paste("chr", 1:22, sep = ""))
names(DMRbackground_auto) <- c("Disc_All", "Disc_Males", "Disc_Females", "Rep_All", "Rep_Males", "Rep_Females")

DMRbackPercent <- mapply(function(x, y){c("Hyper  " = sum(x$width[x$percentDifference > 0]) * 100 / sum(y$width), 
                                          "Hypo" = sum(x$width[x$percentDifference < 0]) * 100 / sum(y$width))},
                         x = DMRs_auto, y = DMRbackground_auto)
colnames(DMRbackPercent) <- c("Discovery All", "Discovery Males", "Discovery Females", 
                              "Replication All", "Replication Males", "Replication Females")
DMRbackPercent <- melt(DMRbackPercent)
colnames(DMRbackPercent) <- c("Direction", "Comparison", "Percent")
DMRbackPercent$Comparison <- factor(DMRbackPercent$Comparison, levels = c("Discovery All", "Discovery Males", "Discovery Females", 
                                                                          "Replication All", "Replication Males", "Replication Females"))
gg <- ggplot(DMRbackPercent, aes(x = Comparison, y = Percent, fill = Direction, color = Direction))
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
        ylab("DMR Width (% of Background)") +
        scale_fill_manual(breaks = c("Hyper  ", "Hypo"), values = c("#FF3366", "#3366CC")) +
        scale_color_manual(breaks = c("Hyper  ", "Hypo"), values = c("#FF3366", "#3366CC")) +
        coord_cartesian(ylim = c(0, 0.85)) +
        scale_y_continuous(expand = c(0.004, 0), breaks = pretty_breaks(n = 4))
ggsave("Figures/DMR Percent of Background by Direction Autosomes Only Stacked Barplot.png", dpi = 600, width = 5, height = 7, 
       units = "in")

# DMR Percent of Background Hyper and Hypomethylated chrX
DMRs_chrX <- lapply(append(DiscDMRs, RepDMRs), subset, chr == "chrX")
names(DMRs_chrX) <- c("Disc_All", "Disc_Males", "Disc_Females", "Rep_All", "Rep_Males", "Rep_Females")
DMRbackground_chrX <- lapply(append(DisDxBackRegions, RepDxBackRegions), subset, chr == "chrX")
names(DMRbackground_chrX) <- c("Disc_All", "Disc_Males", "Disc_Females", "Rep_All", "Rep_Males", "Rep_Females")

DMRbackPercent <- mapply(function(x, y){c("Hyper  " = sum(x$width[x$percentDifference > 0]) * 100 / sum(y$width), 
                                          "Hypo" = sum(x$width[x$percentDifference < 0]) * 100 / sum(y$width))},
                         x = DMRs_chrX, y = DMRbackground_chrX)
colnames(DMRbackPercent) <- c("Discovery All", "Discovery Males", "Discovery Females", 
                              "Replication All", "Replication Males", "Replication Females")
DMRbackPercent <- melt(DMRbackPercent)
colnames(DMRbackPercent) <- c("Direction", "Comparison", "Percent")
DMRbackPercent$Comparison <- factor(DMRbackPercent$Comparison, levels = c("Discovery All", "Discovery Males", "Discovery Females", 
                                                                          "Replication All", "Replication Males", "Replication Females"))
gg <- ggplot(DMRbackPercent, aes(x = Comparison, y = Percent, fill = Direction, color = Direction))
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
        ylab("DMR Width (% of Background)") +
        scale_fill_manual(breaks = c("Hyper  ", "Hypo"), values = c("#FF3366", "#3366CC")) +
        scale_color_manual(breaks = c("Hyper  ", "Hypo"), values = c("#FF3366", "#3366CC")) +
        coord_cartesian(ylim = c(0, 0.85)) +
        scale_y_continuous(expand = c(0.004, 0), breaks = pretty_breaks(n = 4))
ggsave("Figures/DMR Percent of Background by Direction chrX Only Stacked Barplot.png", dpi = 600, width = 5, height = 7, 
       units = "in")

