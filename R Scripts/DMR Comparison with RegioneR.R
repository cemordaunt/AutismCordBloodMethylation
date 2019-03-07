# Discovery and Replication DMR Comparison --------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 2/19/19

# Packages ####
sapply(c("reshape2", "tidyverse", "regioneR", "LOLA"), require, character.only = TRUE)

# Functions ####
loadRegions <- function(file, chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), sort = TRUE){
        if(grepl("txt", file, fixed = TRUE)){
                regions <- read.delim(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
        }
        else{
                regions <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
        }
        if("seqnames" %in% colnames(regions)){
                colnames(regions)[colnames(regions) == "seqnames"] <- "chr"
        }
        regions <- subset(regions, chr %in% chroms)
        regions$chr <- factor(regions$chr, levels = chroms)
        if(sort){
                regions <- regions[order(regions$chr, regions$start),]
        }
        return(regions)
}

DMRpermTest <- function(A, B, genome, universe, Comparison, file){
        message("[DMRpermTest] Performing permutation test of regions using regioneR.")
        pt <- permTest(A = A, B = B, genome = genome, ntimes = 10000, universe = universe, 
                       evaluate.function = c(numOverlaps, meanDistance), randomize.function = resampleRegions, 
                       mc.set.seed = FALSE, force.parallel = TRUE)
        stats <- data.frame("Comparison" = Comparison, "Overlap_observed" = pt$Function1$observed, 
                            "Overlap_zscore" = pt$Function1$zscore, "Overlap_pvalue" = pt$Function1$pval, 
                            "Distance_observed" = pt$Function2$observed, "Distance_zscore" = pt$Function2$zscore, 
                            "Distance_pvalue" = pt$Function2$pval)
        message("[DMRpermTest] Complete! Writing plot and returning stats.")
        pdf(file = file, width = 10, height = 5)
        plot(x = pt, ncol = 2)
        dev.off()
        return(stats)
}

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

