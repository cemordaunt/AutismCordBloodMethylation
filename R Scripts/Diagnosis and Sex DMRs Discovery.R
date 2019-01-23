# Diagnosis and Sex DMRs Discovery ####
# Autism Cord Blood Methylation
# Charles Mordaunt
# 1/22/19

# cmordaunt@epigenerate:/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Discovery

# Functions ####
#' packageLoad
#' @description Install and load desired packages
#' @param packages Character string of desired packages
#' @export packageLoad
packageLoad <- function(packages = packages){
        print(glue::glue("\n","Checking for BiocManager and helpers..."))
        CRAN <- c("BiocManager", "remotes", "magrittr")
        new.CRAN.packages <- CRAN[!(CRAN %in% installed.packages()[,"Package"])]
        if(length(new.CRAN.packages)>0){
                install.packages(new.CRAN.packages, repos ="https://cloud.r-project.org", quiet = TRUE)
        }
        print(glue::glue("Loading package management..."))
        stopifnot(suppressMessages(sapply(CRAN, require, character.only = TRUE)))
        
        new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
        if(length(new.packages)>0){
                print(glue::glue("\n","Installing missing packages..."))
                new.packages <- packages %>%
                        gsub("ggbiplot", "vqv/ggbiplot", .) %>% 
                        gsub("DMRichR", "ben-laufer/DMRichR", .) %>% 
                        gsub("gt", "rstudio/gt", .)
                BiocManager::install(new.packages, ask = FALSE, quiet = TRUE)
        }
        print(glue::glue("Loading packages..."))
        stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))
        suppressWarnings(BiocManager::valid(fix = TRUE, update = TRUE, ask = FALSE))
}

# Load Packages ####
.libPaths("/share/lasallelab/programs/DMRichR/R_3.5")
packageLoad(c("tidyverse", "annotatr", "rGREAT", "enrichR", "ChIPseeker", "BiocParallel", "ggbiplot",
              "liftOver", "openxlsx", "CMplot", "optparse", "gplots", "RColorBrewer", "broom", "lsmeans", "glue"))
BiocManager::install(c("kdkorthauer/dmrseq", "ben-laufer/DMRichR", "rstudio/gt"))
stopifnot(suppressMessages(sapply(c("dmrseq", "DMRichR", "gt"), require, character.only = TRUE)))

# Global Variables ####
genome <- as.character("hg38")
coverage <- as.numeric(1)
perGroup <- (as.numeric(50)/100)
minCpGs <- as.numeric(3)
maxPerms <- as.numeric(10)
testCovariate <- as.character("Diagnosis")
adjustCovariate <- NULL
matchCovariate <- NULL
cores <- 40

# Load and Process Samples ####
name <- gsub( "_.*$","", list.files(path = getwd(), pattern = "*.txt.gz"))
bs.filtered <- processBismark(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                              meta = read.xlsx("sample_info.xlsx", colNames = TRUE) %>% mutate_if(is.character,as.factor),
                              groups = testCovariate, Cov = coverage, mc.cores = cores, per.Group = perGroup)
saveRDS(bs.filtered, "Filtered_BSseq_Discovery50.rds")

# Background Regions ####
background <- getBackground(bs.filtered, minNumRegion = minCpGs, maxGap = 1000)
write.table(background, file = "bsseq_background_Discovery50.csv", sep = ",", quote = FALSE, row.names = FALSE)

# Smoothed Methylation ####
bs.filtered.bsseq <- BSmooth(bs.filtered, BPPARAM = MulticoreParam(workers = cores, progressbar = TRUE))

pData <- pData(bs.filtered.bsseq)
pData$col <- NULL
pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "#3366CC"
pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <-  "#FF3366"
pData$label <- NULL
pData$label[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "TD"
pData$label[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <-  "ASD"
pData(bs.filtered.bsseq) <- pData

saveRDS(bs.filtered.bsseq, "Filtered_Smoothed_BSseq_Discovery50.rds")

# Annotation Databases ####
packageLoad(c("BSgenome.Hsapiens.UCSC.hg38", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"))
goi <- BSgenome.Hsapiens.UCSC.hg38
TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene
annoDb <- "org.Hs.eg.db"
annoTrack <- getAnnot(genome)
saveRDS(annoTrack, "hg38_annoTrack.rds")

# New R session setup ####
.libPaths("/share/lasallelab/programs/DMRichR/R_3.5")
library(bsseq)
library(dmrseq)
library(DMRichR)

bs.filtered <- readRDS("Filtered_BSseq_Discovery50.rds")
set.seed(5)
register(MulticoreParam(1))
minCpGs <- 3
maxPerms <- 10
testCovariate <- "Diagnosis"
cores <- 40
bs.filtered.bsseq <- readRDS("Filtered_Smoothed_BSseq_Discovery50.rds")
annoTrack <- readRDS("hg38_annoTrack.rds")

# Meth ~ Diagnosis DMRs ####
# DMRs (Running)
regionsDx <- dmrseq(bs = bs.filtered, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms,
                    testCovariate = testCovariate, adjustCovariate = NULL, matchCovariate = NULL)
regionsDx$percentDifference <- round(regionsDx$beta/pi * 100)
sigRegionsDx <- regionsDx[regionsDx$pval < 0.05,]
gr2csv(regionsDx, "CandidateRegions_Dx_Discovery50.csv")
gr2csv(sigRegionsDx, "DMRs_Dx_Discovery50.csv")

# Smoothed Methylation
smoothed <- getSmooth(bsseq = bs.filtered.bsseq, regions = sigRegionsDx, 
                      out = "DMR_smoothed_methylation_Dx_Discovery50.txt")

# Plots
pdf("DMRs_Dx_Discovery50.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegionsDx, testCovariate = testCovariate, extend = 5000,
          addRegions = sigRegionsDx, annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE,
          horizLegend = TRUE)
dev.off()

# Meth ~ Diagnosis + AdjSex DMRs ####
# DMRs (Running)
regions <- dmrseq(bs = bs.filtered, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms,
                  testCovariate = testCovariate, adjustCovariate = "Sex", matchCovariate = NULL)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "CandidateRegions_DxAdjSex_Discovery50.csv")
gr2csv(sigRegions, "DMRs_DxAdjSex_Discovery50.csv")

# Smoothed Methylation
smoothed <- getSmooth(bsseq = bs.filtered.bsseq, regions = sigRegions, 
                      out = "DMR_smoothed_methylation_DxAdjSex_Discovery50.txt")

# Plots
pdf("DMRs_DxAdjSex_Discovery50.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate, extend = 5000,
          addRegions = sigRegions, annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE,
          horizLegend = TRUE)
dev.off()

# Meth ~ Diagnosis + MatchSex DMRs ####
# DMRs
regions <- dmrseq(bs = bs.filtered, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms,
                  testCovariate = testCovariate, adjustCovariate = NULL, matchCovariate = "Sex")
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "CandidateRegions_DxMatchSex_Discovery50.csv")
gr2csv(sigRegions, "DMRs_DxMatchSex_Discovery50.csv")

# Smoothed Methylation
smoothed <- getSmooth(bsseq = bs.filtered.bsseq, regions = sigRegions, 
                      out = "DMR_smoothed_methylation_DxMatchSex_Discovery50.txt")

# Plots
pdf("DMRs_DxMatchSex_Discovery50.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate, extend = 5000,
          addRegions = sigRegions, annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE,
          horizLegend = TRUE)
dev.off()

# Meth ~ Diagnosis, Males Only ####
# DMRs
bs.filtered.males <- bs.filtered[,pData(bs.filtered)$Sex == "M"]
regions <- dmrseq(bs = bs.filtered.males, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms,
                  testCovariate = testCovariate, adjustCovariate = NULL, matchCovariate = NULL)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "CandidateRegions_DxMales_Discovery50.csv")
gr2csv(sigRegions, "DMRs_DxMales_Discovery50.csv")

# Smoothed Methylation
bs.filtered.bsseq.males <- bs.filtered.bsseq[,pData(bs.filtered.bsseq)$Sex == "M"]
smoothed <- getSmooth(bsseq = bs.filtered.bsseq.males, regions = sigRegions, 
                      out = "DMR_smoothed_methylation_DxMales_Discovery50.txt")

# Plots
pdf("DMRs_DxMales_Discovery50.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq.males, regions = sigRegions, testCovariate = testCovariate, extend = 5000,
          addRegions = sigRegions, annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE,
          horizLegend = TRUE)
dev.off()

# Meth ~ Diagnosis, Females Only ####
# DMRs
bs.filtered.females <- bs.filtered[,pData(bs.filtered)$Sex == "F"]
regions <- dmrseq(bs = bs.filtered.females, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms,
                  testCovariate = testCovariate, adjustCovariate = NULL, matchCovariate = NULL)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "CandidateRegions_DxFemales_Discovery50.csv")
gr2csv(sigRegions, "DMRs_DxFemales_Discovery50.csv")

# Smoothed Methylation
bs.filtered.bsseq.females <- bs.filtered.bsseq[,pData(bs.filtered.bsseq)$Sex == "F"]
smoothed <- getSmooth(bsseq = bs.filtered.bsseq.females, regions = sigRegions, 
                      out = "DMR_smoothed_methylation_DxFemales_Discovery50.txt")

# Plots
pdf("DMRs_DxFemales_Discovery50.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq.females, regions = sigRegions, testCovariate = testCovariate, extend = 5000,
          addRegions = sigRegions, annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE,
          horizLegend = TRUE)
dev.off()

