# DMR Power Analysis ------------------------------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 3/6/20

# Load Packages ####
.libPaths("/share/lasallelab/Charles/R")
sapply(c("tidyverse", "DropletUtils", "devtools","bsseq", "MethylSeqDesign"), require, character.only = TRUE)

# Functions ####
loadRegions <- function(file, chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), sort = TRUE, 
                        DMRid = FALSE, DMBid = FALSE){
        if(grepl("txt", file, fixed = TRUE)){
                regions <- read.delim(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
        }
        else{
                regions <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
        }
        if(nrow(regions) == 0){
                stop(c(paste("There are no regions in", file, sep = " ")))
        }
        if("seqnames" %in% colnames(regions)){
                colnames(regions)[colnames(regions) == "seqnames"] <- "chr"
        }
        regions <- subset(regions, chr %in% chroms)
        regions$chr <- factor(regions$chr, levels = chroms)
        if(sort){
                regions <- regions[order(regions$chr, regions$start),]
        }
        if(DMRid){
                regions$DMRid <- paste("DMR", 1:nrow(regions), sep = "_")
        }
        if(DMBid){
                regions$DMBid <- paste("DMB", 1:nrow(regions), sep = "_")
        }
        return(regions)
}

makeGRange <- function(DMRs, direction = c("all", "hyper", "hypo")){
        direction <- match.arg(direction)
        if(direction == "hyper"){
                DMRs <- subset(DMRs, percentDifference > 0)
        }
        if(direction == "hypo"){
                DMRs <- subset(DMRs, percentDifference < 0)
        }
        if("chr" %in% colnames(DMRs)){
                GR <- GRanges(seqnames = DMRs$chr, ranges = IRanges(start = DMRs$start, end = DMRs$end))
        } else {
                GR <- GRanges(seqnames = DMRs$seqnames, ranges = IRanges(start = DMRs$start, end = DMRs$end))
        }
        return(GR)
}

# Get Sample Numbers and Coverage ------------------------------------------------
# Sample Numbers ####
samples <- read.csv("Tables/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo and Cell Type.csv", 
                    stringsAsFactors = FALSE)

sampleN <- with(samples, table(Diagnosis_Alg, Platform, Sex))
# Males
# Diagnosis_Alg HiSeq4000 HiSeqX10
#           ASD        21       35
#           TD         17       39

# Females
# Diagnosis_Alg HiSeq4000 HiSeqX10
#           ASD         5       15
#           TD          3       17 

# Coverage ####
meanCov <- aggregate(CG_coverage ~ Platform + Sex, data = samples, FUN = mean)
#  Platform Sex CG_coverage
# HiSeq4000   M    2.431053
#  HiSeqX10   M    6.927919
# HiSeq4000   F    2.349250
#  HiSeqX10   F    6.275188

# Power Analysis for Discovery Males -----------------------------------
# Get Data ####
BSmalesDisc <- readRDS("Discovery/Dx_Males/Filtered_BSseq_Discovery50_males.rds")
pData(BSmalesDisc) # TD first
chroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")
background <- loadRegions("Discovery/Dx_Males/bsseq_background_Discovery50_males.csv", chroms = chroms) %>%
        makeGRange()
meth <- getCoverage(BSmalesDisc, regions = background, type = "M", what = "perRegionTotal")
cov <- getCoverage(BSmalesDisc, regions = background, type = "Cov", what = "perRegionTotal")

# Power Analysis ####
# Started 2020-03-08 18:15:47 PDT
DMRresult <- DMR.analysis(N0 = c(39,35), cov.matrix = cov, methyl.matrix = meth, R = c(1,3,5,6.93), pilot.R = 6.93)
predictResult <- Estimate.EDR.from.pilot(DMRresult, N0 = c(39,35), target.N = seq(5,100,5))
plot <- plotContour(predictResult)
