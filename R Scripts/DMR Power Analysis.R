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
source("R Scripts/DMR Analysis Functions.R")

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
# HiSeq4000   F    2.349250
#  HiSeqX10   F    6.275188
# HiSeq4000   M    2.431053
#  HiSeqX10   M    6.927919

minCov <- aggregate(CG_coverage ~ Platform + Sex, data = samples, FUN = min)
#  Platform Sex CG_coverage
# HiSeq4000   F       1.216
#  HiSeqX10   F       2.270
# HiSeq4000   M       1.109
#  HiSeqX10   M       2.342

maxCov <- aggregate(CG_coverage ~ Platform + Sex, data = samples, FUN = max)
#  Platform Sex CG_coverage
# HiSeq4000   F       3.324
#  HiSeqX10   F       8.662
# HiSeq4000   M       3.837
#  HiSeqX10   M       9.215

# Power Analysis for Discovery Males -----------------------------------
# Get Data (Cluster) ####
BSmalesDisc <- readRDS("Discovery/Dx_Males/Filtered_BSseq_Discovery50_males.rds")
pData(BSmalesDisc) # TD first
chroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")
background <- loadRegions("Discovery/Dx_Males/bsseq_background_Discovery50_males.csv", chroms = chroms) %>%
        makeGRange()
meth <- getCoverage(BSmalesDisc, regions = background, type = "M", what = "perRegionTotal")
cov <- getCoverage(BSmalesDisc, regions = background, type = "Cov", what = "perRegionTotal")

# Power Analysis (Cluster) ####
DMRresult <- DMR.analysis(N0 = c(39,35), cov.matrix = cov, methyl.matrix = meth, R = c(1,3,5,6.93), pilot.R = 6.93)
saveRDS(DMRresult, "Discovery Males Power Analysis DMR Results.rds")
predictResult <- Estimate.EDR.from.pilot(DMRresult, N0 = c(39,35), target.N = c(5,10,20,40,80,160))
saveRDS(predictResult, "Discovery Males Power Analysis Prediction Results.rds")
# Finished

# Power Analysis for Discovery Females -----------------------------------
# Get Data (Cluster) ####
BSfemalesDisc <- readRDS("Discovery/Dx_Females/Filtered_BSseq_Discovery50_females.rds")
pData(BSfemalesDisc) # TD first
chroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrM")
background <- loadRegions("Discovery/Dx_Females/bsseq_background_Discovery50_females.csv", chroms = chroms) %>%
        makeGRange()
meth <- getCoverage(BSfemalesDisc, regions = background, type = "M", what = "perRegionTotal")
cov <- getCoverage(BSfemalesDisc, regions = background, type = "Cov", what = "perRegionTotal")

# Power Analysis (Cluster) ####
DMRresult <- DMR.analysis(N0 = c(17,15), cov.matrix = cov, methyl.matrix = meth, R = c(1,3,5,6.28), pilot.R = 6.28)
saveRDS(DMRresult, "Discovery Females Power Analysis DMR Results.rds")
predictResult <- Estimate.EDR.from.pilot(DMRresult, N0 = c(17,15), target.N = c(5,10,20,40,80,160))
saveRDS(predictResult, "Discovery Females Power Analysis Prediction Results.rds")
# Finished

# Power Analysis for Replication Males (Not Run) -----------------------------------
# Get Data ####
BSmalesRep <- readRDS("Replication/Dx_Males/Filtered_BSseq_Replication50_males.rds")
pData(BSmalesRep) # TD first
chroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")
background <- loadRegions("Replication/Dx_Males/bsseq_background_Replication50_males.csv", chroms = chroms) %>%
        makeGRange()
meth <- getCoverage(BSmalesRep, regions = background, type = "M", what = "perRegionTotal")
cov <- getCoverage(BSmalesRep, regions = background, type = "Cov", what = "perRegionTotal")

# Power Analysis ####
Sys.time()
DMRresult <- DMR.analysis(N0 = c(17,21), cov.matrix = cov, methyl.matrix = meth, R = c(0.5,1,1.5,2,2.43), pilot.R = 2.43)
saveRDS(DMRresult, "Replication Males Power Analysis DMR Results.rds")

predictResult <- Estimate.EDR.from.pilot(DMRresult, N0 = c(17,21), target.N = c(5,10,20,40,80,160))
saveRDS(predictResult, "Replication Males Power Analysis Prediction Results.rds")
Sys.time()

# Power Analysis for Replication Females (Not Run) -----------------------------------
# Get Data ####
BSfemalesRep <- readRDS("Replication/Dx_Females_100/Filtered_BSseq_Replication100_females.rds")
pData(BSfemalesRep) # TD first
chroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrM")
background <- loadRegions("Replication/Dx_Females_100/bsseq_background_Replication100_females.csv", chroms = chroms) %>%
        makeGRange()
meth <- getCoverage(BSfemalesRep, regions = background, type = "M", what = "perRegionTotal")
cov <- getCoverage(BSfemalesRep, regions = background, type = "Cov", what = "perRegionTotal")

# Power Analysis ####
Sys.time()
DMRresult <- DMR.analysis(N0 = c(3,5), cov.matrix = cov, methyl.matrix = meth, R = c(0.5,1,1.5,2,2.35), pilot.R = 2.35)
saveRDS(DMRresult, "Replication Females Power Analysis DMR Results.rds")

predictResult <- Estimate.EDR.from.pilot(DMRresult, N0 = c(3,5), target.N = c(5,10,20,40,80,160))
saveRDS(predictResult, "Replication Females Power Analysis Prediction Results.rds")
Sys.time()

# Power Analysis Results --------------------------------------------------
# Discovery Males ####
DMRresult_DiscMales <- readRDS("R Objects/Discovery Males Power Analysis DMR Results.rds")
DMRresult_DiscMales$q.values <- p.adjust(DMRresult_DiscMales$p.values, method = "fdr")
table(DMRresult_DiscMales$q.values < 0.05) # 66
predictResult_DiscMales <- readRDS("R Objects/Discovery Males Power Analysis Prediction Results.rds")
# $EDR
#      5 vs 5 10 vs 10 20 vs 20 40 vs 40 80 vs 80 160 vs 160
# 1         0    0e+00   0.0010   0.0015   0.0065     0.0575
# 3         0    9e-04   0.0011   0.0066   0.0796     0.2763
# 5         0    1e-03   0.0022   0.0224   0.1827     0.3938
# 6.93      0    1e-03   0.0039   0.0589   0.2504     0.4498
# 
# $DeclareDMR
#      5 vs 5 10 vs 10 20 vs 20 40 vs 40 80 vs 80 160 vs 160
# 1         0        0        1        2       12        114
# 3         0        1        1       12      157        547
# 5         0        1        4       44      363        784
# 6.93      0        1        7      115      490        882
contourPlot <- plotContour(predictResult_DiscMales)

# Discovery Females ####
DMRresult_DiscFemales <- readRDS("R Objects/Discovery Females Power Analysis DMR Results.rds")
DMRresult_DiscFemales$q.values <- p.adjust(DMRresult_DiscFemales$p.values, method = "fdr")
table(DMRresult_DiscFemales$q.values < 0.05) # 31
predictResult_DiscFemales <- readRDS("R Objects/Discovery Females Power Analysis Prediction Results.rds")
# $EDR
#      5 vs 5 10 vs 10 20 vs 20 40 vs 40 80 vs 80 160 vs 160
# 1         0        0        0        0    0e+00     0.0001
# 3         0        0        0        0    0e+00     0.0006
# 5         0        0        0        0    0e+00     0.0033
# 6.28      0        0        0        0    3e-04     0.0057
# 
# $DeclareDMR
#      5 vs 5 10 vs 10 20 vs 20 40 vs 40 80 vs 80 160 vs 160
# 1         0        0        0        0        0          1
# 3         0        0        0        0        0          2
# 5         0        0        0        0        0          7
# 6.28      0        0        0        0        1         12
contourPlot <- plotContour(predictResult_DiscFemales)
rm(DMRresult_DiscMales, DMRresult_DiscFemales, predictResult_DiscMales, predictResult_DiscFemales)
