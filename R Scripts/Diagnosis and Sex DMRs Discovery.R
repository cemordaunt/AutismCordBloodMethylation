# Diagnosis and Sex DMRs Discovery ####
# Autism Cord Blood Methylation
# Charles Mordaunt
# 6/17/19
# Excluded JLCM032B and JLCM050B
# Updated filtering for DxAdjSex

# cmordaunt@epigenerate:/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Discovery
# THREADS=${SLURM_NTASKS}
# MEM=$(expr ${SLURM_MEM_PER_CPU} / 1024)
# echo "Allocated threads: " $THREADS
# echo "Allocated memory: " $MEM

# Load Packages ####
.libPaths("/share/lasallelab/programs/DMRichR/R_3.6")
sapply(c("tidyverse", "openxlsx", "bsseq", "dmrseq", "DMRichR"), require, character.only = TRUE)

# Global Variables ####
genome <- as.character("hg38")
coverage <- as.numeric(1)
perGroup <- (as.numeric(50)/100)
minCpGs <- as.numeric(3)
maxPerms <- as.numeric(10)
testCovariate <- as.character("Diagnosis")
cores <- 5
set.seed(5)
register(MulticoreParam(1))

# Annotation Databases (Done) ####
sapply(c("BSgenome.Hsapiens.UCSC.hg38", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"), require, 
       character.only = TRUE)
goi <- BSgenome.Hsapiens.UCSC.hg38
TxDb <- TxDb.Hsapiens.UCSC.hg38.knownGene
annoDb <- "org.Hs.eg.db"
annoTrack <- getAnnot(genome)
saveRDS(annoTrack, "hg38_annoTrack.rds")

# Meth ~ Diagnosis, exclude chrX and Y, DMRs ####
# Load and Process Samples (Rerun without JLCM032B and JLCM050B, Complete)
bs.filtered <- processBismark(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                              meta = read.xlsx("sample_info.xlsx", colNames = TRUE) %>% mutate_if(is.character,as.factor),
                              groups = testCovariate, Cov = coverage, mc.cores = cores, per.Group = perGroup)
bs.filtered <- chrSelectBSseq(bs.filtered, seqnames = c(paste("chr", 1:22, sep = ""), "chrM")) # Remove chrX, chrY
saveRDS(bs.filtered, "Dx_All/Filtered_BSseq_Discovery50.rds")

# Background Regions (Rerun without JLCM032B and JLCM050B, Complete)
background <- getBackground(bs.filtered, minNumRegion = minCpGs, maxGap = 1000)
write.table(background, file = "Dx_All/bsseq_background_Discovery50.csv", sep = ",", quote = FALSE, row.names = FALSE)

# DMRs and Raw Methylation (Rerun without JLCM032B and JLCM050B, Complete)
regions <- dmrseq(bs = bs.filtered, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms, testCovariate = testCovariate)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "Dx_All/CandidateRegions_DxNoXY_Discovery50.csv")
gr2csv(sigRegions, "Dx_All/DMRs_DxNoXY_Discovery50.csv")
raw <- as.data.frame(getMeth(BSseq = bs.filtered, regions = sigRegions, type = "raw", what = "perRegion"))
raw <- cbind(sigRegions, raw)
write.table(raw, "Dx_All/DMR_raw_methylation_DxNoXY_Discovery50.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Smoothed Methylation (Rerun without JLCM032B and JLCM050B, Complete)
bs.filtered.bsseq <- BSmooth(bs.filtered, BPPARAM = MulticoreParam(workers = cores, progressbar = TRUE))
pData <- pData(bs.filtered.bsseq)
pData$col <- NULL
pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "#3366CC"
pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <-  "#FF3366"
pData$label <- NULL
pData$label[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "TD"
pData$label[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <-  "ASD"
pData(bs.filtered.bsseq) <- pData
saveRDS(bs.filtered.bsseq, "Dx_All/Filtered_Smoothed_BSseq_Discovery50.rds")
smoothed <- getSmooth(bsseq = bs.filtered.bsseq, regions = sigRegions)
write.table(smoothed, "Dx_All/DMR_smoothed_methylation_DxNoXY_Discovery50.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Plots (Rerun without JLCM032B and JLCM050B, Complete)
annoTrack <- readRDS("hg38_annoTrack.rds")
pdf("Dx_All/DMRs_DxNoXY_Discovery50.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate, 
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2,
          addRegions = sigRegions, annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE,
          horizLegend = TRUE)
dev.off()

pdf("Dx_All/DMRs_DxNoXY_Discovery50_noPoints.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate, 
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2,
          addRegions = sigRegions, annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE,
          horizLegend = TRUE, addPoints = FALSE)
dev.off()

# Meth ~ Diagnosis + AdjSex DMRs and Plots ####
# New processBismark
processBismark <- function(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                           meta = openxlsx::read.xlsx("sample_info.xlsx", colNames = TRUE) %>% dplyr::mutate_if(is.character, as.factor),
                           testCovar = testCovariate,
                           adjustCovar = adjustCovariate,
                           matchCovar = NULL,
                           Cov = coverage,
                           mc.cores = cores,
                           per.Group = perGroup){
        
        cat("\n[DMRichR] Processing Bismark cytosine reports \t\t", format(Sys.time(), "%d-%m-%Y %X"), "\n")
        start_time <- Sys.time()
        print(glue::glue("Selecting files..."))
        files.idx <- pmatch(meta$Name, files)
        files <- files[files.idx]
        #names <- as.data.frame(gsub( "_.*$","", files[files.idx])) # For colData, but jumbles file order with parallel processing
        #colnames(names) <- "Name"
        #rownames(names) <- names[,1]
        #names[,1] <- NULL
        
        # glue::glue("Determining parallelization...") # Does not work on some clusters due to use of BiocParallel, but speeds up desktops 
        # if(mc.cores >= 4){
        #  BPPARAM <- BiocParallel::MulticoreParam(workers = floor(mc.cores/4), progressbar = TRUE)
        #  nThread <- as.integer(floor(mc.cores/floor(mc.cores/4)))
        #  glue::glue("Parallel processing will be used with {floor(mc.cores/4)} cores consisting of {nThread} threads each")
        # }else if(mc.cores < 4){
        #  BPPARAM <- BiocParallel::MulticoreParam(workers = 1, progressbar = TRUE)
        #  nThread <- as.integer(1)
        #  glue::glue("Parallel processing will not be used")
        # }
        
        print(glue::glue("Reading cytosine reports..."))
        bs <- read.bismark(files = files,
                           #colData = names,
                           rmZeroCov = FALSE,
                           strandCollapse = TRUE,
                           verbose = TRUE,
                           BPPARAM = MulticoreParam(workers = mc.cores, progressbar = FALSE), # BPPARAM # bpparam() # MulticoreParam(workers = mc.cores, progressbar = TRUE)
                           nThread = 1) # 1L # nThread
        
        print(glue::glue("Assigning sample metadata with {testCovar} as factor of interest..."))
        sampleNames(bs) <- gsub( "_.*$","", sampleNames(bs))
        meta <- meta[order(match(meta[,1],sampleNames(bs))),]
        stopifnot(sampleNames(bs) == as.character(meta$Name))
        pData(bs) <- cbind(pData(bs), meta[2:length(meta)])
        print(pData(bs))
        
        print(glue::glue("Filtering CpGs..."))
        bs <- GenomeInfoDb::keepStandardChromosomes(bs, pruning.mode = "coarse")
        pData(bs)[[testCovar]] <- as.factor(pData(bs)[[testCovar]])
        loci.cov <- getCoverage(bs, type = "Cov")
        
        if(!is.null(adjustCovar)){
                excludeCovar <- NULL
                for(i in 1:length(adjustCovar)){
                        if(is.numeric(pData(bs)[, adjustCovar[i]]) | is.integer(pData(bs)[, adjustCovar[i]])){
                                print(glue::glue("Assuming adjustment covariate {adjustCovar[i]} is continuous and excluding it from filtering..."))
                                excludeCovar <- c(excludeCovar, adjustCovar[i])
                                
                        }else{
                                print(glue::glue("Assuming adjustment covariate {adjustCovar[i]} is discrete and including it for filtering..."))
                        }
                }
                adjustCovar <- adjustCovar[!adjustCovar %in% excludeCovar]
        }
        
        if(!is.null(matchCovar)){
                if(length(matchCovar) > 1){
                        stop(print(glue::glue("Only one matching covariate can be used")))
                        
                }else if(is.numeric(pData(bs)[, matchCovar]) | is.integer(pData(bs)[, matchCovar])){
                        stop(print(glue::glue("Matching covariate {matchCovar} must be discrete")))
                        
                }else{
                        print(glue::glue("Assuming matching covariate {matchCovar} is discrete and including it for filtering..."))
                }
        }
        
        covar.groups <- apply(pData(bs)[, as.character(c(testCovar, adjustCovar, matchCovar))] %>% as.data.frame(), 
                              MARGIN = 1, FUN = paste, collapse = "_") %>% 
                as.factor() # Covariate combination groups
        group.samples <- split(t(loci.cov >= Cov) %>% as.data.frame(), f = covar.groups) %>% 
                mclapply(FUN = as.matrix, mc.cores = mc.cores) %>% 
                mclapply(FUN = DelayedMatrixStats::colSums2, mc.cores = mc.cores) %>%
                simplify2array() %>%
                as.data.frame() # Samples in each cov.group meeting coverage threshold by CpG (slow)
        
        print(glue::glue("Making coverage filter table..."))
        per.Group.seq <- seq(0,1,0.05)
        covFilter <- NULL
        for(i in 1:length(per.Group.seq)){
                groups.n <- (table(covar.groups) * per.Group.seq[i]) %>% ceiling() %>% as.integer()
                per.Group.seq.test <- mapply(function(x, y){x >= y}, 
                                             x = group.samples, 
                                             y = (table(covar.groups) * per.Group.seq[i]) %>% ceiling() %>% as.integer()) # Test if enough samples are in each group by CpG
                CpGs <- sum(DelayedMatrixStats::rowSums2(per.Group.seq.test) >= length(unique(covar.groups))) # Total CpGs meeting coverage threshold in at least per.Group of all covariate combos
                temp <- c(per.Group.seq[i] * 100, groups.n, CpGs, round(CpGs * 100 / length(bs), 2))
                covFilter <- rbind(covFilter, temp)
        }
        covFilter <- as.data.frame(covFilter, row.names = 1:nrow(covFilter))
        colnames(covFilter) <- c("perGroup", paste("n", unique(covar.groups), sep = "_"), "nCpG", "perCpG")
        print(covFilter)
        
        if(per.Group <= 1){
                print(glue::glue("Filtering for {Cov}x coverage in at least {per.Group*100}% of samples for \\
                                 all combinations of covariates..."))
                sample.idx <- which(pData(bs)[[testCovar]] %in% levels(pData(bs)[[testCovar]]))
                per.Group.test <- mapply(function(x, y){x >= y}, 
                                         x = group.samples, 
                                         y = (table(covar.groups) * per.Group) %>% ceiling() %>% as.integer()) # Test if enough samples are in each group by CpG
                loci.idx <- which(DelayedMatrixStats::rowSums2(per.Group.test) >= length(unique(covar.groups))) # Which CpGs meet coverage threshold in at least per.Group of all covariate combos
                bs.filtered <- bs[loci.idx, sample.idx]
                
        }else if(per.Group > 1){
                stop(print(glue::glue("perGroup is {per.Group} and cannot be greater than 1, which is 100% of samples")))
                
        }else{
                stop(print(glue::glue("processBismark arguments")))
        } 
        
        print(glue::glue("processBismark timing..."))
        end_time <- Sys.time()
        print(end_time - start_time)
        
        print(glue::glue("Before filtering for {Cov}x coverage there were {nrow(bs)} CpGs, \\
                         after filtering there are {nrow(bs.filtered)} CpGs, \\
                         which is {round(nrow(bs.filtered)/nrow(bs)*100,1)}% of all CpGs."))
        
        return(bs.filtered)
}

# Load and Process Samples, using new processBismark with covariate filtering (Running on Epigenerate 6/17, Complete)
adjustCovariate <- "Sex"
bs.filtered <- processBismark(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                              meta = read.xlsx("sample_info.xlsx", colNames = TRUE) %>% mutate_if(is.character,as.factor),
                              testCovar = testCovariate, adjustCovar = adjustCovariate, Cov = coverage, mc.cores = cores, per.Group = perGroup)
bs.filtered <- chrSelectBSseq(bs.filtered, seqnames = c(paste("chr", 1:22, sep = ""), "chrX", "chrM")) # Remove chrY
saveRDS(bs.filtered, "Dx_Sex_All/Filtered_BSseq_Discovery50_DxAdjSex.rds")

# Background Regions (Running on Epigenerate 6/17, Complete)
background <- getBackground(bs.filtered, minNumRegion = minCpGs, maxGap = 1000)
write.table(background, file = "Dx_Sex_All/bsseq_background_Discovery50_DxAdjSex.csv", sep = ",", quote = FALSE, row.names = FALSE)

# DMRs and Raw Methylation (Running on Epigenerate 6/17, need to rerun, Rerunning on epigenerate on 6/23, complete)
bs.filtered <- readRDS("Dx_Sex_All/Filtered_BSseq_Discovery50_DxAdjSex.rds")
regions <- dmrseq(bs = bs.filtered, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms, testCovariate = testCovariate,
                  adjustCovariate = adjustCovariate)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "Dx_Sex_All/CandidateRegions_Discovery50_DxAdjSex.csv")
gr2csv(sigRegions, "Dx_Sex_All/DMRs_Discovery50_DxAdjSex.csv")
raw <- as.data.frame(getMeth(BSseq = bs.filtered, regions = sigRegions, type = "raw", what = "perRegion"))
raw <- cbind(sigRegions, raw)
write.table(raw, "Dx_Sex_All/DMR_raw_methylation_Discovery50_DxAdjSex.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Smoothed Methylation (Running on Epigenerate 6/17, need to rerun, Rerunning on epigenerate on 6/23, complete)
# bs.filtered.bsseq <- BSmooth(bs.filtered, BPPARAM = MulticoreParam(workers = cores, progressbar = TRUE))
# pData <- pData(bs.filtered.bsseq)
# pData$col <- NULL
# pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "#3366CC"
# pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <-  "#FF3366"
# pData$label <- NULL
# pData$label[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "TD"
# pData$label[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <-  "ASD"
# pData(bs.filtered.bsseq) <- pData
# saveRDS(bs.filtered.bsseq, "Dx_Sex_All/Filtered_Smoothed_BSseq_Discovery50_DxAdjSex.rds")

bs.filtered.bsseq <- readRDS("Dx_Sex_All/Filtered_Smoothed_BSseq_Discovery50_DxAdjSex.rds")
smoothed <- getSmooth(bsseq = bs.filtered.bsseq, regions = sigRegions)
write.table(smoothed, "Dx_Sex_All/DMR_smoothed_methylation_Discovery50_DxAdjSex.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Plots (Running on Epigenerate 6/17, need to rerun, rerunning on epigenerate on 6/23, complete)
annoTrack <- readRDS("hg38_annoTrack.rds")
pdf("Dx_Sex_All/DMRs_Discovery50_DxAdjSex.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate, 
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2,
          addRegions = sigRegions, annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE,
          horizLegend = TRUE)
dev.off()

pdf("Dx_Sex_All/DMRs_Discovery50_noPoints_DxAdjSex.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate, 
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2,
          addRegions = sigRegions, annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE,
          horizLegend = TRUE, addPoints = FALSE)
dev.off()

# DMRs with Males Only (Rerun without JLCM032B and JLCM050B, Complete) ####
# New R session
# Load and Process Samples (Rerun without JLCM032B and JLCM050B, Complete)
bs.filtered <- processBismark(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                              meta = read.xlsx("sample_info_males.xlsx", colNames = TRUE) %>% mutate_if(is.character,as.factor),
                              groups = testCovariate, Cov = coverage, mc.cores = cores, per.Group = perGroup)
saveRDS(bs.filtered, "Dx_Males/Filtered_BSseq_Discovery50_males.rds")

# Background Regions (Rerun without JLCM032B and JLCM050B, Complete)
background <- getBackground(bs.filtered, minNumRegion = minCpGs, maxGap = 1000)
write.table(background, file = "Dx_Males/bsseq_background_Discovery50_males.csv", sep = ",", quote = FALSE, row.names = FALSE)

# Meth ~ Diagnosis DMRs and Raw Methylation (Rerun without JLCM032B and JLCM050B, Complete)
regions <- dmrseq(bs = bs.filtered, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms,
                  testCovariate = testCovariate, adjustCovariate = NULL, matchCovariate = NULL)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "Dx_Males/CandidateRegions_Dx_Discovery50_males.csv")
gr2csv(sigRegions, "Dx_Males/DMRs_Dx_Discovery50_males.csv")
raw <- as.data.frame(getMeth(BSseq = bs.filtered, regions = sigRegions, type = "raw", what = "perRegion"))
raw <- cbind(sigRegions, raw)
write.table(raw, "Dx_Males/DMR_raw_methylation_Dx_Discovery50_males.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Smoothed Methylation (Rerun without JLCM032B and JLCM050B, Complete)
bs.filtered.bsseq <- BSmooth(bs.filtered, BPPARAM = MulticoreParam(workers = cores, progressbar = TRUE))
pData <- pData(bs.filtered.bsseq)
pData$col <- NULL
pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "#3366CC"
pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <-  "#FF3366"
pData$label <- NULL
pData$label[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "TD"
pData$label[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <-  "ASD"
pData(bs.filtered.bsseq) <- pData
saveRDS(bs.filtered.bsseq, "Dx_Males/Filtered_Smoothed_BSseq_Discovery50_males.rds")
smoothed <- getSmooth(bsseq = bs.filtered.bsseq, regions = sigRegions)
write.table(smoothed, "Dx_Males/DMR_smoothed_methylation_Dx_Discovery50_males.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# Plots (Rerun without JLCM032B and JLCM050B, Complete)
annoTrack <- readRDS("hg38_annoTrack.rds")
pdf("Dx_Males/DMRs_Dx_Discovery50_males_noPoints.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate,
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2, addRegions = sigRegions, 
          annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE, horizLegend = TRUE, addPoints = FALSE)
dev.off()

pdf("Dx_Males/DMRs_Dx_Discovery50_males.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate,
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2, addRegions = sigRegions, 
          annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE, horizLegend = TRUE)
dev.off()

# DMRs with Females Only ####
# New R session
# Load and Process Samples (Done)
name <- gsub( "_.*$","", list.files(path = getwd(), pattern = "*.txt.gz"))
bs.filtered <- processBismark(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                              meta = read.xlsx("sample_info_females.xlsx", colNames = TRUE) %>% mutate_if(is.character,as.factor),
                              groups = testCovariate, Cov = coverage, mc.cores = cores, per.Group = perGroup)
saveRDS(bs.filtered, "Filtered_BSseq_Discovery50_females.rds")

# Background Regions (Done)
background <- getBackground(bs.filtered, minNumRegion = minCpGs, maxGap = 1000)
write.table(background, file = "bsseq_background_Discovery50_females.csv", sep = ",", quote = FALSE, row.names = FALSE)

# Meth ~ Diagnosis DMRs (Done)
regions <- dmrseq(bs = bs.filtered, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms,
                  testCovariate = testCovariate, adjustCovariate = NULL, matchCovariate = NULL)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "CandidateRegions_Dx_Discovery50_females.csv")
gr2csv(sigRegions, "DMRs_Dx_Discovery50_females.csv")

# Smoothed Methylation and Plots (Done)
bs.filtered.bsseq <- BSmooth(bs.filtered, BPPARAM = MulticoreParam(workers = cores, progressbar = TRUE))
pData <- pData(bs.filtered.bsseq)
pData$col <- NULL
pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "#3366CC"
pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <-  "#FF3366"
pData$label <- NULL
pData$label[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "TD"
pData$label[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <-  "ASD"
pData(bs.filtered.bsseq) <- pData
saveRDS(bs.filtered.bsseq, "Filtered_Smoothed_BSseq_Discovery50_females.rds")

smoothed <- getSmooth(bsseq = bs.filtered.bsseq, regions = sigRegions, 
                      out = "DMR_smoothed_methylation_Dx_Discovery50_females.txt")

bs.filtered.bsseq <- readRDS("Filtered_Smoothed_BSseq_Discovery50_females.rds")
testCovariate <- "Diagnosis"
sigRegions <- read.csv("DMRs_Dx_Discovery50_females.csv", header = TRUE, stringsAsFactors = FALSE)
sigRegions <- data.frame2GRanges(sigRegions)
annoTrack <- readRDS("hg38_annoTrack.rds")
# pasted in plotFunctions.R
pdf("DMRs_Dx_Discovery50_females_noPoints.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate, 
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2, addRegions = sigRegions, 
          annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE, horizLegend = TRUE, 
          addPoints = FALSE)
dev.off()

pdf("DMRs_Dx_Discovery50_females.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate, 
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2, addRegions = sigRegions, 
          annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE, horizLegend = TRUE, 
          addPoints = TRUE)
dev.off()

# Get raw DMR methylation (Done)
bs.filtered <- readRDS("Filtered_BSseq_Discovery50_females.rds")
sigRegions <- read.csv("DMRs_Dx_Discovery50_females.csv", header = TRUE, stringsAsFactors = FALSE)
raw <- as.data.frame(getMeth(BSseq = bs.filtered, regions = data.frame2GRanges(sigRegions), type = "raw", what = "perRegion"))
raw <- cbind(sigRegions, raw)
write.table(raw, "DMR_raw_methylation_Dx_Discovery50_females.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# DMRs with EARLI Samples Only ####
# Includes JLCM032B
# Load and Process Samples (Done)
cores <- 5
name <- gsub( "_.*$","", list.files(path = getwd(), pattern = "*.txt.gz"))
bs.filtered <- processBismark(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                              meta = read.xlsx("sample_info_EARLI.xlsx", colNames = TRUE) %>% mutate_if(is.character,as.factor),
                              groups = testCovariate, Cov = coverage, mc.cores = cores, per.Group = perGroup)
saveRDS(bs.filtered, "Filtered_BSseq_Discovery50_EARLI.rds") # Includes chrX, chrY

# Background Regions (Done)
background <- getBackground(bs.filtered, minNumRegion = minCpGs, maxGap = 1000) # Includes chrX, chrY
write.table(background, file = "bsseq_background_Discovery50_EARLI.csv", sep = ",", quote = FALSE, row.names = FALSE)

# Meth ~ Diagnosis DMRs and Plots (Done)
bs.filtered <- chrSelectBSseq(bs.filtered, seqnames = c(paste("chr", 1:22, sep = ""), "chrM")) # Remove chrX, chrY
regions <- dmrseq(bs = bs.filtered, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms,
                  testCovariate = testCovariate, adjustCovariate = NULL, matchCovariate = NULL)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "CandidateRegions_DxNoXY_Discovery50_EARLI.csv")
gr2csv(sigRegions, "DMRs_DxNoXY_Discovery50_EARLI.csv")

# Smoothed Methylation and Plots (Done)
bs.filtered.bsseq <- BSmooth(bs.filtered, BPPARAM = MulticoreParam(workers = cores, progressbar = TRUE))
pData <- pData(bs.filtered.bsseq)
pData$col <- NULL
pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "#3366CC"
pData$col[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <-  "#FF3366"
pData$label <- NULL
pData$label[pData[,testCovariate] == levels(pData[,testCovariate])[1]] <- "TD"
pData$label[pData[,testCovariate] == levels(pData[,testCovariate])[2]] <-  "ASD"
pData(bs.filtered.bsseq) <- pData
saveRDS(bs.filtered.bsseq, "Filtered_Smoothed_BSseq_Discovery50_EARLI_noXY.rds")

annoTrack <- readRDS("hg38_annoTrack.rds")

smoothed <- getSmooth(bsseq = bs.filtered.bsseq, regions = sigRegions, 
                      out = "DMR_smoothed_methylation_DxNoXY_Discovery50_EARLI.txt")

pdf("DMRs_DxNoXY_Discovery50_EARLI_noPoints.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate,
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2, addRegions = sigRegions, 
          annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE, horizLegend = TRUE, addPoints = FALSE)
dev.off()

pdf("DMRs_DxNoXY_Discovery50_EARLI.pdf", height = 4, width = 8)
plotDMRs2(bs.filtered.bsseq, regions = sigRegions, testCovariate = testCovariate,
          extend = 5000 - (end(sigRegions) - start(sigRegions) + 1)/2, addRegions = sigRegions, 
          annoTrack = annoTrack, lwd = 2, qval = FALSE, stat = FALSE, horizLegend = TRUE)
dev.off()

# Get raw DMR methylation (Done)
bs.filtered <- readRDS("Filtered_BSseq_Discovery50_EARLI.rds")
sigRegions <- read.csv("DMRs_DxNoXY_Discovery50_EARLI.csv", header = TRUE, stringsAsFactors = FALSE)
raw <- as.data.frame(getMeth(BSseq = bs.filtered, regions = data.frame2GRanges(sigRegions), type = "raw", what = "perRegion"))
raw <- cbind(sigRegions, raw)
write.table(raw, "DMR_raw_methylation_DxNoXY_Discovery50_EARLI.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# DMRs with Reduced Number of Males (Rerun without JLCM032B and JLCM050B, Complete) ####
# Reduce number of males (laptop)
set.seed(5)
meta <- read.xlsx("DMRs/Discovery/Diagnosis Males 50/sample_info_males.xlsx", colNames = TRUE) 
# Ctrl_TD Exp_ASD 
#      39      37 
metaRep <- read.xlsx("DMRs/Replication/Diagnosis Males 50/sample_info_males.xlsx", colNames = TRUE)
# Ctrl_TD Exp_ASD 
#      17      21
reducedTD <- sample(meta$Name[meta$Diagnosis == "Ctrl_TD"], size = table(metaRep$Diagnosis)["Ctrl_TD"]) %>% sort
reducedASD <- sample(meta$Name[meta$Diagnosis == "Exp_ASD"], size = table(metaRep$Diagnosis)["Exp_ASD"]) %>% sort
metaReduced <- subset(meta, Name %in% reducedTD | Name %in% reducedASD)
write.xlsx(metaReduced, file = "DMRs/Discovery/Diagnosis Males 50/sample_info_reduced_males.xlsx") # Excluded JLCM032B and JLCM050B from this

# Load and Process Samples (Rerun without JLCM032B and JLCM050B, Complete)
bs.filtered <- processBismark(files = list.files(path = getwd(), pattern = "*.txt.gz"),
                              meta = read.xlsx("sample_info_reduced_males.xlsx", colNames = TRUE) %>% 
                                      mutate_if(is.character, as.factor), groups = testCovariate, Cov = coverage, 
                              mc.cores = cores, per.Group = perGroup)
saveRDS(bs.filtered, "Dx_Reduced_Males/Filtered_BSseq_Discovery50_reduced_males.rds")

# Background Regions (Rerun without JLCM032B and JLCM050B, Complete)
background <- getBackground(bs.filtered, minNumRegion = minCpGs)
write.csv(background, file = "Dx_Reduced_Males/bsseq_background_Discovery50_reduced_males.csv", quote = FALSE, row.names = FALSE)

# Meth ~ Diagnosis DMRs (Rerun without JLCM032B and JLCM050B, Complete)
regions <- dmrseq(bs = bs.filtered, cutoff = 0.05, minNumRegion = minCpGs, maxPerms = maxPerms, testCovariate = testCovariate)
regions$percentDifference <- round(regions$beta/pi * 100)
sigRegions <- regions[regions$pval < 0.05,]
gr2csv(regions, "Dx_Reduced_Males/CandidateRegions_Dx_Discovery50_reduced_males.csv")
gr2csv(sigRegions, "Dx_Reduced_Males/DMRs_Dx_Discovery50_reduced_males.csv")

# DMR Comparison (Rerun without JLCM032B and JLCM050B, Complete) ####
library(ChIPpeakAnno)
source("R Scripts/DMR Analysis Functions.R")
samples <- read.xlsx("DMRs/Discovery/Diagnosis 50/sample_info.xlsx", colNames = TRUE) # JLCM032B and JLCM050B excluded
table(samples$Diagnosis, samples$Sex)
#          F  M
# Ctrl_TD 17 39
# Exp_ASD 15 35

# Load DMRs
regions <- list(All = loadRegions("DMRs/Discovery/Diagnosis 50/DMRs_DxNoXY_Discovery50.csv", 
                                  chroms = c(paste("chr", 1:22, sep = ""), "chrM"), DMRid = TRUE),
                All_Sex = loadRegions("DMRs/Discovery/Diagnosis and Sex 50/DMRs_Discovery50_DxAdjSex.csv", 
                                      chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), DMRid = TRUE),
                Males = loadRegions("DMRs/Discovery/Diagnosis Males 50/DMRs_Dx_Discovery50_males.csv", 
                                    chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), DMRid = TRUE),
                Females = loadRegions("DMRs/Discovery/Diagnosis Females 50/DMRs_Dx_Discovery50_females.csv",
                                      chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), DMRid = TRUE))

# Load Background, all background is without coordinate shifting              
background <- list(All = loadRegions("DMRs/Discovery/Diagnosis 50/bsseq_background_Discovery50.csv", 
                                     chroms = c(paste("chr", 1:22, sep = ""), "chrM")),
                   All_Sex = loadRegions("DMRs/Discovery/Diagnosis and Sex 50/bsseq_background_Discovery50_DxAdjSex.csv", 
                                         chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM")),
                   Males = loadRegions("DMRs/Discovery/Diagnosis Males 50/bsseq_background_Discovery50_males.csv", 
                                       chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")),
                   Females = loadRegions("DMRs/Discovery/Diagnosis Females 50/bsseq_background_Discovery50_females.csv",
                                         chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM")))

# Region Stats
regionStats <- getRegionStats(regions, background = background, n = c(106, 106, 74, 32))
write.table(regionStats, "Tables/DMR Region Stats All Chr Dx Discovery 50.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Write BED files
mapply(function(x, y){writeBED(regions = x, file = y)}, x = regions, 
       y = c("UCSC Tracks/Discovery Diagnosis DMRs.bed", "UCSC Tracks/Discovery AdjSex Diagnosis DMRs.bed",
             "UCSC Tracks/Discovery Diagnosis Males DMRs.bed", "UCSC Tracks/Discovery Diagnosis Females DMRs.bed"))

# Region Overlap Venn Diagrams
GR_regions_hyper <- sapply(regions, makeGRange, direction = "hyper")
GR_regions_hypo <- sapply(regions, makeGRange, direction = "hypo")

DMRoverlapVenn(GR_regions_hyper, NameOfPeaks = c("All", "AllSex", "Males", "Females"), 
               file = "Figures/Hyper DMR Overlap Dx Discovery 50.pdf",
               cat.pos = c(0, 0, 0, 0), cat.dist = c(0.1, 0.1, 0.1, 0.09), 
               fill = c("lightblue", "lightpink", "lightgreen", "orange"))
DMRoverlapVenn(GR_regions_hypo, NameOfPeaks = c("All", "AllSex", "Males", "Females"), 
               file = "Figures/Hypo DMR Overlap Dx Discovery 50.pdf",
               cat.pos = c(0, 0, 0, 0), cat.dist = c(0.1, 0.1, 0.1, 0.09), 
               fill = c("lightblue", "lightpink", "lightgreen", "orange"))

# All vs All Adj Sex Overlap
DMRoverlapVenn(GR_regions_hyper[c("All", "All_Sex")], NameOfPeaks = c("All", "AllSex"), 
               file = "Figures/Hyper DMR Overlap Dx Discovery 50 All vs All AdjSex.pdf", rotation.degree = 0,
               cat.pos = c(0, 15), cat.dist = c(0.02, 0.02), fill = c("lightblue", "lightpink"))
DMRoverlapVenn(GR_regions_hypo[c("All", "All_Sex")], NameOfPeaks = c("All", "AllSex"), 
               file = "Figures/Hypo DMR Overlap Dx Discovery 50 All vs All AdjSex.pdf", rotation.degree = 0,
               cat.pos = c(0, 5), cat.dist = c(0.02, 0.02), fill = c("lightblue", "lightpink"))

# Male Female Overlap
DMRoverlapVenn(GR_regions_hyper[c("Males", "Females")], NameOfPeaks = c("Males", "Females"), 
               file = "Figures/Hyper DMR Overlap Dx Discovery 50 Males vs Females.pdf", rotation.degree = 180,
               cat.pos = c(180, 180), cat.dist = c(0.02, 0.02), fill = c("lightblue", "lightpink"))
DMRoverlapVenn(GR_regions_hypo[c("Males", "Females")], NameOfPeaks = c("Males", "Females"), 
               file = "Figures/Hypo DMR Overlap Dx Discovery 50 Males vs Females.pdf", rotation.degree = 180,
               cat.pos = c(180, 180), cat.dist = c(0.02, 0.02), fill = c("lightblue", "lightpink"))

HyperMF <- GenomicRanges::intersect(x = GR_regions_hyper$Males, y = GR_regions_hyper$Females) %>% as.data.frame
HypoMF <- GenomicRanges::intersect(x = GR_regions_hypo$Males, y = GR_regions_hypo$Females) %>% as.data.frame
MF <- rbind(HyperMF, HypoMF)

HyperMalesDMRsInF <- subset(regions$Males[regions$Males$percentDifference > 0,], GR_regions_hyper$Males %over% GR_regions_hyper$Females)
HypoMalesDMRsInF <- subset(regions$Males[regions$Males$percentDifference < 0,], GR_regions_hypo$Males %over% GR_regions_hypo$Females)
MalesDMRsInF <- rbind(HyperMalesDMRsInF, HypoMalesDMRsInF)
write.table(MalesDMRsInF, "Tables/Dx Discovery Males DMRs Overlapping with Females.txt", sep = "\t", quote = FALSE, row.names = FALSE)

HyperFemalesDMRsInM <- subset(regions$Females[regions$Females$percentDifference > 0,], 
                              GR_regions_hyper$Females %over% GR_regions_hyper$Males)
HypoFemalesDMRsInM <- subset(regions$Females[regions$Females$percentDifference < 0,], 
                             GR_regions_hypo$Females %over% GR_regions_hypo$Males)
FemalesDMRsInM <- rbind(HyperFemalesDMRsInM, HypoFemalesDMRsInM)
write.table(FemalesDMRsInM, "Tables/Dx Discovery Females DMRs Overlapping with Males.txt", sep = "\t", quote = FALSE, row.names = FALSE)

rm(HyperMF, HypoMF, HyperMalesDMRsInF, HypoMalesDMRsInF, HyperFemalesDMRsInM, HypoFemalesDMRsInM, venn)









