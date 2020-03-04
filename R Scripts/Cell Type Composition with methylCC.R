# Cell Type Composition with methylCC -------------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 2/25/20

# Previous Workflow ####
#' 1. Get 450K array probes and liftover locations to hg38. Exclude probes mapping to multiple locations
#' 2. Get methylated and total reads at all 450k probes for all samples
#' 3. Subset methylation data for probes included in FlowSorted.CordBlood.450k.ModelPars
#' 4. Subset methylation data for probes covered in > 90% of samples
#' 5. Replace missing values with mean value for that probe
#' 6. Subset Model Pars to match covered probes
#' 7. Estimate cell composition with minfi::projectCellType()
#' 8. Scale cell composition to 100% for each sample

# Packages ####
.libPaths("/share/lasallelab/Charles/R")
sapply(c("tidyverse", "wesanderson", "openxlsx", "rlist", "FlowSorted.CordBlood.450k", "IlluminaHumanMethylation450kmanifest", 
         "IlluminaHumanMethylation450kanno.ilmn12.hg19", "genefilter", "rtracklayer", "AnnotationHub", "minfi", 
         "bsseq", "methylCC"), require, character.only = TRUE)

# Functions ####
find_dmrs <- function(verbose = TRUE, gr_target = NULL, include_cpgs = FALSE, include_dmrs = TRUE, num_cpgs = 50, 
                      num_regions = 50, bumphunter_beta_cutoff = 0.2, dmr_up_cutoff = 0.5, dmr_down_cutoff = 0.4, 
                      dmr_pval_cutoff = 1e-11, cpg_pval_cutoff = 1e-08, cpg_up_dm_cutoff = 0, cpg_down_dm_cutoff = 0, 
                      pairwise_comparison = FALSE, mset_train_flow_sort = NULL){
        # Based on methylCC:::.find_dmrs
        if (is.null(mset_train_flow_sort)) {
                FlowSorted.Blood.450k <- updateObject(FlowSorted.Blood.450k)
                mset_train_flow_sort <- preprocessIllumina(FlowSorted.Blood.450k)
                mset_train_flow_sort <- mapToGenome(mset_train_flow_sort, mergeManifest = FALSE)
                rm(FlowSorted.Blood.450k)
        }
        IDs = c("Gran", "CD4T", "CD8T", "Bcell", "Mono", "NK", "nRBC") # Added nRBCs, removed specific sample exclusion
        mset_train_flow_sort <- mset_train_flow_sort[, (pData(mset_train_flow_sort)$CellType %in% IDs)]
        if (!is.null(gr_target)) {
                zz <- findOverlaps(granges(mset_train_flow_sort), gr_target)
                mset_train_flow_sort <- mset_train_flow_sort[queryHits(zz), ]
                if (verbose) {
                        mes <- "[estimatecc] gr_target is not null. Using %s overlapping CpGs."
                        message(sprintf(mes, nrow(mset_train_flow_sort)))
                }
        }
        pd <- as.data.frame(pData(mset_train_flow_sort))
        gr <- granges(mset_train_flow_sort)
        p_beta <- getBeta(mset_train_flow_sort, type = "Illumina")
        colnames(p_beta) = pd$Sample_Name = rownames(pd) = gsub("\\+", "", pd$Sample_Name)
        cell <- factor(pd$CellType, levels = IDs)
        cell_levels <- levels(cell)
        chr <- as.character(seqnames(gr))
        pos <- start(gr)
        cl <- clusterMaker(chr, pos)
        xmat = cbind(rep(1, length(cell)), model.matrix(~cell - 1))
        colnames(xmat) = c("Intercept", cell_levels)
        if (pairwise_comparison) {
                all_poss = as.matrix(expand.grid(c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1)))
                all_poss = all_poss[2:32, ]
                all_poss <- (all_poss == TRUE)
                colnames(all_poss) <- cell_levels
        }
        else {
                all_poss = diag(length(cell_levels))
                all_poss <- (all_poss == TRUE)
                colnames(all_poss) <- cell_levels
        }
        regions_all <- GRanges()
        zmat <- c()
        for (ind in seq_len(nrow(all_poss))) {
                if (verbose) {
                        if (include_dmrs & include_cpgs) {
                                mes <- "[estimatecc] Searching for %s cell type-specific regions and CpGs."
                                message(sprintf(mes, paste(cell_levels[all_poss[ind, ]], collapse = ",")))
                        }
                        if (include_dmrs & !include_cpgs) {
                                mes <- "[estimatecc] Searching for %s cell type-specific regions."
                                message(sprintf(mes, paste(cell_levels[all_poss[ind, ]], collapse = ",")))
                        }
                        if (!include_dmrs & include_cpgs) {
                                mes <- "[estimatecc] Searching for %s cell type-specific CpGs"
                                message(sprintf(mes, paste(cell_levels[all_poss[ind, ]], collapse = ",")))
                        }
                }
                x_ind = cbind(Intercept = xmat[, "Intercept"], 
                              cellTypes = rowSums(as.matrix(xmat[, cell_levels[all_poss[ind, ]]], 
                                                            ncols = length(cell_levels[all_poss[ind, ]]))))
                if (!include_dmrs) {
                        gr_regions_up <- GRanges()
                        gr_regions_down <- GRanges()
                }
                if (include_dmrs) {
                        bumps = bumphunter(object = p_beta, design = x_ind, chr = chr, pos = pos, cluster = cl, 
                                           cutoff = bumphunter_beta_cutoff, B = 0, smooth = FALSE, 
                                           smoothFunction = loessByCluster)
                        y_regions <- t(apply(bumps$table[, 7:8], 1, function(z) {
                                colMeans(p_beta[(z[1]):(z[2]), , drop = FALSE])
                        }))
                        tmp <- rowttests(y_regions, factor(x_ind[, "cellTypes"]))
                        bumps$table$p.value <- tmp$p.value
                        bumps$table$dm <- tmp$dm
                        bumps$table$dmr_up_max_diff <- apply(abs(sweep(y_regions, 2, x_ind[, "cellTypes"], FUN = "-")), 
                                                             1, max)
                        bumps$table$dmr_down_max_diff <- apply(abs(sweep(y_regions, 2, (1 - x_ind[, "cellTypes"]), 
                                                                         FUN = "-")), 1, max)
                        L = dm <- NULL
                        keep_ind_regions <- (bumps$table$L > 1 | (bumps$table$L == 1 & bumps$table$clusterL == 1)) & 
                                (bumps$table$p.value < dmr_pval_cutoff)
                        bump_mat_up <- bumps$table[keep_ind_regions & bumps$table$dm < 0 & 
                                                           bumps$table$dmr_up_max_diff < dmr_up_cutoff, ]
                        bump_mat_up <- bump_mat_up[order(-bump_mat_up$L, bump_mat_up$dm), ]
                        if (nrow(bump_mat_up) > 0) {
                                gr_regions_up <- makeGRangesFromDataFrame(bump_mat_up, keep.extra.columns = TRUE)
                                mcols(gr_regions_up)$dmr_status <- rep("DMR", length(gr_regions_up))
                                gr_regions_up <- gr_regions_up[, names(mcols(gr_regions_up)) %in% 
                                                                       c("indexStart", "indexEnd", "L", "dm", "p.value", 
                                                                         "dmr_status", "dmr_up_max_diff")]
                                names(mcols(gr_regions_up))[names(mcols(gr_regions_up)) == "dmr_up_max_diff"] <- "dmr_max_diff"
                                gr_regions_up <- gr_regions_up %>% arrange(-L, dm) %>% head(num_regions)
                        }
                        else {
                                gr_regions_up <- GRanges()
                        }
                        bump_mat_down <- bumps$table[(keep_ind_regions) & bumps$table$dm > 0 & 
                                                             bumps$table$dmr_down_max_diff < dmr_down_cutoff, ]
                        bump_mat_down <- bump_mat_down[order(-bump_mat_down$L, -bump_mat_down$dm), ]
                        if (nrow(bump_mat_down) > 0) {
                                gr_regions_down <- makeGRangesFromDataFrame(bump_mat_down, keep.extra.columns = TRUE)
                                mcols(gr_regions_down)$dmr_status <- rep("DMR", length(gr_regions_down))
                                gr_regions_down <- gr_regions_down[, names(mcols(gr_regions_down)) %in% 
                                                                           c("indexStart", "indexEnd", "L", "dm", "p.value", 
                                                                             "dmr_status", "dmr_down_max_diff")]
                                names(mcols(gr_regions_down))[names(mcols(gr_regions_down)) == "dmr_down_max_diff"] <- "dmr_max_diff"
                                gr_regions_down <- gr_regions_down %>% arrange(-L, -dm) %>% head(num_regions)
                        }
                        else {
                                gr_regions_down <- GRanges()
                        }
                }
                if (include_cpgs) {
                        tstats <- rowttests(p_beta, factor(x_ind[, "cellTypes"]))
                        tstats <- tstats[(tstats[, "p.value"] < cpg_pval_cutoff), ]
                        tstats_up <- tstats[order(tstats[, "dm"], decreasing = FALSE), ]
                        tstats_up <- tstats_up[tstats_up$dm < cpg_up_dm_cutoff, ]
                        probe_keep <- rownames(tstats_up)[seq_len(min(nrow(tstats_up), num_cpgs))]
                        if (length(probe_keep) > 0) {
                                gr_probe <- granges(mset_train_flow_sort[probe_keep, ])
                                mcols(gr_probe) <- tstats[probe_keep, c("dm", "p.value")]
                                mcols(gr_probe)$L <- rep(1, length(probe_keep))
                                mcols(gr_probe)$indexStart <- match(probe_keep, rownames(mset_train_flow_sort))
                                mcols(gr_probe)$indexEnd <- match(probe_keep, rownames(mset_train_flow_sort))
                                mcols(gr_probe)$dmr_status <- rep("CpG", length(gr_probe))
                                gr_regions_up <- unique(c(gr_regions_up, 
                                                          gr_probe[, c("indexStart", "indexEnd", "L", "p.value", "dm", 
                                                                       "dmr_status")]))
                                gr_regions_up <- gr_regions_up %>% arrange(-L, dm) # Removed extra num_regions subset
                        }
                        tstats_down <- tstats[order(tstats[, "dm"], decreasing = TRUE), ]
                        tstats_down <- tstats_down[tstats_down$dm > cpg_down_dm_cutoff, ]
                        probe_keep <- rownames(tstats_down)[seq_len(min(nrow(tstats_down), num_cpgs))]
                        if (length(probe_keep) > 0) {
                                gr_probe <- granges(mset_train_flow_sort[probe_keep, ])
                                mcols(gr_probe) <- tstats[probe_keep, c("dm", "p.value")]
                                mcols(gr_probe)$L <- rep(1, length(probe_keep))
                                mcols(gr_probe)$indexStart <- match(probe_keep, rownames(mset_train_flow_sort))
                                mcols(gr_probe)$indexEnd <- match(probe_keep, rownames(mset_train_flow_sort))
                                mcols(gr_probe)$dmr_status <- rep("CpG", length(gr_probe))
                                gr_regions_down <- unique(c(gr_regions_down, 
                                                            gr_probe[, c("indexStart", "indexEnd", "L", "p.value", "dm", 
                                                                         "dmr_status")]))
                                gr_regions_down <- gr_regions_down %>% arrange(-L, -dm) # Removed extra num_regions subset
                        }
                }
                mcols(gr_regions_up)$status <- rep("Up", length(gr_regions_up))
                mcols(gr_regions_down)$status <- rep("Down", length(gr_regions_down))
                bump_mat_all <- c(gr_regions_up, gr_regions_down)
                mcols(bump_mat_all)$cellType <- rep(paste(cell_levels[all_poss[ind, ]], collapse = ","), 
                                                    length(bump_mat_all))
                if (verbose) {
                        if (include_dmrs & include_cpgs) {
                                mes <- "[estimatecc] Found %s %s cell type-specific regions and CpGs."
                                message(sprintf(mes, length(bump_mat_all), paste(cell_levels[all_poss[ind, ]], 
                                                                                 collapse = ",")))
                        }
                        if (include_dmrs & !include_cpgs) {
                                mes <- "[estimatecc] Found %s %s cell type-specific regions."
                                message(sprintf(mes, length(bump_mat_all), paste(cell_levels[all_poss[ind, ]], 
                                                                                 collapse = ",")))
                        }
                        if (!include_dmrs & include_cpgs) {
                                mes <- "[estimatecc] Found %s %s cell type-specific CpGs."
                                message(sprintf(mes, length(bump_mat_all), paste(cell_levels[all_poss[ind, ]], 
                                                                                 collapse = ",")))
                        }
                }
                if (length(bump_mat_all) > 0) {
                        regions_all <- c(regions_all, bump_mat_all)
                }
                if (length(gr_regions_up) > 0) {
                        zmat <- rbind(zmat, t(replicate(length(gr_regions_up), as.numeric(all_poss[ind, ])))) # Removed extra num_regions subset
                }
                if (length(gr_regions_down) > 0) {
                        zmat <- rbind(zmat, t(replicate(length(gr_regions_down), as.numeric(!all_poss[ind, ])))) # Removed extra num_regions subset
                }
        }
        colnames(zmat) <- cell_levels
        y_regions <- t(apply(as.data.frame(mcols(regions_all))[, seq_len(2)], 1, function(ind) {
                colMeans(p_beta[(ind[1]):(ind[2]), , drop = FALSE])
        }))
        profiles <- vapply(methylCC:::.splitit(cell), FUN = function(ind) {
                rowMeans(y_regions[, ind])
        }, FUN.VALUE = numeric(nrow(y_regions)))
        removeMe <- duplicated(regions_all)
        list(regions_all = regions_all[!removeMe, ], zmat = zmat[!removeMe, ], y_regions = y_regions[!removeMe, ], 
             profiles = profiles[!removeMe, ], cell = cell, cell_mat = all_poss, cell_levels = cell_levels, pd = pd)
}

# Identify Covered CpGs ---------------------------------------------------
# Load BSseq Objects ####
BSmalesDisc <- readRDS("Discovery/Dx_Males/Filtered_BSseq_Discovery50_males.rds")
BSmalesRep <- readRDS("Replication/Dx_Males/Filtered_BSseq_Replication50_males.rds")
BSfemalesDisc <- readRDS("Discovery/Dx_Females/Filtered_BSseq_Discovery50_females.rds")
BSfemalesRep <- readRDS("Replication/Dx_Females_100/Filtered_BSseq_Replication100_females.rds")

# Intersect CpGs ####
GRmalesDisc <- granges(BSmalesDisc)
GRmalesRep <- granges(BSmalesRep)
GRfemalesDisc <- granges(BSfemalesDisc)
GRfemalesRep <- granges(BSfemalesRep)
GRall <- intersect(GRmalesDisc, GRmalesRep) %>% intersect(GRfemalesDisc) %>% intersect(GRfemalesRep) # 10882594 CpGs
saveRDS(GRall, file = "Intersected_GRanges_Males_Females_Discovery_Replication.rds")

# Create Combined BSseq Object ####
BSall <- combine(BSmalesDisc, BSmalesRep, BSfemalesDisc, BSfemalesRep) %>% subsetByOverlaps(ranges = GRall)
saveRDS(BSall, file = "Intersected_BSseq_Males_Females_Discovery_Replication.rds")

# Define Cord Blood Cell-Type DMRs and Test on Sorted Reference Data ------
# Create MethylSet Object ####
FlowSorted.CordBlood.450k <- updateObject(FlowSorted.CordBlood.450k)
cordMethylSet <- preprocessIllumina(FlowSorted.CordBlood.450k) %>% mapToGenome(mergeManifest = FALSE) # hg19
seqlevelsStyle(cordMethylSet) <- "UCSC"
names(cordMethylSet@colData@listData)[names(cordMethylSet@colData@listData) == "X"] <- "Sample_Name"
pData <- pData(cordMethylSet) %>% as.data.frame()

# LiftOver MethylSet to hg38 ####
AnnotationHub <- AnnotationHub()
query(AnnotationHub, pattern = c("chain", "hg38", "hg19"))
# AH14108 | hg38ToHg19.over.chain.gz
# AH14150 | hg19ToHg38.over.chain.gz
liftOver <- liftOver(cordMethylSet, chain = AnnotationHub[["AH14150"]]) %>% unlist()
table(duplicated(names(liftOver))) # FALSE, all unique
cordMethylSet <- cordMethylSet[names(liftOver),]
table(names(liftOver) == names(cordMethylSet)) # All TRUE
table(names(cordMethylSet) == rownames(minfi::getMeth(cordMethylSet))) # All TRUE
rowRanges(cordMethylSet) <- liftOver # Replace hg19 coordinates with hg38
seqnames(cordMethylSet)
# factor-Rle of length 485344 with 33 runs
# Lengths: 28119     7 13977     8  4699 13685     6 21115 25158 ...  5922 12045     2 13466 10379  4243  8516 11229   416
# Values :  chr1 chr21  chr1  chr9  chr1  chr2  chr1  chr2  chr3 ... chr18 chr19  chr3 chr19 chr20 chr21 chr22  chrX  chrY
# Levels(24): chr1 chr21 chr9 chr2 chr3 chr4 chr5 chr6 chr7 chr8 ... chr13 chr14 chr16 chr17 chr18 chr19 chr20 chr22 chrX chrY
is.unsorted(cordMethylSet) # TRUE
cordMethylSet <- sortSeqlevels(cordMethylSet) %>% sort()
saveRDS(cordMethylSet, "R Objects/Cord Cell Data.rds")
rm(AnnotationHub, FlowSorted.CordBlood.450k, liftOver)

# Call DMRs ####
# Intersect array with CpGs used in all 4 DMR comparisons
GRall <- readRDS("R Objects/Intersected_GRanges_Males_Females_Discovery_Replication.rds")
cutoffs <- seq(0.5,0.8,0.05)
cordDMRs <- lapply(cutoffs, function(x){
        find_dmrs(mset_train_flow_sort = cordMethylSet, gr_target = GRall, num_regions = 100, dmr_pval_cutoff = 1e-8, 
                  dmr_up_cutoff = x, dmr_down_cutoff = x)
})
names(cordDMRs) <- paste("cutoff", cutoffs, sep = "_")
saveRDS(cordDMRs, "R Objects/Cord Cell DMRs.rds")

# Estimate Cell Counts ####
cordDMRs <- cordDMRs[!names(cordDMRs) == "cutoff_0.5"] # No regions in NK and CD8T cells
cordCounts <- lapply(cordDMRs, function(x){
        estimatecc(cordMethylSet, find_dmrs_object = x) %>% cell_counts() %>% cbind(pData, .)
})
saveRDS(cordCounts, "R Objects/Cord Cell Counts.rds")

# Proportion of Sorted Cell Type ####
cellTypes <- unique(pData$CellType) %>% as.character() %>% .[!. == "WholeBlood"]
proportion <- lapply(cordCounts, function(z){
        sapply(cellTypes, function(x) mean(z[z$CellType == x,x]))
}) %>% list.rbind() %>% as.data.frame()
proportion$Mean <- rowMeans(proportion)
proportion$Cutoff <- rownames(proportion)
proportion <- proportion[,c("Cutoff", colnames(proportion)[!colnames(proportion) == "Cutoff"])]
proportion$Mean # 0.5081934 0.5297888 0.5364002 0.5179606 0.4933898 0.4847378 Best at 0.65

# Error in Whole Blood Samples ####
expected <- c(0.08,0.15,0.1,0.45,0.1,0.02,0.1)
error <- lapply(cordCounts, function(z){
        mapply(function(x,y) mean(abs(z[z$CellType == "WholeBlood",x] - y)), x = cellTypes, y = expected)
}) %>% list.rbind() %>% as.data.frame()
error$Mean <- rowMeans(error)
error$Cutoff <- rownames(error)
error <- error[,c("Cutoff", colnames(error)[!colnames(error) == "Cutoff"])]
error$Mean # 0.04840655 0.04857092 0.04680500 0.05310193 0.06742275 0.07305515 Best at 0.65

# Define DMRs with num_regions = 50 and Test ------------------------------
# Call DMRs ####
cutoffs <- seq(0.55,0.8,0.05)
cordDMRs <- lapply(cutoffs, function(x){
        find_dmrs(mset_train_flow_sort = cordMethylSet, gr_target = GRall, num_regions = 50, dmr_pval_cutoff = 1e-8, 
                  dmr_up_cutoff = x, dmr_down_cutoff = x)
})
names(cordDMRs) <- paste("cutoff", cutoffs, sep = "_")
saveRDS(cordDMRs, "R Objects/Cord Cell DMRs num_regions 50.rds")

# Estimate Cell Counts ####
cordCounts <- lapply(cordDMRs, function(x){
        estimatecc(cordMethylSet, find_dmrs_object = x) %>% cell_counts() %>% cbind(pData, .)
})
saveRDS(cordCounts, "R Objects/Cord Cell Counts num_regions 50.rds")

# Proportion of Sorted Cell Type ####
proportion <- lapply(cordCounts, function(z){
        sapply(cellTypes, function(x) mean(z[z$CellType == x,x]))
}) %>% list.rbind() %>% as.data.frame()
proportion$Mean <- rowMeans(proportion)
proportion$Cutoff <- rownames(proportion)
proportion <- proportion[,c("Cutoff", colnames(proportion)[!colnames(proportion) == "Cutoff"])]
proportion$Mean   # 0.4969695 0.5553030 0.5258447 0.5052229 0.4845632 0.4750341 Best at 0.6

# Error in Whole Blood Samples ####
error <- lapply(cordCounts, function(z){
        mapply(function(x,y) mean(abs(z[z$CellType == "WholeBlood",x] - y)), x = cellTypes, y = expected)
}) %>% list.rbind() %>% as.data.frame()
error$Mean <- rowMeans(error)
error$Cutoff <- rownames(error)
error <- error[,c("Cutoff", colnames(error)[!colnames(error) == "Cutoff"])]
error$Mean        # 0.05267417 0.04691637 0.05375292 0.05895134 0.06834236 0.07013429 Best at 0.6

# Define CpGs and Test ---------------------------------------
# Call CpGs ####
cutoffs <- seq(0,0.5,0.1) # Default is 0
cordCpGs <- lapply(cutoffs, function(x){
        find_dmrs(mset_train_flow_sort = cordMethylSet, include_cpgs = TRUE, include_dmrs = FALSE, gr_target = GRall, 
                  num_cpgs = 100, cpg_pval_cutoff = 1e-8, cpg_up_dm_cutoff = -x, cpg_down_dm_cutoff = x)
})
names(cordCpGs) <- paste("cutoff", cutoffs, sep = "_")
saveRDS(cordCpGs, "R Objects/Cord Cell CpGs.rds")

# Estimate Cell Counts ####
cordCpGcounts <- lapply(cordCpGs, function(x){
        estimatecc(cordMethylSet, find_dmrs_object = x) %>% cell_counts() %>% cbind(pData, .)
})
saveRDS(cordCpGcounts, "R Objects/Cord Cell Counts CpGs.rds")

# Proportion of Sorted Cell Type ####
proportion <- lapply(cordCpGcounts, function(z){
        sapply(cellTypes, function(x) mean(z[z$CellType == x,x]))
}) %>% list.rbind() %>% as.data.frame()
proportion$Mean <- rowMeans(proportion)
proportion$Cutoff <- rownames(proportion)
proportion <- proportion[,c("Cutoff", colnames(proportion)[!colnames(proportion) == "Cutoff"])]
proportion$Mean # 0.4983005 0.4982870 0.4983988 0.5000257 0.5090718 0.5046370

# Error in Whole Blood Samples ####
error <- lapply(cordCpGcounts, function(z){
        mapply(function(x,y) mean(abs(z[z$CellType == "WholeBlood",x] - y)), x = cellTypes, y = expected)
}) %>% list.rbind() %>% as.data.frame()
error$Mean <- rowMeans(error)
error$Cutoff <- rownames(error)
error <- error[,c("Cutoff", colnames(error)[!colnames(error) == "Cutoff"])]
error$Mean # 0.08585276 0.08586336 0.08585141 0.08936568 0.08831503 0.07339794

# Define CpGs with num_cpgs = 50 and Test ------------------------------
# Call CpGs ####
cutoffs <- seq(0,0.5,0.1) # Default is 0
cordCpGs <- lapply(cutoffs, function(x){
        find_dmrs(mset_train_flow_sort = cordMethylSet, include_cpgs = TRUE, include_dmrs = FALSE, gr_target = GRall, 
                  num_cpgs = 50, cpg_pval_cutoff = 1e-8, cpg_up_dm_cutoff = -x, cpg_down_dm_cutoff = x)
})
names(cordCpGs) <- paste("cutoff", cutoffs, sep = "_")
saveRDS(cordCpGs, "R Objects/Cord Cell CpGs.rds")

# Estimate Cell Counts ####
cordCpGcounts <- lapply(cordCpGs, function(x){
        estimatecc(cordMethylSet, find_dmrs_object = x) %>% cell_counts() %>% cbind(pData, .)
})
saveRDS(cordCpGcounts, "R Objects/Cord Cell Counts CpGs.rds")

# Proportion of Sorted Cell Type ####
proportion <- lapply(cordCpGcounts, function(z){
        sapply(cellTypes, function(x) mean(z[z$CellType == x,x]))
}) %>% list.rbind() %>% as.data.frame()
proportion$Mean <- rowMeans(proportion)
proportion$Cutoff <- rownames(proportion)
proportion <- proportion[,c("Cutoff", colnames(proportion)[!colnames(proportion) == "Cutoff"])]
proportion$Mean # 0.5089637 0.5089487 0.5089725 0.5115660 0.5138068 0.5436219

# Error in Whole Blood Samples ####
error <- lapply(cordCpGcounts, function(z){
        mapply(function(x,y) mean(abs(z[z$CellType == "WholeBlood",x] - y)), x = cellTypes, y = expected)
}) %>% list.rbind() %>% as.data.frame()
error$Mean <- rowMeans(error)
error$Cutoff <- rownames(error)
error <- error[,c("Cutoff", colnames(error)[!colnames(error) == "Cutoff"])]
error$Mean # 0.08081603 0.08079256 0.08082589 0.08368527 0.08590094 0.06929703

# Best Settings: DMRs num_regions = 50, cutoff = 0.6 ####
cordDMRs <- cordDMRs$cutoff_0.6
saveRDS(cordDMRs, file = "R Objects/Cord Cell Type DMRs num_regions 50 cutoff 06.rds")
cordDMRs_df <- as.data.frame(cordDMRs$regions_all, row.names = 1:length(cordDMRs$regions_all))
write.table(cordDMRs_df, file = "Tables/Cord Cell Type DMRs num_regions 50 cutoff 06.txt", sep = "\t", quote = FALSE,
            row.names = FALSE)
cordCounts <- cordCounts$cutoff_0.6
write.table(cordCounts, file = "Tables/Cord Cell Type Reference Counts num_regions 50 cutoff 06.txt", sep = "\t", quote = FALSE,
            row.names = FALSE)

# Cell Type DMRs
with(cordDMRs_df, table(status, cellType))
# status Bcell CD4T CD8T Gran Mono NK nRBC
#   Down    50   50   13   50   50 50   50
#   Up      50   10    1   25   12 36   25
nrow(cordDMRs_df) # 472 DMRs
sum(cordDMRs_df$L) # 1046 CpGs

# Proportion Detail
sapply(cellTypes, function(x) mean(cordCounts[cordCounts$CellType == x,x]))
#     Bcell      CD4T      CD8T      Gran      Mono        NK      nRBC 
# 0.6758082 0.5480523 0.4174173 0.6531876 0.6992372 0.4905554 0.4028630

# Whole Blood Proportion
sapply(cellTypes, function(x) mean(cordCounts[cordCounts$CellType == "WholeBlood",x]))
#      Bcell       CD4T       CD8T       Gran       Mono         NK       nRBC 
# 0.07885333 0.13431421 0.08646268 0.36817930 0.17163929 0.07528425 0.08526694 

# Error Detail
mapply(function(x,y) mean(abs(cordCounts[cordCounts$CellType == "WholeBlood",x] - y)), x = cellTypes, y = expected)
#      Bcell       CD4T       CD8T       Gran       Mono         NK       nRBC 
# 0.02124538 0.02645216 0.03825444 0.08182070 0.07163929 0.05528425 0.03371837 

# Average Error Detail
sapply(cellTypes, function(x) mean(cordCounts[cordCounts$CellType == "WholeBlood",x])) - expected
# -0.001146671 -0.015685793 -0.013537321 -0.081820700  0.071639290  0.055284250 -0.014733055
rm(cordCpGcounts, cordCpGs, cordCpGsDefault, error, GRall, proportion, specificity, cutoffs, expected, find_dmrs)

# Compare with Previous ####
# Overlap CpGs
minfi_probes <- read.delim("Tables/Cord Cell Type Probes for Minfi.txt", stringsAsFactors = FALSE)
minfi_probes <- GRanges(seqnames = minfi_probes$Chromosome, ranges = IRanges(start = minfi_probes$Start, end = minfi_probes$End))
table(overlapsAny(query = minfi_probes, subject = cordDMRs$regions_all))
# FALSE  TRUE 
#   495    67 minfi CpG probes overlap methylCC DMRs

# Compare Cell Counts
cordCounts_m <- cordCounts[,c("Sample_Name", cellTypes)] %>% reshape2::melt(id.vars = "Sample_Name")
cordCounts_m$Sample_Name <- as.character(cordCounts_m$Sample_Name)
colnames(cordCounts_m)[colnames(cordCounts_m) == "variable"] <- "cellType"
colnames(cordCounts_m)[colnames(cordCounts_m) == "value"] <- "methylCC_count"
cordCounts_m$methylCC_count <- cordCounts_m$methylCC_count * 100

minfi_counts <- read.delim("Tables/Cord Cell Type Reference Counts minfi.txt", stringsAsFactors = FALSE)
minfi_counts <- minfi_counts[,c("X", cellTypes)] %>% reshape2::melt(id.vars = "X")
colnames(minfi_counts)[colnames(minfi_counts) == "X"] <- "Sample_Name"
colnames(minfi_counts)[colnames(minfi_counts) == "variable"] <- "cellType"
colnames(minfi_counts)[colnames(minfi_counts) == "value"] <- "minfi_count"
minfi_counts$minfi_count <- minfi_counts$minfi_count * 100
table(cordCounts_m$Sample_Name == minfi_counts$Sample_Name) # TRUE
table(cordCounts_m$cellType == minfi_counts$cellType) # TRUE
cordCounts_m$minfi_count <- minfi_counts$minfi_count

# Correlation Stats
corAllRef <- cor.test(x = cordCounts_m$minfi_count, y = cordCounts_m$methylCC_count) %>% .[c("estimate", "p.value")] %>%
        simplify2array()
# estimate.cor       p.value 
# 9.149201e-01 2.869317e-288 

corByCellTypeRef <- split(cordCounts_m, f = cordCounts_m$cellType) %>% 
        lapply(function(x){
                cor.test(x = x$minfi_count, y = x$methylCC_count) %>% .[c("estimate", "p.value")] %>% simplify2array()
        }) %>%
        list.rbind()
#       estimate.cor      p.value
# Bcell    0.9920507 1.205217e-93
# CD4T     0.9583246 2.598329e-57
# CD8T     0.6472908 1.130850e-13
# Gran     0.9481247 1.421349e-52
# Mono     0.9688664 1.174227e-63
# NK       0.9138145 1.044536e-41
# nRBC     0.9811353 1.274988e-74

# Scatterplot
gg <- ggplot(data = cordCounts_m, aes(x = minfi_count, y = methylCC_count)) +
        geom_smooth(method = "lm") +
        geom_point(aes(color = cellType), size = 2) +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'vertical', legend.position = c(1.13,0.86), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              plot.margin = unit(c(1,7,1,1), "lines")) +
        scale_x_continuous(breaks = breaks_pretty(n = 5)) +
        scale_y_continuous(breaks = breaks_pretty(n = 5)) +
        coord_cartesian(xlim = c(0,100), ylim = c(0,100)) +
        xlab("minfi Estimate") +
        ylab("methylCC Estimate") +
        scale_color_manual(values = wes_palette(n = 7, name = "FantasticFox1", type = "continuous"))
ggsave("Figures/minfi vs methylCC Reference Cell Composition.pdf", plot = gg, dpi = 600, width = 8, height = 7, units = "in")
rm(corByCellTypeRef, cordCounts_m, cordCounts, cordMethylSet, gg, minfi_counts, minfi_probes, refCellCounts, corAllRef, expected, pData)

# Estimate Cell Composition in Cord Blood WGBS and Compare with Previous ----------------------------
# Estimate from WGBS Data ####
sapply(c("tidyverse", "bsseq", "methylCC"), require, character.only = TRUE)
cordDMRs <- readRDS("Cord Cell Type DMRs num_regions 50 cutoff 06.rds")
BSall <- readRDS("Intersected_BSseq_Males_Females_Discovery_Replication.rds")
pData <- pData(BSall) %>% as.data.frame()
cordWGBScounts <- estimatecc(BSall, find_dmrs_object = cordDMRs) %>% cell_counts() %>% cbind(pData, .)
cordWGBScounts$Sequencing_ID <- rownames(cordWGBScounts)
write.table(cordWGBScounts, file = "Cord Cell Type WGBS Counts num_regions 50 cutoff 06.txt", sep = "\t", quote = FALSE,
            row.names = FALSE)

# Compare to Previous Method ####
# Data
cordWGBScounts <- read.csv("Tables/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo and Cell Type.csv")
cordWGBScounts <- cordWGBScounts[,c("Sequencing_ID", cellTypes)] %>% reshape2::melt(id.vars = "Sequencing_ID")
colnames(cordWGBScounts)[colnames(cordWGBScounts) == "variable"] <- "cellType"
colnames(cordWGBScounts)[colnames(cordWGBScounts) == "value"] <- "minfi_count"

methylCC_counts <- read.delim("Tables/Cord Cell Type WGBS Counts num_regions 50 cutoff 06.txt", stringsAsFactors = FALSE)
methylCC_counts <- methylCC_counts[,c("Sequencing_ID", cellTypes)] %>% reshape2::melt(id.vars = "Sequencing_ID")
colnames(methylCC_counts)[colnames(methylCC_counts) == "variable"] <- "cellType"
colnames(methylCC_counts)[colnames(methylCC_counts) == "value"] <- "methylCC_count"
methylCC_counts$methylCC_count <- methylCC_counts$methylCC_count * 100
cordWGBScounts <- merge(x = cordWGBScounts, y = methylCC_counts, by = c("Sequencing_ID", "cellType"), all = TRUE)

# Correlation Stats
corAll <- cor.test(x = cordWGBScounts$minfi_count, y = cordWGBScounts$methylCC_count) %>% .[c("estimate", "p.value")] %>%
        simplify2array()
# estimate.cor       p.value 
#    0.8462167 1.276502e-292 

corByCellType <- split(cordWGBScounts, f = cordWGBScounts$cellType) %>% 
        lapply(function(x){
                cor.test(x = x$minfi_count, y = x$methylCC_count) %>% .[c("estimate", "p.value")] %>% simplify2array()
        }) %>%
        list.rbind()
#       estimate.cor      p.value
# Bcell    0.5974167 4.491960e-16
# CD4T     0.5917150 9.955419e-16
# CD8T     0.1647048 4.258843e-02
# Gran     0.9328199 2.240642e-68
# Mono     0.3366080 2.233911e-05
# NK       0.2615187 1.136193e-03
# nRBC     0.5404471 6.643822e-13

# Scatterplot
gg <- ggplot(data = cordWGBScounts, aes(x = minfi_count, y = methylCC_count)) +
        geom_smooth(method = "lm") +
        geom_point(aes(color = cellType), size = 2) +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'vertical', legend.position = c(1.13,0.86), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              plot.margin = unit(c(1,7,1,1), "lines")) +
        scale_x_continuous(breaks = breaks_pretty(n = 5)) +
        scale_y_continuous(breaks = breaks_pretty(n = 5)) +
        coord_cartesian(xlim = c(0,100), ylim = c(0,100)) +
        xlab("minfi Estimate") +
        ylab("methylCC Estimate") +
        scale_color_manual(values = wes_palette(n = 7, name = "FantasticFox1", type = "continuous"))
ggsave("Figures/minfi vs methylCC WGBS Cell Composition.pdf", plot = gg, dpi = 600, width = 8, height = 7, units = "in")
