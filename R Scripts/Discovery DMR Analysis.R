# Discovery DMR Analysis --------------------------------------------------
# Diagnosis and Sex
# Autism Cord Blood Methylation
# Charles Mordaunt
# Excluded JLCM032B and JLCM050B
# 6/4/19

# Packages ####
sapply(c("tidyverse", "ggdendro", "scales", "ggplot2", "ggbiplot", "reshape", "grid", "RColorBrewer", "CMplot", "rlist",
         "annotatr", "GenomicRanges", "LOLA", "rtracklayer", "R.utils", "rGREAT", "limma", "DMRichR"), require, character.only = TRUE)

# Functions ####
source("R Scripts/DMR Analysis Functions.R")

# Diagnosis DMRs, no chrXY ------------------------------------------------
# Data ####
# DMRs
DMRs <- loadRegions("DMRs/Discovery/Diagnosis 50/DMR_smoothed_methylation_DxNoXY_Discovery50.txt",
                    chroms = c(paste("chr", 1:22, sep = ""), "chrM"), sort = TRUE)
DMRs$DMRid <- paste("DMR", 1:nrow(DMRs), sep = "_")
meth <- DMRs[,c("chr", "start", "end", "DMRid", colnames(DMRs)[grepl("JLCM", colnames(DMRs), fixed = TRUE)])]
DMRs <- DMRs[,c("chr", "start", "end", "DMRid", "width", "L", "percentDifference", "stat", "pval", "qval")]
raw <- loadRegions("DMRs/Discovery/Diagnosis 50/DMR_raw_methylation_DxNoXY_Discovery50.txt", 
                   chroms = c(paste("chr", 1:22, sep = ""), "chrM"), sort = TRUE)
raw$DMRid <- paste("DMR", 1:nrow(raw), sep = "_")
raw <- cbind(raw[,c("chr", "start", "end", "DMRid", colnames(raw)[grepl("JLCM", colnames(raw), fixed = TRUE)])])

# Candidates
candidates <- loadRegions("DMRs/Discovery/Diagnosis 50/CandidateRegions_DxNoXY_Discovery50.csv",
                          chroms = c(paste("chr", 1:22, sep = ""), "chrM"), sort = TRUE)
candidates <- candidates[,c("chr", "start", "end", "width", "L", "percentDifference", "stat", "pval", "qval")]

# Background
background <- loadRegions("DMRs/Discovery/Diagnosis 50/bsseq_background_Discovery50.csv",
                          chroms = c(paste("chr", 1:22, sep = ""), "chrM"), sort = TRUE)

# Samples
samples <- read.delim(file = "Merged Database/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo.txt", 
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE)
samples <- subset(samples, Sequencing_ID %in% colnames(meth))

# Manhattan and QQ plots ####
# Manhattan plot
CMplotDMR(candidates, prefix = "Figures/Diagnosis DMRs NoXY", plot.type = "m")

# QQ plot
CMplotDMR(candidates, prefix = "Figures/Diagnosis DMRs NoXY", plot.type = "q")
rm(candidates)

# Raw Meth Heatmap ####
# Meth Data
meth_heat <- raw[,grepl("JLCM", colnames(raw), fixed = TRUE)]
methdiff <- (meth_heat - rowMeans(meth_heat, na.rm = TRUE)) %>% as.matrix %>% "*" (100)
hm.lim <- quantile(abs(methdiff), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
table(is.na(methdiff))

# Pheno Data
phenoData <- samples[,c("Sequencing_ID", "Diagnosis_Alg", "Study", "Sex")]
colnames(phenoData) <- c("Sample", "Diagnosis", "Study", "Sex")
phenoData <- phenoData[match(colnames(methdiff), phenoData$Sample),]
table(colnames(methdiff) == phenoData$Sample) # All TRUE
phenoData$Sample <- factor(phenoData$Sample, levels = unique(phenoData$Sample), ordered = TRUE)
phenoData$Diagnosis <- factor(phenoData$Diagnosis, levels = c("TD", "ASD"), ordered = TRUE)
phenoData$Study <- factor(phenoData$Study, levels = c("MARBLES", "EARLI"), ordered = TRUE)
phenoData$Sex <- factor(phenoData$Sex, levels = c("M", "F"))

# Plot Heatmap
pdf(file="Figures/Diagnosis DMRs NoXY Raw Methylation Heatmap.pdf", width = 10, height = 8, onefile = FALSE)
buildHeatmap(x = methdiff, phenoData = phenoData, hm.low = - hm.lim, hm.high = hm.lim, 
             pheno.breaks = c("TD", "ASD", "MARBLES", "EARLI", "M", "F"), 
             pheno.values = c("TD" = "#3366CC", "ASD" = "#FF3366", "MARBLES" = "#FFFF33", "EARLI" = "#FF6633",
                              "M" = "#009933", "F" = "#9900CC")) %>%
        printHeatmap(widths = c(0.02, 0.76, 0.16, 0.06), heights = c(0.02, 0.15, 0.09, 0.72, 0.02),
                     heatmap.legend.position = c(0.69, -0.81), pheno.legend.position = c(0.925, 0.885))
dev.off()
rm(hm.lim, methdiff)

# Raw Meth PCA ####
meth_heat <- raw[,grepl("JLCM", colnames(raw), fixed = TRUE)] %>% as.matrix
table(colnames(meth_heat) == phenoData$Sample) # All TRUE
table(is.na(meth_heat))
for(i in 1:nrow(meth_heat)){
        temp <- as.numeric(meth_heat[i,])
        temp[is.na(temp)] <- mean(temp, na.rm = TRUE) # Replace missing values with mean meth of that DMR
        meth_heat[i,] <- temp
}
data.pca <- prcomp(t(meth_heat), center = TRUE, scale. = TRUE)
ggbiplotPCA(data.pca = data.pca, groups = phenoData$Diagnosis, pc = c(1,2), 
            file = "Figures/Diagnosis DMRs NoXY Raw Meth PCA by Diagnosis.png", xlim = c(-15, 15), ylim = c(-15, 15))
ggbiplotPCA(data.pca = data.pca, groups = phenoData$Study, pc = c(1,2), 
            file = "Figures/Diagnosis DMRs NoXY Raw Meth PCA by Study.png", xlim = c(-15, 15), ylim = c(-15, 15),
            breaks = c("MARBLES", "EARLI"), values = c("MARBLES" = "#3366CC", "EARLI" = "#FF3366"),
            legend.position = c(0.77, 1.03))
ggbiplotPCA(data.pca = data.pca, groups = phenoData$Sex, pc = c(1,2), 
            file = "Figures/Diagnosis DMRs NoXY Raw Meth PCA by Sex.png", xlim = c(-15, 15), ylim = c(-15, 15),
            breaks = c("M", "F"), values = c("M" = "#3366CC", "F" = "#FF3366"),
            legend.position = c(0.9, 1.03))
rm(data.pca, i, temp)

# Raw Meth Covariate Association ####
# Prep Data
meth_cov <- raw[,grepl("JLCM", colnames(raw), fixed = TRUE)] %>% t %>% as.data.frame
colnames(meth_cov) <- meth$DMRid
meth_cov$Sequencing_ID <- rownames(meth_cov)
samples_cov <- merge(x = samples, y = meth_cov, by = "Sequencing_ID", all = FALSE, sort = FALSE)
catVars <- c("Study", "Platform", "Sex", "Site", "Diagnosis_Alg", "MomEdu_detail", "DM1or2", "GDM", "PE", 
             "home_ownership", "marital_status", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3", "AllEQ_PV_YN_Mo_2", 
             "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", "AllEQ_PV_YN_Mo4", 
             "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", "AllEQ_PV_YN_Mo9")
catVars <- catVars[!catVars %in% c("Platform")] # Exclude catVars with only 1 level
samples_cov$Study <- factor(samples_cov$Study, levels = c("MARBLES", "EARLI"))
samples_cov$Platform <- factor(samples_cov$Platform, levels = c("HiSeqX10", "HiSeq4000"))
samples_cov$Sex <- factor(samples_cov$Sex, levels = c("M", "F"))
samples_cov$Site <- factor(samples_cov$Site, levels = c("UC Davis", "Kaiser Permanente", "Drexel", "Johns Hopkins University"))
samples_cov$Diagnosis_Alg <- factor(samples_cov$Diagnosis_Alg, levels = c("TD", "ASD"))
samples_cov$MomEdu_detail <- factor(samples_cov$MomEdu_detail, levels = c(6, 1:5, 7,8))
samples_cov$home_ownership[samples_cov$home_ownership == 99] <- NA
samples_cov$marital_status[samples_cov$marital_status == 99] <- NA
factorCols <- c("DM1or2", "GDM", "PE", "marital_status", "home_ownership", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3", 
                "AllEQ_PV_YN_Mo_2", "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", 
                "AllEQ_PV_YN_Mo4", "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", 
                "AllEQ_PV_YN_Mo9")
samples_cov[,factorCols] <- lapply(samples_cov[,factorCols], as.factor)
contVars <- colnames(samples_cov)[!colnames(samples_cov) %in% catVars & !colnames(samples_cov) %in% meth$DMRid &
                                          !colnames(samples_cov) %in% c("COI_ID", "Sequencing_ID", "Cord_Blood_IBC", 
                                                                        "Platform", "Sex")]
# Get Stats
covStats <- DMRmethLm(DMRs = meth$DMRid, catVars = catVars, contVars = contVars, sampleData = samples_cov, 
                      file = "Tables/Diagnosis DMRs NoXY Raw Methylation by Covariate Stats.txt")
# Missing Data
# Error with DMR_19 and DM1or2
# Error with DMR_25 and PE
# Error with DMR_90 and DM1or2
# Error with DMR_222 and PE
covSum <- DMRmethLmSum(covStats, file = "Tables/Diagnosis DMRs NoXY Raw Methylation by Covariate Summary.txt")

# Plot Heatmap
variables<-c("Diagnosis_Alg", "Sex", "Study", "Site_Drexel", "Site_Johns_Hopkins_University", "Site_Kaiser_Permanente", 
             "ADOScs", "MSLelcStandard36", "MSLelTscore36", "MSLfmTscore36", "MSLrlTscore36", "MSLvrTscore36", "ga_w", 
             "bw_g", "percent_trimmed", "percent_aligned", "percent_duplicate", "dedup_reads_M", "C_coverage", 
             "CG_coverage", "percent_cpg_meth", "percent_chg_meth", "percent_chh_meth", "MomEdu_detail_1", 
             "MomEdu_detail_2", "MomEdu_detail_3", "MomEdu_detail_4", "MomEdu_detail_5", "MomEdu_detail_7", 
             "MomEdu_detail_8", "home_ownership", "marital_status", "MomAgeYr", "Mat_Height_cm", "Mat_Weight_kg_PrePreg",
             "Mat_BMI_PrePreg", "DM1or2", "GDM", "PE", "parity", "dad_age", "SmokeYN_Pregnancy","cotinine_urine_ngml", 
             "final_creatinine_mgdl", "AllEQ_PV_YN_Mo_3", "AllEQ_PV_YN_Mo_2", "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", 
             "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", "AllEQ_PV_YN_Mo4", "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", 
             "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", "AllEQ_PV_YN_Mo9", "AllEQ_tot_All_FA_mcg_Mo_3", 
             "AllEQ_tot_All_FA_mcg_Mo_2", "AllEQ_tot_All_FA_mcg_Mo_1", "AllEQ_tot_All_FA_mcg_Mo1", 
             "AllEQ_tot_All_FA_mcg_Mo2", "AllEQ_tot_All_FA_mcg_Mo3", "AllEQ_tot_All_FA_mcg_Mo4", 
             "AllEQ_tot_All_FA_mcg_Mo5", "AllEQ_tot_All_FA_mcg_Mo6", "AllEQ_tot_All_FA_mcg_Mo7", 
             "AllEQ_tot_All_FA_mcg_Mo8","AllEQ_tot_All_FA_mcg_Mo9")
covHeatmap(covStats, variableOrdering = "manual", regionOrdering = "variable", variables = variables,
           sortVariable = "Diagnosis_Alg", 
           file = "Figures/Diagnosis DMRs NoXY Raw Meth Covariate Heatmap Sorted by Diagnosis.png")
covHeatmap(covStats, variableOrdering = "hierarchical", regionOrdering = "hierarchical", 
           file = "Figures/Diagnosis DMRs NoXY Raw Meth Covariate Heatmap Clustered.png")

# DMR Annotation ####
regDomains <- read.delim("Tables/Regulatory domains hg38.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
DMRs_anno <- getDMRanno(DMRstats = DMRs, regDomains = regDomains, file = "Tables/Diagnosis DMRs NoXY Annotation.txt")
DMRs_Dx_geneList <- list("All" = getDMRgeneList(DMRstats = DMRs, regDomains = regDomains, direction = "all", type = "gene_name"),
                         "Hyper" = getDMRgeneList(DMRstats = DMRs, regDomains = regDomains, direction = "hyper", type = "gene_name"),
                         "Hypo" = getDMRgeneList(DMRstats = DMRs, regDomains = regDomains, direction = "hypo", type = "gene_name"),
                         "Background" = getDMRgeneList(DMRstats = background, regDomains = regDomains, direction = "all", type = "gene_name"))
rm(regDomains, DMRs_anno)

# GREAT Analysis ####
# Make Files for GREAT (hg19, DMRs redefined to match background, background < 1M regions)
prepGREAT(DMRs = makeGRange(DMRs, direction = "all"), Background = makeGRange(background, direction = "all"), 
          fileName = "UCSC Tracks/Diagnosis NoXY DMRs Discovery hg19 for GREAT.bed", writeBack = TRUE,
          backName = "UCSC Tracks/Diagnosis NoXY Background Discovery hg19 for GREAT.bed", 
          chroms = c(paste("chr", 1:22, sep = ""), "chrM"))
prepGREAT(DMRs = makeGRange(DMRs, direction = "hyper"), Background = makeGRange(background, direction = "all"), 
          fileName = "UCSC Tracks/Diagnosis NoXY Hyper DMRs Discovery hg19 for GREAT.bed", writeBack = FALSE,
          chroms = c(paste("chr", 1:22, sep = ""), "chrM"))
prepGREAT(DMRs = makeGRange(DMRs, direction = "hypo"), Background = makeGRange(background, direction = "all"), 
          fileName = "UCSC Tracks/Diagnosis NoXY Hypo DMRs Discovery hg19 for GREAT.bed", writeBack = FALSE,
          chroms = c(paste("chr", 1:22, sep = ""), "chrM"))

# Get Enrichments from GREAT
BEDfile_Back <- "UCSC Tracks/Diagnosis NoXY Background Discovery hg19 for GREAT.bed.gz"
greatAll <- getGREATenrichments(BEDfile_DMR = "UCSC Tracks/Diagnosis NoXY DMRs Discovery hg19 for GREAT.bed.gz", 
                                BEDfile_Back = BEDfile_Back)
greatHyper <- getGREATenrichments(BEDfile_DMR = "UCSC Tracks/Diagnosis NoXY Hyper DMRs Discovery hg19 for GREAT.bed.gz", 
                                  BEDfile_Back = BEDfile_Back)
greatHypo <- getGREATenrichments(BEDfile_DMR = "UCSC Tracks/Diagnosis NoXY Hypo DMRs Discovery hg19 for GREAT.bed.gz", 
                                 BEDfile_Back = BEDfile_Back)
greatCombined <- rbind(greatAll, greatHyper, greatHypo)
greatCombined$Direction <- factor(c(rep("All", nrow(greatAll)), rep("Hyper", nrow(greatHyper)), 
                                    rep("Hypo", nrow(greatHypo))), levels = c("All", "Hyper", "Hypo"))
#write.table(greatCombined, "Tables/Diagnosis NoXY DMRs Discovery GREAT Combined Results.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
# No enriched terms

rm(greatAll, greatCombined, greatHyper, greatHypo, BEDfile_Back)

# Diagnosis and Sex DMRs --------------------------------------------------
# Includes JLCM032B and JLCM050B
# Data ####
# DMRs
DMRs <- loadRegions("DMRs/Discovery/Diagnosis and Sex 50/DMR_smoothed_methylation_DxAdjSex_Discovery50.txt",
                    chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)
DMRs$DMRid <- paste("DMR", 1:nrow(DMRs), sep = "_")
meth <- DMRs[,c("chr", "start", "end", "DMRid", colnames(DMRs)[grepl("JLCM", colnames(DMRs), fixed = TRUE)])]
DMRs <- DMRs[,c("chr", "start", "end", "DMRid", "width", "L", "percentDifference", "stat", "pval", "qval")]
raw <- loadRegions("DMRs/Discovery/Diagnosis and Sex 50/DMR_raw_methylation_DxAdjSex_Discovery50.txt", 
                   chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)
raw$DMRid <- paste("DMR", 1:nrow(raw), sep = "_")
raw <- cbind(raw[,c("chr", "start", "end", "DMRid", colnames(raw)[grepl("JLCM", colnames(raw), fixed = TRUE)])])

# Candidates
candidates <- loadRegions("DMRs/Discovery/Diagnosis and Sex 50/CandidateRegions_DxAdjSex_Discovery50.csv",
                          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)
candidates <- candidates[,c("chr", "start", "end", "width", "L", "percentDifference", "stat", "pval", "qval")]

# Background
background <- loadRegions("DMRs/Discovery/Diagnosis and Sex 50/bsseq_background_Discovery50.csv",
                          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)

# Samples
samples <- read.delim(file = "Merged Database/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo.txt", 
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE)
samples <- subset(samples, Sequencing_ID %in% colnames(meth))

# Manhattan and QQ plots ####
# Manhattan plot
CMplotDMR(candidates, prefix = "Figures/Diagnosis and Sex DMRs", plot.type = "m")

# QQ plot
CMplotDMR(candidates, prefix = "Figures/Diagnosis and Sex DMRs", plot.type = "q")
rm(candidates)

# Heatmap ####
# Meth Data
meth_heat <- meth[,grepl("JLCM", colnames(meth), fixed = TRUE)]

# Pheno Data
phenoData <- samples[,c("Sequencing_ID", "Diagnosis_Alg", "Study", "Sex")]
colnames(phenoData) <- c("Sample", "Diagnosis", "Study", "Sex")
phenoData <- phenoData[match(colnames(meth_heat), phenoData$Sample),]
table(colnames(meth_heat) == phenoData$Sample) # All TRUE
phenoData$Sample <- factor(phenoData$Sample, levels = unique(phenoData$Sample))
phenoData$Diagnosis <- factor(phenoData$Diagnosis, levels = c("TD", "ASD"))
phenoData$Study <- factor(phenoData$Study, levels = c("MARBLES", "EARLI"))
phenoData$Sex <- factor(phenoData$Sex, levels = c("M", "F"))

# Plot Heatmap
methdiff <- (meth_heat - rowMeans(meth_heat, na.rm = TRUE)) %>% as.matrix %>% "*" (100)
hm.lim <- quantile(abs(methdiff), probs = 0.99, names = FALSE) %>% ceiling # 25

pdf(file="Figures/Diagnosis and Sex DMRs Methylation Heatmap.pdf", width = 10, height = 8, onefile = FALSE)
buildHeatmap(x = methdiff, phenoData = phenoData[c("Sample", "Diagnosis", "Study")], hm.low = - hm.lim, hm.high = hm.lim, 
             pheno.breaks = c("TD", "ASD", "MARBLES", "EARLI", "M", "F"), 
             pheno.values = c("TD" = "#3366CC", "ASD" = "#FF3366", "MARBLES" = "#FFFF33", "EARLI" = "#FF6633")) %>%
        printHeatmap
dev.off()
rm(hm.lim, methdiff)

# Plot Heatmap with Meth Data Adjusted for Sex
mod <- model.matrix(~ Diagnosis, data = phenoData)
meth_adj <- removeBatchEffect(meth_heat, batch = phenoData$Sex, design = mod)
methdiff <- (meth_adj - rowMeans(meth_adj, na.rm = TRUE)) %>% as.matrix %>% "*" (100)
hm.lim <- quantile(abs(methdiff), probs = 0.99, names = FALSE) %>% ceiling # 25

pdf(file="Figures/Diagnosis and Sex DMRs Methylation AdjSex Heatmap.pdf", width = 10, height = 8, onefile = FALSE)
buildHeatmap(x = methdiff, phenoData = phenoData[c("Sample", "Diagnosis", "Study")], hm.low = - hm.lim, hm.high = hm.lim, 
             pheno.breaks = c("TD", "ASD", "MARBLES", "EARLI", "M", "F"), 
             pheno.values = c("TD" = "#3366CC", "ASD" = "#FF3366", "MARBLES" = "#FFFF33", "EARLI" = "#FF6633")) %>%
        printHeatmap
dev.off()
rm(hm.lim, methdiff)

# PCA ####
# Smoothed Methylation
table(colnames(meth_heat) == phenoData$Sample) # All TRUE
data.pca <- prcomp(meth_heat %>% as.matrix %>% t, center = TRUE, scale. = TRUE)
ggbiplotPCA(data.pca = data.pca, groups = phenoData$Diagnosis, pc = c(1,2), 
            file = "Figures/Diagnosis and Sex DMRs PCA by Diagnosis.png", xlim = c(-10, 10), ylim = c(-10, 10))
ggbiplotPCA(data.pca = data.pca, groups = phenoData$Study, pc = c(1,2), 
            file = "Figures/Diagnosis and Sex DMRs PCA by Study.png", xlim = c(-10, 10), ylim = c(-10, 10),
            breaks = c("MARBLES", "EARLI"), values = c("MARBLES" = "#3366CC", "EARLI" = "#FF3366"),
            legend.position = c(0.77, 1.03))
ggbiplotPCA(data.pca = data.pca, groups = phenoData$Sex, pc = c(1,2), 
            file = "Figures/Diagnosis and Sex DMRs PCA by Sex.png", xlim = c(-10, 10), ylim = c(-10, 10),
            breaks = c("M", "F"), values = c("M" = "#3366CC", "F" = "#FF3366"),
            legend.position = c(0.77, 1.03))
rm(data.pca)

# Smoothed Methylation Adjusted for Sex
table(colnames(meth_adj) == phenoData$Sample) # All TRUE
data.pca <- prcomp(meth_adj %>% t, center = TRUE, scale. = TRUE)
ggbiplotPCA(data.pca = data.pca, groups = phenoData$Diagnosis, pc = c(1,2), 
            file = "Figures/Diagnosis and Sex DMRs Meth AdjSex PCA by Diagnosis.png", xlim = c(-10, 10), ylim = c(-10, 10))
ggbiplotPCA(data.pca = data.pca, groups = phenoData$Study, pc = c(1,2), 
            file = "Figures/Diagnosis and Sex DMRs Meth AdjSex PCA by Study.png", xlim = c(-10, 10), ylim = c(-10, 10),
            breaks = c("MARBLES", "EARLI"), values = c("MARBLES" = "#3366CC", "EARLI" = "#FF3366"),
            legend.position = c(0.77, 1.03))
ggbiplotPCA(data.pca = data.pca, groups = phenoData$Sex, pc = c(1,2), 
            file = "Figures/Diagnosis and Sex DMRs Meth AdjSex PCA by Sex.png", xlim = c(-10, 10), ylim = c(-10, 10),
            breaks = c("M", "F"), values = c("M" = "#3366CC", "F" = "#FF3366"),
            legend.position = c(0.77, 1.03))
rm(data.pca)

# Covariate Association ####
# Prep Data
meth_cov <- meth_heat %>% t %>% as.data.frame
colnames(meth_cov) <- meth$DMRid
meth_cov$Sequencing_ID <- rownames(meth_cov)
samples_cov <- merge(x = samples, y = meth_cov, by = "Sequencing_ID", all = FALSE, sort = FALSE)
catVars <- c("Study", "Platform", "Sex", "Site", "Diagnosis_Alg", "MomEdu_detail", "DM1or2", "GDM", "PE", 
             "home_ownership", "marital_status", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3", "AllEQ_PV_YN_Mo_2", 
             "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", "AllEQ_PV_YN_Mo4", 
             "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", "AllEQ_PV_YN_Mo9")
catVars <- catVars[!catVars %in% c("Platform")] # Exclude catVars with only 1 level
samples_cov$Study <- factor(samples_cov$Study, levels = c("MARBLES", "EARLI"))
samples_cov$Platform <- factor(samples_cov$Platform, levels = c("HiSeqX10", "HiSeq4000"))
samples_cov$Sex <- factor(samples_cov$Sex, levels = c("M", "F"))
samples_cov$Site <- factor(samples_cov$Site, levels = c("UC Davis", "Kaiser Permanente", "Drexel", "Johns Hopkins University"))
samples_cov$Diagnosis_Alg <- factor(samples_cov$Diagnosis_Alg, levels = c("TD", "ASD"))
samples_cov$MomEdu_detail <- factor(samples_cov$MomEdu_detail, levels = c(6, 1:5, 7,8))
samples_cov$home_ownership[samples_cov$home_ownership == 99] <- NA
samples_cov$marital_status[samples_cov$marital_status == 99] <- NA
factorCols <- c("DM1or2", "GDM", "PE", "marital_status", "home_ownership", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3", 
                "AllEQ_PV_YN_Mo_2", "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", 
                "AllEQ_PV_YN_Mo4", "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", 
                "AllEQ_PV_YN_Mo9")
samples_cov[,factorCols] <- lapply(samples_cov[,factorCols], as.factor)
contVars <- colnames(samples_cov)[!colnames(samples_cov) %in% catVars & !colnames(samples_cov) %in% meth$DMRid &
                                          !colnames(samples_cov) %in% c("COI_ID", "Sequencing_ID", "Cord_Blood_IBC", 
                                                                        "Platform", "Sex")]
variables<-c("Diagnosis_Alg", "Sex", "Study", "Site_Drexel", "Site_Johns_Hopkins_University", "Site_Kaiser_Permanente", 
             "ADOScs", "MSLelcStandard36", "MSLelTscore36", "MSLfmTscore36", "MSLrlTscore36", "MSLvrTscore36", "ga_w", 
             "bw_g", "percent_trimmed", "percent_aligned", "percent_duplicate", "dedup_reads_M", "C_coverage", 
             "CG_coverage", "percent_cpg_meth", "percent_chg_meth", "percent_chh_meth", "MomEdu_detail_1", 
             "MomEdu_detail_2", "MomEdu_detail_3", "MomEdu_detail_4", "MomEdu_detail_5", "MomEdu_detail_7", 
             "MomEdu_detail_8", "home_ownership", "marital_status", "MomAgeYr", "Mat_Height_cm", "Mat_Weight_kg_PrePreg",
             "Mat_BMI_PrePreg", "DM1or2", "GDM", "PE", "parity", "dad_age", "SmokeYN_Pregnancy","cotinine_urine_ngml", 
             "final_creatinine_mgdl", "AllEQ_PV_YN_Mo_3", "AllEQ_PV_YN_Mo_2", "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", 
             "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", "AllEQ_PV_YN_Mo4", "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", 
             "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", "AllEQ_PV_YN_Mo9", "AllEQ_tot_All_FA_mcg_Mo_3", 
             "AllEQ_tot_All_FA_mcg_Mo_2", "AllEQ_tot_All_FA_mcg_Mo_1", "AllEQ_tot_All_FA_mcg_Mo1", 
             "AllEQ_tot_All_FA_mcg_Mo2", "AllEQ_tot_All_FA_mcg_Mo3", "AllEQ_tot_All_FA_mcg_Mo4", 
             "AllEQ_tot_All_FA_mcg_Mo5", "AllEQ_tot_All_FA_mcg_Mo6", "AllEQ_tot_All_FA_mcg_Mo7", 
             "AllEQ_tot_All_FA_mcg_Mo8","AllEQ_tot_All_FA_mcg_Mo9")

# Get Stats
covStats <- DMRmethLm(DMRs = meth$DMRid, catVars = catVars, contVars = contVars, sampleData = samples_cov, 
                      file = "Tables/Diagnosis and Sex DMRs Methylation by Covariate Stats.txt")
covSum <- DMRmethLmSum(covStats, file = "Tables/Diagnosis and Sex DMRs Methylation by Covariate Summary.txt")

# Plot Heatmap
covHeatmap(covStats, variableOrdering = "manual", regionOrdering = "variable", variables = variables,
           sortVariable = "Diagnosis_Alg", 
           file = "Figures/Diagnosis and Sex DMRs Covariate Heatmap Sorted by Diagnosis.png")
covHeatmap(covStats, variableOrdering = "hierarchical", regionOrdering = "hierarchical", 
           file = "Figures/Diagnosis and Sex DMRs Covariate Heatmap Clustered.png")

# Get Stats Adjusted for Sex
catVars <- catVars[!catVars == "Sex"]
covStats_adj <- DMRmethLm(DMRs = meth$DMRid, catVars = catVars, contVars = contVars, sampleData = samples_cov, 
                      file = "Tables/Diagnosis and Sex DMRs Methylation by Covariate with AdjSex Stats.txt", adj = "Sex")
covSum_adj <- DMRmethLmSum(covStats_adj, file = "Tables/Diagnosis and Sex DMRs Methylation by Covariate with AdjSex Summary.txt")

# Plot Heatmap Adjusted for Sex
variables <- variables[!variables == "Sex"]
covHeatmap(covStats_adj, variableOrdering = "manual", regionOrdering = "variable", variables = variables,
           sortVariable = "Diagnosis_Alg", 
           file = "Figures/Diagnosis and Sex DMRs Covariate with AdjSex Heatmap Sorted by Diagnosis.png")
covHeatmap(covStats_adj, variableOrdering = "hierarchical", regionOrdering = "hierarchical", 
           file = "Figures/Diagnosis and Sex DMRs Covariate with AdjSex Heatmap Clustered.png")

# DMR Annotation ####
regDomains <- read.delim("Tables/Regulatory domains hg38.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
DMRs_anno <- getDMRanno(DMRstats = DMRs, regDomains = regDomains, file = "Tables/Diagnosis and Sex DMRs Annotation.txt")
DMRs_DxAdjSex_geneList <- list("All" = getDMRgeneList(DMRstats = DMRs, regDomains = regDomains, direction = "all", type = "gene_name"),
                      "Hyper" = getDMRgeneList(DMRstats = DMRs, regDomains = regDomains, direction = "hyper", type = "gene_name"),
                      "Hypo" = getDMRgeneList(DMRstats = DMRs, regDomains = regDomains, direction = "hypo", type = "gene_name"),
                      "Background" = getDMRgeneList(DMRstats = background, regDomains = regDomains, direction = "all", type = "gene_name"))
rm(regDomains, DMRs_anno)

# GREAT Analysis ####
# Make Files for GREAT (hg19, DMRs redefined to match background, background < 1M regions)
prepGREAT(DMRs = makeGRange(DMRs, direction = "all"), Background = makeGRange(background, direction = "all"), 
          fileName = "UCSC Tracks/Diagnosis and Sex DMRs Discovery hg19 for GREAT.bed", writeBack = TRUE,
          backName = "UCSC Tracks/Diagnosis and Sex Background Discovery hg19 for GREAT.bed", 
          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"))
prepGREAT(DMRs = makeGRange(DMRs, direction = "hyper"), Background = makeGRange(background, direction = "all"), 
          fileName = "UCSC Tracks/Diagnosis and Sex Hyper DMRs Discovery hg19 for GREAT.bed", writeBack = FALSE,
          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"))
prepGREAT(DMRs = makeGRange(DMRs, direction = "hypo"), Background = makeGRange(background, direction = "all"), 
          fileName = "UCSC Tracks/Diagnosis and Sex Hypo DMRs Discovery hg19 for GREAT.bed", writeBack = FALSE,
          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"))

# Get Enrichments from GREAT
BEDfile_Back <- "UCSC Tracks/Diagnosis and Sex Background Discovery hg19 for GREAT.bed.gz"
greatAll <- getGREATenrichments(BEDfile_DMR = "UCSC Tracks/Diagnosis and Sex DMRs Discovery hg19 for GREAT.bed.gz", 
                                BEDfile_Back = BEDfile_Back)
greatHyper <- getGREATenrichments(BEDfile_DMR = "UCSC Tracks/Diagnosis and Sex Hyper DMRs Discovery hg19 for GREAT.bed.gz", 
                                  BEDfile_Back = BEDfile_Back)
greatHypo <- getGREATenrichments(BEDfile_DMR = "UCSC Tracks/Diagnosis and Sex Hypo DMRs Discovery hg19 for GREAT.bed.gz", 
                                 BEDfile_Back = BEDfile_Back)
greatCombined <- rbind(greatAll, greatHyper, greatHypo)
greatCombined$Direction <- factor(c(rep("All", nrow(greatAll)), rep("Hyper", nrow(greatHyper)), 
                                    rep("Hypo", nrow(greatHypo))), levels = c("All", "Hyper", "Hypo"))
write.table(greatCombined, "Tables/Diagnosis and Sex DMRs Discovery GREAT Combined Results.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
# No enriched terms

rm(greatAll, greatCombined, greatHyper, greatHypo, background, DMRs, BEDfile_Back, samples)


# Raw Meth Heatmap ####
# Meth Data
meth_heat <- raw[,grepl("JLCM", colnames(raw), fixed = TRUE)]
methdiff <- (meth_heat - rowMeans(meth_heat, na.rm = TRUE)) %>% as.matrix %>% "*" (100)
hm.lim <- quantile(abs(methdiff), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
table(is.na(methdiff))

# Pheno Data
phenoData <- samples[,c("Sequencing_ID", "Diagnosis_Alg", "Study", "Sex")]
colnames(phenoData) <- c("Sample", "Diagnosis", "Study", "Sex")
phenoData <- phenoData[match(colnames(methdiff), phenoData$Sample),]
table(colnames(methdiff) == phenoData$Sample) # All TRUE
phenoData$Sample <- factor(phenoData$Sample, levels = unique(phenoData$Sample), ordered = TRUE)
phenoData$Diagnosis <- factor(phenoData$Diagnosis, levels = c("TD", "ASD"), ordered = TRUE)
phenoData$Study <- factor(phenoData$Study, levels = c("MARBLES", "EARLI"), ordered = TRUE)
phenoData$Sex <- factor(phenoData$Sex, levels = c("M", "F"))

# Plot Heatmap
pdf(file="Figures/Diagnosis and Sex DMRs Raw Methylation Heatmap.pdf", width = 10, height = 8, onefile = FALSE)
buildHeatmap(x = methdiff, phenoData = phenoData[,c("Sample", "Diagnosis", "Study")], hm.low = - hm.lim, hm.high = hm.lim, 
             pheno.breaks = c("TD", "ASD", "MARBLES", "EARLI"), 
             pheno.values = c("TD" = "#3366CC", "ASD" = "#FF3366", "MARBLES" = "#FFFF33", "EARLI" = "#FF6633")) %>%
        printHeatmap
dev.off()
rm(hm.lim, methdiff)

# Plot Heatmap with Meth Data Adjusted for Sex
mod <- model.matrix(~ Diagnosis, data = phenoData)
table(is.na(meth_heat))
meth_adj <- removeBatchEffect(meth_heat, batch = phenoData$Sex, design = mod)
methdiff <- (meth_adj - rowMeans(meth_adj, na.rm = TRUE)) %>% as.matrix %>% "*" (100)
hm.lim <- quantile(abs(methdiff), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling

pdf(file="Figures/Diagnosis and Sex DMRs Raw Methylation AdjSex Heatmap.pdf", width = 10, height = 8, onefile = FALSE)
buildHeatmap(x = methdiff, phenoData = phenoData[c("Sample", "Diagnosis", "Study")], hm.low = - hm.lim, hm.high = hm.lim, 
             pheno.breaks = c("TD", "ASD", "MARBLES", "EARLI", "M", "F"), 
             pheno.values = c("TD" = "#3366CC", "ASD" = "#FF3366", "MARBLES" = "#FFFF33", "EARLI" = "#FF6633")) %>%
        printHeatmap
dev.off()

# Raw Meth PCA ####
meth_heat <- raw[,grepl("JLCM", colnames(raw), fixed = TRUE)] %>% as.matrix
table(colnames(meth_heat) == phenoData$Sample) # All TRUE
table(is.na(meth_heat))
for(i in 1:nrow(meth_heat)){
        temp <- as.numeric(meth_heat[i,])
        temp[is.na(temp)] <- mean(temp, na.rm = TRUE) # Replace missing values with mean meth of that DMR
        meth_heat[i,] <- temp
}
table(is.na(meth_heat))
data.pca <- prcomp(t(meth_heat), center = TRUE, scale. = TRUE)
ggbiplotPCA(data.pca = data.pca, groups = phenoData$Diagnosis, pc = c(1,2), 
            file = "Figures/Diagnosis and Sex DMRs Raw Meth PCA by Diagnosis.png", xlim = c(-10, 10), ylim = c(-10, 10))
ggbiplotPCA(data.pca = data.pca, groups = phenoData$Study, pc = c(1,2), 
            file = "Figures/Diagnosis and Sex DMRs Raw Meth PCA by Study.png", xlim = c(-10, 10), ylim = c(-10, 10),
            breaks = c("MARBLES", "EARLI"), values = c("MARBLES" = "#3366CC", "EARLI" = "#FF3366"),
            legend.position = c(0.77, 1.03))
ggbiplotPCA(data.pca = data.pca, groups = phenoData$Sex, pc = c(1,2), 
            file = "Figures/Diagnosis and Sex DMRs Raw Meth PCA by Sex.png", xlim = c(-10, 10), ylim = c(-10, 10),
            breaks = c("M", "F"), values = c("M" = "#3366CC", "F" = "#FF3366"),
            legend.position = c(0.77, 1.03))
rm(data.pca)

# Raw Methylation Adjusted for Sex
table(colnames(meth_adj) == phenoData$Sample) # All TRUE
table(is.na(meth_adj))
for(i in 1:nrow(meth_adj)){
        temp <- as.numeric(meth_adj[i,])
        temp[is.na(temp)] <- mean(temp, na.rm = TRUE) # Replace missing values with mean meth of that DMR
        meth_adj[i,] <- temp
}
table(is.na(meth_adj))
data.pca <- prcomp(meth_adj %>% t, center = TRUE, scale. = TRUE)
ggbiplotPCA(data.pca = data.pca, groups = phenoData$Diagnosis, pc = c(1,2), 
            file = "Figures/Diagnosis and Sex DMRs Raw Meth AdjSex PCA by Diagnosis.png", xlim = c(-10, 10), ylim = c(-10, 10))
ggbiplotPCA(data.pca = data.pca, groups = phenoData$Study, pc = c(1,2), 
            file = "Figures/Diagnosis and Sex DMRs Raw Meth AdjSex PCA by Study.png", xlim = c(-10, 10), ylim = c(-10, 10),
            breaks = c("MARBLES", "EARLI"), values = c("MARBLES" = "#3366CC", "EARLI" = "#FF3366"),
            legend.position = c(0.77, 1.03))
ggbiplotPCA(data.pca = data.pca, groups = phenoData$Sex, pc = c(1,2), 
            file = "Figures/Diagnosis and Sex DMRs Raw Meth AdjSex PCA by Sex.png", xlim = c(-10, 10), ylim = c(-10, 10),
            breaks = c("M", "F"), values = c("M" = "#3366CC", "F" = "#FF3366"),
            legend.position = c(0.77, 1.03))
rm(data.pca)

# Raw Meth Covariate Association ####
# Prep Data
meth_cov <- raw[,grepl("JLCM", colnames(raw), fixed = TRUE)] %>% t %>% as.data.frame
colnames(meth_cov) <- meth$DMRid
meth_cov$Sequencing_ID <- rownames(meth_cov)
samples_cov <- merge(x = samples, y = meth_cov, by = "Sequencing_ID", all = FALSE, sort = FALSE)
catVars <- c("Study", "Platform", "Sex", "Site", "Diagnosis_Alg", "MomEdu_detail", "DM1or2", "GDM", "PE", 
             "home_ownership", "marital_status", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3", "AllEQ_PV_YN_Mo_2", 
             "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", "AllEQ_PV_YN_Mo4", 
             "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", "AllEQ_PV_YN_Mo9")
catVars <- catVars[!catVars %in% c("Platform")] # Exclude catVars with only 1 level
samples_cov$Study <- factor(samples_cov$Study, levels = c("MARBLES", "EARLI"))
samples_cov$Platform <- factor(samples_cov$Platform, levels = c("HiSeqX10", "HiSeq4000"))
samples_cov$Sex <- factor(samples_cov$Sex, levels = c("M", "F"))
samples_cov$Site <- factor(samples_cov$Site, levels = c("UC Davis", "Kaiser Permanente", "Drexel", "Johns Hopkins University"))
samples_cov$Diagnosis_Alg <- factor(samples_cov$Diagnosis_Alg, levels = c("TD", "ASD"))
samples_cov$MomEdu_detail <- factor(samples_cov$MomEdu_detail, levels = c(6, 1:5, 7,8))
samples_cov$home_ownership[samples_cov$home_ownership == 99] <- NA
samples_cov$marital_status[samples_cov$marital_status == 99] <- NA
factorCols <- c("DM1or2", "GDM", "PE", "marital_status", "home_ownership", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3", 
                "AllEQ_PV_YN_Mo_2", "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", 
                "AllEQ_PV_YN_Mo4", "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", 
                "AllEQ_PV_YN_Mo9")
samples_cov[,factorCols] <- lapply(samples_cov[,factorCols], as.factor)
contVars <- colnames(samples_cov)[!colnames(samples_cov) %in% catVars & !colnames(samples_cov) %in% meth$DMRid &
                                          !colnames(samples_cov) %in% c("COI_ID", "Sequencing_ID", "Cord_Blood_IBC", 
                                                                        "Platform", "Sex")]
# Get Stats
covStats <- DMRmethLm(DMRs = meth$DMRid, catVars = catVars, contVars = contVars, sampleData = samples_cov, 
                      file = "Tables/Diagnosis and Sex DMRs Raw Methylation by Covariate Stats.txt")

covSum <- DMRmethLmSum(covStats, file = "Tables/Diagnosis and Sex DMRs Raw Methylation by Covariate Summary.txt")

# Plot Heatmap
variables<-c("Diagnosis_Alg", "Sex", "Study", "Site_Drexel", "Site_Johns_Hopkins_University", "Site_Kaiser_Permanente", 
             "ADOScs", "MSLelcStandard36", "MSLelTscore36", "MSLfmTscore36", "MSLrlTscore36", "MSLvrTscore36", "ga_w", 
             "bw_g", "percent_trimmed", "percent_aligned", "percent_duplicate", "dedup_reads_M", "C_coverage", 
             "CG_coverage", "percent_cpg_meth", "percent_chg_meth", "percent_chh_meth", "MomEdu_detail_1", 
             "MomEdu_detail_2", "MomEdu_detail_3", "MomEdu_detail_4", "MomEdu_detail_5", "MomEdu_detail_7", 
             "MomEdu_detail_8", "home_ownership", "marital_status", "MomAgeYr", "Mat_Height_cm", "Mat_Weight_kg_PrePreg",
             "Mat_BMI_PrePreg", "DM1or2", "GDM", "PE", "parity", "dad_age", "SmokeYN_Pregnancy","cotinine_urine_ngml", 
             "final_creatinine_mgdl", "AllEQ_PV_YN_Mo_3", "AllEQ_PV_YN_Mo_2", "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", 
             "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", "AllEQ_PV_YN_Mo4", "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", 
             "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", "AllEQ_PV_YN_Mo9", "AllEQ_tot_All_FA_mcg_Mo_3", 
             "AllEQ_tot_All_FA_mcg_Mo_2", "AllEQ_tot_All_FA_mcg_Mo_1", "AllEQ_tot_All_FA_mcg_Mo1", 
             "AllEQ_tot_All_FA_mcg_Mo2", "AllEQ_tot_All_FA_mcg_Mo3", "AllEQ_tot_All_FA_mcg_Mo4", 
             "AllEQ_tot_All_FA_mcg_Mo5", "AllEQ_tot_All_FA_mcg_Mo6", "AllEQ_tot_All_FA_mcg_Mo7", 
             "AllEQ_tot_All_FA_mcg_Mo8","AllEQ_tot_All_FA_mcg_Mo9")
covHeatmap(covStats, variableOrdering = "manual", regionOrdering = "variable", variables = variables,
           sortVariable = "Diagnosis_Alg", 
           file = "Figures/Diagnosis and Sex DMRs Raw Meth Covariate Heatmap Sorted by Diagnosis.png")
covHeatmap(covStats, variableOrdering = "hierarchical", regionOrdering = "hierarchical", 
           file = "Figures/Diagnosis and Sex DMRs Raw Meth Covariate Heatmap Clustered.png")

# Get Stats Adjusted for Sex
catVars <- catVars[!catVars == "Sex"]
covStats_adj <- DMRmethLm(DMRs = meth$DMRid, catVars = catVars, contVars = contVars, sampleData = samples_cov, 
                          file = "Tables/Diagnosis and Sex DMRs Raw Methylation by Covariate with AdjSex Stats.txt", adj = "Sex")
covSum_adj <- DMRmethLmSum(covStats_adj, file = "Tables/Diagnosis and Sex DMRs Raw Methylation by Covariate with AdjSex Summary.txt")

# Plot Heatmap Adjusted for Sex
variables <- variables[!variables == "Sex"]
covHeatmap(covStats_adj, variableOrdering = "manual", regionOrdering = "variable", variables = variables,
           sortVariable = "Diagnosis_Alg", 
           file = "Figures/Diagnosis and Sex DMRs Raw Meth Covariate with AdjSex Heatmap Sorted by Diagnosis.png")
covHeatmap(covStats_adj, variableOrdering = "hierarchical", regionOrdering = "hierarchical", 
           file = "Figures/Diagnosis and Sex DMRs Raw Meth Covariate with AdjSex Heatmap Clustered.png")

# Males Diagnosis DMRs ----------------------------------------------------
# Data #### 
# DMRs
DMRs <- loadRegions("DMRs/Discovery/Diagnosis Males 50/DMRs_Dx_Discovery50_males.csv",
                    chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), sort = TRUE)
DMRs$DMRid <- paste("DMR", 1:nrow(DMRs), sep = "_")
DMRs <- DMRs[,c("chr", "start", "end", "DMRid", "width", "L", "percentDifference", "stat", "pval", "qval")]

meth <- loadRegions("DMRs/Discovery/Diagnosis Males 50/DMR_smoothed_methylation_Dx_Discovery50_males.txt",
                     chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), sort = TRUE)
meth$DMRid <- paste("DMR", 1:nrow(meth), sep = "_")
meth <- meth[,c("chr", "start", "end", "DMRid", colnames(meth)[grepl("JLCM", colnames(meth), fixed = TRUE)])]

raw <- loadRegions("DMRs/Discovery/Diagnosis Males 50/DMR_raw_methylation_Dx_Discovery50_males.txt", 
                   chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), sort = TRUE)
raw$DMRid <- paste("DMR", 1:nrow(raw), sep = "_")
raw <- cbind(raw[,c("chr", "start", "end", "DMRid", colnames(raw)[grepl("JLCM", colnames(raw), fixed = TRUE)])])

# Candidates
candidates <- loadRegions("DMRs/Discovery/Diagnosis Males 50/CandidateRegions_Dx_Discovery50_males.csv",
                          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), sort = TRUE)
candidates <- candidates[,c("chr", "start", "end", "width", "L", "percentDifference", "stat", "pval", "qval")]

# Background
background <- loadRegions("DMRs/Discovery/Diagnosis Males 50/bsseq_background_Discovery50_males.csv",
                          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), sort = TRUE)

# Samples
samples <- read.delim(file = "Merged Database/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo.txt", 
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE)
samples <- subset(samples, Sequencing_ID %in% colnames(meth))

# Manhattan and QQ plots ####
# Manhattan plot
CMplotDMR(candidates, prefix = "Figures/Diagnosis Males DMRs", plot.type = "m")

# QQ plot
CMplotDMR(candidates, prefix = "Figures/Diagnosis Males DMRs", plot.type = "q")
rm(candidates)

# Heatmap ####
# Meth Data
meth_heat <- meth[,grepl("JLCM", colnames(meth), fixed = TRUE)]
methdiff <- (meth_heat - rowMeans(meth_heat, na.rm = TRUE)) %>% as.matrix %>% "*" (100)
hm.lim <- quantile(abs(methdiff), probs = 0.99, names = FALSE) %>% ceiling # 25

# Pheno Data
phenoData <- samples[,c("Sequencing_ID", "Diagnosis_Alg", "Study")]
colnames(phenoData) <- c("Sample", "Diagnosis", "Study")
phenoData <- phenoData[match(colnames(methdiff), phenoData$Sample),]
table(colnames(methdiff) == phenoData$Sample) # All TRUE
phenoData$Sample <- factor(phenoData$Sample, levels = unique(phenoData$Sample), ordered = TRUE)
phenoData$Diagnosis <- factor(phenoData$Diagnosis, levels = c("TD", "ASD"), ordered = TRUE)
phenoData$Study <- factor(phenoData$Study, levels = c("MARBLES", "EARLI"), ordered = TRUE)

# Plot Heatmap
pdf(file="Figures/Diagnosis Males DMRs Methylation Heatmap.pdf", width = 10, height = 8, onefile = FALSE)
buildHeatmap(x = methdiff, phenoData = phenoData, hm.low = - hm.lim, hm.high = hm.lim, 
             pheno.breaks = c("TD", "ASD", "MARBLES", "EARLI"), 
             pheno.values = c("TD" = "#3366CC", "ASD" = "#FF3366", "MARBLES" = "#FFFF33", "EARLI" = "#FF6633")) %>%
        printHeatmap
dev.off()
rm(hm.lim, methdiff)

# PCA ####
table(colnames(meth_heat) == phenoData$Sample) # All TRUE
data.pca <- prcomp(meth_heat %>% as.matrix %>% t, center = TRUE, scale. = TRUE)
ggbiplotPCA(data.pca = data.pca, groups = phenoData$Diagnosis, pc = c(1,2), 
            file = "Figures/Males Diagnosis DMRs PCA by Diagnosis.png", xlim = c(-30, 30), ylim = c(-30, 30))
ggbiplotPCA(data.pca = data.pca, groups = phenoData$Study, pc = c(1,2), 
            file = "Figures/Males Diagnosis DMRs PCA by Study.png", xlim = c(-30, 30), ylim = c(-30, 30),
            breaks = c("MARBLES", "EARLI"), values = c("MARBLES" = "#3366CC", "EARLI" = "#FF3366"),
            legend.position = c(0.77, 1.03))
rm(data.pca)

# Covariate Association ####
# Prep Data
meth_cov <- meth_heat %>% t %>% as.data.frame
colnames(meth_cov) <- meth$DMRid
meth_cov$Sequencing_ID <- rownames(meth_cov)
samples_cov <- merge(x = samples, y = meth_cov, by = "Sequencing_ID", all = FALSE, sort = FALSE)
catVars <- c("Study", "Platform", "Sex", "Site", "Diagnosis_Alg", "MomEdu_detail", "DM1or2", "GDM", "PE", 
             "home_ownership", "marital_status", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3", "AllEQ_PV_YN_Mo_2", 
             "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", "AllEQ_PV_YN_Mo4", 
             "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", "AllEQ_PV_YN_Mo9")
catVars <- catVars[!catVars %in% c("Platform", "Sex", "PE")] # Exclude catVars with only 1 level
samples_cov$Study <- factor(samples_cov$Study, levels = c("MARBLES", "EARLI"))
samples_cov$Platform <- factor(samples_cov$Platform, levels = c("HiSeqX10", "HiSeq4000"))
samples_cov$Sex <- factor(samples_cov$Sex, levels = c("M", "F"))
samples_cov$Site <- factor(samples_cov$Site, levels = c("UC Davis", "Kaiser Permanente", "Drexel", "Johns Hopkins University"))
samples_cov$Diagnosis_Alg <- factor(samples_cov$Diagnosis_Alg, levels = c("TD", "ASD"))
samples_cov$MomEdu_detail <- factor(samples_cov$MomEdu_detail, levels = c(6, 1:5, 7,8))
samples_cov$home_ownership[samples_cov$home_ownership == 99] <- NA
samples_cov$marital_status[samples_cov$marital_status == 99] <- NA
factorCols <- c("DM1or2", "GDM", "PE", "marital_status", "home_ownership", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3", 
                "AllEQ_PV_YN_Mo_2", "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", 
                "AllEQ_PV_YN_Mo4", "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", 
                "AllEQ_PV_YN_Mo9")
samples_cov[,factorCols] <- lapply(samples_cov[,factorCols], as.factor)
contVars <- colnames(samples_cov)[!colnames(samples_cov) %in% catVars & !colnames(samples_cov) %in% meth$DMRid &
                                          !colnames(samples_cov) %in% c("COI_ID", "Sequencing_ID", "Cord_Blood_IBC", 
                                                                        "Platform", "Sex", "PE")]

# Get Stats
covStats <- DMRmethLm(DMRs = meth$DMRid, catVars = catVars, contVars = contVars, sampleData = samples_cov, 
                      file = "Tables/Males Diagnosis DMRs Methylation by Covariate Stats.txt")
covSum <- DMRmethLmSum(covStats, file = "Tables/Males Diagnosis DMRs Methylation by Covariate Summary.txt")

# Plot Heatmap
variables<-c("Diagnosis_Alg", "Study", "Site_Drexel", "Site_Johns_Hopkins_University", "Site_Kaiser_Permanente", 
             "ADOScs", "MSLelcStandard36", "MSLelTscore36", "MSLfmTscore36", "MSLrlTscore36", "MSLvrTscore36", "ga_w", 
             "bw_g", "percent_trimmed", "percent_aligned", "percent_duplicate", "dedup_reads_M", "C_coverage", 
             "CG_coverage", "percent_cpg_meth", "percent_chg_meth", "percent_chh_meth", "MomEdu_detail_1", 
             "MomEdu_detail_2", "MomEdu_detail_3", "MomEdu_detail_4", "MomEdu_detail_5", "MomEdu_detail_7", 
             "MomEdu_detail_8", "home_ownership", "marital_status", "MomAgeYr", "Mat_Height_cm", "Mat_Weight_kg_PrePreg",
             "Mat_BMI_PrePreg", "DM1or2", "GDM", "parity", "dad_age", "SmokeYN_Pregnancy","cotinine_urine_ngml", 
             "final_creatinine_mgdl", "AllEQ_PV_YN_Mo_3", "AllEQ_PV_YN_Mo_2", "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", 
             "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", "AllEQ_PV_YN_Mo4", "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", 
             "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", "AllEQ_PV_YN_Mo9", "AllEQ_tot_All_FA_mcg_Mo_3", 
             "AllEQ_tot_All_FA_mcg_Mo_2", "AllEQ_tot_All_FA_mcg_Mo_1", "AllEQ_tot_All_FA_mcg_Mo1", 
             "AllEQ_tot_All_FA_mcg_Mo2", "AllEQ_tot_All_FA_mcg_Mo3", "AllEQ_tot_All_FA_mcg_Mo4", 
             "AllEQ_tot_All_FA_mcg_Mo5", "AllEQ_tot_All_FA_mcg_Mo6", "AllEQ_tot_All_FA_mcg_Mo7", 
             "AllEQ_tot_All_FA_mcg_Mo8","AllEQ_tot_All_FA_mcg_Mo9")

covHeatmap(covStats, variableOrdering = "manual", regionOrdering = "variable", variables = variables,
           sortVariable = "Diagnosis_Alg", 
           file = "Figures/Males Diagnosis DMRs Covariate Heatmap Sorted by Diagnosis.png")
covHeatmap(covStats, variableOrdering = "hierarchical", regionOrdering = "hierarchical", 
           file = "Figures/Males Diagnosis DMRs Covariate Heatmap Clustered.png")

# DMR Annotation ####
regDomains <- read.delim("Tables/Regulatory domains hg38.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
DMRs_anno <- getDMRanno(DMRstats = DMRs, regDomains = regDomains, file = "Tables/Males Diagnosis DMRs Annotation.txt")
DMRs_males_geneList <- list("All" = getDMRgeneList(DMRstats = DMRs, regDomains = regDomains, direction = "all", type = "gene_name"),
                              "Hyper" = getDMRgeneList(DMRstats = DMRs, regDomains = regDomains, direction = "hyper", type = "gene_name"),
                              "Hypo" = getDMRgeneList(DMRstats = DMRs, regDomains = regDomains, direction = "hypo", type = "gene_name"),
                              "Background" = getDMRgeneList(DMRstats = background, regDomains = regDomains, direction = "all", type = "gene_name"))
rm(regDomains, DMRs_anno)

# GREAT Analysis ####
# Make Files for GREAT (hg19, DMRs redefined to match background, background < 1M regions)
prepGREAT(DMRs = makeGRange(DMRs, direction = "all"), Background = makeGRange(background, direction = "all"), 
          fileName = "UCSC Tracks/Males Diagnosis DMRs Discovery hg19 for GREAT.bed", writeBack = TRUE,
          backName = "UCSC Tracks/Males Diagnosis Background Discovery hg19 for GREAT.bed", 
          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"))
prepGREAT(DMRs = makeGRange(DMRs, direction = "hyper"), Background = makeGRange(background, direction = "all"), 
          fileName = "UCSC Tracks/Males Diagnosis Hyper DMRs Discovery hg19 for GREAT.bed", writeBack = FALSE,
          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"))
prepGREAT(DMRs = makeGRange(DMRs, direction = "hypo"), Background = makeGRange(background, direction = "all"), 
          fileName = "UCSC Tracks/Males Diagnosis Hypo DMRs Discovery hg19 for GREAT.bed", writeBack = FALSE,
          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"))

# Get Enrichments from GREAT
BEDfile_Back <- "UCSC Tracks/Males Diagnosis Background Discovery hg19 for GREAT.bed.gz"
greatAll <- getGREATenrichments(BEDfile_DMR = "UCSC Tracks/Males Diagnosis DMRs Discovery hg19 for GREAT.bed.gz", 
                                BEDfile_Back = BEDfile_Back)
greatHyper <- getGREATenrichments(BEDfile_DMR = "UCSC Tracks/Males Diagnosis Hyper DMRs Discovery hg19 for GREAT.bed.gz", 
                                  BEDfile_Back = BEDfile_Back)
greatHypo <- getGREATenrichments(BEDfile_DMR = "UCSC Tracks/Males Diagnosis Hypo DMRs Discovery hg19 for GREAT.bed.gz", 
                                 BEDfile_Back = BEDfile_Back)
greatCombined <- rbind(greatAll, greatHyper, greatHypo)
greatCombined$Direction <- factor(c(rep("All", nrow(greatAll)), rep("Hyper", nrow(greatHyper)), 
                                    rep("Hypo", nrow(greatHypo))), levels = c("All", "Hyper", "Hypo"))
#write.table(greatCombined, "Tables/Males Diagnosis DMRs Discovery GREAT Combined Results.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
# No enriched terms

rm(greatAll, greatCombined, greatHyper, greatHypo, background, DMRs, BEDfile_Back, samples)

# Raw Meth Heatmap ####
# Meth Data
meth_heat <- raw[,grepl("JLCM", colnames(raw), fixed = TRUE)]
methdiff <- (meth_heat - rowMeans(meth_heat, na.rm = TRUE)) %>% as.matrix %>% "*" (100)
hm.lim <- quantile(abs(methdiff), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
table(is.na(methdiff))

# Pheno Data
phenoData <- samples[,c("Sequencing_ID", "Diagnosis_Alg", "Study")]
colnames(phenoData) <- c("Sample", "Diagnosis", "Study")
phenoData <- phenoData[match(colnames(methdiff), phenoData$Sample),]
table(colnames(methdiff) == phenoData$Sample) # All TRUE
phenoData$Sample <- factor(phenoData$Sample, levels = unique(phenoData$Sample), ordered = TRUE)
phenoData$Diagnosis <- factor(phenoData$Diagnosis, levels = c("TD", "ASD"), ordered = TRUE)
phenoData$Study <- factor(phenoData$Study, levels = c("MARBLES", "EARLI"), ordered = TRUE)

# Plot Heatmap
pdf(file="Figures/Males Diagnosis DMRs Raw Methylation Heatmap.pdf", width = 10, height = 8, onefile = FALSE)
buildHeatmap(x = methdiff, phenoData = phenoData[,c("Sample", "Diagnosis", "Study")], hm.low = - hm.lim, hm.high = hm.lim, 
             pheno.breaks = c("TD", "ASD", "MARBLES", "EARLI"), 
             pheno.values = c("TD" = "#3366CC", "ASD" = "#FF3366", "MARBLES" = "#FFFF33", "EARLI" = "#FF6633")) %>%
        printHeatmap
dev.off()
rm(hm.lim, methdiff)

# Raw Meth PCA ####
meth_heat <- raw[,grepl("JLCM", colnames(raw), fixed = TRUE)] %>% as.matrix
table(colnames(meth_heat) == phenoData$Sample) # All TRUE
table(is.na(meth_heat))
for(i in 1:nrow(meth_heat)){
        temp <- as.numeric(meth_heat[i,])
        temp[is.na(temp)] <- mean(temp, na.rm = TRUE) # Replace missing values with mean meth of that DMR
        meth_heat[i,] <- temp
}
table(is.na(meth_heat))
data.pca <- prcomp(t(meth_heat), center = TRUE, scale. = TRUE)
ggbiplotPCA(data.pca = data.pca, groups = phenoData$Diagnosis, pc = c(1,2), 
            file = "Figures/Males Diagnosis DMRs Raw Meth PCA by Diagnosis.png", xlim = c(-30, 30), ylim = c(-30, 30))
ggbiplotPCA(data.pca = data.pca, groups = phenoData$Study, pc = c(1,2), 
            file = "Figures/Males Diagnosis DMRs Raw Meth PCA by Study.png", xlim = c(-30, 30), ylim = c(-30, 30),
            breaks = c("MARBLES", "EARLI"), values = c("MARBLES" = "#3366CC", "EARLI" = "#FF3366"),
            legend.position = c(0.77, 1.03))
rm(data.pca)

# Raw Meth Covariate Association ####
# Prep Data
meth_cov <- raw[,grepl("JLCM", colnames(raw), fixed = TRUE)] %>% t %>% as.data.frame
colnames(meth_cov) <- meth$DMRid
meth_cov$Sequencing_ID <- rownames(meth_cov)
samples_cov <- merge(x = samples, y = meth_cov, by = "Sequencing_ID", all = FALSE, sort = FALSE)
catVars <- c("Study", "Platform", "Sex", "Site", "Diagnosis_Alg", "MomEdu_detail", "DM1or2", "GDM", "PE", 
             "home_ownership", "marital_status", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3", "AllEQ_PV_YN_Mo_2", 
             "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", "AllEQ_PV_YN_Mo4", 
             "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", "AllEQ_PV_YN_Mo9")
catVars <- catVars[!catVars %in% c("Platform", "Sex", "PE")] # Exclude catVars with only 1 level
samples_cov$Study <- factor(samples_cov$Study, levels = c("MARBLES", "EARLI"))
samples_cov$Platform <- factor(samples_cov$Platform, levels = c("HiSeqX10", "HiSeq4000"))
samples_cov$Sex <- factor(samples_cov$Sex, levels = c("M", "F"))
samples_cov$Site <- factor(samples_cov$Site, levels = c("UC Davis", "Kaiser Permanente", "Drexel", "Johns Hopkins University"))
samples_cov$Diagnosis_Alg <- factor(samples_cov$Diagnosis_Alg, levels = c("TD", "ASD"))
samples_cov$MomEdu_detail <- factor(samples_cov$MomEdu_detail, levels = c(6, 1:5, 7,8))
samples_cov$home_ownership[samples_cov$home_ownership == 99] <- NA
samples_cov$marital_status[samples_cov$marital_status == 99] <- NA
factorCols <- c("DM1or2", "GDM", "PE", "marital_status", "home_ownership", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3", 
                "AllEQ_PV_YN_Mo_2", "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", 
                "AllEQ_PV_YN_Mo4", "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", 
                "AllEQ_PV_YN_Mo9")
samples_cov[,factorCols] <- lapply(samples_cov[,factorCols], as.factor)
contVars <- colnames(samples_cov)[!colnames(samples_cov) %in% catVars & !colnames(samples_cov) %in% meth$DMRid &
                                          !colnames(samples_cov) %in% c("COI_ID", "Sequencing_ID", "Cord_Blood_IBC", 
                                                                        "Platform", "Sex", "PE")]
# Get Stats
covStats <- DMRmethLm(DMRs = meth$DMRid, catVars = catVars, contVars = contVars, sampleData = samples_cov, 
                      file = "Tables/Males Diagnosis DMRs Raw Methylation by Covariate Stats.txt")
# Error with DMR_19 and DM1or2
# Error with DMR_40 and DM1or2
# Error with DMR_64 and SmokeYN_Pregnancy
# Error with DMR_88 and DM1or2
# Error with DMR_111 and DM1or2
# Error with DMR_138 and DM1or2
# Error with DMR_187 and DM1or2
# Error with DMR_257 and DM1or2
# Error with DMR_286 and DM1or2
# Error with DMR_537 and DM1or2
covSum <- DMRmethLmSum(covStats, file = "Tables/Males Diagnosis DMRs Raw Methylation by Covariate Summary.txt")

# Plot Heatmap
variables<-c("Diagnosis_Alg", "Study", "Site_Drexel", "Site_Johns_Hopkins_University", "Site_Kaiser_Permanente", 
             "ADOScs", "MSLelcStandard36", "MSLelTscore36", "MSLfmTscore36", "MSLrlTscore36", "MSLvrTscore36", "ga_w", 
             "bw_g", "percent_trimmed", "percent_aligned", "percent_duplicate", "dedup_reads_M", "C_coverage", 
             "CG_coverage", "percent_cpg_meth", "percent_chg_meth", "percent_chh_meth", "MomEdu_detail_1", 
             "MomEdu_detail_2", "MomEdu_detail_3", "MomEdu_detail_4", "MomEdu_detail_5", "MomEdu_detail_7", 
             "MomEdu_detail_8", "home_ownership", "marital_status", "MomAgeYr", "Mat_Height_cm", "Mat_Weight_kg_PrePreg",
             "Mat_BMI_PrePreg", "DM1or2", "GDM", "parity", "dad_age", "SmokeYN_Pregnancy","cotinine_urine_ngml", 
             "final_creatinine_mgdl", "AllEQ_PV_YN_Mo_3", "AllEQ_PV_YN_Mo_2", "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", 
             "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", "AllEQ_PV_YN_Mo4", "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", 
             "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", "AllEQ_PV_YN_Mo9", "AllEQ_tot_All_FA_mcg_Mo_3", 
             "AllEQ_tot_All_FA_mcg_Mo_2", "AllEQ_tot_All_FA_mcg_Mo_1", "AllEQ_tot_All_FA_mcg_Mo1", 
             "AllEQ_tot_All_FA_mcg_Mo2", "AllEQ_tot_All_FA_mcg_Mo3", "AllEQ_tot_All_FA_mcg_Mo4", 
             "AllEQ_tot_All_FA_mcg_Mo5", "AllEQ_tot_All_FA_mcg_Mo6", "AllEQ_tot_All_FA_mcg_Mo7", 
             "AllEQ_tot_All_FA_mcg_Mo8","AllEQ_tot_All_FA_mcg_Mo9") # Excluded Sex and PE
covHeatmap(covStats, variableOrdering = "manual", regionOrdering = "variable", variables = variables,
           sortVariable = "Diagnosis_Alg", 
           file = "Figures/Males Diagnosis DMRs Raw Meth Covariate Heatmap Sorted by Diagnosis.png")
covHeatmap(covStats, variableOrdering = "hierarchical", regionOrdering = "hierarchical", 
           file = "Figures/Males Diagnosis DMRs Raw Meth Covariate Heatmap Clustered.png")

# Females Diagnosis DMRs --------------------------------------------------
# Data ####
# DMRs
DMRs <- loadRegions("DMRs/Discovery/Diagnosis Females 50/DMR_smoothed_methylation_Dx_Discovery50_females.txt",
                    chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)
DMRs$DMRid <- paste("DMR", 1:nrow(DMRs), sep = "_")
meth <- DMRs[,c("chr", "start", "end", "DMRid", colnames(DMRs)[grepl("JLCM", colnames(DMRs), fixed = TRUE)])]
DMRs <- DMRs[,c("chr", "start", "end", "DMRid", "width", "L", "percentDifference", "stat", "pval", "qval")]

raw <- loadRegions("DMRs/Discovery/Diagnosis Females 50/DMR_raw_methylation_Dx_Discovery50_females.txt", 
                   chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)
raw$DMRid <- paste("DMR", 1:nrow(raw), sep = "_")
raw <- cbind(raw[,c("chr", "start", "end", "DMRid", colnames(raw)[grepl("JLCM", colnames(raw), fixed = TRUE)])])

# Candidates
candidates <- loadRegions("DMRs/Discovery/Diagnosis Females 50/CandidateRegions_Dx_Discovery50_females.csv",
                          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)
candidates <- candidates[,c("chr", "start", "end", "width", "L", "percentDifference", "stat", "pval", "qval")]

# Background
background <- loadRegions("DMRs/Discovery/Diagnosis Females 50/bsseq_background_Discovery50_females.csv",
                          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)

# Samples
samples <- read.delim(file = "Merged Database/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo.txt", 
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE)
samples <- subset(samples, Sequencing_ID %in% colnames(meth))

# Manhattan and QQ plots ####
# Manhattan plot
CMplotDMR(candidates, prefix = "Figures/Diagnosis Females DMRs", plot.type = "m")

# QQ plot
CMplotDMR(candidates, prefix = "Figures/Diagnosis Females DMRs", plot.type = "q")
rm(candidates)
# Heatmap ####
# Meth Data
meth_heat <- meth[,grepl("JLCM", colnames(meth), fixed = TRUE)]
methdiff <- (meth_heat - rowMeans(meth_heat, na.rm = TRUE)) %>% as.matrix %>% "*" (100)
hm.lim <- quantile(abs(methdiff), probs = 0.99, names = FALSE) %>% ceiling # 25

# Pheno Data
phenoData <- samples[,c("Sequencing_ID", "Diagnosis_Alg", "Study")]
colnames(phenoData) <- c("Sample", "Diagnosis", "Study")
phenoData <- phenoData[match(colnames(methdiff), phenoData$Sample),]
table(colnames(methdiff) == phenoData$Sample) # All TRUE
phenoData$Sample <- factor(phenoData$Sample, levels = unique(phenoData$Sample), ordered = TRUE)
phenoData$Diagnosis <- factor(phenoData$Diagnosis, levels = c("TD", "ASD"), ordered = TRUE)
phenoData$Study <- factor(phenoData$Study, levels = c("MARBLES", "EARLI"), ordered = TRUE)

# Plot Heatmap
pdf(file="Figures/Diagnosis Females DMRs Methylation Heatmap.pdf", width = 10, height = 8, onefile = FALSE)
buildHeatmap(x = methdiff, phenoData = phenoData, hm.low = - hm.lim, hm.high = hm.lim, 
             pheno.breaks = c("TD", "ASD", "MARBLES", "EARLI"), 
             pheno.values = c("TD" = "#3366CC", "ASD" = "#FF3366", "MARBLES" = "#FFFF33", "EARLI" = "#FF6633")) %>%
        printHeatmap
dev.off()
rm(hm.lim, methdiff)

# PCA ####
table(colnames(meth_heat) == phenoData$Sample) # All TRUE
data.pca <- prcomp(meth_heat %>% as.matrix %>% t, center = TRUE, scale. = TRUE)
ggbiplotPCA(data.pca = data.pca, groups = phenoData$Diagnosis, pc = c(1,2), 
            file = "Figures/Females Diagnosis DMRs PCA by Diagnosis.png", xlim = c(-40, 40), ylim = c(-40, 40))
ggbiplotPCA(data.pca = data.pca, groups = phenoData$Study, pc = c(1,2), 
            file = "Figures/Females Diagnosis DMRs PCA by Study.png", xlim = c(-40, 40), ylim = c(-40, 40),
            breaks = c("MARBLES", "EARLI"), values = c("MARBLES" = "#3366CC", "EARLI" = "#FF3366"),
            legend.position = c(0.77, 1.03))
rm(data.pca)

# Covariate Association ####
# Prep Data
meth_cov <- meth_heat %>% t %>% as.data.frame
colnames(meth_cov) <- meth$DMRid
meth_cov$Sequencing_ID <- rownames(meth_cov)
samples_cov <- merge(x = samples, y = meth_cov, by = "Sequencing_ID", all = FALSE, sort = FALSE)
catVars <- c("Study", "Platform", "Sex", "Site", "Diagnosis_Alg", "MomEdu_detail", "DM1or2", "GDM", "PE", 
             "home_ownership", "marital_status", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3", "AllEQ_PV_YN_Mo_2", 
             "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", "AllEQ_PV_YN_Mo4", 
             "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", "AllEQ_PV_YN_Mo9")
catVars <- catVars[!catVars %in% c("Platform", "Sex")] # Exclude catVars with only 1 level
samples_cov$Study <- factor(samples_cov$Study, levels = c("MARBLES", "EARLI"))
samples_cov$Platform <- factor(samples_cov$Platform, levels = c("HiSeqX10", "HiSeq4000"))
samples_cov$Sex <- factor(samples_cov$Sex, levels = c("M", "F"))
samples_cov$Site <- factor(samples_cov$Site, levels = c("UC Davis", "Kaiser Permanente", "Drexel", "Johns Hopkins University"))
samples_cov$Diagnosis_Alg <- factor(samples_cov$Diagnosis_Alg, levels = c("TD", "ASD"))
samples_cov$MomEdu_detail <- factor(samples_cov$MomEdu_detail, levels = c(6, 1:5, 7,8))
samples_cov$home_ownership[samples_cov$home_ownership == 99] <- NA
samples_cov$marital_status[samples_cov$marital_status == 99] <- NA
factorCols <- c("DM1or2", "GDM", "PE", "marital_status", "home_ownership", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3", 
                "AllEQ_PV_YN_Mo_2", "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", 
                "AllEQ_PV_YN_Mo4", "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", 
                "AllEQ_PV_YN_Mo9")
samples_cov[,factorCols] <- lapply(samples_cov[,factorCols], as.factor)
contVars <- colnames(samples_cov)[!colnames(samples_cov) %in% catVars & !colnames(samples_cov) %in% meth$DMRid &
                                          !colnames(samples_cov) %in% c("COI_ID", "Sequencing_ID", "Cord_Blood_IBC", 
                                                                        "Platform", "Sex")]

# Get Stats
covStats <- DMRmethLm(DMRs = meth$DMRid, catVars = catVars, contVars = contVars, sampleData = samples_cov, 
                      file = "Tables/Females Diagnosis DMRs Methylation by Covariate Stats.txt")
covSum <- DMRmethLmSum(covStats, file = "Tables/Females Diagnosis DMRs Methylation by Covariate Summary.txt")
covSum$Variable[covSum$FDRsig > 0]
# [1] "percent_cpg_meth"  "percent_chh_meth"  "GDM"               "C_coverage"        "dedup_reads_M"    
# [6] "CG_coverage"       "percent_duplicate" "Site" 

# Plot Heatmap
variables<-c("Diagnosis_Alg", "Study", "Site_Drexel", "Site_Johns_Hopkins_University", "Site_Kaiser_Permanente", 
             "ADOScs", "MSLelcStandard36", "MSLelTscore36", "MSLfmTscore36", "MSLrlTscore36", "MSLvrTscore36", "ga_w", 
             "bw_g", "percent_trimmed", "percent_aligned", "percent_duplicate", "dedup_reads_M", "C_coverage", 
             "CG_coverage", "percent_cpg_meth", "percent_chg_meth", "percent_chh_meth", "MomEdu_detail_1", 
             "MomEdu_detail_2", "MomEdu_detail_3", "MomEdu_detail_4", "MomEdu_detail_5", "MomEdu_detail_7", 
             "MomEdu_detail_8", "home_ownership", "marital_status", "MomAgeYr", "Mat_Height_cm", "Mat_Weight_kg_PrePreg",
             "Mat_BMI_PrePreg", "DM1or2", "GDM", "PE", "parity", "dad_age", "SmokeYN_Pregnancy","cotinine_urine_ngml", 
             "final_creatinine_mgdl", "AllEQ_PV_YN_Mo_3", "AllEQ_PV_YN_Mo_2", "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", 
             "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", "AllEQ_PV_YN_Mo4", "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", 
             "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", "AllEQ_PV_YN_Mo9", "AllEQ_tot_All_FA_mcg_Mo_3", 
             "AllEQ_tot_All_FA_mcg_Mo_2", "AllEQ_tot_All_FA_mcg_Mo_1", "AllEQ_tot_All_FA_mcg_Mo1", 
             "AllEQ_tot_All_FA_mcg_Mo2", "AllEQ_tot_All_FA_mcg_Mo3", "AllEQ_tot_All_FA_mcg_Mo4", 
             "AllEQ_tot_All_FA_mcg_Mo5", "AllEQ_tot_All_FA_mcg_Mo6", "AllEQ_tot_All_FA_mcg_Mo7", 
             "AllEQ_tot_All_FA_mcg_Mo8","AllEQ_tot_All_FA_mcg_Mo9")

covHeatmap(covStats, variableOrdering = "manual", regionOrdering = "variable", variables = variables,
           sortVariable = "Diagnosis_Alg", 
           file = "Figures/Females Diagnosis DMRs Covariate Heatmap Sorted by Diagnosis.png")
covHeatmap(covStats, variableOrdering = "hierarchical", regionOrdering = "hierarchical", 
           file = "Figures/Females Diagnosis DMRs Covariate Heatmap Clustered.png")

# DMR Annotation ####
regDomains <- read.delim("Tables/Regulatory domains hg38.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
DMRs_anno <- getDMRanno(DMRstats = DMRs, regDomains = regDomains, file = "Tables/Females Diagnosis DMRs Annotation.txt")
DMRs_females_geneList <- list("All" = getDMRgeneList(DMRstats = DMRs, regDomains = regDomains, direction = "all", type = "gene_name"),
                      "Hyper" = getDMRgeneList(DMRstats = DMRs, regDomains = regDomains, direction = "hyper", type = "gene_name"),
                      "Hypo" = getDMRgeneList(DMRstats = DMRs, regDomains = regDomains, direction = "hypo", type = "gene_name"),
                      "Background" = getDMRgeneList(DMRstats = background, regDomains = regDomains, direction = "all", type = "gene_name"))
rm(regDomains, DMRs_anno)

# GREAT Analysis ####
# Make Files for GREAT (hg19, DMRs redefined to match background, background < 1M regions)
prepGREAT(DMRs = makeGRange(DMRs, direction = "all"), Background = makeGRange(background, direction = "all"), 
          fileName = "UCSC Tracks/Females Diagnosis DMRs Discovery hg19 for GREAT.bed", writeBack = TRUE,
          backName = "UCSC Tracks/Females Diagnosis Background Discovery hg19 for GREAT.bed", 
          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"))
prepGREAT(DMRs = makeGRange(DMRs, direction = "hyper"), Background = makeGRange(background, direction = "all"), 
          fileName = "UCSC Tracks/Females Diagnosis Hyper DMRs Discovery hg19 for GREAT.bed", writeBack = FALSE,
          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"))
prepGREAT(DMRs = makeGRange(DMRs, direction = "hypo"), Background = makeGRange(background, direction = "all"), 
          fileName = "UCSC Tracks/Females Diagnosis Hypo DMRs Discovery hg19 for GREAT.bed", writeBack = FALSE,
          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"))

# Get Enrichments from GREAT
BEDfile_Back <- "UCSC Tracks/Females Diagnosis Background Discovery hg19 for GREAT.bed.gz"
greatAll <- getGREATenrichments(BEDfile_DMR = "UCSC Tracks/Females Diagnosis DMRs Discovery hg19 for GREAT.bed.gz", 
                                BEDfile_Back = BEDfile_Back)
greatHyper <- getGREATenrichments(BEDfile_DMR = "UCSC Tracks/Females Diagnosis Hyper DMRs Discovery hg19 for GREAT.bed.gz", 
                                BEDfile_Back = BEDfile_Back)
greatHypo <- getGREATenrichments(BEDfile_DMR = "UCSC Tracks/Females Diagnosis Hypo DMRs Discovery hg19 for GREAT.bed.gz", 
                                BEDfile_Back = BEDfile_Back)
greatCombined <- rbind(greatAll, greatHyper, greatHypo)
greatCombined$Direction <- factor(c(rep("All", nrow(greatAll)), rep("Hyper", nrow(greatHyper)), 
                                    rep("Hypo", nrow(greatHypo))), levels = c("All", "Hyper", "Hypo"))
write.table(greatCombined, "Tables/Females Diagnosis DMRs Discovery GREAT Combined Results.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# Plot All Results
greatPlotData <- subset(greatCombined, ID %in% # Get top 10 from each direction
                                unique(c(greatCombined$ID[greatCombined$Direction == "All"][1:10], 
                                         greatCombined$ID[greatCombined$Direction == "Hyper"][1:10],
                                         greatCombined$ID[greatCombined$Direction == "Hypo"][1:10])))
greatPlotData$name <- factor(greatPlotData$name, levels = rev(unique(greatPlotData$name)))
plotGREAT(greatPlotData, file = "Figures/Females Diagnosis Discovery GREAT Plot.png", axis.text.y.size = 9, 
          axis.text.y.width = 35, legend.position = c(1.28, 0.87))

# Plot GO Terms
greatPlotData <- subset(greatCombined, Ontology %in% c("GO Biological Process", "GO Cellular Component", "GO Molecular Function"))
greatPlotData <- subset(greatPlotData, ID %in% # Get top 10 from each direction
                                unique(c(greatPlotData$ID[greatPlotData$Direction == "All"][1:10], 
                                         greatPlotData$ID[greatPlotData$Direction == "Hyper"][1:10],
                                         greatPlotData$ID[greatPlotData$Direction == "Hypo"][1:10])))
greatPlotData$name <- factor(greatPlotData$name, levels = rev(unique(greatPlotData$name)))
plotGREAT(greatPlotData, file = "Figures/Females Diagnosis Discovery GREAT Plot GO Terms.png", axis.text.y.size = 10, 
          axis.text.y.width = 40, legend.position = c(1.31, 0.87))

# Plot Pathways
greatPlotData <- subset(greatCombined, Ontology %in% c("GO Biological Process", "PANTHER Pathway", "BioCyc Pathway",
                                                       "MSigDB Pathway"))
greatPlotData <- subset(greatPlotData, ID %in% # Get top 10 from each direction
                                unique(c(greatPlotData$ID[greatPlotData$Direction == "All"][1:10], 
                                         greatPlotData$ID[greatPlotData$Direction == "Hyper"][1:10],
                                         greatPlotData$ID[greatPlotData$Direction == "Hypo"][1:10])))
greatPlotData$name <- factor(greatPlotData$name, levels = rev(unique(greatPlotData$name)))
plotGREAT(greatPlotData, file = "Figures/Females Diagnosis Discovery GREAT Plot Pathways.png", axis.text.y.size = 9, 
          axis.text.y.width = 40, legend.position = c(1.31, 0.865))

# Plot Immunologic Signatures
greatPlotData <- subset(greatCombined, Ontology %in% c("MSigDB Immunologic Signatures"))
greatPlotData <- subset(greatPlotData, ID %in% # Get top 10 from each direction
                                unique(c(greatPlotData$ID[greatPlotData$Direction == "All"][1:5], 
                                         greatPlotData$ID[greatPlotData$Direction == "Hyper"][1:5],
                                         greatPlotData$ID[greatPlotData$Direction == "Hypo"][1:5])))
greatPlotData$name <- factor(greatPlotData$name, levels = rev(unique(greatPlotData$name)))
plotGREAT(greatPlotData, file = "Figures/Females Diagnosis Discovery GREAT Plot Immune.png", axis.text.y.size = 10, 
          axis.text.y.width = 60, legend.position = c(1.28, 0.87), wrap = TRUE)

rm(greatAll, greatCombined, greatHyper, greatHypo, background, DMRs, BEDfile_Back, samples)

# Raw Meth Heatmap ####
# Meth Data
meth_heat <- raw[,grepl("JLCM", colnames(raw), fixed = TRUE)]
methdiff <- (meth_heat - rowMeans(meth_heat, na.rm = TRUE)) %>% as.matrix %>% "*" (100)
hm.lim <- quantile(abs(methdiff), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
table(is.na(methdiff))

# Pheno Data
phenoData <- samples[,c("Sequencing_ID", "Diagnosis_Alg", "Study")]
colnames(phenoData) <- c("Sample", "Diagnosis", "Study")
phenoData <- phenoData[match(colnames(methdiff), phenoData$Sample),]
table(colnames(methdiff) == phenoData$Sample) # All TRUE
phenoData$Sample <- factor(phenoData$Sample, levels = unique(phenoData$Sample), ordered = TRUE)
phenoData$Diagnosis <- factor(phenoData$Diagnosis, levels = c("TD", "ASD"), ordered = TRUE)
phenoData$Study <- factor(phenoData$Study, levels = c("MARBLES", "EARLI"), ordered = TRUE)

# Plot Heatmap
pdf(file="Figures/Females Diagnosis DMRs Raw Methylation Heatmap.pdf", width = 10, height = 8, onefile = FALSE)
buildHeatmap(x = methdiff, phenoData = phenoData[,c("Sample", "Diagnosis", "Study")], hm.low = - hm.lim, hm.high = hm.lim, 
             pheno.breaks = c("TD", "ASD", "MARBLES", "EARLI"), 
             pheno.values = c("TD" = "#3366CC", "ASD" = "#FF3366", "MARBLES" = "#FFFF33", "EARLI" = "#FF6633")) %>%
        printHeatmap
dev.off()
rm(hm.lim, methdiff)

# Raw Meth PCA ####
meth_heat <- raw[,grepl("JLCM", colnames(raw), fixed = TRUE)] %>% as.matrix
table(colnames(meth_heat) == phenoData$Sample) # All TRUE
table(is.na(meth_heat))
for(i in 1:nrow(meth_heat)){
        temp <- as.numeric(meth_heat[i,])
        temp[is.na(temp)] <- mean(temp, na.rm = TRUE) # Replace missing values with mean meth of that DMR
        meth_heat[i,] <- temp
}
table(is.na(meth_heat))
data.pca <- prcomp(t(meth_heat), center = TRUE, scale. = TRUE)
ggbiplotPCA(data.pca = data.pca, groups = phenoData$Diagnosis, pc = c(1,2), 
            file = "Figures/Females Diagnosis DMRs Raw Meth PCA by Diagnosis.png", xlim = c(-60, 60), ylim = c(-60, 60))
ggbiplotPCA(data.pca = data.pca, groups = phenoData$Study, pc = c(1,2), 
            file = "Figures/Females Diagnosis DMRs Raw Meth PCA by Study.png", xlim = c(-60, 60), ylim = c(-60, 60),
            breaks = c("MARBLES", "EARLI"), values = c("MARBLES" = "#3366CC", "EARLI" = "#FF3366"),
            legend.position = c(0.77, 1.03))
rm(data.pca)

# Raw Meth Covariate Association ####
# Prep Data
meth_cov <- raw[,grepl("JLCM", colnames(raw), fixed = TRUE)] %>% t %>% as.data.frame
colnames(meth_cov) <- meth$DMRid
meth_cov$Sequencing_ID <- rownames(meth_cov)
samples_cov <- merge(x = samples, y = meth_cov, by = "Sequencing_ID", all = FALSE, sort = FALSE)
catVars <- c("Study", "Platform", "Sex", "Site", "Diagnosis_Alg", "MomEdu_detail", "DM1or2", "GDM", "PE", 
             "home_ownership", "marital_status", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3", "AllEQ_PV_YN_Mo_2", 
             "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", "AllEQ_PV_YN_Mo4", 
             "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", "AllEQ_PV_YN_Mo9")
catVars <- catVars[!catVars %in% c("Platform", "Sex")] # Exclude catVars with only 1 level
samples_cov$Study <- factor(samples_cov$Study, levels = c("MARBLES", "EARLI"))
samples_cov$Platform <- factor(samples_cov$Platform, levels = c("HiSeqX10", "HiSeq4000"))
samples_cov$Sex <- factor(samples_cov$Sex, levels = c("M", "F"))
samples_cov$Site <- factor(samples_cov$Site, levels = c("UC Davis", "Kaiser Permanente", "Drexel", "Johns Hopkins University"))
samples_cov$Diagnosis_Alg <- factor(samples_cov$Diagnosis_Alg, levels = c("TD", "ASD"))
samples_cov$MomEdu_detail <- factor(samples_cov$MomEdu_detail, levels = c(6, 1:5, 7,8))
samples_cov$home_ownership[samples_cov$home_ownership == 99] <- NA
samples_cov$marital_status[samples_cov$marital_status == 99] <- NA
factorCols <- c("DM1or2", "GDM", "PE", "marital_status", "home_ownership", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3", 
                "AllEQ_PV_YN_Mo_2", "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", 
                "AllEQ_PV_YN_Mo4", "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", 
                "AllEQ_PV_YN_Mo9")
samples_cov[,factorCols] <- lapply(samples_cov[,factorCols], as.factor)
contVars <- colnames(samples_cov)[!colnames(samples_cov) %in% catVars & !colnames(samples_cov) %in% meth$DMRid &
                                          !colnames(samples_cov) %in% c("COI_ID", "Sequencing_ID", "Cord_Blood_IBC", 
                                                                        "Platform", "Sex")]
# Get Stats
covStats <- DMRmethLm(DMRs = meth$DMRid, catVars = catVars, contVars = contVars, sampleData = samples_cov, 
                      file = "Tables/Females Diagnosis DMRs Raw Methylation by Covariate Stats.txt")
# Some DMR + cov combinations had too much missing data
# Error with DMR_238 and DM1or2
# Error with DMR_743 and DM1or2
# Error with DMR_752 and AllEQ_PV_YN_Mo9
# Error with DMR_861 and AllEQ_PV_YN_Mo4
# Error with DMR_861 and AllEQ_PV_YN_Mo5
# Error with DMR_872 and DM1or2
# Error with DMR_1355 and DM1or2
# Error with DMR_1390 and DM1or2
# Error with DMR_1399 and DM1or2
# Error with DMR_1420 and DM1or2
# Error with DMR_1499 and AllEQ_PV_YN_Mo9
# Error with DMR_1640 and DM1or2
# Error with DMR_1713 and AllEQ_PV_YN_Mo4
# Error with DMR_1713 and AllEQ_PV_YN_Mo5
# Error with DMR_1713 and AllEQ_PV_YN_Mo6
# Error with DMR_1713 and AllEQ_PV_YN_Mo8
# Error with DMR_1713 and AllEQ_PV_YN_Mo9
covSum <- DMRmethLmSum(covStats, file = "Tables/Females Diagnosis DMRs Raw Methylation by Covariate Summary.txt")

# Plot Heatmap
variables<-c("Diagnosis_Alg", "Study", "Site_Drexel", "Site_Johns_Hopkins_University", "Site_Kaiser_Permanente", 
             "ADOScs", "MSLelcStandard36", "MSLelTscore36", "MSLfmTscore36", "MSLrlTscore36", "MSLvrTscore36", "ga_w", 
             "bw_g", "percent_trimmed", "percent_aligned", "percent_duplicate", "dedup_reads_M", "C_coverage", 
             "CG_coverage", "percent_cpg_meth", "percent_chg_meth", "percent_chh_meth", "MomEdu_detail_1", 
             "MomEdu_detail_2", "MomEdu_detail_3", "MomEdu_detail_4", "MomEdu_detail_5", "MomEdu_detail_7", 
             "MomEdu_detail_8", "home_ownership", "marital_status", "MomAgeYr", "Mat_Height_cm", "Mat_Weight_kg_PrePreg",
             "Mat_BMI_PrePreg", "DM1or2", "GDM", "PE", "parity", "dad_age", "SmokeYN_Pregnancy","cotinine_urine_ngml", 
             "final_creatinine_mgdl", "AllEQ_PV_YN_Mo_3", "AllEQ_PV_YN_Mo_2", "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", 
             "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", "AllEQ_PV_YN_Mo4", "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", 
             "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", "AllEQ_PV_YN_Mo9", "AllEQ_tot_All_FA_mcg_Mo_3", 
             "AllEQ_tot_All_FA_mcg_Mo_2", "AllEQ_tot_All_FA_mcg_Mo_1", "AllEQ_tot_All_FA_mcg_Mo1", 
             "AllEQ_tot_All_FA_mcg_Mo2", "AllEQ_tot_All_FA_mcg_Mo3", "AllEQ_tot_All_FA_mcg_Mo4", 
             "AllEQ_tot_All_FA_mcg_Mo5", "AllEQ_tot_All_FA_mcg_Mo6", "AllEQ_tot_All_FA_mcg_Mo7", 
             "AllEQ_tot_All_FA_mcg_Mo8","AllEQ_tot_All_FA_mcg_Mo9")
covHeatmap(covStats, variableOrdering = "manual", regionOrdering = "variable", variables = variables,
           sortVariable = "Diagnosis_Alg", 
           file = "Figures/Females Diagnosis DMRs Raw Meth Covariate Heatmap Sorted by Diagnosis.png")
covHeatmap(covStats, variableOrdering = "hierarchical", regionOrdering = "hierarchical", 
           file = "Figures/Females Diagnosis DMRs Raw Meth Covariate Heatmap Clustered.png")
