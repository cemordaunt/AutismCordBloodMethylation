# Replication DMR Analysis --------------------------------------------------
# Diagnosis and Sex
# Autism Cord Blood Methylation
# Charles Mordaunt
# 2/14/19

# Packages ####
sapply(c("tidyverse", "ggdendro", "scales", "ggplot2", "ggbiplot", "reshape", "grid", "RColorBrewer", "CMplot", "rlist",
         "annotatr", "GenomicRanges", "LOLA", "rtracklayer", "R.utils", "rGREAT", "limma", "DMRichR"), require, character.only = TRUE)

# Functions ####
source("R Scripts/DMR Analysis Functions.R")

# Diagnosis DMRs, no chrXY ------------------------------------------------
# Data ####
# DMRs
DMRs <- loadRegions("DMRs/Replication/Diagnosis 50/DMRs_DxNoXY_Replication50.csv",
                    chroms = c(paste("chr", 1:22, sep = ""), "chrM"), sort = TRUE)
DMRs$DMRid <- paste("DMR", 1:nrow(DMRs), sep = "_")
DMRs <- DMRs[,c("chr", "start", "end", "DMRid", "width", "L", "percentDifference", "stat", "pval", "qval")]
raw <- loadRegions("DMRs/Replication/Diagnosis 50/DMR_raw_methylation_DxNoXY_Replication50.txt", 
                   chroms = c(paste("chr", 1:22, sep = ""), "chrM"), sort = TRUE)
raw$DMRid <- paste("DMR", 1:nrow(raw), sep = "_")
raw <- cbind(raw[,c("chr", "start", "end", "DMRid", colnames(raw)[grepl("JLCM", colnames(raw), fixed = TRUE)])])

# Candidates
candidates <- loadRegions("DMRs/Replication/Diagnosis 50/CandidateRegions_DxNoXY_Replication50.csv",
                          chroms = c(paste("chr", 1:22, sep = ""), "chrM"), sort = TRUE)
candidates <- candidates[,c("chr", "start", "end", "width", "L", "percentDifference", "stat", "pval", "qval")]

# Background
background <- loadRegions("DMRs/Replication/Diagnosis 50/bsseq_background_Replication50.csv",
                          chroms = c(paste("chr", 1:22, sep = ""), "chrM"), sort = TRUE)

# Samples
samples <- read.delim(file = "Merged Database/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo.txt", 
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE)
samples <- subset(samples, Sequencing_ID %in% colnames(raw))

# Manhattan and QQ plots ####
# Manhattan plot
CMplotDMR(candidates, prefix = "Figures/Replication Diagnosis DMRs NoXY", plot.type = "m")

# QQ plot
CMplotDMR(candidates, prefix = "Figures/Replication Diagnosis DMRs NoXY", plot.type = "q")
rm(candidates)

# Raw Meth Heatmap ####
# Meth Data
meth_heat <- raw[,grepl("JLCM", colnames(raw), fixed = TRUE)]
methdiff <- (meth_heat - rowMeans(meth_heat, na.rm = TRUE)) %>% as.matrix %>% "*" (100)
hm.lim <- quantile(abs(methdiff), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
table(is.na(methdiff))

# Pheno Data
phenoData <- samples[,c("Sequencing_ID", "Diagnosis_Alg", "Sex")]
colnames(phenoData) <- c("Sample", "Diagnosis", "Sex")
phenoData <- phenoData[match(colnames(methdiff), phenoData$Sample),]
table(colnames(methdiff) == phenoData$Sample) # All TRUE
phenoData$Sample <- factor(phenoData$Sample, levels = unique(phenoData$Sample), ordered = TRUE)
phenoData$Diagnosis <- factor(phenoData$Diagnosis, levels = c("TD", "ASD"), ordered = TRUE)
phenoData$Sex <- factor(phenoData$Sex, levels = c("M", "F"))

# Plot Heatmap
pdf(file = "Figures/Replication Diagnosis DMRs NoXY Raw Methylation Heatmap.pdf", width = 10, height = 8, onefile = FALSE)
buildHeatmap(x = methdiff, phenoData = phenoData[,c("Sample", "Diagnosis", "Sex")], hm.low = - hm.lim, hm.high = hm.lim, 
             pheno.breaks = c("TD", "ASD", "M", "F"), 
             pheno.values = c("TD" = "#3366CC", "ASD" = "#FF3366", "M" = "#FFFF33", "F" = "#FF6633")) %>%
        printHeatmap(pheno.legend.position = c(0.897, 0.915))
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
            file = "Figures/Replication Diagnosis DMRs NoXY Raw Meth PCA by Diagnosis.png", xlim = c(-80, 80), 
            ylim = c(-80, 80))
ggbiplotPCA(data.pca = data.pca, groups = phenoData$Sex, pc = c(1,2), 
            file = "Figures/Replication Diagnosis DMRs NoXY Raw Meth PCA by Sex.png", xlim = c(-80, 80), 
            ylim = c(-80, 80), breaks = c("M", "F"), values = c("M" = "#3366CC", "F" = "#FF3366"),
            legend.position = c(0.9, 1.03))
rm(data.pca, i, temp, phenoData, meth_heat)

# Raw Meth Covariate Association ####
# Prep Data
meth_cov <- raw[,grepl("JLCM", colnames(raw), fixed = TRUE)] %>% t %>% as.data.frame
colnames(meth_cov) <- raw$DMRid
meth_cov$Sequencing_ID <- rownames(meth_cov)
samples_cov <- merge(x = samples, y = meth_cov, by = "Sequencing_ID", all = FALSE, sort = FALSE)
catVars <- c("Study", "Platform", "Sex", "Site", "Diagnosis_Alg", "MomEdu_detail", "DM1or2", "GDM", "PE", 
             "home_ownership", "marital_status", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3", "AllEQ_PV_YN_Mo_2", 
             "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", "AllEQ_PV_YN_Mo4", 
             "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", "AllEQ_PV_YN_Mo9")
catVars <- catVars[!catVars %in% c("Platform", "Study", "Site")] # Exclude catVars with only 1 level
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
contVars <- colnames(samples_cov)[!colnames(samples_cov) %in% catVars & !colnames(samples_cov) %in% raw$DMRid &
                                          !colnames(samples_cov) %in% c("COI_ID", "Sequencing_ID", "Cord_Blood_IBC", 
                                                                        "Platform", "Study", "Site")]
# Get Stats
covStats <- DMRmethLm(DMRs = raw$DMRid, catVars = catVars, contVars = contVars, sampleData = samples_cov, 
                      file = "Tables/Replication Diagnosis DMRs NoXY Raw Methylation by Covariate Stats.txt")
# Missing Data
# Error with DMR_568 and DM1or2
# Error with DMR_601 and DM1or2
# Error with DMR_1574 and DM1or2
# Error with DMR_2724 and DM1or2
# Error with DMR_3821 and DM1or2
# Error with DMR_4096 and DM1or2
covSum <- DMRmethLmSum(covStats, file = "Tables/Replication Diagnosis DMRs NoXY Raw Methylation by Covariate Summary.txt")

# Plot Heatmap
variables <- c("Diagnosis_Alg", "Sex", "Study", "Site_Drexel", "Site_Johns_Hopkins_University", "Site_Kaiser_Permanente", 
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
           file = "Figures/Replication Diagnosis DMRs NoXY Raw Meth Covariate Heatmap Sorted by Diagnosis.png")
covHeatmap(covStats, variableOrdering = "hierarchical", regionOrdering = "hierarchical", 
           file = "Figures/Replication Diagnosis DMRs NoXY Raw Meth Covariate Heatmap Clustered.png")
rm(meth_cov, samples_cov, catVars, contVars, factorCols, variables, covStats, covSum)

# DMR Annotation ####
regDomains <- read.delim("Tables/Regulatory domains hg38.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
DMRs_anno <- getDMRanno(DMRstats = DMRs, regDomains = regDomains, file = "Tables/Replication Diagnosis DMRs NoXY Annotation.txt")
rm(regDomains, DMRs_anno)

# GREAT Analysis ####
# Make Files for GREAT (hg19, DMRs redefined to match background, background < 1M regions)
prepGREAT(DMRs = makeGRange(DMRs, direction = "all"), Background = makeGRange(background, direction = "all"), 
          fileName = "UCSC Tracks/Diagnosis NoXY DMRs Replication hg19 for GREAT.bed", writeBack = TRUE,
          backName = "UCSC Tracks/Diagnosis NoXY Background Replication hg19 for GREAT.bed", 
          chroms = c(paste("chr", 1:22, sep = ""), "chrM"))
prepGREAT(DMRs = makeGRange(DMRs, direction = "hyper"), Background = makeGRange(background, direction = "all"), 
          fileName = "UCSC Tracks/Diagnosis NoXY Hyper DMRs Replication hg19 for GREAT.bed", writeBack = FALSE,
          chroms = c(paste("chr", 1:22, sep = ""), "chrM"))
prepGREAT(DMRs = makeGRange(DMRs, direction = "hypo"), Background = makeGRange(background, direction = "all"), 
          fileName = "UCSC Tracks/Diagnosis NoXY Hypo DMRs Replication hg19 for GREAT.bed", writeBack = FALSE,
          chroms = c(paste("chr", 1:22, sep = ""), "chrM"))

# Get Enrichments from GREAT
BEDfile_Back <- "UCSC Tracks/Diagnosis NoXY Background Replication hg19 for GREAT.bed.gz"
greatAll <- getGREATenrichments(BEDfile_DMR = "UCSC Tracks/Diagnosis NoXY DMRs Replication hg19 for GREAT.bed.gz", 
                                BEDfile_Back = BEDfile_Back)
greatHyper <- getGREATenrichments(BEDfile_DMR = "UCSC Tracks/Diagnosis NoXY Hyper DMRs Replication hg19 for GREAT.bed.gz", 
                                  BEDfile_Back = BEDfile_Back)
greatHypo <- getGREATenrichments(BEDfile_DMR = "UCSC Tracks/Diagnosis NoXY Hypo DMRs Replication hg19 for GREAT.bed.gz", 
                                 BEDfile_Back = BEDfile_Back)
greatCombined <- rbind(greatAll, greatHyper, greatHypo)
greatCombined$Direction <- factor(c(rep("All", nrow(greatAll)), rep("Hyper", nrow(greatHyper)), 
                                    rep("Hypo", nrow(greatHypo))), levels = c("All", "Hyper", "Hypo"))
write.table(greatCombined, "Tables/Replication Diagnosis NoXY DMRs GREAT Combined Results.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# Plot All Enrichments
greatPlotData <- greatCombined
greatPlotData$name <- factor(greatPlotData$name, levels = rev(unique(greatPlotData$name)))
plotGREAT(greatPlotData, file = "Figures/Replication Diagnosis NoXY DMRs GREAT Plot.png", axis.text.y.size = 9, 
          axis.text.y.width = 35, legend.position = c(1.28, 0.87))

# Plot Immune Enrichments
greatPlotData <- subset(greatCombined, Ontology == "MSigDB Immunologic Signatures")
greatPlotData$name <- factor(greatPlotData$name, levels = rev(unique(greatPlotData$name)))
plotGREAT(greatPlotData, file = "Figures/Replication Diagnosis NoXY DMRs GREAT Immune Plot.png", axis.text.y.size = 9, 
          axis.text.y.width = 60, legend.position = c(1.28, 0.87)) # wrapping is on

# Plot Not Immune Enrichments
greatPlotData <- subset(greatCombined, !Ontology == "MSigDB Immunologic Signatures")
greatPlotData$name <- factor(greatPlotData$name, levels = rev(unique(greatPlotData$name)))
plotGREAT(greatPlotData, file = "Figures/Replication Diagnosis NoXY DMRs GREAT Not Immune Plot.png", 
          axis.text.y.size = 10, axis.text.y.width = 35, legend.position = c(1.24, 0.87)) # wrapping is off

rm(greatAll, greatCombined, greatHyper, greatHypo, background, DMRs, BEDfile_Back, samples, greatPlotData)

# Diagnosis and Sex DMRs --------------------------------------------------
# Data ####
# DMRs
DMRs <- loadRegions("DMRs/Replication/Diagnosis and Sex 50/DMRs_DxAdjSex_Replication50.csv",
                    chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)
DMRs$DMRid <- paste("DMR", 1:nrow(DMRs), sep = "_")
DMRs <- DMRs[,c("chr", "start", "end", "DMRid", "width", "L", "percentDifference", "stat", "pval", "qval")]
raw <- loadRegions("DMRs/Replication/Diagnosis and Sex 50/DMR_raw_methylation_DxAdjSex_Replication50.txt", 
                   chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)
raw$DMRid <- paste("DMR", 1:nrow(raw), sep = "_")
raw <- cbind(raw[,c("chr", "start", "end", "DMRid", colnames(raw)[grepl("JLCM", colnames(raw), fixed = TRUE)])])

# Candidates
candidates <- loadRegions("DMRs/Replication/Diagnosis and Sex 50/backgroundRegions50.csv",
                          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)
candidates <- candidates[,c("chr", "start", "end", "width", "L", "percentDifference", "stat", "pval", "qval")]

# Background
background <- loadRegions("DMRs/Replication/Diagnosis and Sex 50/bsseq_background_Replication50.csv",
                          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)

# Samples
samples <- read.delim(file = "Merged Database/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo.txt", 
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE)
samples <- subset(samples, Sequencing_ID %in% colnames(raw))

# Manhattan and QQ plots ####
# Manhattan plot
CMplotDMR(candidates, prefix = "Figures/Replication Diagnosis and Sex DMRs", plot.type = "m")

# QQ plot
CMplotDMR(candidates, prefix = "Figures/Replication Diagnosis and Sex DMRs", plot.type = "q")
rm(candidates)

# Raw Meth Heatmap ####
# Meth Data
meth_heat <- raw[,grepl("JLCM", colnames(raw), fixed = TRUE)]
methdiff <- (meth_heat - rowMeans(meth_heat, na.rm = TRUE)) %>% as.matrix %>% "*" (100)
hm.lim <- quantile(abs(methdiff), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
table(is.na(methdiff))

# Pheno Data
phenoData <- samples[,c("Sequencing_ID", "Diagnosis_Alg", "Sex")]
colnames(phenoData) <- c("Sample", "Diagnosis", "Sex")
phenoData <- phenoData[match(colnames(methdiff), phenoData$Sample),]
table(colnames(methdiff) == phenoData$Sample) # All TRUE
phenoData$Sample <- factor(phenoData$Sample, levels = unique(phenoData$Sample), ordered = TRUE)
phenoData$Diagnosis <- factor(phenoData$Diagnosis, levels = c("TD", "ASD"), ordered = TRUE)
phenoData$Sex <- factor(phenoData$Sex, levels = c("M", "F"))

# Plot Heatmap
pdf(file="Figures/Replication Diagnosis and Sex DMRs Raw Methylation Heatmap.pdf", width = 10, height = 8, onefile = FALSE)
buildHeatmap(x = methdiff, phenoData = phenoData[,c("Sample", "Diagnosis", "Sex")], hm.low = - hm.lim, hm.high = hm.lim, 
             pheno.breaks = c("TD", "ASD", "M", "F"), 
             pheno.values = c("TD" = "#3366CC", "ASD" = "#FF3366", "M" = "#FFFF33", "F" = "#FF6633")) %>%
        printHeatmap(pheno.legend.position = c(0.897, 0.915))
dev.off()
rm(hm.lim, methdiff)

# Plot Heatmap with Meth Data Adjusted for Sex
mod <- model.matrix(~ Diagnosis, data = phenoData)
table(is.na(meth_heat))
meth_adj <- removeBatchEffect(meth_heat, batch = phenoData$Sex, design = mod)
methdiff <- (meth_adj - rowMeans(meth_adj, na.rm = TRUE)) %>% as.matrix %>% "*" (100)
hm.lim <- quantile(abs(methdiff), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling

pdf(file="Figures/Replication Diagnosis and Sex DMRs Raw Methylation AdjSex Heatmap.pdf", width = 10, height = 8, onefile = FALSE)
buildHeatmap(x = methdiff, phenoData = phenoData[c("Sample", "Diagnosis", "Sex")], hm.low = - hm.lim, hm.high = hm.lim, 
             pheno.breaks = c("TD", "ASD", "MARBLES", "EARLI", "M", "F"), 
             pheno.values = c("TD" = "#3366CC", "ASD" = "#FF3366", "M" = "#FFFF33", "F" = "#FF6633")) %>%
        printHeatmap(pheno.legend.position = c(0.897, 0.915))
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
            file = "Figures/Replication Diagnosis and Sex DMRs Raw Meth PCA by Diagnosis.png", 
            xlim = c(-80, 80), ylim = c(-80, 80))
ggbiplotPCA(data.pca = data.pca, groups = phenoData$Sex, pc = c(1,2), 
            file = "Figures/Replication Diagnosis and Sex DMRs Raw Meth PCA by Sex.png", 
            xlim = c(-80, 80), ylim = c(-80, 80),
            breaks = c("M", "F"), values = c("M" = "#3366CC", "F" = "#FF3366"),
            legend.position = c(0.9, 1.03))
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
            file = "Figures/Replication Diagnosis and Sex DMRs Raw Meth AdjSex PCA by Diagnosis.png", 
            xlim = c(-80, 80), ylim = c(-80, 80))
ggbiplotPCA(data.pca = data.pca, groups = phenoData$Sex, pc = c(1,2), 
            file = "Figures/Replication Diagnosis and Sex DMRs Raw Meth AdjSex PCA by Sex.png", 
            xlim = c(-80, 80), ylim = c(-80, 80),
            breaks = c("M", "F"), values = c("M" = "#3366CC", "F" = "#FF3366"),
            legend.position = c(0.9, 1.03))
rm(data.pca)

# Raw Meth Covariate Association ####
# Prep Data
meth_cov <- raw[,grepl("JLCM", colnames(raw), fixed = TRUE)] %>% t %>% as.data.frame
colnames(meth_cov) <- raw$DMRid
meth_cov$Sequencing_ID <- rownames(meth_cov)
samples_cov <- merge(x = samples, y = meth_cov, by = "Sequencing_ID", all = FALSE, sort = FALSE)
catVars <- c("Study", "Platform", "Sex", "Site", "Diagnosis_Alg", "MomEdu_detail", "DM1or2", "GDM", "PE", 
             "home_ownership", "marital_status", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3", "AllEQ_PV_YN_Mo_2", 
             "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", "AllEQ_PV_YN_Mo4", 
             "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", "AllEQ_PV_YN_Mo9")
catVars <- catVars[!catVars %in% c("Platform", "Study", "Site")] # Exclude catVars with only 1 level
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
contVars <- colnames(samples_cov)[!colnames(samples_cov) %in% catVars & !colnames(samples_cov) %in% raw$DMRid &
                                          !colnames(samples_cov) %in% c("COI_ID", "Sequencing_ID", "Cord_Blood_IBC", 
                                                                        "Platform", "Study", "Site")]
# Get Stats
covStats <- DMRmethLm(DMRs = raw$DMRid, catVars = catVars, contVars = contVars, sampleData = samples_cov, 
                      file = "Tables/Replication Diagnosis and Sex DMRs Raw Methylation by Covariate Stats.txt")
# Missing Data
# Error with DMR_900 and DM1or2
# Error with DMR_1587 and DM1or2
# Error with DMR_2455 and DM1or2
covSum <- DMRmethLmSum(covStats, file = "Tables/Replication Diagnosis and Sex DMRs Raw Methylation by Covariate Summary.txt")

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
           file = "Figures/Replication Diagnosis and Sex DMRs Raw Meth Covariate Heatmap Sorted by Diagnosis.png")
covHeatmap(covStats, variableOrdering = "hierarchical", regionOrdering = "hierarchical", 
           file = "Figures/Replication Diagnosis and Sex DMRs Raw Meth Covariate Heatmap Clustered.png")

# Get Stats Adjusted for Sex
catVars <- catVars[!catVars == "Sex"]
covStats_adj <- DMRmethLm(DMRs = raw$DMRid, catVars = catVars, contVars = contVars, sampleData = samples_cov, 
                          file = "Tables/Replication Diagnosis and Sex DMRs Raw Methylation by Covariate with AdjSex Stats.txt", adj = "Sex")
# Missing Data
# Error with DMR_900 and DM1or2
# Error with DMR_1587 and DM1or2
# Error with DMR_2455 and DM1or2
covSum_adj <- DMRmethLmSum(covStats_adj, file = "Tables/Replication Diagnosis and Sex DMRs Raw Methylation by Covariate with AdjSex Summary.txt")

# Plot Heatmap Adjusted for Sex
variables <- variables[!variables == "Sex"]
covHeatmap(covStats_adj, variableOrdering = "manual", regionOrdering = "variable", variables = variables,
           sortVariable = "Diagnosis_Alg", 
           file = "Figures/Replication Diagnosis and Sex DMRs Raw Meth Covariate with AdjSex Heatmap Sorted by Diagnosis.png")
covHeatmap(covStats_adj, variableOrdering = "hierarchical", regionOrdering = "hierarchical", 
           file = "Figures/Replication Diagnosis and Sex DMRs Raw Meth Covariate with AdjSex Heatmap Clustered.png")
rm(covStats, covStats_adj, covSum, covSum_adj, meth_adj, meth_cov, meth_heat, methdiff, mod, phenoData, raw, 
   samples, samples_cov, catVars, contVars, factorCols, hm.lim, i, temp, variables)

# DMR Annotation ####
regDomains <- read.delim("Tables/Regulatory domains hg38.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
DMRs_anno <- getDMRanno(DMRstats = DMRs, regDomains = regDomains, file = "Tables/Replication Diagnosis and Sex DMRs Annotation.txt")
rm(regDomains, DMRs_anno)

# GREAT Analysis ####
# Make Files for GREAT (hg19, DMRs redefined to match background, background < 1M regions)
prepGREAT(DMRs = makeGRange(DMRs, direction = "all"), Background = makeGRange(background, direction = "all"), 
          fileName = "UCSC Tracks/Diagnosis and Sex DMRs Replication hg19 for GREAT.bed", writeBack = TRUE,
          backName = "UCSC Tracks/Diagnosis and Sex Background Replication hg19 for GREAT.bed", 
          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"))
prepGREAT(DMRs = makeGRange(DMRs, direction = "hyper"), Background = makeGRange(background, direction = "all"), 
          fileName = "UCSC Tracks/Diagnosis and Sex Hyper DMRs Replication hg19 for GREAT.bed", writeBack = FALSE,
          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"))
prepGREAT(DMRs = makeGRange(DMRs, direction = "hypo"), Background = makeGRange(background, direction = "all"), 
          fileName = "UCSC Tracks/Diagnosis and Sex Hypo DMRs Replication hg19 for GREAT.bed", writeBack = FALSE,
          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"))

# Get Enrichments from GREAT
BEDfile_Back <- "UCSC Tracks/Diagnosis and Sex Background Replication hg19 for GREAT.bed.gz"
greatAll <- getGREATenrichments(BEDfile_DMR = "UCSC Tracks/Diagnosis and Sex DMRs Replication hg19 for GREAT.bed.gz", 
                                BEDfile_Back = BEDfile_Back)
greatHyper <- getGREATenrichments(BEDfile_DMR = "UCSC Tracks/Diagnosis and Sex Hyper DMRs Replication hg19 for GREAT.bed.gz", 
                                  BEDfile_Back = BEDfile_Back)
greatHypo <- getGREATenrichments(BEDfile_DMR = "UCSC Tracks/Diagnosis and Sex Hypo DMRs Replication hg19 for GREAT.bed.gz", 
                                 BEDfile_Back = BEDfile_Back)
greatCombined <- rbind(greatAll, greatHyper, greatHypo)
greatCombined$Direction <- factor(c(rep("All", nrow(greatAll)), rep("Hyper", nrow(greatHyper)), 
                                    rep("Hypo", nrow(greatHypo))), levels = c("All", "Hyper", "Hypo"))
write.table(greatCombined, "Tables/Replication Diagnosis and Sex DMRs GREAT Combined Results.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# Plot All Enrichments
greatPlotData <- subset(greatCombined, ID %in% # Get top 20 from each direction
                                unique(c(greatCombined$ID[greatCombined$Direction == "All"][1:20], 
                                         greatCombined$ID[greatCombined$Direction == "Hyper"][1:20],
                                         greatCombined$ID[greatCombined$Direction == "Hypo"][1:20])))
greatPlotData$name <- factor(greatPlotData$name, levels = rev(unique(greatPlotData$name)))
plotGREAT(greatPlotData, file = "Figures/Replication Diagnosis and Sex DMRs GREAT Plot.png", axis.text.y.size = 9, 
          axis.text.y.width = 35, legend.position = c(1.28, 0.87))

# Plot Immune Enrichments
greatPlotData <- subset(greatCombined, Ontology == "MSigDB Immunologic Signatures")
greatPlotData <- subset(greatPlotData, ID %in% # Get top 10 from each direction
                                unique(c(greatPlotData$ID[greatPlotData$Direction == "All"][1:10], 
                                         greatPlotData$ID[greatPlotData$Direction == "Hyper"][1:10],
                                         greatPlotData$ID[greatPlotData$Direction == "Hypo"][1:10])))
greatPlotData$name <- factor(greatPlotData$name, levels = rev(unique(greatPlotData$name)))
plotGREAT(greatPlotData, file = "Figures/Replication Diagnosis and Sex DMRs GREAT Immune Plot.png", axis.text.y.size = 9, 
          axis.text.y.width = 60, legend.position = c(1.28, 0.87)) # wrapping is on

# Plot Not Immune Enrichments
greatPlotData <- subset(greatCombined, !Ontology == "MSigDB Immunologic Signatures")
greatPlotData <- subset(greatPlotData, ID %in% # Get top 20 from each direction
                                unique(c(greatPlotData$ID[greatPlotData$Direction == "All"][1:20], 
                                         greatPlotData$ID[greatPlotData$Direction == "Hyper"][1:20],
                                         greatPlotData$ID[greatPlotData$Direction == "Hypo"][1:20])))
greatPlotData$name <- factor(greatPlotData$name, levels = rev(unique(greatPlotData$name)))
plotGREAT(greatPlotData, file = "Figures/Replication Diagnosis and Sex DMRs GREAT Not Immune Plot.png", 
          axis.text.y.size = 10, axis.text.y.width = 35, legend.position = c(1.2, 0.87)) # wrapping is off

rm(greatAll, greatCombined, greatHyper, greatHypo, background, DMRs, BEDfile_Back, samples)

# Males Diagnosis DMRs ------------------------------------------------
# Data ####
# DMRs
DMRs <- loadRegions("DMRs/Replication/Diagnosis Males 50/DMRs_Dx_Replication50_males.csv",
                    chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), sort = TRUE)
DMRs$DMRid <- paste("DMR", 1:nrow(DMRs), sep = "_")
DMRs <- DMRs[,c("chr", "start", "end", "DMRid", "width", "L", "percentDifference", "stat", "pval", "qval")]
raw <- loadRegions("DMRs/Replication/Diagnosis Males 50/DMR_raw_methylation_Dx_Replication50_males.txt", 
                   chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), sort = TRUE)
raw$DMRid <- paste("DMR", 1:nrow(raw), sep = "_")
raw <- cbind(raw[,c("chr", "start", "end", "DMRid", colnames(raw)[grepl("JLCM", colnames(raw), fixed = TRUE)])])

# Candidates
candidates <- loadRegions("DMRs/Replication/Diagnosis Males 50/CandidateRegions_Dx_Replication50_males.csv",
                          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), sort = TRUE)
candidates <- candidates[,c("chr", "start", "end", "width", "L", "percentDifference", "stat", "pval", "qval")]

# Background
background <- loadRegions("DMRs/Replication/Diagnosis Males 50/bsseq_background_Replication50_males.csv",
                          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), sort = TRUE)

# Samples
samples <- read.delim(file = "Merged Database/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo.txt", 
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE)
samples <- subset(samples, Sequencing_ID %in% colnames(raw))

# Manhattan and QQ plots ####
# Manhattan plot
CMplotDMR(candidates, prefix = "Figures/Replication Males Diagnosis DMRs", plot.type = "m")

# QQ plot
CMplotDMR(candidates, prefix = "Figures/Replication Males Diagnosis DMRs", plot.type = "q")
rm(candidates)

# Raw Meth Heatmap ####
# Meth Data
meth_heat <- raw[,grepl("JLCM", colnames(raw), fixed = TRUE)]
methdiff <- (meth_heat - rowMeans(meth_heat, na.rm = TRUE)) %>% as.matrix %>% "*" (100)
hm.lim <- quantile(abs(methdiff), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
table(is.na(methdiff))

# Pheno Data
phenoData <- samples[,c("Sequencing_ID", "Diagnosis_Alg", "Sex")]
colnames(phenoData) <- c("Sample", "Diagnosis", "Sex")
phenoData <- phenoData[match(colnames(methdiff), phenoData$Sample),]
table(colnames(methdiff) == phenoData$Sample) # All TRUE
phenoData$Sample <- factor(phenoData$Sample, levels = unique(phenoData$Sample), ordered = TRUE)
phenoData$Diagnosis <- factor(phenoData$Diagnosis, levels = c("TD", "ASD"), ordered = TRUE)
phenoData$Sex <- factor(phenoData$Sex, levels = c("M", "F"))

# Plot Heatmap
pdf(file = "Figures/Replication Males Diagnosis DMRs Raw Methylation Heatmap.pdf", width = 10, height = 8, onefile = FALSE)
buildHeatmap(x = methdiff, phenoData = phenoData[,c("Sample", "Diagnosis", "Sex")], hm.low = - hm.lim, hm.high = hm.lim, 
             pheno.breaks = c("TD", "ASD", "M", "F"), 
             pheno.values = c("TD" = "#3366CC", "ASD" = "#FF3366", "M" = "#FFFF33", "F" = "#FF6633")) %>%
        printHeatmap(pheno.legend.position = c(0.897, 0.915))
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
            file = "Figures/Replication Males Diagnosis DMRs Raw Meth PCA by Diagnosis.png", xlim = c(-80, 80), 
            ylim = c(-80, 80))
rm(data.pca, i, temp, phenoData, meth_heat)

# Raw Meth Covariate Association ####
# Prep Data
meth_cov <- raw[,grepl("JLCM", colnames(raw), fixed = TRUE)] %>% t %>% as.data.frame
colnames(meth_cov) <- raw$DMRid
meth_cov$Sequencing_ID <- rownames(meth_cov)
samples_cov <- merge(x = samples, y = meth_cov, by = "Sequencing_ID", all = FALSE, sort = FALSE)
catVars <- c("Study", "Platform", "Sex", "Site", "Diagnosis_Alg", "MomEdu_detail", "DM1or2", "GDM", "PE", 
             "home_ownership", "marital_status", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3", "AllEQ_PV_YN_Mo_2", 
             "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", "AllEQ_PV_YN_Mo4", 
             "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", "AllEQ_PV_YN_Mo9")
catVars <- catVars[!catVars %in% c("Platform", "Study", "Site", "Sex", "DM1or2")] # Exclude catVars with only 1 level
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
contVars <- colnames(samples_cov)[!colnames(samples_cov) %in% catVars & !colnames(samples_cov) %in% raw$DMRid &
                                          !colnames(samples_cov) %in% c("COI_ID", "Sequencing_ID", "Cord_Blood_IBC", 
                                                                        "Platform", "Study", "Site", "Sex", "DM1or2")]
# Get Stats
covStats <- DMRmethLm(DMRs = raw$DMRid, catVars = catVars, contVars = contVars, sampleData = samples_cov, 
                      file = "Tables/Replication Males Diagnosis DMRs Raw Methylation by Covariate Stats.txt")
# Missing Data
# Error with DMR_3747 and PE

covSum <- DMRmethLmSum(covStats, file = "Tables/Replication Males Diagnosis DMRs Raw Methylation by Covariate Summary.txt")

# Plot Heatmap
variables <- c("Diagnosis_Alg", "Sex", "Study", "Site_Drexel", "Site_Johns_Hopkins_University", "Site_Kaiser_Permanente", 
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
           file = "Figures/Replication Males Diagnosis DMRs Raw Meth Covariate Heatmap Sorted by Diagnosis.png")
covHeatmap(covStats, variableOrdering = "hierarchical", regionOrdering = "hierarchical", 
           file = "Figures/Replication Males Diagnosis DMRs Raw Meth Covariate Heatmap Clustered.png")
rm(meth_cov, samples_cov, catVars, contVars, factorCols, variables, covStats, covSum)

# DMR Annotation ####
regDomains <- read.delim("Tables/Regulatory domains hg38.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
DMRs_anno <- getDMRanno(DMRstats = DMRs, regDomains = regDomains, file = "Tables/Replication Males Diagnosis DMRs Annotation.txt")
rm(regDomains, DMRs_anno)

# GREAT Analysis ####
# Make Files for GREAT (hg19, DMRs redefined to match background, background < 1M regions)
prepGREAT(DMRs = makeGRange(DMRs, direction = "all"), Background = makeGRange(background, direction = "all"), 
          fileName = "UCSC Tracks/Males Diagnosis DMRs Replication hg19 for GREAT.bed", writeBack = TRUE,
          backName = "UCSC Tracks/Males Diagnosis Background Replication hg19 for GREAT.bed", 
          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"))
prepGREAT(DMRs = makeGRange(DMRs, direction = "hyper"), Background = makeGRange(background, direction = "all"), 
          fileName = "UCSC Tracks/Males Diagnosis Hyper DMRs Replication hg19 for GREAT.bed", writeBack = FALSE,
          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"))
prepGREAT(DMRs = makeGRange(DMRs, direction = "hypo"), Background = makeGRange(background, direction = "all"), 
          fileName = "UCSC Tracks/Males Diagnosis Hypo DMRs Replication hg19 for GREAT.bed", writeBack = FALSE,
          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"))
rm(background, DMRs, raw, samples)

# Get Enrichments from GREAT
BEDfile_Back <- "UCSC Tracks/Males Diagnosis Background Replication hg19 for GREAT.bed.gz"
greatAll <- getGREATenrichments(BEDfile_DMR = "UCSC Tracks/Males Diagnosis DMRs Replication hg19 for GREAT.bed.gz", 
                                BEDfile_Back = BEDfile_Back)
greatHyper <- getGREATenrichments(BEDfile_DMR = "UCSC Tracks/Males Diagnosis Hyper DMRs Replication hg19 for GREAT.bed.gz", 
                                  BEDfile_Back = BEDfile_Back)
greatHypo <- getGREATenrichments(BEDfile_DMR = "UCSC Tracks/Males Diagnosis Hypo DMRs Replication hg19 for GREAT.bed.gz", 
                                 BEDfile_Back = BEDfile_Back)
greatCombined <- rbind(greatAll, greatHyper, greatHypo)
greatCombined$Direction <- factor(c(rep("All", nrow(greatAll)), rep("Hyper", nrow(greatHyper)), 
                                    rep("Hypo", nrow(greatHypo))), levels = c("All", "Hyper", "Hypo"))
write.table(greatCombined, "Tables/Replication Males Diagnosis DMRs GREAT Combined Results.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# Plot All Enrichments
greatPlotData <- subset(greatCombined, ID %in% # Get top 10 from each direction
                                unique(c(greatCombined$ID[greatCombined$Direction == "All"][1:10], 
                                greatCombined$ID[greatCombined$Direction == "Hyper"][1:10],
                                greatCombined$ID[greatCombined$Direction == "Hypo"][1:10])))
greatPlotData$name <- factor(greatPlotData$name, levels = rev(unique(greatPlotData$name)))
plotGREAT(greatPlotData, file = "Figures/Replication Males Diagnosis DMRs GREAT Plot.png", axis.text.y.size = 9, 
          axis.text.y.width = 35, legend.position = c(1.28, 0.87))

# Plot Immune Enrichments
greatPlotData <- subset(greatCombined, Ontology == "MSigDB Immunologic Signatures")
greatPlotData <- subset(greatPlotData, ID %in% # Get top 10 from each direction
                                unique(c(greatPlotData$ID[greatPlotData$Direction == "All"][1:10], 
                                         greatPlotData$ID[greatPlotData$Direction == "Hyper"][1:10],
                                         greatPlotData$ID[greatPlotData$Direction == "Hypo"][1:10])))
greatPlotData$name <- factor(greatPlotData$name, levels = rev(unique(greatPlotData$name)))
plotGREAT(greatPlotData, file = "Figures/Replication Males Diagnosis DMRs GREAT Immune Plot.png", axis.text.y.size = 9, 
          axis.text.y.width = 60, legend.position = c(1.28, 0.87)) # wrapping is on

# Plot GO Terms
greatPlotData <- subset(greatCombined, Ontology %in% c("GO Biological Process", "GO Cellular Component", "GO Molecular Function"))
greatPlotData <- subset(greatPlotData, ID %in% # Get top 20 from each direction
                                unique(c(greatPlotData$ID[greatPlotData$Direction == "All"][1:20], 
                                         greatPlotData$ID[greatPlotData$Direction == "Hyper"][1:20],
                                         greatPlotData$ID[greatPlotData$Direction == "Hypo"][1:20])))
greatPlotData$name <- factor(greatPlotData$name, levels = rev(unique(greatPlotData$name)))
plotGREAT(greatPlotData, file = "Figures/Replication Males Diagnosis DMRs GREAT GO Term Plot.png", axis.text.y.size = 9, 
          axis.text.y.width = 35, legend.position = c(1.28, 0.87)) # wrapping is off

# Plot Mouse Phenotype
greatPlotData <- subset(greatCombined, Ontology == "Mouse Phenotype")
greatPlotData <- subset(greatPlotData, ID %in% # Get top 20 from each direction
                                unique(c(greatPlotData$ID[greatPlotData$Direction == "All"][1:20], 
                                         greatPlotData$ID[greatPlotData$Direction == "Hyper"][1:20],
                                         greatPlotData$ID[greatPlotData$Direction == "Hypo"][1:20])))
greatPlotData$name <- factor(greatPlotData$name, levels = rev(unique(greatPlotData$name)))
plotGREAT(greatPlotData, file = "Figures/Replication Males Diagnosis DMRs GREAT Mouse Pheno Plot.png", axis.text.y.size = 9, 
          axis.text.y.width = 35, legend.position = c(1.2, 0.87)) # wrapping is off

rm(greatAll, greatCombined, greatHyper, greatHypo, BEDfile_Back, greatPlotData)

# Females Diagnosis DMRs ------------------------------------------------
# Data ####
# DMRs
DMRs <- loadRegions("DMRs/Replication/Diagnosis Females 50/DMRs_Dx_Replication50_females.csv",
                    chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)
DMRs$DMRid <- paste("DMR", 1:nrow(DMRs), sep = "_")
DMRs <- DMRs[,c("chr", "start", "end", "DMRid", "width", "L", "percentDifference", "stat", "pval", "qval")]
raw <- loadRegions("DMRs/Replication/Diagnosis Females 50/DMR_raw_methylation_Dx_Replication50_females.txt", 
                   chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)
raw$DMRid <- paste("DMR", 1:nrow(raw), sep = "_")
raw <- cbind(raw[,c("chr", "start", "end", "DMRid", colnames(raw)[grepl("JLCM", colnames(raw), fixed = TRUE)])])

# Candidates
candidates <- loadRegions("DMRs/Replication/Diagnosis Females 50/CandidateRegions_Dx_Replication50_females.csv",
                          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)
candidates <- candidates[,c("chr", "start", "end", "width", "L", "percentDifference", "stat", "pval", "qval")]

# Background
background <- loadRegions("DMRs/Replication/Diagnosis Females 50/bsseq_background_Replication50_females.csv",
                          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)

# Samples
samples <- read.delim(file = "Merged Database/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo.txt", 
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE)
samples <- subset(samples, Sequencing_ID %in% colnames(raw))

# Manhattan and QQ plots ####
# Manhattan plot
CMplotDMR(candidates, prefix = "Figures/Replication Females Diagnosis DMRs", plot.type = "m", bin.max = 1800)

# QQ plot
CMplotDMR(candidates, prefix = "Figures/Replication Females Diagnosis DMRs", plot.type = "q")
rm(candidates)

# Raw Meth Heatmap ####
# Meth Data
meth_heat <- raw[,grepl("JLCM", colnames(raw), fixed = TRUE)]
methdiff <- (meth_heat - rowMeans(meth_heat, na.rm = TRUE)) %>% as.matrix %>% "*" (100)
hm.lim <- quantile(abs(methdiff), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
table(is.na(methdiff))

# Pheno Data
phenoData <- samples[,c("Sequencing_ID", "Diagnosis_Alg", "Sex")]
colnames(phenoData) <- c("Sample", "Diagnosis", "Sex")
phenoData <- phenoData[match(colnames(methdiff), phenoData$Sample),]
table(colnames(methdiff) == phenoData$Sample) # All TRUE
phenoData$Sample <- factor(phenoData$Sample, levels = unique(phenoData$Sample), ordered = TRUE)
phenoData$Diagnosis <- factor(phenoData$Diagnosis, levels = c("TD", "ASD"), ordered = TRUE)
phenoData$Sex <- factor(phenoData$Sex, levels = c("M", "F"))

# Plot Heatmap
pdf(file = "Figures/Replication Females Diagnosis DMRs Raw Methylation Heatmap.pdf", width = 10, height = 8, onefile = FALSE)
buildHeatmap(x = methdiff, phenoData = phenoData[,c("Sample", "Diagnosis", "Sex")], hm.low = - hm.lim, hm.high = hm.lim, 
             pheno.breaks = c("TD", "ASD", "M", "F"), 
             pheno.values = c("TD" = "#3366CC", "ASD" = "#FF3366", "M" = "#FFFF33", "F" = "#FF6633")) %>%
        printHeatmap(pheno.legend.position = c(0.897, 0.915))
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
            file = "Figures/Replication Females Diagnosis DMRs Raw Meth PCA by Diagnosis.png", xlim = c(-200, 200), 
            ylim = c(-200, 200))
rm(data.pca, i, temp, phenoData, meth_heat)

# Raw Meth Covariate Association ####
# Prep Data
meth_cov <- raw[,grepl("JLCM", colnames(raw), fixed = TRUE)] %>% t %>% as.data.frame
colnames(meth_cov) <- raw$DMRid
meth_cov$Sequencing_ID <- rownames(meth_cov)
samples_cov <- merge(x = samples, y = meth_cov, by = "Sequencing_ID", all = FALSE, sort = FALSE)
catVars <- c("Study", "Platform", "Sex", "Site", "Diagnosis_Alg", "MomEdu_detail", "DM1or2", "GDM", "PE", 
             "home_ownership", "marital_status", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3", "AllEQ_PV_YN_Mo_2", 
             "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", "AllEQ_PV_YN_Mo4", 
             "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", "AllEQ_PV_YN_Mo9")
catVars <- catVars[!catVars %in% c("Platform", "Study", "Site", "Sex", "PE", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3",
                                   "AllEQ_PV_YN_Mo_2", "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1")] # Exclude catVars with only 1 level
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
contVars <- colnames(samples_cov)[!colnames(samples_cov) %in% catVars & !colnames(samples_cov) %in% raw$DMRid &
                                          !colnames(samples_cov) %in% c("COI_ID", "Sequencing_ID", "Cord_Blood_IBC", 
                                                                        "Platform", "Study", "Site", "Sex", "PE", 
                                                                        "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3",
                                                                        "AllEQ_PV_YN_Mo_2", "AllEQ_PV_YN_Mo_1", 
                                                                        "AllEQ_PV_YN_Mo1")]
oneLevel <- NULL
for(i in 1:length(catVars)){
        if(length(levels(samples_cov[,catVars[i]])) == 1){oneLevel <- c(oneLevel, catVars[i])}
}
# "PE" "SmokeYN_Pregnancy" "AllEQ_PV_YN_Mo_3" "AllEQ_PV_YN_Mo_2" "AllEQ_PV_YN_Mo_1" "AllEQ_PV_YN_Mo1"

# Get Stats
covStats <- DMRmethLm(DMRs = raw$DMRid, catVars = catVars, contVars = contVars, sampleData = samples_cov, 
                      file = "Tables/Replication Females Diagnosis DMRs Raw Methylation by Covariate Stats.txt")
# Missing Data
# Error with DMR_793 and DM1or2
# Error with DMR_898 and GDM
# Error with DMR_898 and home_ownership
# Error with DMR_1136 and GDM
# Error with DMR_1136 and home_ownership
# Error with DMR_2924 and DM1or2
# Error with DMR_4767 and GDM
# Error with DMR_4767 and home_ownership
# Error with DMR_4767 and AllEQ_PV_YN_Mo4
# Error with DMR_5778 and DM1or2
# Error with DMR_6094 and DM1or2
# Error with DMR_6484 and DM1or2
# Error with DMR_12508 and GDM
# Error with DMR_12508 and home_ownership
# Error with DMR_12508 and marital_status
# Error with DMR_12508 and AllEQ_PV_YN_Mo3
# Error with DMR_12508 and AllEQ_PV_YN_Mo5
# Error with DMR_15295 and DM1or2
# Error with DMR_15675 and DM1or2
# Error with DMR_16278 and DM1or2
# Error with DMR_16293 and GDM
# Error with DMR_16293 and home_ownership
# Error with DMR_16725 and DM1or2

covSum <- DMRmethLmSum(covStats, file = "Tables/Replication Females Diagnosis DMRs Raw Methylation by Covariate Summary.txt")

# Plot Heatmap
variables <- c("Diagnosis_Alg", "Sex", "Study", "Site_Drexel", "Site_Johns_Hopkins_University", "Site_Kaiser_Permanente", 
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
           file = "Figures/Replication Females Diagnosis DMRs Raw Meth Covariate Heatmap Sorted by Diagnosis.png")
covHeatmap(covStats, variableOrdering = "hierarchical", regionOrdering = "hierarchical", 
           file = "Figures/Replication Females Diagnosis DMRs Raw Meth Covariate Heatmap Clustered.png")

rm(meth_cov, samples_cov, catVars, contVars, factorCols, variables, covStats, covSum)

# DMR Annotation ####
regDomains <- read.delim("Tables/Regulatory domains hg38.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
DMRs_anno <- getDMRanno(DMRstats = DMRs, regDomains = regDomains, file = "Tables/Replication Females Diagnosis DMRs Annotation.txt")
rm(regDomains, DMRs_anno)

# GREAT Analysis ####
# Make Files for GREAT (hg19, DMRs redefined to match background, background < 1M regions)
prepGREAT(DMRs = makeGRange(DMRs, direction = "all"), Background = makeGRange(background, direction = "all"), 
          fileName = "UCSC Tracks/Females Diagnosis DMRs Replication hg19 for GREAT.bed", writeBack = TRUE,
          backName = "UCSC Tracks/Females Diagnosis Background Replication hg19 for GREAT.bed", 
          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"))
prepGREAT(DMRs = makeGRange(DMRs, direction = "hyper"), Background = makeGRange(background, direction = "all"), 
          fileName = "UCSC Tracks/Females Diagnosis Hyper DMRs Replication hg19 for GREAT.bed", writeBack = FALSE,
          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"))
prepGREAT(DMRs = makeGRange(DMRs, direction = "hypo"), Background = makeGRange(background, direction = "all"), 
          fileName = "UCSC Tracks/Females Diagnosis Hypo DMRs Replication hg19 for GREAT.bed", writeBack = FALSE,
          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"))
rm(background, DMRs, raw, samples)

# Get Enrichments from GREAT
BEDfile_Back <- "UCSC Tracks/Females Diagnosis Background Replication hg19 for GREAT.bed.gz"
greatAll <- getGREATenrichments(BEDfile_DMR = "UCSC Tracks/Females Diagnosis DMRs Replication hg19 for GREAT.bed.gz", 
                                BEDfile_Back = BEDfile_Back)
greatHyper <- getGREATenrichments(BEDfile_DMR = "UCSC Tracks/Females Diagnosis Hyper DMRs Replication hg19 for GREAT.bed.gz", 
                                  BEDfile_Back = BEDfile_Back)
greatHypo <- getGREATenrichments(BEDfile_DMR = "UCSC Tracks/Females Diagnosis Hypo DMRs Replication hg19 for GREAT.bed.gz", 
                                 BEDfile_Back = BEDfile_Back)
greatCombined <- rbind(greatAll, greatHyper, greatHypo)
greatCombined$Direction <- factor(c(rep("All", nrow(greatAll)), rep("Hyper", nrow(greatHyper)), 
                                    rep("Hypo", nrow(greatHypo))), levels = c("All", "Hyper", "Hypo"))
write.table(greatCombined, "Tables/Replication Females Diagnosis DMRs GREAT Combined Results.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
table(greatCombined$Direction, greatCombined$Ontology)

# Plot All Enrichments
greatPlotData <- subset(greatCombined, ID %in% # Get top 10 from each direction
                                unique(c(greatCombined$ID[greatCombined$Direction == "All"][1:10], 
                                         greatCombined$ID[greatCombined$Direction == "Hyper"][1:10],
                                         greatCombined$ID[greatCombined$Direction == "Hypo"][1:10])))
greatPlotData$name <- factor(greatPlotData$name, levels = rev(unique(greatPlotData$name)))
plotGREAT(greatPlotData, file = "Figures/Replication Females Diagnosis DMRs GREAT Plot.png", axis.text.y.size = 9, 
          axis.text.y.width = 35, legend.position = c(1.25, 0.87), wrap = FALSE)

# Plot Immune Enrichments
greatPlotData <- subset(greatCombined, Ontology == "MSigDB Immunologic Signatures")
greatPlotData <- subset(greatPlotData, ID %in% # Get top 10 from each direction
                                unique(c(greatPlotData$ID[greatPlotData$Direction == "All"][1:10], 
                                         greatPlotData$ID[greatPlotData$Direction == "Hyper"][1:10],
                                         greatPlotData$ID[greatPlotData$Direction == "Hypo"][1:10])))
greatPlotData$name <- factor(greatPlotData$name, levels = rev(unique(greatPlotData$name)))
plotGREAT(greatPlotData, file = "Figures/Replication Females Diagnosis DMRs GREAT Immune Plot.png", axis.text.y.size = 9, 
          axis.text.y.width = 60, legend.position = c(1.28, 0.87), wrap = TRUE)

# Plot GO Terms
greatPlotData <- subset(greatCombined, Ontology %in% c("GO Biological Process", "GO Cellular Component", "GO Molecular Function"))
greatPlotData <- subset(greatPlotData, ID %in% # Get top 20 from each direction
                                unique(c(greatPlotData$ID[greatPlotData$Direction == "All"][1:20], 
                                         greatPlotData$ID[greatPlotData$Direction == "Hyper"][1:20],
                                         greatPlotData$ID[greatPlotData$Direction == "Hypo"][1:20])))
greatPlotData$name <- factor(greatPlotData$name, levels = rev(unique(greatPlotData$name)))
plotGREAT(greatPlotData, file = "Figures/Replication Females Diagnosis DMRs GREAT GO Term Plot.png", axis.text.y.size = 9, 
          axis.text.y.width = 35, legend.position = c(1.28, 0.87), wrap = FALSE) 

# Plot Mouse Phenotype
greatPlotData <- subset(greatCombined, Ontology == "Mouse Phenotype")
greatPlotData <- subset(greatPlotData, ID %in% # Get top 20 from each direction
                                unique(c(greatPlotData$ID[greatPlotData$Direction == "All"][1:20], 
                                         greatPlotData$ID[greatPlotData$Direction == "Hyper"][1:20],
                                         greatPlotData$ID[greatPlotData$Direction == "Hypo"][1:20])))
greatPlotData$name <- factor(greatPlotData$name, levels = rev(unique(greatPlotData$name)))
plotGREAT(greatPlotData, file = "Figures/Replication Females Diagnosis DMRs GREAT Mouse Pheno Plot.png", axis.text.y.size = 9, 
          axis.text.y.width = 35, legend.position = c(1.2, 0.87), wrap = FALSE) 

# Plot MSigDB Pathway
greatPlotData <- subset(greatCombined, Ontology == "MSigDB Pathway")
greatPlotData <- subset(greatPlotData, ID %in% # Get top 15 from each direction
                                unique(c(greatPlotData$ID[greatPlotData$Direction == "All"][1:15], 
                                         greatPlotData$ID[greatPlotData$Direction == "Hyper"][1:15],
                                         greatPlotData$ID[greatPlotData$Direction == "Hypo"][1:15])))
greatPlotData$name <- factor(greatPlotData$name, levels = rev(unique(greatPlotData$name)))
plotGREAT(greatPlotData, file = "Figures/Replication Females Diagnosis DMRs GREAT MSigDB Pathway Plot.png", axis.text.y.size = 9, 
          axis.text.y.width = 35, legend.position = c(1.25, 0.87), wrap = FALSE) 

rm(greatAll, greatCombined, greatHyper, greatHypo, BEDfile_Back, greatPlotData)

# Females Diagnosis 100 DMRs ------------------------------------------------
# Data ####
# DMRs
DMRs <- loadRegions("DMRs/Replication/Diagnosis Females 100/DMRs_Dx_Replication100_females.csv",
                    chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)
DMRs$DMRid <- paste("DMR", 1:nrow(DMRs), sep = "_")
DMRs <- DMRs[,c("chr", "start", "end", "DMRid", "width", "L", "percentDifference", "stat", "pval", "qval")]
raw <- loadRegions("DMRs/Replication/Diagnosis Females 100/DMR_raw_methylation_Dx_Replication100_females.txt", 
                   chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)
raw$DMRid <- paste("DMR", 1:nrow(raw), sep = "_")
raw <- cbind(raw[,c("chr", "start", "end", "DMRid", colnames(raw)[grepl("JLCM", colnames(raw), fixed = TRUE)])])

# Candidates
candidates <- loadRegions("DMRs/Replication/Diagnosis Females 100/CandidateRegions_Dx_Replication100_females.csv",
                          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)
candidates <- candidates[,c("chr", "start", "end", "width", "L", "percentDifference", "stat", "pval", "qval")]

# Background
background <- loadRegions("DMRs/Replication/Diagnosis Females 100/bsseq_background_Replication100_females.csv",
                          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"), sort = TRUE)

# Samples
samples <- read.delim(file = "Merged Database/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo.txt", 
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE)
samples <- subset(samples, Sequencing_ID %in% colnames(raw))

# Manhattan and QQ plots ####
# Manhattan plot
CMplotDMR(candidates, prefix = "Figures/Replication Females Diagnosis 100 DMRs", plot.type = "m", bin.max = 900)

# QQ plot
CMplotDMR(candidates, prefix = "Figures/Replication Females Diagnosis 100 DMRs", plot.type = "q")
rm(candidates)

# Raw Meth Heatmap ####
# Meth Data
meth_heat <- raw[,grepl("JLCM", colnames(raw), fixed = TRUE)]
methdiff <- (meth_heat - rowMeans(meth_heat, na.rm = TRUE)) %>% as.matrix %>% "*" (100)
hm.lim <- quantile(abs(methdiff), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
table(is.na(methdiff))

# Pheno Data
phenoData <- samples[,c("Sequencing_ID", "Diagnosis_Alg", "Sex")]
colnames(phenoData) <- c("Sample", "Diagnosis", "Sex")
phenoData <- phenoData[match(colnames(methdiff), phenoData$Sample),]
table(colnames(methdiff) == phenoData$Sample) # All TRUE
phenoData$Sample <- factor(phenoData$Sample, levels = unique(phenoData$Sample), ordered = TRUE)
phenoData$Diagnosis <- factor(phenoData$Diagnosis, levels = c("TD", "ASD"), ordered = TRUE)
phenoData$Sex <- factor(phenoData$Sex, levels = c("M", "F"))

# Plot Heatmap
pdf(file = "Figures/Replication Females Diagnosis 100 DMRs Raw Methylation Heatmap.pdf", width = 10, height = 8, onefile = FALSE)
buildHeatmap(x = methdiff, phenoData = phenoData[,c("Sample", "Diagnosis", "Sex")], hm.low = - hm.lim, hm.high = hm.lim, 
             pheno.breaks = c("TD", "ASD", "M", "F"), 
             pheno.values = c("TD" = "#3366CC", "ASD" = "#FF3366", "M" = "#FFFF33", "F" = "#FF6633")) %>%
        printHeatmap(pheno.legend.position = c(0.897, 0.915))
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
            file = "Figures/Replication Females Diagnosis 100 DMRs Raw Meth PCA by Diagnosis.png", xlim = c(-100, 100), 
            ylim = c(-100, 100))
rm(data.pca, i, temp, phenoData, meth_heat)

# Raw Meth Covariate Association ####
# Prep Data
meth_cov <- raw[,grepl("JLCM", colnames(raw), fixed = TRUE)] %>% t %>% as.data.frame
colnames(meth_cov) <- raw$DMRid
meth_cov$Sequencing_ID <- rownames(meth_cov)
samples_cov <- merge(x = samples, y = meth_cov, by = "Sequencing_ID", all = FALSE, sort = FALSE)
catVars <- c("Study", "Platform", "Sex", "Site", "Diagnosis_Alg", "MomEdu_detail", "DM1or2", "GDM", "PE", 
             "home_ownership", "marital_status", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3", "AllEQ_PV_YN_Mo_2", 
             "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", "AllEQ_PV_YN_Mo4", 
             "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", "AllEQ_PV_YN_Mo9")
catVars <- catVars[!catVars %in% c("Platform", "Study", "Site", "Sex", "PE", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3",
                                   "AllEQ_PV_YN_Mo_2", "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1")] # Exclude catVars with only 1 level
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
contVars <- colnames(samples_cov)[!colnames(samples_cov) %in% catVars & !colnames(samples_cov) %in% raw$DMRid &
                                          !colnames(samples_cov) %in% c("COI_ID", "Sequencing_ID", "Cord_Blood_IBC", 
                                                                        "Platform", "Study", "Site", "Sex", "PE", 
                                                                        "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3",
                                                                        "AllEQ_PV_YN_Mo_2", "AllEQ_PV_YN_Mo_1", 
                                                                        "AllEQ_PV_YN_Mo1")]

# Get Stats
covStats <- DMRmethLm(DMRs = raw$DMRid, catVars = catVars, contVars = contVars, sampleData = samples_cov, 
                      file = "Tables/Replication Females Diagnosis 100 DMRs Raw Methylation by Covariate Stats.txt")
# Missing Data

covSum <- DMRmethLmSum(covStats, file = "Tables/Replication Females Diagnosis 100 DMRs Raw Methylation by Covariate Summary.txt")

# Plot Heatmap
variables <- c("Diagnosis_Alg", "Sex", "Study", "Site_Drexel", "Site_Johns_Hopkins_University", "Site_Kaiser_Permanente", 
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
           file = "Figures/Replication Females Diagnosis 100 DMRs Raw Meth Covariate Heatmap Sorted by Diagnosis.png")
covHeatmap(covStats, variableOrdering = "hierarchical", regionOrdering = "hierarchical", 
           file = "Figures/Replication Females Diagnosis 100 DMRs Raw Meth Covariate Heatmap Clustered.png")

rm(meth_cov, samples_cov, catVars, contVars, factorCols, variables, covStats, covSum)

# DMR Annotation ####
regDomains <- read.delim("Tables/Regulatory domains hg38.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
DMRs_anno <- getDMRanno(DMRstats = DMRs, regDomains = regDomains, file = "Tables/Replication Females Diagnosis 100 DMRs Annotation.txt")
rm(regDomains, DMRs_anno)

# GREAT Analysis ####
# Make Files for GREAT (hg19, DMRs redefined to match background, background < 1M regions)
prepGREAT(DMRs = makeGRange(DMRs, direction = "all"), Background = makeGRange(background, direction = "all"), 
          fileName = "UCSC Tracks/Females Diagnosis 100 DMRs Replication hg19 for GREAT.bed", writeBack = TRUE,
          backName = "UCSC Tracks/Females Diagnosis 100 Background Replication hg19 for GREAT.bed", 
          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"))
prepGREAT(DMRs = makeGRange(DMRs, direction = "hyper"), Background = makeGRange(background, direction = "all"), 
          fileName = "UCSC Tracks/Females Diagnosis 100 Hyper DMRs Replication hg19 for GREAT.bed", writeBack = FALSE,
          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"))
prepGREAT(DMRs = makeGRange(DMRs, direction = "hypo"), Background = makeGRange(background, direction = "all"), 
          fileName = "UCSC Tracks/Females Diagnosis 100 Hypo DMRs Replication hg19 for GREAT.bed", writeBack = FALSE,
          chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrM"))
rm(background, DMRs, raw, samples)

# Get Enrichments from GREAT
BEDfile_Back <- "UCSC Tracks/Females Diagnosis 100 Background Replication hg19 for GREAT.bed.gz"
greatAll <- getGREATenrichments(BEDfile_DMR = "UCSC Tracks/Females Diagnosis 100 DMRs Replication hg19 for GREAT.bed.gz", 
                                BEDfile_Back = BEDfile_Back)
greatHyper <- getGREATenrichments(BEDfile_DMR = "UCSC Tracks/Females Diagnosis 100 Hyper DMRs Replication hg19 for GREAT.bed.gz", 
                                  BEDfile_Back = BEDfile_Back)
greatHypo <- getGREATenrichments(BEDfile_DMR = "UCSC Tracks/Females Diagnosis 100 Hypo DMRs Replication hg19 for GREAT.bed.gz", 
                                 BEDfile_Back = BEDfile_Back)
greatCombined <- rbind(greatAll, greatHyper, greatHypo)
greatCombined$Direction <- factor(c(rep("All", nrow(greatAll)), rep("Hyper", nrow(greatHyper)), 
                                    rep("Hypo", nrow(greatHypo))), levels = c("All", "Hyper", "Hypo"))
write.table(greatCombined, "Tables/Replication Females Diagnosis 100 DMRs GREAT Combined Results.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
(table(greatCombined$Direction, greatCombined$Ontology) * 100 / nrow(greatCombined)) %>% round(1)

# Plot All Enrichments
greatPlotData <- subset(greatCombined, ID %in% # Get top 5 from each direction
                                unique(c(greatCombined$ID[greatCombined$Direction == "All"][1:5], 
                                         greatCombined$ID[greatCombined$Direction == "Hyper"][1:5],
                                         greatCombined$ID[greatCombined$Direction == "Hypo"][1:5])))
greatPlotData$name <- factor(greatPlotData$name, levels = rev(unique(greatPlotData$name)))
plotGREAT(greatPlotData, file = "Figures/Replication Females Diagnosis 100 DMRs GREAT Plot.png", axis.text.y.size = 9, 
          axis.text.y.width = 60, legend.position = c(1.24, 0.87), wrap = TRUE)

# Plot Immune Enrichments (Same as all enrichments)
greatPlotData <- subset(greatCombined, Ontology == "MSigDB Immunologic Signatures")
greatPlotData <- subset(greatPlotData, ID %in% # Get top 5 from each direction
                                unique(c(greatPlotData$ID[greatPlotData$Direction == "All"][1:5], 
                                         greatPlotData$ID[greatPlotData$Direction == "Hyper"][1:5],
                                         greatPlotData$ID[greatPlotData$Direction == "Hypo"][1:5])))
greatPlotData$name <- factor(greatPlotData$name, levels = rev(unique(greatPlotData$name)))
plotGREAT(greatPlotData, file = "Figures/Replication Females Diagnosis 100 DMRs GREAT Immune Plot.png", 
          axis.text.y.size = 9, axis.text.y.width = 60, legend.position = c(1.24, 0.87), wrap = TRUE)

# Plot GO Terms
greatPlotData <- subset(greatCombined, Ontology %in% c("GO Biological Process", "GO Cellular Component", "GO Molecular Function"))
greatPlotData <- subset(greatPlotData, ID %in% # Get top 10 from each direction
                                unique(c(greatPlotData$ID[greatPlotData$Direction == "All"][1:10], 
                                         greatPlotData$ID[greatPlotData$Direction == "Hyper"][1:10],
                                         greatPlotData$ID[greatPlotData$Direction == "Hypo"][1:10])))
greatPlotData$name <- factor(greatPlotData$name, levels = rev(unique(greatPlotData$name)))
plotGREAT(greatPlotData, file = "Figures/Replication Females Diagnosis 100 DMRs GREAT GO Term Plot.png", 
          axis.text.y.size = 9, axis.text.y.width = 35, legend.position = c(1.2, 0.87), wrap = FALSE) 

# Plot Mouse Phenotype
greatPlotData <- subset(greatCombined, Ontology == "Mouse Phenotype")
greatPlotData <- subset(greatPlotData, ID %in% # Get top 15 from each direction
                                unique(c(greatPlotData$ID[greatPlotData$Direction == "All"][1:15], 
                                         greatPlotData$ID[greatPlotData$Direction == "Hyper"][1:15],
                                         greatPlotData$ID[greatPlotData$Direction == "Hypo"][1:15])))
greatPlotData$name <- factor(greatPlotData$name, levels = rev(unique(greatPlotData$name)))
plotGREAT(greatPlotData, file = "Figures/Replication Females Diagnosis 100 DMRs GREAT Mouse Pheno Plot.png", 
          axis.text.y.size = 9, axis.text.y.width = 35, legend.position = c(1.18, 0.87), wrap = FALSE) 

# Plot MSigDB Pathway
greatPlotData <- subset(greatCombined, Ontology == "MSigDB Pathway")
greatPlotData <- subset(greatPlotData, ID %in% # Get top 15 from each direction
                                unique(c(greatPlotData$ID[greatPlotData$Direction == "All"][1:15], 
                                         greatPlotData$ID[greatPlotData$Direction == "Hyper"][1:15],
                                         greatPlotData$ID[greatPlotData$Direction == "Hypo"][1:15])))
greatPlotData$name <- factor(greatPlotData$name, levels = rev(unique(greatPlotData$name)))
plotGREAT(greatPlotData, file = "Figures/Replication Females Diagnosis 100 DMRs GREAT MSigDB Pathway Plot.png", 
          axis.text.y.size = 9, axis.text.y.width = 35, legend.position = c(1.23, 0.87), wrap = FALSE) 

# Plot Human Phenotype
greatPlotData <- subset(greatCombined, Ontology == "Human Phenotype")
greatPlotData <- subset(greatPlotData, ID %in% # Get top 15 from each direction
                                unique(c(greatPlotData$ID[greatPlotData$Direction == "All"][1:15], 
                                         greatPlotData$ID[greatPlotData$Direction == "Hyper"][1:15],
                                         greatPlotData$ID[greatPlotData$Direction == "Hypo"][1:15])))
greatPlotData$name <- factor(greatPlotData$name, levels = rev(unique(greatPlotData$name)))
plotGREAT(greatPlotData, file = "Figures/Replication Females Diagnosis 100 DMRs GREAT Human Pheno Plot.png", 
          axis.text.y.size = 10, axis.text.y.width = 35, legend.position = c(1.16, 0.87), wrap = FALSE) 

rm(greatAll, greatCombined, greatHyper, greatHypo, BEDfile_Back, greatPlotData)

