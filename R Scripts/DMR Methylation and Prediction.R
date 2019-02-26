# DMR Methylation and Prediction ------------------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 2/23/19

# Packages ####
.libPaths("/share/lasallelab/programs/DMRichR/R_3.5")
sapply(c("reshape", "grid", "scales", "tidyverse", "ggdendro", "rlist", "openxlsx", "bsseq", "annotatr", "dmrseq", "DMRichR",
         "randomForest", "randomForestExplainer"), require, character.only = TRUE)

# Functions ####
source("R Scripts/DMR Analysis Functions.R")
regDomains <- read.delim("Tables/Regulatory domains hg38.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Get Discovery DMR Methylation in Replication Samples --------------------
# Discovery Diagnosis DMRs All Samples
bs.filtered <- readRDS("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Replication/BSseq_05Group.rds") # Replication samples
sigRegions <- read.csv("DMRs_DxNoXY_Discovery50.csv", header = TRUE, stringsAsFactors = FALSE)
raw <- as.data.frame(getMeth(BSseq = bs.filtered, regions = data.frame2GRanges(sigRegions), type = "raw", what = "perRegion"))
raw <- cbind(sigRegions, raw)
write.table(raw, "DxNoXY_Discovery50_DMR_raw_methylation_in_Replication.txt", sep = "\t", quote = FALSE, row.names = FALSE)
rm(bs.filtered, sigRegions, raw)

# Discovery Diagnosis and Sex DMRs All Samples
bs.filtered <- readRDS("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Replication/BSseq_05Group.rds") # Replication samples
sigRegions <- read.csv("DMRs_DxAdjSex_Discovery50.csv", header = TRUE, stringsAsFactors = FALSE)
raw <- as.data.frame(getMeth(BSseq = bs.filtered, regions = data.frame2GRanges(sigRegions), type = "raw", what = "perRegion"))
raw <- cbind(sigRegions, raw)
write.table(raw, "DxAdjSex_Discovery50_DMR_raw_methylation_in_Replication.txt", sep = "\t", quote = FALSE, row.names = FALSE)
rm(bs.filtered, sigRegions, raw)

# Discovery Diagnosis DMRs Males
bs.filtered <- readRDS("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Replication/Filtered_BSseq_Replication50_males.rds") # Replication samples
sigRegions <- read.csv("DMRs_Dx_Discovery50_males.csv", header = TRUE, stringsAsFactors = FALSE)
raw <- as.data.frame(getMeth(BSseq = bs.filtered, regions = data.frame2GRanges(sigRegions), type = "raw", what = "perRegion"))
raw <- cbind(sigRegions, raw)
write.table(raw, "Dx_Discovery50_males_DMR_raw_methylation_in_Replication.txt", sep = "\t", quote = FALSE, row.names = FALSE)
rm(bs.filtered, sigRegions, raw)

# Discovery Diagnosis DMRs Females
bs.filtered <- readRDS("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports/Replication/Filtered_BSseq_Replication50_females.rds") # Replication samples
sigRegions <- read.csv("DMRs_Dx_Discovery50_females.csv", header = TRUE, stringsAsFactors = FALSE)
raw <- as.data.frame(getMeth(BSseq = bs.filtered, regions = data.frame2GRanges(sigRegions), type = "raw", what = "perRegion"))
raw <- cbind(sigRegions, raw)
write.table(raw, "Dx_Discovery50_females_DMR_raw_methylation_in_Replication.txt", sep = "\t", quote = FALSE, row.names = FALSE)
rm(bs.filtered, sigRegions, raw)
# Run again with perGroup = 100 raw methylation?

# Differential Methylation Analysis ---------------------------------------
# Discovery Diagnosis DMRs All Samples ####
# Data
DMRmeth <- loadRegions("Tables/DxNoXY_Discovery50_DMR_raw_methylation_in_Replication.txt", 
                   chroms = c(paste("chr", 1:22, sep = ""), "chrM"), sort = TRUE)
DMRmeth$DMRid <- paste("DMR", 1:nrow(DMRmeth), sep = "_")
meth <- DMRmeth[,grepl("JLCM", colnames(DMRmeth), fixed = TRUE)]
phenoData <- read.delim(file = "Merged Database/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo.txt", 
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
        subset(Sequencing_ID %in% colnames(meth), select = c("Sequencing_ID", "Diagnosis_Alg", "Sex"))
colnames(phenoData) <- c("Sample", "Diagnosis", "Sex")
phenoData <- phenoData[match(colnames(meth), phenoData$Sample),]
table(colnames(meth) == phenoData$Sample) # All TRUE
phenoData$Sample <- factor(phenoData$Sample, levels = unique(phenoData$Sample), ordered = TRUE)
phenoData$Diagnosis <- factor(phenoData$Diagnosis, levels = c("TD", "ASD"), ordered = TRUE)
phenoData$Sex <- factor(phenoData$Sex, levels = c("M", "F"), ordered = TRUE)

# Remove DMRs missing data in all samples
missing <- apply(meth, 1, function(x) table(is.na(x))["TRUE"]) 
meth <- subset(meth, !missing == nrow(phenoData) | is.na(missing)) %>% t %>% as.data.frame
DMRmeth <- subset(DMRmeth, !missing == nrow(phenoData) | is.na(missing))
colnames(meth) <- DMRmeth$DMRid
meth$Sample <- rownames(meth)

# Differential Methylation
phenoData <- merge(x = phenoData, y = meth, by = "Sample", all = FALSE, sort = FALSE)
stats <- DMRmethLm(DMRmeth$DMRid, catVars = "Diagnosis", contVars = NULL, sampleData = phenoData, 
                   file = "Tables/Differential Methylation Discovery Diagnosis DMRs All in Replication Samples.txt")
stats <- cbind(DMRmeth[,c("chr", "start", "end", "width", "L", "area", "beta", "stat", "pval", "qval", 
                          "percentDifference")], stats)
sigDMRs <- subset(stats, pvalue < 0.05 & percentDifference * Estimate > 0) # same direction
colnames(sigDMRs)[colnames(sigDMRs) == "Region"] <- "DMRid"
sigDMRs <- getDMRanno(DMRstats = sigDMRs, regDomains = regDomains, 
                           file = "Tables/Significant Differential Methylation Discovery Diagnosis DMRs All in Replication Samples.txt")
nrow(sigDMRs) # 5

# Heatmap
meth_heat <- subset(DMRmeth, DMRid %in% sigDMRs$DMRid, select = colnames(DMRmeth)[grepl("JLCM", colnames(DMRmeth), fixed = TRUE)])
methdiff <- (meth_heat - rowMeans(meth_heat, na.rm = TRUE)) %>% as.matrix %>% "*" (100)
hm.lim <- quantile(abs(methdiff), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
table(colnames(methdiff) == phenoData$Sample) # All TRUE

pdf(file = "Figures/Differential Methylation Discovery Diagnosis DMRs All in Replication Samples Heatmap.pdf", width = 10, height = 6, onefile = FALSE)
buildHeatmap(x = methdiff, phenoData = phenoData[,c("Sample", "Diagnosis", "Sex")], hm.low = - hm.lim, hm.high = hm.lim, 
             pheno.breaks = c("TD", "ASD", "M", "F"), 
             pheno.values = c("TD" = "#3366CC", "ASD" = "#FF3366", "M" = "#FFFF33", "F" = "#FF6633")) %>%
        printHeatmap(widths = c(0.02, 0.8, 0.16, 0.02), heights = c(0.02, 0.18, 0.1, 0.68, 0.02),
                     heatmap.legend.position = c(0.65, -0.47), pheno.legend.position = c(0.93, 0.915))
dev.off()
rm(sigDMRs, stats, missing, meth_heat, methdiff, hm.lim, phenoData, meth)

# Random Forest Classification
samples <- read.delim(file = "Merged Database/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo.txt", 
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE)
discDMRmeth <- loadRegions("DMRs/Discovery/Diagnosis 50/DMR_raw_methylation_DxNoXY_Discovery50.txt", 
                           chroms = c(paste("chr", 1:22, sep = ""), "chrM"), sort = TRUE)
discDMRmeth$DMRid <- paste("DMR", 1:nrow(discDMRmeth), sep = "_")
discDMRmeth <- subset(discDMRmeth, DMRid %in% DMRmeth$DMRid)
table(DMRmeth$start == discDMRmeth$start) # All TRUE
set.seed(5)

DMRforest <- DMRrandomForest(discDMRmeth = discDMRmeth, repDMRmeth = DMRmeth, samples = samples)

pdf(file = "Figures/Discovery Diagnosis DMRs Random Forest Initial Model Plot.pdf", width = 6, height = 6)
plot(DMRforest$InitialModel)
dev.off()
pdf(file = "Figures/Discovery Diagnosis DMRs Random Forest Optimized Model Plot.pdf", width = 6, height = 6)
plot(DMRforest$OptimizedModel)
dev.off()

explain_forest(DMRforest$InitialModel)
file.copy(from = "Your_forest_explained.html", 
          to = "Tables/Discovery Diagnosis DMRs Random Forest Initial Model Explained.html")
explain_forest(DMRforest$OptimizedModel)
file.copy(from = "Your_forest_explained.html", 
          to = "Tables/Discovery Diagnosis DMRs Random Forest Optimized Model Explained.html")

rm(discDMRmeth, DMRforest, DMRmeth, samples)

# Discovery Diagnosis and Sex DMRs All Samples ####
# Data
DMRmeth <- loadRegions("Tables/DxAdjSex_Discovery50_DMR_raw_methylation_in_Replication.txt", 
                       chroms = c(paste("chr", 1:22, sep = ""), "chrM"), sort = TRUE)
DMRmeth$DMRid <- paste("DMR", 1:nrow(DMRmeth), sep = "_")
meth <- DMRmeth[,grepl("JLCM", colnames(DMRmeth), fixed = TRUE)]
phenoData <- read.delim(file = "Merged Database/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo.txt", 
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
        subset(Sequencing_ID %in% colnames(meth), select = c("Sequencing_ID", "Diagnosis_Alg", "Sex"))
colnames(phenoData) <- c("Sample", "Diagnosis", "Sex")
phenoData <- phenoData[match(colnames(meth), phenoData$Sample),]
table(colnames(meth) == phenoData$Sample) # All TRUE
phenoData$Sample <- factor(phenoData$Sample, levels = unique(phenoData$Sample), ordered = TRUE)
phenoData$Diagnosis <- factor(phenoData$Diagnosis, levels = c("TD", "ASD"), ordered = TRUE)
phenoData$Sex <- factor(phenoData$Sex, levels = c("M", "F"), ordered = TRUE)

# Remove DMRs missing data in all samples
missing <- apply(meth, 1, function(x) table(is.na(x))["TRUE"]) 
meth <- subset(meth, !missing == nrow(phenoData) | is.na(missing)) %>% t %>% as.data.frame
DMRmeth <- subset(DMRmeth, !missing == nrow(phenoData) | is.na(missing))
colnames(meth) <- DMRmeth$DMRid
meth$Sample <- rownames(meth)

# Differential Methylation
phenoData <- merge(x = phenoData, y = meth, by = "Sample", all = FALSE, sort = FALSE)
stats <- DMRmethLm(DMRmeth$DMRid, catVars = "Diagnosis", contVars = NULL, sampleData = phenoData, 
                   file = "Tables/Differential Methylation Discovery Diagnosis AdjSex DMRs All in Replication Samples.txt")
stats <- cbind(DMRmeth[,c("chr", "start", "end", "width", "L", "area", "beta", "stat", "pval", "qval", 
                          "percentDifference")], stats)
sigDMRs <- subset(stats, pvalue < 0.05 & percentDifference * Estimate > 0) # same direction
colnames(sigDMRs)[colnames(sigDMRs) == "Region"] <- "DMRid"
sigDMRs <- getDMRanno(DMRstats = sigDMRs, regDomains = regDomains, 
                      file = "Tables/Significant Differential Methylation Discovery Diagnosis AdjSex DMRs All in Replication Samples.txt")
nrow(sigDMRs) # 1
rm(meth, phenoData, sigDMRs, stats, missing)

# Random Forest Classification
samples <- read.delim(file = "Merged Database/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo.txt", 
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE)
discDMRmeth <- loadRegions("DMRs/Discovery/Diagnosis and Sex 50/DMR_raw_methylation_DxAdjSex_Discovery50.txt", 
                           chroms = c(paste("chr", 1:22, sep = ""), "chrM"), sort = TRUE)
discDMRmeth$DMRid <- paste("DMR", 1:nrow(discDMRmeth), sep = "_")
discDMRmeth <- subset(discDMRmeth, DMRid %in% DMRmeth$DMRid)
table(DMRmeth$start == discDMRmeth$start) # All TRUE
set.seed(5)

DMRforest <- DMRrandomForest(discDMRmeth = discDMRmeth, repDMRmeth = DMRmeth, samples = samples)

pdf(file = "Figures/Discovery Diagnosis and Sex DMRs Random Forest Initial Model Plot.pdf", width = 6, height = 6)
plot(DMRforest$InitialModel)
dev.off()
pdf(file = "Figures/Discovery Diagnosis and Sex DMRs Random Forest Optimized Model Plot.pdf", width = 6, height = 6)
plot(DMRforest$OptimizedModel)
dev.off()

explain_forest(DMRforest$InitialModel)
file.copy(from = "Your_forest_explained.html", 
          to = "Tables/Discovery Diagnosis and Sex DMRs Random Forest Initial Model Explained.html")
explain_forest(DMRforest$OptimizedModel)
file.copy(from = "Your_forest_explained.html", 
          to = "Tables/Discovery Diagnosis and Sex DMRs Random Forest Optimized Model Explained.html")

rm(discDMRmeth, DMRforest, DMRmeth, samples)

# Discovery Diagnosis DMRs Males ####
# Data
DMRmeth <- loadRegions("Tables/Dx_Discovery50_males_DMR_raw_methylation_in_Replication.txt", 
                       chroms = c(paste("chr", 1:22, sep = ""), "chrM"), sort = TRUE)
DMRmeth$DMRid <- paste("DMR", 1:nrow(DMRmeth), sep = "_")
meth <- DMRmeth[,grepl("JLCM", colnames(DMRmeth), fixed = TRUE)]
phenoData <- read.delim(file = "Merged Database/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo.txt", 
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
        subset(Sequencing_ID %in% colnames(meth), select = c("Sequencing_ID", "Diagnosis_Alg", "Sex"))
colnames(phenoData) <- c("Sample", "Diagnosis", "Sex")
phenoData <- phenoData[match(colnames(meth), phenoData$Sample),]
table(colnames(meth) == phenoData$Sample) # All TRUE
phenoData$Sample <- factor(phenoData$Sample, levels = unique(phenoData$Sample), ordered = TRUE)
phenoData$Diagnosis <- factor(phenoData$Diagnosis, levels = c("TD", "ASD"), ordered = TRUE)
phenoData$Sex <- factor(phenoData$Sex, levels = c("M", "F"), ordered = TRUE)

# Remove DMRs missing data in all samples
missing <- apply(meth, 1, function(x) table(is.na(x))["TRUE"]) 
meth <- subset(meth, !missing == nrow(phenoData) | is.na(missing)) %>% t %>% as.data.frame
DMRmeth <- subset(DMRmeth, !missing == nrow(phenoData) | is.na(missing))
colnames(meth) <- DMRmeth$DMRid
meth$Sample <- rownames(meth)

# Differential Methylation
phenoData <- merge(x = phenoData, y = meth, by = "Sample", all = FALSE, sort = FALSE)
stats <- DMRmethLm(DMRmeth$DMRid, catVars = "Diagnosis", contVars = NULL, sampleData = phenoData, 
                   file = "Tables/Differential Methylation Discovery Diagnosis Males DMRs All in Replication Samples.txt")
stats <- cbind(DMRmeth[,c("chr", "start", "end", "width", "L", "area", "beta", "stat", "pval", "qval", 
                          "percentDifference")], stats)
sigDMRs <- subset(stats, pvalue < 0.05 & percentDifference * Estimate > 0) # same direction
colnames(sigDMRs)[colnames(sigDMRs) == "Region"] <- "DMRid"
sigDMRs <- getDMRanno(DMRstats = sigDMRs, regDomains = regDomains, 
                      file = "Tables/Significant Differential Methylation Discovery Diagnosis Males DMRs All in Replication Samples.txt")
nrow(sigDMRs) # 11

# Heatmap
meth_heat <- subset(DMRmeth, DMRid %in% sigDMRs$DMRid, select = colnames(DMRmeth)[grepl("JLCM", colnames(DMRmeth), fixed = TRUE)])
methdiff <- (meth_heat - rowMeans(meth_heat, na.rm = TRUE)) %>% as.matrix %>% "*" (100)
hm.lim <- quantile(abs(methdiff), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
table(colnames(methdiff) == phenoData$Sample) # All TRUE

pdf(file = "Figures/Differential Methylation Discovery Diagnosis DMRs Males in Replication Samples Heatmap.pdf", width = 10, height = 6, onefile = FALSE)
buildHeatmap(x = methdiff, phenoData = phenoData[,c("Sample", "Diagnosis", "Sex")], hm.low = - hm.lim, hm.high = hm.lim, 
             pheno.breaks = c("TD", "ASD", "M", "F"), 
             pheno.values = c("TD" = "#3366CC", "ASD" = "#FF3366", "M" = "#FFFF33", "F" = "#FF6633")) %>%
        printHeatmap(widths = c(0.02, 0.8, 0.16, 0.02), heights = c(0.02, 0.18, 0.1, 0.68, 0.02),
                     heatmap.legend.position = c(0.65, -0.47), pheno.legend.position = c(0.93, 0.915))
dev.off()

rm(meth, phenoData, sigDMRs, stats, missing, meth_heat, methdiff, hm.lim)

# Random Forest Classification
samples <- read.delim(file = "Merged Database/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo.txt", 
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE)
discDMRmeth <- loadRegions("DMRs/Discovery/Diagnosis Males 50/DMR_raw_methylation_Dx_Discovery50_males.txt", 
                           chroms = c(paste("chr", 1:22, sep = ""), "chrM"), sort = TRUE)
discDMRmeth$DMRid <- paste("DMR", 1:nrow(discDMRmeth), sep = "_")
discDMRmeth <- subset(discDMRmeth, DMRid %in% DMRmeth$DMRid)
table(DMRmeth$start == discDMRmeth$start) # All TRUE
set.seed(5)

DMRforest <- DMRrandomForest(discDMRmeth = discDMRmeth, repDMRmeth = DMRmeth, samples = samples)

pdf(file = "Figures/Discovery Diagnosis Males DMRs Random Forest Initial Model Plot.pdf", width = 6, height = 6)
plot(DMRforest$InitialModel)
dev.off()
pdf(file = "Figures/Discovery Diagnosis Males DMRs Random Forest Optimized Model Plot.pdf", width = 6, height = 6)
plot(DMRforest$OptimizedModel)
dev.off()

explain_forest(DMRforest$InitialModel)
file.copy(from = "Your_forest_explained.html", 
          to = "Tables/Discovery Diagnosis Males DMRs Random Forest Initial Model Explained.html")
explain_forest(DMRforest$OptimizedModel)
file.copy(from = "Your_forest_explained.html", 
          to = "Tables/Discovery Diagnosis Males DMRs Random Forest Optimized Model Explained.html")

rm(discDMRmeth, DMRforest, DMRmeth, samples)

# Discovery Diagnosis DMRs Females ####
# Data
DMRmeth <- loadRegions("Tables/Dx_Discovery50_females_DMR_raw_methylation_in_Replication.txt", 
                       chroms = c(paste("chr", 1:22, sep = ""), "chrM"), sort = TRUE)
DMRmeth$DMRid <- paste("DMR", 1:nrow(DMRmeth), sep = "_")
meth <- DMRmeth[,grepl("JLCM", colnames(DMRmeth), fixed = TRUE)]
phenoData <- read.delim(file = "Merged Database/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo.txt", 
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
        subset(Sequencing_ID %in% colnames(meth), select = c("Sequencing_ID", "Diagnosis_Alg", "Sex"))
colnames(phenoData) <- c("Sample", "Diagnosis", "Sex")
phenoData <- phenoData[match(colnames(meth), phenoData$Sample),]
table(colnames(meth) == phenoData$Sample) # All TRUE
phenoData$Sample <- factor(phenoData$Sample, levels = unique(phenoData$Sample), ordered = TRUE)
phenoData$Diagnosis <- factor(phenoData$Diagnosis, levels = c("TD", "ASD"), ordered = TRUE)
phenoData$Sex <- factor(phenoData$Sex, levels = c("M", "F"), ordered = TRUE)

# Remove DMRs missing data in all samples
missing <- apply(meth, 1, function(x) table(is.na(x))["TRUE"]) 
meth <- subset(meth, !missing == nrow(phenoData) | is.na(missing)) %>% t %>% as.data.frame
DMRmeth <- subset(DMRmeth, !missing == nrow(phenoData) | is.na(missing))
colnames(meth) <- DMRmeth$DMRid
meth$Sample <- rownames(meth)

# Differential Methylation
phenoData <- merge(x = phenoData, y = meth, by = "Sample", all = FALSE, sort = FALSE)
stats <- DMRmethLm(DMRmeth$DMRid, catVars = "Diagnosis", contVars = NULL, sampleData = phenoData, 
                   file = "Tables/Differential Methylation Discovery Diagnosis Females DMRs All in Replication Samples.txt")
stats <- cbind(DMRmeth[,c("chr", "start", "end", "width", "L", "area", "beta", "stat", "pval", "qval", 
                          "percentDifference")], stats)
sigDMRs <- subset(stats, pvalue < 0.05 & percentDifference * Estimate > 0) # same direction
colnames(sigDMRs)[colnames(sigDMRs) == "Region"] <- "DMRid"
sigDMRs <- getDMRanno(DMRstats = sigDMRs, regDomains = regDomains, 
                      file = "Tables/Significant Differential Methylation Discovery Diagnosis Females DMRs All in Replication Samples.txt")
nrow(sigDMRs) # 38

# Heatmap
meth_heat <- subset(DMRmeth, DMRid %in% sigDMRs$DMRid, select = colnames(DMRmeth)[grepl("JLCM", colnames(DMRmeth), fixed = TRUE)])
methdiff <- (meth_heat - rowMeans(meth_heat, na.rm = TRUE)) %>% as.matrix %>% "*" (100)
hm.lim <- quantile(abs(methdiff), probs = 0.975, names = FALSE, na.rm = TRUE) %>% ceiling
table(colnames(methdiff) == phenoData$Sample) # All TRUE

pdf(file = "Figures/Differential Methylation Discovery Diagnosis DMRs Females in Replication Samples Heatmap.pdf", width = 10, height = 8, onefile = FALSE)
buildHeatmap(x = methdiff, phenoData = phenoData[,c("Sample", "Diagnosis", "Sex")], hm.low = - hm.lim, hm.high = hm.lim, 
             pheno.breaks = c("TD", "ASD", "M", "F"), 
             pheno.values = c("TD" = "#3366CC", "ASD" = "#FF3366", "M" = "#FFFF33", "F" = "#FF6633")) %>%
        printHeatmap(widths = c(0.02, 0.8, 0.16, 0.02), heights = c(0.02, 0.18, 0.1, 0.68, 0.02),
                     heatmap.legend.position = c(0.56, -0.42), pheno.legend.position = c(0.915, 0.87))
dev.off()

rm(meth, phenoData, sigDMRs, stats, missing, meth_heat, methdiff, hm.lim)

# Random Forest Classification
samples <- read.delim(file = "Merged Database/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo.txt", 
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE)
discDMRmeth <- loadRegions("DMRs/Discovery/Diagnosis Females 50/DMR_raw_methylation_Dx_Discovery50_females.txt", 
                           chroms = c(paste("chr", 1:22, sep = ""), "chrM"), sort = TRUE)
discDMRmeth$DMRid <- paste("DMR", 1:nrow(discDMRmeth), sep = "_")
discDMRmeth <- subset(discDMRmeth, DMRid %in% DMRmeth$DMRid)
table(DMRmeth$start == discDMRmeth$start) # All TRUE
set.seed(5)

DMRforest <- DMRrandomForest(discDMRmeth = discDMRmeth, repDMRmeth = DMRmeth, samples = samples)

pdf(file = "Figures/Discovery Diagnosis Females DMRs Random Forest Initial Model Plot.pdf", width = 6, height = 6)
plot(DMRforest$InitialModel)
dev.off()
pdf(file = "Figures/Discovery Diagnosis Females DMRs Random Forest Optimized Model Plot.pdf", width = 6, height = 6)
plot(DMRforest$OptimizedModel)
dev.off()

explain_forest(DMRforest$InitialModel)
file.copy(from = "Your_forest_explained.html", 
          to = "Tables/Discovery Diagnosis Females DMRs Random Forest Initial Model Explained.html")
explain_forest(DMRforest$OptimizedModel)
file.copy(from = "Your_forest_explained.html", 
          to = "Tables/Discovery Diagnosis Females DMRs Random Forest Optimized Model Explained.html")

rm(discDMRmeth, DMRforest, DMRmeth, samples, regDomains)
