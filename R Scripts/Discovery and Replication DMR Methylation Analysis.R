# Discovery and Replication DMR Methylation Analysis ----------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 9/16/19

# Packages ####
.libPaths("/share/lasallelab/programs/DMRichR/R_3.6")
sapply(c("reshape", "scales", "tidyverse", "bsseq", "annotatr", "dmrseq", "rlist"), require, character.only = TRUE)

# Functions ####
# Laptop
source("R Scripts/DMR Analysis Functions.R")

# Data ####
# Laptop
regDomains <- read.delim("Tables/Regulatory domains hg38.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
samples <- read.csv("Tables/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo and Cell Type.csv", header = TRUE,
                    stringsAsFactors = FALSE)
samples$Sequencing_ID <- factor(samples$Sequencing_ID, levels = unique(samples$Sequencing_ID))
samples$Diagnosis_Alg <- factor(samples$Diagnosis_Alg, levels = c("TD", "ASD"))
samples$Sex <- factor(samples$Sex, levels = c("M", "F"))

# Get DMR Methylation in Other Set --------------------
# Discovery DMR Methylation in Replication Samples ####
# Cluster
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports")

# Males
BSobj <- readRDS("Replication/Dx_Males/Filtered_BSseq_Replication50_males.rds")
DMRs <- read.csv("Discovery/Dx_Males/DMRs_Dx_Discovery50_males.csv", header = TRUE, stringsAsFactors = FALSE)
meth <- getMeth(BSseq = BSobj, regions = data.frame2GRanges(DMRs), type = "raw", what = "perRegion") %>% 
        as.data.frame() %>% cbind(DMRs, .)
write.csv(meth, "Discovery/Dx_Males/Dx_Discovery50_males_DMR_raw_methylation_in_Replication.csv", quote = FALSE, row.names = FALSE)

# Females
BSobj <- readRDS("Replication/Dx_Females_100/Filtered_BSseq_Replication100_females.rds")
DMRs <- read.csv("Discovery/Dx_Females/DMRs_Dx_Discovery50_females.csv", header = TRUE, stringsAsFactors = FALSE)
meth <- getMeth(BSseq = BSobj, regions = data.frame2GRanges(DMRs), type = "raw", what = "perRegion") %>% 
        as.data.frame() %>% cbind(DMRs, .)
write.csv(meth, "Discovery/Dx_Females/Dx_Discovery50_females_DMR_raw_methylation_in_Replication.csv", quote = FALSE, row.names = FALSE)

# Replication DMR Methylation in Discovery Samples ####
# Cluster
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports")

# Males
BSobj <- readRDS("Discovery/Dx_Males/Filtered_BSseq_Discovery50_males.rds")
DMRs <- read.csv("Replication/Dx_Males/DMRs_Dx_Replication50_males.csv", header = TRUE, stringsAsFactors = FALSE)
meth <- getMeth(BSseq = BSobj, regions = data.frame2GRanges(DMRs), type = "raw", what = "perRegion") %>% 
        as.data.frame() %>% cbind(DMRs, .)
write.csv(meth, "Replication/Dx_Males/Dx_Replication50_males_DMR_raw_methylation_in_Discovery.csv", quote = FALSE, row.names = FALSE)

# Females
BSobj <- readRDS("Discovery/Dx_Females/Filtered_BSseq_Discovery50_females.rds")
DMRs <- read.csv("Replication/Dx_Females_100/DMRs_Dx_Replication100_females.csv", header = TRUE, stringsAsFactors = FALSE)
meth <- getMeth(BSseq = BSobj, regions = data.frame2GRanges(DMRs), type = "raw", what = "perRegion") %>% 
        as.data.frame() %>% cbind(DMRs, .)
write.csv(meth, "Replication/Dx_Females_100/Dx_Replication100_females_DMR_raw_methylation_in_Discovery.csv", quote = FALSE, row.names = FALSE)

# Males Differential Methylation Analysis by Diagnosis ---------------------------------------
# Sample Data ####
malesDisc <- subset(samples, Sex == "M" & Platform == "HiSeqX10", select = c("Sequencing_ID", "Diagnosis_Alg", "Sex"))
malesRep <- subset(samples, Sex == "M" & Platform == "HiSeq4000", select = c("Sequencing_ID", "Diagnosis_Alg", "Sex"))
maleChroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")

# Discovery Males DMRs ####
# Data
methDisc <- loadRegions("DMRs/Discovery/Diagnosis Males 50/DMR_raw_methylation_Dx_Discovery50_males.txt", 
                        chroms = maleChroms, sort = TRUE, DMRid = TRUE)
methRep <- loadRegions("DMRs/Discovery/Diagnosis Males 50/Dx_Discovery50_males_DMR_raw_methylation_in_Replication.csv", 
                       chroms = maleChroms, sort = TRUE, DMRid = TRUE)
table(methDisc$start == methRep$start) # All TRUE

# Remove DMRs Covered in < 50% either TD or ASD in Replication Samples
meth <- methRep[,grepl("JLCM", colnames(methRep), fixed = TRUE)]
missingTD <- apply(meth[,as.character(malesRep$Sequencing_ID[malesRep$Diagnosis_Alg == "TD"])], 
                   1, function(x) table(is.na(x))["TRUE"]) / 
        length(as.character(malesRep$Sequencing_ID[malesRep$Diagnosis_Alg == "TD"])) # percent of TD with no reads
missingTD[is.na(missingTD)] <- 0
missingASD <- apply(meth[,as.character(malesRep$Sequencing_ID[malesRep$Diagnosis_Alg == "ASD"])], 
                    1, function(x) table(is.na(x))["TRUE"]) /
        length(as.character(malesRep$Sequencing_ID[malesRep$Diagnosis_Alg == "ASD"])) # percent of ASD with no reads
missingASD[is.na(missingASD)] <- 0
methRep <- subset(methRep, missingTD <= 0.5 & missingASD <= 0.5) # 584 / 635 DMRs remaining
methDisc <- subset(methDisc, DMRid %in% methRep$DMRid)

# Differential Methylation in Discovery Samples
meth <- t(methDisc[, grepl("JLCM", colnames(methDisc), fixed = TRUE)]) %>% as.data.frame
colnames(meth) <- methDisc$DMRid
meth$Sequencing_ID <- rownames(meth)
meth <- merge(x = malesDisc, y = meth, by = "Sequencing_ID", all = FALSE)
stats_DiscDMRsInDisc <- DMRmethLm(methDisc$DMRid, catVars = "Diagnosis_Alg", contVars = NULL, sampleData = meth,
                                  file = "Tables/Differential Methylation Discovery Males DMRs in Discovery Males.txt")
table(stats_DiscDMRsInDisc$pvalue < 0.05) # 551 / 584 sig DMRs

# Differential Methylation in Replication Samples
meth <- t(methRep[, grepl("JLCM", colnames(methRep), fixed = TRUE)]) %>% as.data.frame
colnames(meth) <- methRep$DMRid
meth$Sequencing_ID <- rownames(meth)
meth <- merge(x = malesRep, y = meth, by = "Sequencing_ID", all = FALSE)
stats_DiscDMRsInRep <- DMRmethLm(methRep$DMRid, catVars = "Diagnosis_Alg", contVars = NULL, sampleData = meth,
                                  file = "Tables/Differential Methylation Discovery Males DMRs in Replication Males.txt")
table(stats_DiscDMRsInRep$pvalue < 0.05) # 27 / 584 sig DMRs

# Subset by Differential Methylation in Both Sets and Annotate
colnames(stats_DiscDMRsInDisc) <- paste("Disc", colnames(stats_DiscDMRsInDisc), sep = "_")
colnames(stats_DiscDMRsInRep) <- paste("Rep", colnames(stats_DiscDMRsInRep), sep = "_")
stats_DiscDMRs <- merge(x = stats_DiscDMRsInDisc, y = stats_DiscDMRsInRep, by.x = "Disc_Region", by.y = "Rep_Region", all = FALSE, 
                        sort = FALSE)
stats_DiscDMRs <- stats_DiscDMRs[, !colnames(stats_DiscDMRs) %in% c("Rep_Variable", "Rep_Term")]
colnames(stats_DiscDMRs) <- str_replace_all(colnames(stats_DiscDMRs), 
                                            pattern = c("Disc_Region" = "Region", "Disc_Variable" = "Variable", 
                                                        "Disc_Term" = "Term"))
stats_DiscDMRsSig <- subset(stats_DiscDMRs, Disc_pvalue < 0.05 & Rep_pvalue < 0.05 & Disc_Estimate * Rep_Estimate > 0) # 15 / 584 replicated DMRs 
stats_DiscDMRsSig <- merge(x = methDisc[,c("DMRid", "chr", "start", "end", "width", "L", "area", "beta", "percentDifference",
                                           "stat", "pval", "qval")], 
                           y = stats_DiscDMRsSig, by.x = "DMRid", by.y = "Region", all = FALSE, sort = FALSE)
stats_DiscDMRsSig <- getDMRanno(DMRstats = stats_DiscDMRsSig, regDomains = regDomains, 
                                file = "Tables/Replicated Differential Methylation Discovery Males Diagnosis DMRs with Annotation.txt")
genes_DiscDMRsSig <- getDMRgeneList(DMRstats = stats_DiscDMRsSig, regDomains = regDomains, direction = "all", type = "gene_name")

# Replication Males DMRs ####
# Data
methDisc <- loadRegions("DMRs/Replication/Diagnosis Males 50/Dx_Replication50_males_DMR_raw_methylation_in_Discovery.csv", 
                        chroms = maleChroms, sort = TRUE, DMRid = TRUE)
methRep <- loadRegions("DMRs/Replication/Diagnosis Males 50/DMR_raw_methylation_Dx_Replication50_males.txt", 
                       chroms = maleChroms, sort = TRUE, DMRid = TRUE)
table(methDisc$start == methRep$start) # All TRUE

# Remove DMRs Covered in < 50% either TD or ASD in Discovery Samples
meth <- methDisc[,grepl("JLCM", colnames(methDisc), fixed = TRUE)]
missingTD <- apply(meth[,as.character(malesDisc$Sequencing_ID[malesDisc$Diagnosis_Alg == "TD"])], 
                   1, function(x) table(is.na(x))["TRUE"]) / 
        length(as.character(malesDisc$Sequencing_ID[malesDisc$Diagnosis_Alg == "TD"])) # percent of TD with no reads
missingTD[is.na(missingTD)] <- 0
missingASD <- apply(meth[,as.character(malesDisc$Sequencing_ID[malesDisc$Diagnosis_Alg == "ASD"])], 
                    1, function(x) table(is.na(x))["TRUE"]) /
        length(as.character(malesDisc$Sequencing_ID[malesDisc$Diagnosis_Alg == "ASD"])) # percent of ASD with no reads
missingASD[is.na(missingASD)] <- 0
methDisc <- subset(methDisc, missingTD <= 0.5 & missingASD <= 0.5) # 4648 / 4650 DMRs remaining
methRep <- subset(methRep, DMRid %in% methDisc$DMRid)

# Differential Methylation in Discovery Samples
meth <- t(methDisc[, grepl("JLCM", colnames(methDisc), fixed = TRUE)]) %>% as.data.frame
colnames(meth) <- methDisc$DMRid
meth$Sequencing_ID <- rownames(meth)
meth <- merge(x = malesDisc, y = meth, by = "Sequencing_ID", all = FALSE)
stats_RepDMRsInDisc <- DMRmethLm(methDisc$DMRid, catVars = "Diagnosis_Alg", contVars = NULL, sampleData = meth,
                                 file = "Tables/Differential Methylation Replication Males DMRs in Discovery Males.txt")
table(stats_RepDMRsInDisc$pvalue < 0.05) # 234 / 4648 sig DMRs

# Differential Methylation in Replication Samples
meth <- t(methRep[, grepl("JLCM", colnames(methRep), fixed = TRUE)]) %>% as.data.frame
colnames(meth) <- methRep$DMRid
meth$Sequencing_ID <- rownames(meth)
meth <- merge(x = malesRep, y = meth, by = "Sequencing_ID", all = FALSE)
stats_RepDMRsInRep <- DMRmethLm(methRep$DMRid, catVars = "Diagnosis_Alg", contVars = NULL, sampleData = meth,
                                file = "Tables/Differential Methylation Replication Males DMRs in Replication Males.txt")
table(stats_RepDMRsInRep$pvalue < 0.05) # 4310 / 4648 sig DMRs

# Subset by Differential Methylation in Both Sets and Annotate
colnames(stats_RepDMRsInDisc) <- paste("Disc", colnames(stats_RepDMRsInDisc), sep = "_")
colnames(stats_RepDMRsInRep) <- paste("Rep", colnames(stats_RepDMRsInRep), sep = "_")
stats_RepDMRs <- merge(x = stats_RepDMRsInDisc, y = stats_RepDMRsInRep, by.x = "Disc_Region", by.y = "Rep_Region", all = FALSE, 
                        sort = FALSE)
stats_RepDMRs <- stats_RepDMRs[, !colnames(stats_RepDMRs) %in% c("Rep_Variable", "Rep_Term")]
colnames(stats_RepDMRs) <- str_replace_all(colnames(stats_RepDMRs), 
                                            pattern = c("Disc_Region" = "Region", "Disc_Variable" = "Variable", 
                                                        "Disc_Term" = "Term"))
stats_RepDMRsSig <- subset(stats_RepDMRs, Disc_pvalue < 0.05 & Rep_pvalue < 0.05 & Disc_Estimate * Rep_Estimate > 0) # 139 / 4648 replicated DMRs
stats_RepDMRsSig <- merge(x = methRep[,c("DMRid", "chr", "start", "end", "width", "L", "area", "beta", "percentDifference",
                                         "stat", "pval", "qval")], 
                          y = stats_RepDMRsSig, by.x = "DMRid", by.y = "Region", all = FALSE, sort = FALSE)
stats_RepDMRsSig <- getDMRanno(DMRstats = stats_RepDMRsSig, regDomains = regDomains, 
                                file = "Tables/Replicated Differential Methylation Replication Males Diagnosis DMRs with Annotation.txt")
genes_RepDMRsSig <- getDMRgeneList(DMRstats = stats_RepDMRsSig, regDomains = regDomains, direction = "all", type = "gene_name")

# Females Differential Methylation Analysis by Diagnosis ---------------------------------------
# Sample Data ####
femalesDisc <- subset(samples, Sex == "F" & Platform == "HiSeqX10", select = c("Sequencing_ID", "Diagnosis_Alg", "Sex"))
femalesRep <- subset(samples, Sex == "F" & Platform == "HiSeq4000", select = c("Sequencing_ID", "Diagnosis_Alg", "Sex"))
femaleChroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrM")

# Discovery Females DMRs ####
# Data
methDisc <- loadRegions("DMRs/Discovery/Diagnosis Females 50/DMR_raw_methylation_Dx_Discovery50_females.txt", 
                        chroms = femaleChroms, sort = TRUE, DMRid = TRUE)
methRep <- loadRegions("DMRs/Discovery/Diagnosis Females 50/Dx_Discovery50_females_DMR_raw_methylation_in_Replication.csv", 
                       chroms = femaleChroms, sort = TRUE, DMRid = TRUE)
table(methDisc$start == methRep$start) # All TRUE

# Remove DMRs Covered in < 50% either TD or ASD in Replication Samples
meth <- methRep[,grepl("JLCM", colnames(methRep), fixed = TRUE)]
missingTD <- apply(meth[,as.character(femalesRep$Sequencing_ID[femalesRep$Diagnosis_Alg == "TD"])], 
                   1, function(x) table(is.na(x))["TRUE"]) / 
        length(as.character(femalesRep$Sequencing_ID[femalesRep$Diagnosis_Alg == "TD"])) # percent of TD with no reads
missingTD[is.na(missingTD)] <- 0
missingASD <- apply(meth[,as.character(femalesRep$Sequencing_ID[femalesRep$Diagnosis_Alg == "ASD"])], 
                    1, function(x) table(is.na(x))["TRUE"]) /
        length(as.character(femalesRep$Sequencing_ID[femalesRep$Diagnosis_Alg == "ASD"])) # percent of ASD with no reads
missingASD[is.na(missingASD)] <- 0
methRep <- subset(methRep, missingTD <= 0.5 & missingASD <= 0.5) # 1420 / 1952 DMRs remaining
methDisc <- subset(methDisc, DMRid %in% methRep$DMRid)

# Differential Methylation in Discovery Samples
meth <- t(methDisc[, grepl("JLCM", colnames(methDisc), fixed = TRUE)]) %>% as.data.frame
colnames(meth) <- methDisc$DMRid
meth$Sequencing_ID <- rownames(meth)
meth <- merge(x = femalesDisc, y = meth, by = "Sequencing_ID", all = FALSE)
stats_DiscDMRsInDisc <- DMRmethLm(methDisc$DMRid, catVars = "Diagnosis_Alg", contVars = NULL, sampleData = meth,
                                  file = "Tables/Differential Methylation Discovery Females DMRs in Discovery Females.txt")
table(stats_DiscDMRsInDisc$pvalue < 0.05) # 1220 / 1420 sig DMRs

# Differential Methylation in Replication Samples
meth <- t(methRep[, grepl("JLCM", colnames(methRep), fixed = TRUE)]) %>% as.data.frame
colnames(meth) <- methRep$DMRid
meth$Sequencing_ID <- rownames(meth)
meth <- merge(x = femalesRep, y = meth, by = "Sequencing_ID", all = FALSE)
stats_DiscDMRsInRep <- DMRmethLm(methRep$DMRid, catVars = "Diagnosis_Alg", contVars = NULL, sampleData = meth,
                                 file = "Tables/Differential Methylation Discovery Females DMRs in Replication Females.txt")
table(stats_DiscDMRsInRep$pvalue < 0.05) # 53 / 1420 sig DMRs

# Subset by Differential Methylation in Both Sets and Annotate
colnames(stats_DiscDMRsInDisc) <- paste("Disc", colnames(stats_DiscDMRsInDisc), sep = "_")
colnames(stats_DiscDMRsInRep) <- paste("Rep", colnames(stats_DiscDMRsInRep), sep = "_")
stats_DiscDMRs <- merge(x = stats_DiscDMRsInDisc, y = stats_DiscDMRsInRep, by.x = "Disc_Region", by.y = "Rep_Region", all = FALSE, 
                        sort = FALSE)
stats_DiscDMRs <- stats_DiscDMRs[, !colnames(stats_DiscDMRs) %in% c("Rep_Variable", "Rep_Term")]
colnames(stats_DiscDMRs) <- str_replace_all(colnames(stats_DiscDMRs), 
                                            pattern = c("Disc_Region" = "Region", "Disc_Variable" = "Variable", 
                                                        "Disc_Term" = "Term"))
stats_DiscDMRsSig <- subset(stats_DiscDMRs, Disc_pvalue < 0.05 & Rep_pvalue < 0.05 & Disc_Estimate * Rep_Estimate > 0) # 23 / 1420 replicated DMRs 
stats_DiscDMRsSig <- merge(x = methDisc[,c("DMRid", "chr", "start", "end", "width", "L", "area", "beta", "percentDifference",
                                           "stat", "pval", "qval")], 
                           y = stats_DiscDMRsSig, by.x = "DMRid", by.y = "Region", all = FALSE, sort = FALSE)
stats_DiscDMRsSig <- getDMRanno(DMRstats = stats_DiscDMRsSig, regDomains = regDomains, 
                                file = "Tables/Replicated Differential Methylation Discovery Females Diagnosis DMRs with Annotation.txt")
genes_DiscDMRsSig <- getDMRgeneList(DMRstats = stats_DiscDMRsSig, regDomains = regDomains, direction = "all", type = "gene_name")

# Replication Females DMRs ####
# Data
methDisc <- loadRegions("DMRs/Replication/Diagnosis Females 100/Dx_Replication100_females_DMR_raw_methylation_in_Discovery.csv", 
                        chroms = femaleChroms, sort = TRUE, DMRid = TRUE)
methRep <- loadRegions("DMRs/Replication/Diagnosis Females 100/DMR_raw_methylation_Dx_Replication100_females.txt", 
                       chroms = femaleChroms, sort = TRUE, DMRid = TRUE)
table(methDisc$start == methRep$start) # All TRUE

# Remove DMRs Covered in < 50% either TD or ASD in Discovery Samples
meth <- methDisc[,grepl("JLCM", colnames(methDisc), fixed = TRUE)]
missingTD <- apply(meth[,as.character(femalesDisc$Sequencing_ID[femalesDisc$Diagnosis_Alg == "TD"])], 
                   1, function(x) table(is.na(x))["TRUE"]) / 
        length(as.character(femalesDisc$Sequencing_ID[femalesDisc$Diagnosis_Alg == "TD"])) # percent of TD with no reads
missingTD[is.na(missingTD)] <- 0
missingASD <- apply(meth[,as.character(femalesDisc$Sequencing_ID[femalesDisc$Diagnosis_Alg == "ASD"])], 
                    1, function(x) table(is.na(x))["TRUE"]) /
        length(as.character(femalesDisc$Sequencing_ID[femalesDisc$Diagnosis_Alg == "ASD"])) # percent of ASD with no reads
missingASD[is.na(missingASD)] <- 0
methDisc <- subset(methDisc, missingTD <= 0.5 & missingASD <= 0.5) # 8726 / 8728 DMRs remaining
methRep <- subset(methRep, DMRid %in% methDisc$DMRid)

# Differential Methylation in Discovery Samples
meth <- t(methDisc[, grepl("JLCM", colnames(methDisc), fixed = TRUE)]) %>% as.data.frame
colnames(meth) <- methDisc$DMRid
meth$Sequencing_ID <- rownames(meth)
meth <- merge(x = femalesDisc, y = meth, by = "Sequencing_ID", all = FALSE)
stats_RepDMRsInDisc <- DMRmethLm(methDisc$DMRid, catVars = "Diagnosis_Alg", contVars = NULL, sampleData = meth,
                                 file = "Tables/Differential Methylation Replication Females DMRs in Discovery Females.txt")
table(stats_RepDMRsInDisc$pvalue < 0.05) # 354 / 8726 sig DMRs

# Differential Methylation in Replication Samples
meth <- t(methRep[, grepl("JLCM", colnames(methRep), fixed = TRUE)]) %>% as.data.frame
colnames(meth) <- methRep$DMRid
meth$Sequencing_ID <- rownames(meth)
meth <- merge(x = femalesRep, y = meth, by = "Sequencing_ID", all = FALSE)
stats_RepDMRsInRep <- DMRmethLm(methRep$DMRid, catVars = "Diagnosis_Alg", contVars = NULL, sampleData = meth,
                                file = "Tables/Differential Methylation Replication Females DMRs in Replication Females.txt")
table(stats_RepDMRsInRep$pvalue < 0.05) # 6309 / 8726 sig DMRs

# Subset by Differential Methylation in Both Sets and Annotate
colnames(stats_RepDMRsInDisc) <- paste("Disc", colnames(stats_RepDMRsInDisc), sep = "_")
colnames(stats_RepDMRsInRep) <- paste("Rep", colnames(stats_RepDMRsInRep), sep = "_")
stats_RepDMRs <- merge(x = stats_RepDMRsInDisc, y = stats_RepDMRsInRep, by.x = "Disc_Region", by.y = "Rep_Region", all = FALSE, 
                       sort = FALSE)
stats_RepDMRs <- stats_RepDMRs[, !colnames(stats_RepDMRs) %in% c("Rep_Variable", "Rep_Term")]
colnames(stats_RepDMRs) <- str_replace_all(colnames(stats_RepDMRs), 
                                           pattern = c("Disc_Region" = "Region", "Disc_Variable" = "Variable", 
                                                       "Disc_Term" = "Term"))
stats_RepDMRsSig <- subset(stats_RepDMRs, Disc_pvalue < 0.05 & Rep_pvalue < 0.05 & Disc_Estimate * Rep_Estimate > 0) # 122 / 8726 replicated DMRs
stats_RepDMRsSig <- merge(x = methRep[,c("DMRid", "chr", "start", "end", "width", "L", "area", "beta", "percentDifference",
                                         "stat", "pval", "qval")], 
                          y = stats_RepDMRsSig, by.x = "DMRid", by.y = "Region", all = FALSE, sort = FALSE)
stats_RepDMRsSig <- getDMRanno(DMRstats = stats_RepDMRsSig, regDomains = regDomains, 
                               file = "Tables/Replicated Differential Methylation Replication Females Diagnosis DMRs with Annotation.txt")
genes_RepDMRsSig <- getDMRgeneList(DMRstats = stats_RepDMRsSig, regDomains = regDomains, direction = "all", type = "gene_name")

# Males Differential Methylation Analysis by ADOS ---------------------------------------
# Sample Data ####
malesDisc <- subset(samples, Sex == "M" & Platform == "HiSeqX10", select = c("Sequencing_ID", "Diagnosis_Alg", "Sex", "ADOScs"))
malesRep <- subset(samples, Sex == "M" & Platform == "HiSeq4000", select = c("Sequencing_ID", "Diagnosis_Alg", "Sex", "ADOScs"))
maleChroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")

# Discovery Males DMRs ####
# Data
methDisc <- loadRegions("DMRs/Discovery/Diagnosis Males 50/DMR_raw_methylation_Dx_Discovery50_males.txt", 
                        chroms = maleChroms, sort = TRUE, DMRid = TRUE)
methRep <- loadRegions("DMRs/Discovery/Diagnosis Males 50/Dx_Discovery50_males_DMR_raw_methylation_in_Replication.csv", 
                       chroms = maleChroms, sort = TRUE, DMRid = TRUE)
table(methDisc$start == methRep$start) # All TRUE

# Remove DMRs Covered in < 50% either TD or ASD in Replication Samples
meth <- methRep[,grepl("JLCM", colnames(methRep), fixed = TRUE)]
missingTD <- apply(meth[,as.character(malesRep$Sequencing_ID[malesRep$Diagnosis_Alg == "TD"])], 
                   1, function(x) table(is.na(x))["TRUE"]) / 
        length(as.character(malesRep$Sequencing_ID[malesRep$Diagnosis_Alg == "TD"])) # percent of TD with no reads
missingTD[is.na(missingTD)] <- 0
missingASD <- apply(meth[,as.character(malesRep$Sequencing_ID[malesRep$Diagnosis_Alg == "ASD"])], 
                    1, function(x) table(is.na(x))["TRUE"]) /
        length(as.character(malesRep$Sequencing_ID[malesRep$Diagnosis_Alg == "ASD"])) # percent of ASD with no reads
missingASD[is.na(missingASD)] <- 0
methRep <- subset(methRep, missingTD <= 0.5 & missingASD <= 0.5) # 584 / 635 DMRs remaining
methDisc <- subset(methDisc, DMRid %in% methRep$DMRid)

# Differential Methylation in Discovery Samples
meth <- t(methDisc[, grepl("JLCM", colnames(methDisc), fixed = TRUE)]) %>% as.data.frame
colnames(meth) <- methDisc$DMRid
meth$Sequencing_ID <- rownames(meth)
meth <- merge(x = malesDisc, y = meth, by = "Sequencing_ID", all = FALSE)
stats_DiscDMRsInDisc <- DMRmethLm(methDisc$DMRid, catVars = NULL, contVars = "ADOScs", sampleData = meth,
                                  file = "Tables/Differential Methylation by ADOS Discovery Males DMRs in Discovery Males.txt")
table(stats_DiscDMRsInDisc$pvalue < 0.05) # 516 / 584 sig DMRs

# Differential Methylation in Replication Samples
meth <- t(methRep[, grepl("JLCM", colnames(methRep), fixed = TRUE)]) %>% as.data.frame
colnames(meth) <- methRep$DMRid
meth$Sequencing_ID <- rownames(meth)
meth <- merge(x = malesRep, y = meth, by = "Sequencing_ID", all = FALSE)
stats_DiscDMRsInRep <- DMRmethLm(methRep$DMRid, catVars = NULL, contVars = "ADOScs", sampleData = meth,
                                 file = "Tables/Differential Methylation by ADOS Discovery Males DMRs in Replication Males.txt")
table(stats_DiscDMRsInRep$pvalue < 0.05) # 26 / 584 sig DMRs

# Subset by Differential Methylation in Both Sets and Annotate
colnames(stats_DiscDMRsInDisc) <- paste("Disc", colnames(stats_DiscDMRsInDisc), sep = "_")
colnames(stats_DiscDMRsInRep) <- paste("Rep", colnames(stats_DiscDMRsInRep), sep = "_")
stats_DiscDMRs <- merge(x = stats_DiscDMRsInDisc, y = stats_DiscDMRsInRep, by.x = "Disc_Region", by.y = "Rep_Region", all = FALSE, 
                        sort = FALSE)
stats_DiscDMRs <- stats_DiscDMRs[, !colnames(stats_DiscDMRs) %in% c("Rep_Variable", "Rep_Term")]
colnames(stats_DiscDMRs) <- str_replace_all(colnames(stats_DiscDMRs), 
                                            pattern = c("Disc_Region" = "Region", "Disc_Variable" = "Variable", 
                                                        "Disc_Term" = "Term"))
stats_DiscDMRsSig <- subset(stats_DiscDMRs, Disc_pvalue < 0.05 & Rep_pvalue < 0.05 & Disc_Estimate * Rep_Estimate > 0) # 15 / 584 replicated DMRs 
stats_DiscDMRsSig <- merge(x = methDisc[,c("DMRid", "chr", "start", "end", "width", "L", "area", "beta", "percentDifference",
                                           "stat", "pval", "qval")], 
                           y = stats_DiscDMRsSig, by.x = "DMRid", by.y = "Region", all = FALSE, sort = FALSE)
stats_DiscDMRsSig <- getDMRanno(DMRstats = stats_DiscDMRsSig, regDomains = regDomains, 
                                file = "Tables/Replicated Differential Methylation by ADOS Discovery Males Diagnosis DMRs with Annotation.txt")
genes_DiscDMRsSig <- getDMRgeneList(DMRstats = stats_DiscDMRsSig, regDomains = regDomains, direction = "all", type = "gene_name")

# Replication Males DMRs ####
# Data
methDisc <- loadRegions("DMRs/Replication/Diagnosis Males 50/Dx_Replication50_males_DMR_raw_methylation_in_Discovery.csv", 
                        chroms = maleChroms, sort = TRUE, DMRid = TRUE)
methRep <- loadRegions("DMRs/Replication/Diagnosis Males 50/DMR_raw_methylation_Dx_Replication50_males.txt", 
                       chroms = maleChroms, sort = TRUE, DMRid = TRUE)
table(methDisc$start == methRep$start) # All TRUE

# Remove DMRs Covered in < 50% either TD or ASD in Discovery Samples
meth <- methDisc[,grepl("JLCM", colnames(methDisc), fixed = TRUE)]
missingTD <- apply(meth[,as.character(malesDisc$Sequencing_ID[malesDisc$Diagnosis_Alg == "TD"])], 
                   1, function(x) table(is.na(x))["TRUE"]) / 
        length(as.character(malesDisc$Sequencing_ID[malesDisc$Diagnosis_Alg == "TD"])) # percent of TD with no reads
missingTD[is.na(missingTD)] <- 0
missingASD <- apply(meth[,as.character(malesDisc$Sequencing_ID[malesDisc$Diagnosis_Alg == "ASD"])], 
                    1, function(x) table(is.na(x))["TRUE"]) /
        length(as.character(malesDisc$Sequencing_ID[malesDisc$Diagnosis_Alg == "ASD"])) # percent of ASD with no reads
missingASD[is.na(missingASD)] <- 0
methDisc <- subset(methDisc, missingTD <= 0.5 & missingASD <= 0.5) # 4648 / 4650 DMRs remaining
methRep <- subset(methRep, DMRid %in% methDisc$DMRid)

# Differential Methylation in Discovery Samples
meth <- t(methDisc[, grepl("JLCM", colnames(methDisc), fixed = TRUE)]) %>% as.data.frame
colnames(meth) <- methDisc$DMRid
meth$Sequencing_ID <- rownames(meth)
meth <- merge(x = malesDisc, y = meth, by = "Sequencing_ID", all = FALSE)
stats_RepDMRsInDisc <- DMRmethLm(methDisc$DMRid, catVars = NULL, contVars = "ADOScs", sampleData = meth,
                                 file = "Tables/Differential Methylation by ADOS Replication Males DMRs in Discovery Males.txt")
table(stats_RepDMRsInDisc$pvalue < 0.05) # 223 / 4648 sig DMRs

# Differential Methylation in Replication Samples
meth <- t(methRep[, grepl("JLCM", colnames(methRep), fixed = TRUE)]) %>% as.data.frame
colnames(meth) <- methRep$DMRid
meth$Sequencing_ID <- rownames(meth)
meth <- merge(x = malesRep, y = meth, by = "Sequencing_ID", all = FALSE)
stats_RepDMRsInRep <- DMRmethLm(methRep$DMRid, catVars = NULL, contVars = "ADOScs", sampleData = meth,
                                file = "Tables/Differential Methylation by ADOS Replication Males DMRs in Replication Males.txt")
table(stats_RepDMRsInRep$pvalue < 0.05) # 3549 / 4648 sig DMRs

# Subset by Differential Methylation in Both Sets and Annotate
colnames(stats_RepDMRsInDisc) <- paste("Disc", colnames(stats_RepDMRsInDisc), sep = "_")
colnames(stats_RepDMRsInRep) <- paste("Rep", colnames(stats_RepDMRsInRep), sep = "_")
stats_RepDMRs <- merge(x = stats_RepDMRsInDisc, y = stats_RepDMRsInRep, by.x = "Disc_Region", by.y = "Rep_Region", all = FALSE, 
                       sort = FALSE)
stats_RepDMRs <- stats_RepDMRs[, !colnames(stats_RepDMRs) %in% c("Rep_Variable", "Rep_Term")]
colnames(stats_RepDMRs) <- str_replace_all(colnames(stats_RepDMRs), 
                                           pattern = c("Disc_Region" = "Region", "Disc_Variable" = "Variable", 
                                                       "Disc_Term" = "Term"))
stats_RepDMRsSig <- subset(stats_RepDMRs, Disc_pvalue < 0.05 & Rep_pvalue < 0.05 & Disc_Estimate * Rep_Estimate > 0) # 119 / 4648 replicated DMRs
stats_RepDMRsSig <- merge(x = methRep[,c("DMRid", "chr", "start", "end", "width", "L", "area", "beta", "percentDifference",
                                         "stat", "pval", "qval")], 
                          y = stats_RepDMRsSig, by.x = "DMRid", by.y = "Region", all = FALSE, sort = FALSE)
stats_RepDMRsSig <- getDMRanno(DMRstats = stats_RepDMRsSig, regDomains = regDomains, 
                               file = "Tables/Replicated Differential Methylation by ADOS Replication Males Diagnosis DMRs with Annotation.txt")
genes_RepDMRsSig <- getDMRgeneList(DMRstats = stats_RepDMRsSig, regDomains = regDomains, direction = "all", type = "gene_name")

# Females Differential Methylation Analysis by ADOS ---------------------------------------
# Sample Data ####
femalesDisc <- subset(samples, Sex == "F" & Platform == "HiSeqX10", select = c("Sequencing_ID", "Diagnosis_Alg", "Sex", "ADOScs"))
femalesRep <- subset(samples, Sex == "F" & Platform == "HiSeq4000", select = c("Sequencing_ID", "Diagnosis_Alg", "Sex", "ADOScs"))
femaleChroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrM")

# Discovery Females DMRs ####
# Data
methDisc <- loadRegions("DMRs/Discovery/Diagnosis Females 50/DMR_raw_methylation_Dx_Discovery50_females.txt", 
                        chroms = femaleChroms, sort = TRUE, DMRid = TRUE)
methRep <- loadRegions("DMRs/Discovery/Diagnosis Females 50/Dx_Discovery50_females_DMR_raw_methylation_in_Replication.csv", 
                       chroms = femaleChroms, sort = TRUE, DMRid = TRUE)
table(methDisc$start == methRep$start) # All TRUE

# Remove DMRs Covered in < 50% either TD or ASD in Replication Samples
meth <- methRep[,grepl("JLCM", colnames(methRep), fixed = TRUE)]
missingTD <- apply(meth[,as.character(femalesRep$Sequencing_ID[femalesRep$Diagnosis_Alg == "TD"])], 
                   1, function(x) table(is.na(x))["TRUE"]) / 
        length(as.character(femalesRep$Sequencing_ID[femalesRep$Diagnosis_Alg == "TD"])) # percent of TD with no reads
missingTD[is.na(missingTD)] <- 0
missingASD <- apply(meth[,as.character(femalesRep$Sequencing_ID[femalesRep$Diagnosis_Alg == "ASD"])], 
                    1, function(x) table(is.na(x))["TRUE"]) /
        length(as.character(femalesRep$Sequencing_ID[femalesRep$Diagnosis_Alg == "ASD"])) # percent of ASD with no reads
missingASD[is.na(missingASD)] <- 0
methRep <- subset(methRep, missingTD <= 0.5 & missingASD <= 0.5) # 1420 / 1952 DMRs remaining
methDisc <- subset(methDisc, DMRid %in% methRep$DMRid)

# Differential Methylation in Discovery Samples
meth <- t(methDisc[, grepl("JLCM", colnames(methDisc), fixed = TRUE)]) %>% as.data.frame
colnames(meth) <- methDisc$DMRid
meth$Sequencing_ID <- rownames(meth)
meth <- merge(x = femalesDisc, y = meth, by = "Sequencing_ID", all = FALSE)
stats_DiscDMRsInDisc <- DMRmethLm(methDisc$DMRid, catVars = NULL, contVars = "ADOScs", sampleData = meth,
                                  file = "Tables/Differential Methylation by ADOS Discovery Females DMRs in Discovery Females.txt")
table(stats_DiscDMRsInDisc$pvalue < 0.05) # 1044 / 1420 sig DMRs

# Differential Methylation in Replication Samples
meth <- t(methRep[, grepl("JLCM", colnames(methRep), fixed = TRUE)]) %>% as.data.frame
colnames(meth) <- methRep$DMRid
meth$Sequencing_ID <- rownames(meth)
meth <- merge(x = femalesRep, y = meth, by = "Sequencing_ID", all = FALSE)
stats_DiscDMRsInRep <- DMRmethLm(methRep$DMRid, catVars = NULL, contVars = "ADOScs", sampleData = meth,
                                 file = "Tables/Differential Methylation by ADOS Discovery Females DMRs in Replication Females.txt")
table(stats_DiscDMRsInRep$pvalue < 0.05) # 78 / 1420 sig DMRs

# Subset by Differential Methylation in Both Sets and Annotate
colnames(stats_DiscDMRsInDisc) <- paste("Disc", colnames(stats_DiscDMRsInDisc), sep = "_")
colnames(stats_DiscDMRsInRep) <- paste("Rep", colnames(stats_DiscDMRsInRep), sep = "_")
stats_DiscDMRs <- merge(x = stats_DiscDMRsInDisc, y = stats_DiscDMRsInRep, by.x = "Disc_Region", by.y = "Rep_Region", all = FALSE, 
                        sort = FALSE)
stats_DiscDMRs <- stats_DiscDMRs[, !colnames(stats_DiscDMRs) %in% c("Rep_Variable", "Rep_Term")]
colnames(stats_DiscDMRs) <- str_replace_all(colnames(stats_DiscDMRs), 
                                            pattern = c("Disc_Region" = "Region", "Disc_Variable" = "Variable", 
                                                        "Disc_Term" = "Term"))
stats_DiscDMRsSig <- subset(stats_DiscDMRs, Disc_pvalue < 0.05 & Rep_pvalue < 0.05 & Disc_Estimate * Rep_Estimate > 0) # 23 / 1420 replicated DMRs 
stats_DiscDMRsSig <- merge(x = methDisc[,c("DMRid", "chr", "start", "end", "width", "L", "area", "beta", "percentDifference",
                                           "stat", "pval", "qval")], 
                           y = stats_DiscDMRsSig, by.x = "DMRid", by.y = "Region", all = FALSE, sort = FALSE)
stats_DiscDMRsSig <- getDMRanno(DMRstats = stats_DiscDMRsSig, regDomains = regDomains, 
                                file = "Tables/Replicated Differential Methylation by ADOS Discovery Females Diagnosis DMRs with Annotation.txt")
genes_DiscDMRsSig <- getDMRgeneList(DMRstats = stats_DiscDMRsSig, regDomains = regDomains, direction = "all", type = "gene_name")

# Replication Females DMRs ####
# Data
methDisc <- loadRegions("DMRs/Replication/Diagnosis Females 100/Dx_Replication100_females_DMR_raw_methylation_in_Discovery.csv", 
                        chroms = femaleChroms, sort = TRUE, DMRid = TRUE)
methRep <- loadRegions("DMRs/Replication/Diagnosis Females 100/DMR_raw_methylation_Dx_Replication100_females.txt", 
                       chroms = femaleChroms, sort = TRUE, DMRid = TRUE)
table(methDisc$start == methRep$start) # All TRUE

# Remove DMRs Covered in < 50% either TD or ASD in Discovery Samples
meth <- methDisc[,grepl("JLCM", colnames(methDisc), fixed = TRUE)]
missingTD <- apply(meth[,as.character(femalesDisc$Sequencing_ID[femalesDisc$Diagnosis_Alg == "TD"])], 
                   1, function(x) table(is.na(x))["TRUE"]) / 
        length(as.character(femalesDisc$Sequencing_ID[femalesDisc$Diagnosis_Alg == "TD"])) # percent of TD with no reads
missingTD[is.na(missingTD)] <- 0
missingASD <- apply(meth[,as.character(femalesDisc$Sequencing_ID[femalesDisc$Diagnosis_Alg == "ASD"])], 
                    1, function(x) table(is.na(x))["TRUE"]) /
        length(as.character(femalesDisc$Sequencing_ID[femalesDisc$Diagnosis_Alg == "ASD"])) # percent of ASD with no reads
missingASD[is.na(missingASD)] <- 0
methDisc <- subset(methDisc, missingTD <= 0.5 & missingASD <= 0.5) # 8726 / 8728 DMRs remaining
methRep <- subset(methRep, DMRid %in% methDisc$DMRid)

# Differential Methylation in Discovery Samples
meth <- t(methDisc[, grepl("JLCM", colnames(methDisc), fixed = TRUE)]) %>% as.data.frame
colnames(meth) <- methDisc$DMRid
meth$Sequencing_ID <- rownames(meth)
meth <- merge(x = femalesDisc, y = meth, by = "Sequencing_ID", all = FALSE)
stats_RepDMRsInDisc <- DMRmethLm(methDisc$DMRid, catVars = NULL, contVars = "ADOScs", sampleData = meth,
                                 file = "Tables/Differential Methylation by ADOS Replication Females DMRs in Discovery Females.txt")
table(stats_RepDMRsInDisc$pvalue < 0.05) # 446 / 8726 sig DMRs

# Differential Methylation in Replication Samples
meth <- t(methRep[, grepl("JLCM", colnames(methRep), fixed = TRUE)]) %>% as.data.frame
colnames(meth) <- methRep$DMRid
meth$Sequencing_ID <- rownames(meth)
meth <- merge(x = femalesRep, y = meth, by = "Sequencing_ID", all = FALSE)
stats_RepDMRsInRep <- DMRmethLm(methRep$DMRid, catVars = NULL, contVars = "ADOScs", sampleData = meth,
                                file = "Tables/Differential Methylation by ADOS Replication Females DMRs in Replication Females.txt")
table(stats_RepDMRsInRep$pvalue < 0.05) # 4463 / 8726 sig DMRs

# Subset by Differential Methylation in Both Sets and Annotate
colnames(stats_RepDMRsInDisc) <- paste("Disc", colnames(stats_RepDMRsInDisc), sep = "_")
colnames(stats_RepDMRsInRep) <- paste("Rep", colnames(stats_RepDMRsInRep), sep = "_")
stats_RepDMRs <- merge(x = stats_RepDMRsInDisc, y = stats_RepDMRsInRep, by.x = "Disc_Region", by.y = "Rep_Region", all = FALSE, 
                       sort = FALSE)
stats_RepDMRs <- stats_RepDMRs[, !colnames(stats_RepDMRs) %in% c("Rep_Variable", "Rep_Term")]
colnames(stats_RepDMRs) <- str_replace_all(colnames(stats_RepDMRs), 
                                           pattern = c("Disc_Region" = "Region", "Disc_Variable" = "Variable", 
                                                       "Disc_Term" = "Term"))
stats_RepDMRsSig <- subset(stats_RepDMRs, Disc_pvalue < 0.05 & Rep_pvalue < 0.05 & Disc_Estimate * Rep_Estimate > 0) # 97 / 8726 replicated DMRs
stats_RepDMRsSig <- merge(x = methRep[,c("DMRid", "chr", "start", "end", "width", "L", "area", "beta", "percentDifference",
                                         "stat", "pval", "qval")], 
                          y = stats_RepDMRsSig, by.x = "DMRid", by.y = "Region", all = FALSE, sort = FALSE)
stats_RepDMRsSig <- getDMRanno(DMRstats = stats_RepDMRsSig, regDomains = regDomains, 
                               file = "Tables/Replicated Differential Methylation by ADOS Replication Females Diagnosis DMRs with Annotation.txt")
genes_RepDMRsSig <- getDMRgeneList(DMRstats = stats_RepDMRsSig, regDomains = regDomains, direction = "all", type = "gene_name")

# ASD Only Males Differential Methylation Analysis by ADOS ---------------------------------------
# Sample Data ####
malesDisc <- subset(samples, Sex == "M" & Platform == "HiSeqX10" & Diagnosis_Alg == "ASD", 
                    select = c("Sequencing_ID", "Diagnosis_Alg", "Sex", "ADOScs"))
malesRep <- subset(samples, Sex == "M" & Platform == "HiSeq4000" & Diagnosis_Alg == "ASD", 
                   select = c("Sequencing_ID", "Diagnosis_Alg", "Sex", "ADOScs"))
maleChroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")

# Discovery Males DMRs ####
# Data
methDisc <- loadRegions("DMRs/Discovery/Diagnosis Males 50/DMR_raw_methylation_Dx_Discovery50_males.txt", 
                        chroms = maleChroms, sort = TRUE, DMRid = TRUE)
methRep <- loadRegions("DMRs/Discovery/Diagnosis Males 50/Dx_Discovery50_males_DMR_raw_methylation_in_Replication.csv", 
                       chroms = maleChroms, sort = TRUE, DMRid = TRUE)
table(methDisc$start == methRep$start) # All TRUE

# Remove DMRs Covered in < 50% ASD in Replication Samples
meth <- methRep[,as.character(malesRep$Sequencing_ID)]
missingASD <- apply(meth, 1, function(x) table(is.na(x))["TRUE"]) / length(as.character(malesRep$Sequencing_ID)) # percent of ASD with no reads
missingASD[is.na(missingASD)] <- 0
methRep <- subset(methRep, missingASD <= 0.5) # 584 / 635 DMRs remaining
methDisc <- subset(methDisc, DMRid %in% methRep$DMRid)

# Differential Methylation in Discovery Samples
meth <- t(methDisc[, as.character(malesDisc$Sequencing_ID)]) %>% as.data.frame
colnames(meth) <- methDisc$DMRid
meth$Sequencing_ID <- rownames(meth)
meth <- merge(x = malesDisc, y = meth, by = "Sequencing_ID", all = FALSE)
stats_DiscDMRsInDisc <- DMRmethLm(methDisc$DMRid, catVars = NULL, contVars = "ADOScs", sampleData = meth,
                                  file = "Tables/Differential Methylation by ADOS Discovery Males DMRs in Discovery Males ASD Only.txt")
table(stats_DiscDMRsInDisc$pvalue < 0.05) # 24 / 584 sig DMRs
sd(meth$ADOScs) #1.635735

# Differential Methylation in Replication Samples
meth <- t(methRep[, as.character(malesRep$Sequencing_ID)]) %>% as.data.frame
colnames(meth) <- methRep$DMRid
meth$Sequencing_ID <- rownames(meth)
meth <- merge(x = malesRep, y = meth, by = "Sequencing_ID", all = FALSE)
stats_DiscDMRsInRep <- DMRmethLm(methRep$DMRid, catVars = NULL, contVars = "ADOScs", sampleData = meth,
                                 file = "Tables/Differential Methylation by ADOS Discovery Males DMRs in Replication Males ASD Only.txt")
table(stats_DiscDMRsInRep$pvalue < 0.05) # 38 / 584 sig DMRs

# Subset by Differential Methylation in Both Sets and Annotate
colnames(stats_DiscDMRsInDisc) <- paste("Disc", colnames(stats_DiscDMRsInDisc), sep = "_")
colnames(stats_DiscDMRsInRep) <- paste("Rep", colnames(stats_DiscDMRsInRep), sep = "_")
stats_DiscDMRs <- merge(x = stats_DiscDMRsInDisc, y = stats_DiscDMRsInRep, by.x = "Disc_Region", by.y = "Rep_Region", all = FALSE, 
                        sort = FALSE)
stats_DiscDMRs <- stats_DiscDMRs[, !colnames(stats_DiscDMRs) %in% c("Rep_Variable", "Rep_Term")]
colnames(stats_DiscDMRs) <- str_replace_all(colnames(stats_DiscDMRs), 
                                            pattern = c("Disc_Region" = "Region", "Disc_Variable" = "Variable", 
                                                        "Disc_Term" = "Term"))
stats_DiscDMRsSig <- subset(stats_DiscDMRs, Disc_pvalue < 0.05 & Rep_pvalue < 0.05 & Disc_Estimate * Rep_Estimate > 0) # 1 / 584 replicated DMRs 
stats_DiscDMRsSig <- merge(x = methDisc[,c("DMRid", "chr", "start", "end", "width", "L", "area", "beta", "percentDifference",
                                           "stat", "pval", "qval")], 
                           y = stats_DiscDMRsSig, by.x = "DMRid", by.y = "Region", all = FALSE, sort = FALSE)
stats_DiscDMRsSig <- getDMRanno(DMRstats = stats_DiscDMRsSig, regDomains = regDomains, 
                                file = "Tables/Replicated Differential Methylation by ADOS in ASD Only Discovery Males Diagnosis DMRs with Annotation.txt")

# Replicated DMRs Scatterplots
stats_DiscDMRsSig$DMRid #DMR_437
DMR <- "DMR_437"
methPlot <- data.frame(Sequencing_ID = c(as.character(malesDisc$Sequencing_ID), as.character(malesRep$Sequencing_ID)),
                       SampleSet = c(rep("Discovery", nrow(malesDisc)), rep("Replication", nrow(malesRep))),
                       ADOScs = c(malesDisc$ADOScs, malesRep$ADOScs), DMRid = DMR, 
                       Methylation = unlist(c(methDisc[methDisc$DMRid == DMR, as.character(malesDisc$Sequencing_ID)],
                                              methRep[methRep$DMRid == DMR, as.character(malesRep$Sequencing_ID)]))*100)
g <- ggplot(data = methPlot)
g + 
        geom_smooth(aes(x = ADOScs, y = Methylation), method = "lm") +
        geom_point(aes(x = ADOScs, y = Methylation), color = "#FF3366", size = 2) +
        facet_grid(cols = vars(SampleSet)) +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = "None", panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black", size = 18), legend.background = element_blank(), 
              plot.margin = unit(c(0.5,1,1,1), "lines"), strip.background = element_blank(), 
              strip.text = element_text(size = 20), axis.title = element_text(size = 20)) +
        coord_cartesian(ylim = c(45, 100)) +
        scale_x_continuous(breaks = pretty_breaks(n = 5)) +
        scale_y_continuous(breaks = pretty_breaks(n = 5)) +
        xlab("ADOS Comparison Score") +
        ylab("DMR_437 Methylation (%)") +
ggsave("Figures/DMR_437 Percent Methylation Plot Replicated ADOS in ASD Only Males.png", dpi = 600, width = 8, height = 5, 
       units = "in")

# DMR Methylation ~ ADOS Volcano Plots
volcanoPlot <- data.frame(DMRid = stats_DiscDMRs$Region,
                          SampleSet = rep(c("Discovery", "Replication"), each = nrow(stats_DiscDMRs)),
                          Estimate = c(stats_DiscDMRs$Disc_Estimate * sd(malesDisc$ADOScs, na.rm = TRUE),
                                       stats_DiscDMRs$Rep_Estimate * sd(malesRep$ADOScs, na.rm = TRUE)) * 100,
                          log_pvalue = c(stats_DiscDMRs$Disc_log_pvalue, stats_DiscDMRs$Rep_log_pvalue))
volcanoPlot$Significant <- (volcanoPlot$log_pvalue > -log10(0.05)) %>% factor(levels = c("TRUE", "FALSE"))
## Estimate is now in units of change in DMR methylation per 1 ADOScs point
g <- ggplot(data = volcanoPlot)
g + 
        geom_point(data = volcanoPlot, aes(x = Estimate, y = log_pvalue, color = Significant), size = 1) +
        facet_grid(cols = vars(SampleSet), scales = "free") +
        theme_bw(base_size = 25) +
        theme(legend.direction = "vertical", legend.position = c(1.12, 0.88), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_text(size = 20),
              axis.text = element_text(color = "black", size = 18), legend.background = element_blank(), 
              plot.margin = unit(c(0.5,7,1,1), "lines"), strip.background = element_blank(),
              strip.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 18)) +
        scale_x_continuous(breaks = pretty_breaks(n = 4)) +
        scale_y_continuous(breaks = pretty_breaks(n = 4), expand = expand_scale(mult = c(0.01, 0.03))) +
        scale_color_manual(name = expression(italic(p)*"< 0.05"), breaks = c("TRUE", "FALSE"), values = c("#FF3366", "#3366CC")) +
        xlab("DMR Methylation / ADOS Comparison Score") +
        ylab(expression("-log"[10]*"("*italic(p)*"-value)")) +
        ggsave("Figures/Discovery DMR Methylation by ADOS in ASD Only Males Volcano Plot.png", dpi = 600, width = 8, 
               height = 5, units = "in")

# Replication Males DMRs ####
# Data
methDisc <- loadRegions("DMRs/Replication/Diagnosis Males 50/Dx_Replication50_males_DMR_raw_methylation_in_Discovery.csv", 
                        chroms = maleChroms, sort = TRUE, DMRid = TRUE)
methRep <- loadRegions("DMRs/Replication/Diagnosis Males 50/DMR_raw_methylation_Dx_Replication50_males.txt", 
                       chroms = maleChroms, sort = TRUE, DMRid = TRUE)
table(methDisc$start == methRep$start) # All TRUE

# Remove DMRs Covered in < 50% ASD in Discovery Samples
meth <- methDisc[,as.character(malesDisc$Sequencing_ID)]
missingASD <- apply(meth, 1, function(x) table(is.na(x))["TRUE"]) / length(as.character(malesDisc$Sequencing_ID)) # percent of ASD with no reads
missingASD[is.na(missingASD)] <- 0
methDisc <- subset(methDisc, missingASD <= 0.5) # 4648 / 4650 DMRs remaining
methRep <- subset(methRep, DMRid %in% methDisc$DMRid)

# Differential Methylation in Discovery Samples
meth <- t(methDisc[, as.character(malesDisc$Sequencing_ID)]) %>% as.data.frame
colnames(meth) <- methDisc$DMRid
meth$Sequencing_ID <- rownames(meth)
meth <- merge(x = malesDisc, y = meth, by = "Sequencing_ID", all = FALSE)
stats_RepDMRsInDisc <- DMRmethLm(methDisc$DMRid, catVars = NULL, contVars = "ADOScs", sampleData = meth,
                                 file = "Tables/Differential Methylation by ADOS Replication Males DMRs in Discovery Males ASD Only.txt")
table(stats_RepDMRsInDisc$pvalue < 0.05) # 234 / 4648 sig DMRs

# Differential Methylation in Replication Samples
meth <- t(methRep[, as.character(malesRep$Sequencing_ID)]) %>% as.data.frame
colnames(meth) <- methRep$DMRid
meth$Sequencing_ID <- rownames(meth)
meth <- merge(x = malesRep, y = meth, by = "Sequencing_ID", all = FALSE)
stats_RepDMRsInRep <- DMRmethLm(methRep$DMRid, catVars = NULL, contVars = "ADOScs", sampleData = meth,
                                file = "Tables/Differential Methylation by ADOS Replication Males DMRs in Replication Males ASD Only.txt")
table(stats_RepDMRsInRep$pvalue < 0.05) # 286 / 4648 sig DMRs

# Subset by Differential Methylation in Both Sets and Annotate
colnames(stats_RepDMRsInDisc) <- paste("Disc", colnames(stats_RepDMRsInDisc), sep = "_")
colnames(stats_RepDMRsInRep) <- paste("Rep", colnames(stats_RepDMRsInRep), sep = "_")
stats_RepDMRs <- merge(x = stats_RepDMRsInDisc, y = stats_RepDMRsInRep, by.x = "Disc_Region", by.y = "Rep_Region", all = FALSE, 
                       sort = FALSE)
stats_RepDMRs <- stats_RepDMRs[, !colnames(stats_RepDMRs) %in% c("Rep_Variable", "Rep_Term")]
colnames(stats_RepDMRs) <- str_replace_all(colnames(stats_RepDMRs), 
                                           pattern = c("Disc_Region" = "Region", "Disc_Variable" = "Variable", 
                                                       "Disc_Term" = "Term"))
stats_RepDMRsSig <- subset(stats_RepDMRs, Disc_pvalue < 0.05 & Rep_pvalue < 0.05 & Disc_Estimate * Rep_Estimate > 0) # 5 / 4648 replicated DMRs
stats_RepDMRsSig <- merge(x = methRep[,c("DMRid", "chr", "start", "end", "width", "L", "area", "beta", "percentDifference",
                                         "stat", "pval", "qval")], 
                          y = stats_RepDMRsSig, by.x = "DMRid", by.y = "Region", all = FALSE, sort = FALSE)
stats_RepDMRsSig <- getDMRanno(DMRstats = stats_RepDMRsSig, regDomains = regDomains, 
                               file = "Tables/Replicated Differential Methylation by ADOS in ASD Only Replication Males Diagnosis DMRs with Annotation.txt")

# ASD Only Females Differential Methylation Analysis by ADOS ---------------------------------------
# Sample Data ####
femalesDisc <- subset(samples, Sex == "F" & Platform == "HiSeqX10" & Diagnosis_Alg == "ASD", 
                    select = c("Sequencing_ID", "Diagnosis_Alg", "Sex", "ADOScs"))
femalesRep <- subset(samples, Sex == "F" & Platform == "HiSeq4000" & Diagnosis_Alg == "ASD", 
                   select = c("Sequencing_ID", "Diagnosis_Alg", "Sex", "ADOScs"))
femaleChroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrM")

# Discovery Females DMRs ####
# Data
methDisc <- loadRegions("DMRs/Discovery/Diagnosis Females 50/DMR_raw_methylation_Dx_Discovery50_females.txt", 
                        chroms = femaleChroms, sort = TRUE, DMRid = TRUE)
methRep <- loadRegions("DMRs/Discovery/Diagnosis Females 50/Dx_Discovery50_females_DMR_raw_methylation_in_Replication.csv", 
                       chroms = femaleChroms, sort = TRUE, DMRid = TRUE)
table(methDisc$start == methRep$start) # All TRUE

# Remove DMRs Covered in < 50% ASD in Replication Samples
meth <- methRep[,as.character(femalesRep$Sequencing_ID)]
missingASD <- apply(meth, 1, function(x) table(is.na(x))["TRUE"]) / length(as.character(femalesRep$Sequencing_ID)) # percent of ASD with no reads
missingASD[is.na(missingASD)] <- 0
methRep <- subset(methRep, missingASD <= 0.5) # 1420 / 1952 DMRs remaining
methDisc <- subset(methDisc, DMRid %in% methRep$DMRid)

# Differential Methylation in Discovery Samples
meth <- t(methDisc[, as.character(femalesDisc$Sequencing_ID)]) %>% as.data.frame
colnames(meth) <- methDisc$DMRid
meth$Sequencing_ID <- rownames(meth)
meth <- merge(x = femalesDisc, y = meth, by = "Sequencing_ID", all = FALSE)
stats_DiscDMRsInDisc <- DMRmethLm(methDisc$DMRid, catVars = NULL, contVars = "ADOScs", sampleData = meth,
                                  file = "Tables/Differential Methylation by ADOS Discovery Females DMRs in Discovery Females ASD Only.txt")
table(stats_DiscDMRsInDisc$pvalue < 0.05) # 92 / 1420 sig DMRs
sd(meth$ADOScs, na.rm = TRUE) #1.489893

# Differential Methylation in Replication Samples
meth <- t(methRep[, as.character(femalesRep$Sequencing_ID)]) %>% as.data.frame
colnames(meth) <- methRep$DMRid
meth$Sequencing_ID <- rownames(meth)
meth <- merge(x = femalesRep, y = meth, by = "Sequencing_ID", all = FALSE)
stats_DiscDMRsInRep <- DMRmethLm(methRep$DMRid, catVars = NULL, contVars = "ADOScs", sampleData = meth,
                                 file = "Tables/Differential Methylation by ADOS Discovery Females DMRs in Replication Females ASD Only.txt")
table(stats_DiscDMRsInRep$pvalue < 0.05) # 81 / 1420 sig DMRs
sd(meth$ADOScs, na.rm = TRUE) #2.302173

# Subset by Differential Methylation in Both Sets and Annotate
colnames(stats_DiscDMRsInDisc) <- paste("Disc", colnames(stats_DiscDMRsInDisc), sep = "_")
colnames(stats_DiscDMRsInRep) <- paste("Rep", colnames(stats_DiscDMRsInRep), sep = "_")
stats_DiscDMRs <- merge(x = stats_DiscDMRsInDisc, y = stats_DiscDMRsInRep, by.x = "Disc_Region", by.y = "Rep_Region", all = FALSE, 
                        sort = FALSE)
stats_DiscDMRs <- stats_DiscDMRs[, !colnames(stats_DiscDMRs) %in% c("Rep_Variable", "Rep_Term")]
colnames(stats_DiscDMRs) <- str_replace_all(colnames(stats_DiscDMRs), 
                                            pattern = c("Disc_Region" = "Region", "Disc_Variable" = "Variable", 
                                                        "Disc_Term" = "Term"))
stats_DiscDMRsSig <- subset(stats_DiscDMRs, Disc_pvalue < 0.05 & Rep_pvalue < 0.05 & Disc_Estimate * Rep_Estimate > 0) # 3 / 1420 replicated DMRs 
stats_DiscDMRsSig <- merge(x = methDisc[,c("DMRid", "chr", "start", "end", "width", "L", "area", "beta", "percentDifference",
                                           "stat", "pval", "qval")], 
                           y = stats_DiscDMRsSig, by.x = "DMRid", by.y = "Region", all = FALSE, sort = FALSE)
stats_DiscDMRsSig <- getDMRanno(DMRstats = stats_DiscDMRsSig, regDomains = regDomains, 
                                file = "Tables/Replicated Differential Methylation by ADOS in ASD Only Discovery Females Diagnosis DMRs with Annotation.txt")

# Replicated DMRs Scatterplots
stats_DiscDMRsSig$DMRid #[1] "DMR_723"  "DMR_957"  "DMR_1053"
DMR <- "DMR_723"
methPlot <- data.frame(Sequencing_ID = c(as.character(femalesDisc$Sequencing_ID), as.character(femalesRep$Sequencing_ID)),
                       SampleSet = c(rep("Discovery", nrow(femalesDisc)), rep("Replication", nrow(femalesRep))),
                       ADOScs = c(femalesDisc$ADOScs, femalesRep$ADOScs), DMRid = DMR, 
                       Methylation = unlist(c(methDisc[methDisc$DMRid == DMR, as.character(femalesDisc$Sequencing_ID)],
                                              methRep[methRep$DMRid == DMR, as.character(femalesRep$Sequencing_ID)]))*100)
g <- ggplot(data = methPlot)
g + 
        geom_smooth(aes(x = ADOScs, y = Methylation), method = "lm") +
        geom_point(aes(x = ADOScs, y = Methylation), color = "#FF3366", size = 2) +
        facet_grid(cols = vars(SampleSet)) +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = "None", panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black", size = 18), legend.background = element_blank(), 
              plot.margin = unit(c(0.5,1,1,1), "lines"), strip.background = element_blank(), 
              strip.text = element_text(size = 20), axis.title = element_text(size = 20)) +
        coord_cartesian(ylim = c(40, 100)) +
        scale_x_continuous(breaks = pretty_breaks(n = 5)) +
        scale_y_continuous(breaks = pretty_breaks(n = 5)) +
        xlab("ADOS Comparison Score") +
        ylab("DMR_723 Methylation (%)") +
        ggsave("Figures/DMR_723 Percent Methylation Plot Replicated ADOS in ASD Only Females.png", dpi = 600, width = 8, height = 5, 
               units = "in")

DMR <- "DMR_957"
methPlot <- data.frame(Sequencing_ID = c(as.character(femalesDisc$Sequencing_ID), as.character(femalesRep$Sequencing_ID)),
                       SampleSet = c(rep("Discovery", nrow(femalesDisc)), rep("Replication", nrow(femalesRep))),
                       ADOScs = c(femalesDisc$ADOScs, femalesRep$ADOScs), DMRid = DMR, 
                       Methylation = unlist(c(methDisc[methDisc$DMRid == DMR, as.character(femalesDisc$Sequencing_ID)],
                                              methRep[methRep$DMRid == DMR, as.character(femalesRep$Sequencing_ID)]))*100)
g <- ggplot(data = methPlot)
g + 
        geom_smooth(aes(x = ADOScs, y = Methylation), method = "lm") +
        geom_point(aes(x = ADOScs, y = Methylation), color = "#FF3366", size = 2) +
        facet_grid(cols = vars(SampleSet)) +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = "None", panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black", size = 18), legend.background = element_blank(), 
              plot.margin = unit(c(0.5,1,1,1), "lines"), strip.background = element_blank(), 
              strip.text = element_text(size = 20), axis.title = element_text(size = 20)) +
        coord_cartesian(ylim = c(0, 100)) +
        scale_x_continuous(breaks = pretty_breaks(n = 5)) +
        scale_y_continuous(breaks = pretty_breaks(n = 5)) +
        xlab("ADOS Comparison Score") +
        ylab("DMR_957 Methylation (%)") +
        ggsave("Figures/DMR_957 Percent Methylation Plot Replicated ADOS in ASD Only Females.png", dpi = 600, width = 8, height = 5, 
               units = "in")

DMR <- "DMR_1053"
methPlot <- data.frame(Sequencing_ID = c(as.character(femalesDisc$Sequencing_ID), as.character(femalesRep$Sequencing_ID)),
                       SampleSet = c(rep("Discovery", nrow(femalesDisc)), rep("Replication", nrow(femalesRep))),
                       ADOScs = c(femalesDisc$ADOScs, femalesRep$ADOScs), DMRid = DMR, 
                       Methylation = unlist(c(methDisc[methDisc$DMRid == DMR, as.character(femalesDisc$Sequencing_ID)],
                                              methRep[methRep$DMRid == DMR, as.character(femalesRep$Sequencing_ID)]))*100)
g <- ggplot(data = methPlot)
g + 
        geom_smooth(aes(x = ADOScs, y = Methylation), method = "lm") +
        geom_point(aes(x = ADOScs, y = Methylation), color = "#FF3366", size = 2) +
        facet_grid(cols = vars(SampleSet)) +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = "None", panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black", size = 18), legend.background = element_blank(), 
              plot.margin = unit(c(0.5,1,1,1), "lines"), strip.background = element_blank(), 
              strip.text = element_text(size = 20), axis.title = element_text(size = 20)) +
        coord_cartesian(ylim = c(15, 100)) +
        scale_x_continuous(breaks = pretty_breaks(n = 5)) +
        scale_y_continuous(breaks = pretty_breaks(n = 5)) +
        xlab("ADOS Comparison Score") +
        ylab("DMR_1053 Methylation (%)") +
        ggsave("Figures/DMR_1053 Percent Methylation Plot Replicated ADOS in ASD Only Females.png", dpi = 600, width = 8, height = 5, 
               units = "in")

# DMR Methylation ~ ADOS Volcano Plots
volcanoPlot <- data.frame(DMRid = stats_DiscDMRs$Region,
                          SampleSet = rep(c("Discovery", "Replication"), each = nrow(stats_DiscDMRs)),
                          Estimate = c(stats_DiscDMRs$Disc_Estimate * sd(femalesDisc$ADOScs, na.rm = TRUE),
                                       stats_DiscDMRs$Rep_Estimate * sd(femalesRep$ADOScs, na.rm = TRUE)) * 100,
                          log_pvalue = c(stats_DiscDMRs$Disc_log_pvalue, stats_DiscDMRs$Rep_log_pvalue))
volcanoPlot$Significant <- (volcanoPlot$log_pvalue > -log10(0.05)) %>% factor(levels = c("TRUE", "FALSE"))
## Estimate is now in units of change in DMR methylation per 1 ADOScs point
g <- ggplot(data = volcanoPlot)
g + 
        geom_point(data = volcanoPlot, aes(x = Estimate, y = log_pvalue, color = Significant), size = 1) +
        facet_grid(cols = vars(SampleSet), scales = "free") +
        theme_bw(base_size = 25) +
        theme(legend.direction = "vertical", legend.position = c(1.12, 0.88), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_text(size = 20),
              axis.text = element_text(color = "black", size = 18), legend.background = element_blank(), 
              plot.margin = unit(c(0.5,7,1,1), "lines"), strip.background = element_blank(),
              strip.text = element_text(size = 20), axis.title = element_text(size = 20), legend.text = element_text(size = 18)) +
        scale_x_continuous(breaks = pretty_breaks(n = 5)) +
        scale_y_continuous(breaks = pretty_breaks(n = 4), expand = expand_scale(mult = c(0.01, 0.03))) +
        scale_color_manual(name = expression(italic(p)*"< 0.05"), breaks = c("TRUE", "FALSE"), values = c("#FF3366", "#3366CC")) +
        xlab("DMR Methylation / ADOS Comparison Score") +
        ylab(expression("-log"[10]*"("*italic(p)*"-value)")) +
        ggsave("Figures/Discovery DMR Methylation by ADOS in ASD Only Females Volcano Plot.png", dpi = 600, width = 8, 
               height = 5, units = "in")

# Replication Females DMRs ####
# Data
methDisc <- loadRegions("DMRs/Replication/Diagnosis Females 100/Dx_Replication100_females_DMR_raw_methylation_in_Discovery.csv", 
                        chroms = femaleChroms, sort = TRUE, DMRid = TRUE)
methRep <- loadRegions("DMRs/Replication/Diagnosis Females 100/DMR_raw_methylation_Dx_Replication100_females.txt", 
                       chroms = femaleChroms, sort = TRUE, DMRid = TRUE)
table(methDisc$start == methRep$start) # All TRUE

# Remove DMRs Covered in < 50% ASD in Discovery Samples
meth <- methDisc[,as.character(femalesDisc$Sequencing_ID)]
missingASD <- apply(meth, 1, function(x) table(is.na(x))["TRUE"]) / length(as.character(femalesDisc$Sequencing_ID)) # percent of ASD with no reads
missingASD[is.na(missingASD)] <- 0
methDisc <- subset(methDisc, missingASD <= 0.5) # 8726 / 8728 DMRs remaining
methRep <- subset(methRep, DMRid %in% methDisc$DMRid)

# Differential Methylation in Discovery Samples
meth <- t(methDisc[, as.character(femalesDisc$Sequencing_ID)]) %>% as.data.frame
colnames(meth) <- methDisc$DMRid
meth$Sequencing_ID <- rownames(meth)
meth <- merge(x = femalesDisc, y = meth, by = "Sequencing_ID", all = FALSE)
stats_RepDMRsInDisc <- DMRmethLm(methDisc$DMRid, catVars = NULL, contVars = "ADOScs", sampleData = meth,
                                 file = "Tables/Differential Methylation by ADOS Replication Females DMRs in Discovery Females ASD Only.txt")
table(stats_RepDMRsInDisc$pvalue < 0.05) # 599 / 8726 sig DMRs
sd(meth$ADOScs, na.rm = TRUE) #1.489893

# Differential Methylation in Replication Samples
meth <- t(methRep[, as.character(femalesRep$Sequencing_ID)]) %>% as.data.frame
colnames(meth) <- methRep$DMRid
meth$Sequencing_ID <- rownames(meth)
meth <- merge(x = femalesRep, y = meth, by = "Sequencing_ID", all = FALSE)
stats_RepDMRsInRep <- DMRmethLm(methRep$DMRid, catVars = NULL, contVars = "ADOScs", sampleData = meth,
                                file = "Tables/Differential Methylation by ADOS Replication Females DMRs in Replication Females ASD Only.txt")
table(stats_RepDMRsInRep$pvalue < 0.05) # 584 / 8726 sig DMRs
sd(meth$ADOScs, na.rm = TRUE) #2.302173

# Subset by Differential Methylation in Both Sets and Annotate
colnames(stats_RepDMRsInDisc) <- paste("Disc", colnames(stats_RepDMRsInDisc), sep = "_")
colnames(stats_RepDMRsInRep) <- paste("Rep", colnames(stats_RepDMRsInRep), sep = "_")
stats_RepDMRs <- merge(x = stats_RepDMRsInDisc, y = stats_RepDMRsInRep, by.x = "Disc_Region", by.y = "Rep_Region", all = FALSE, 
                       sort = FALSE)
stats_RepDMRs <- stats_RepDMRs[, !colnames(stats_RepDMRs) %in% c("Rep_Variable", "Rep_Term")]
colnames(stats_RepDMRs) <- str_replace_all(colnames(stats_RepDMRs), 
                                           pattern = c("Disc_Region" = "Region", "Disc_Variable" = "Variable", 
                                                       "Disc_Term" = "Term"))
stats_RepDMRsSig <- subset(stats_RepDMRs, Disc_pvalue < 0.05 & Rep_pvalue < 0.05 & Disc_Estimate * Rep_Estimate > 0) # 22 / 8726 replicated DMRs
stats_RepDMRsSig <- merge(x = methRep[,c("DMRid", "chr", "start", "end", "width", "L", "area", "beta", "percentDifference",
                                         "stat", "pval", "qval")], 
                          y = stats_RepDMRsSig, by.x = "DMRid", by.y = "Region", all = FALSE, sort = FALSE)
stats_RepDMRsSig <- getDMRanno(DMRstats = stats_RepDMRsSig, regDomains = regDomains, 
                               file = "Tables/Replicated Differential Methylation by ADOS in ASD Only Replication Females Diagnosis DMRs with Annotation.txt")




