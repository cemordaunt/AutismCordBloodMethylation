# DMR Methylation by Diagnosis and Cell Proportions -----------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 5/7/20

# Setup -------------------------------------------------------------------
# Load Packages ####
sapply(c("reshape2", "scales", "tidyverse", "bsseq", "annotatr", "dmrseq", "rlist"), require, character.only = TRUE)

# Load Functions ####
source("R Scripts/DMR Analysis Functions.R")

# Modified DMRmethLm and .methLm to adjust for all 7 cord blood cell types, returns cell type stats too, all rows in FDR 
.methLmCell <- function(catVars, contVars, sampleData, meth, adj = TRUE){
        # Analyzes categorical and continuous variables for association with methylation, 
        # using linear regression with adjustment by cell type
        stats <- NULL
        methDataCol <- as.numeric(sampleData[,meth])
        if(!is.null(catVars)){
                for(i in 1:length(catVars)){
                        sampleDataCol <- sampleData[,catVars[i]]
                        if(adj){
                                temp <- tryCatch(summary(lm(methDataCol ~ sampleDataCol + sampleData$Bcell + sampleData$CD4T + 
                                                                    sampleData$CD8T + sampleData$Gran + sampleData$Mono + 
                                                                    sampleData$NK + sampleData$nRBC))$coefficients[-1,],
                                                 error = function(x){
                                                         message("Error with ", meth, " and ", catVars[i])
                                                         rep(NA, 4) # Returns NA if error in lm
                                                 })
                        }
                        else {
                                temp <- tryCatch(summary(lm(methDataCol ~ sampleDataCol))$coefficients[-1,],
                                                 error = function(x){
                                                         message("Error with ", meth, " and ", catVars[i])
                                                         rep(NA, 4) # Returns NA if error in lm
                                                 })
                        }
                        if(length(temp) == 4){
                                temp <- c(catVars[i], levels(sampleDataCol)[2], temp)
                        } 
                        else {
                                temp <- cbind(rep(catVars[i], nrow(temp)), 
                                              gsub("sampleDataCol", replacement = "", x = rownames(temp), fixed = TRUE), 
                                              temp)
                        }
                        stats <- rbind(stats, temp)
                }
        }
        if(!is.null(contVars)){
                for(i in 1:length(contVars)){
                        sampleDataCol <- as.numeric(sampleData[,contVars[i]])
                        sampleDataCol <- sampleDataCol/sd(sampleDataCol, na.rm = TRUE)
                        if(adj){
                                temp <- summary(lm(methDataCol ~ sampleDataCol + sampleData$Bcell + sampleData$CD4T + 
                                                           sampleData$CD8T + sampleData$Gran + sampleData$Mono + 
                                                           sampleData$NK + sampleData$nRBC))$coefficients[-1,]
                        }
                        else {
                                temp <- summary(lm(methDataCol ~ sampleDataCol))$coefficients[-1,]
                        }
                        temp <- c(contVars[i], contVars[i], temp)
                        stats <- rbind(stats, temp)
                }
        }
        rownames(stats) <- 1:nrow(stats)
        colnames(stats) <- c("Variable", "Term", "Estimate", "StdError", "tvalue", "pvalue")
        stats <- as.data.frame(stats, stringsAsFactors = FALSE)
        stats$Region <- meth
        return(stats)
}

DMRmethLmCell <- function(DMRs, catVars, contVars, sampleData, file, adj = TRUE){
        if(!is.null(adj)){
                message("[DMRmethLm] Getting DMR methylation by covariate stats using adjustment variable")
        }
        else {
                message("[DMRmethLm] Getting DMR methylation by covariate stats")
        }
        covStats <- lapply(DMRs, function(x){
                .methLmCell(catVars = catVars, contVars = contVars, sampleData = sampleData, meth = x, adj = adj) 
        }) %>% list.rbind()
        covStats$Variable <- as.character(covStats$Variable)
        covStats$Variable[covStats$Variable %in% c("Site", "MomEdu_detail")] <- paste(covStats$Variable[covStats$Variable %in% c("Site", "MomEdu_detail")],
                                                                                      covStats$Term[covStats$Variable %in% c("Site", "MomEdu_detail")], sep = "_")
        covStats$Variable <- gsub(" ", "_", covStats$Variable)
        covStats$Variable <- factor(covStats$Variable, levels = unique(covStats$Variable))
        covStats$Region <- factor(covStats$Region, levels = unique(covStats$Region))
        covStats[,c("Estimate", "StdError", "tvalue", "pvalue")] <- lapply(covStats[,c("Estimate", "StdError", "tvalue", 
                                                                                       "pvalue")], as.numeric)
        covStats$log_pvalue <- -log10(covStats$pvalue)
        covStats$qvalue <- p.adjust(covStats$pvalue, method = "fdr")
        covStats <- covStats[,c("Region", "Variable", "Term", "Estimate", "StdError", "tvalue", "pvalue", "log_pvalue", 
                                "qvalue")]
        message("[DMRmethLm] Complete, writing file")
        write.table(covStats, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
        return(covStats)
}

# Load Sample Data ####
samples <- read.csv("Tables/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo and Cell Type.csv", 
                    stringsAsFactors = FALSE)
samples$Diagnosis_Alg <- factor(samples$Diagnosis_Alg, levels = c("TD", "ASD"))

# Males Discovery Analysis ------------------------------------------------
# Load Methylation Data ####
subSamples <- subset(samples, Sex == "M" & Platform == "HiSeqX10", 
                     select = c("Sequencing_ID", "Diagnosis_Alg", "Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC"))
maleChroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")
DMRmeth <- loadRegions("DMRs/Discovery/Diagnosis Males 50/DMR_raw_methylation_Dx_Discovery50_males.txt", 
                       chroms = maleChroms, sort = TRUE, DMRid = TRUE)
meth <- DMRmeth[,grepl("JLCM", colnames(DMRmeth))] %>% t() %>% as.data.frame()
colnames(meth) <- DMRmeth$DMRid
meth$Sequencing_ID <- rownames(meth)
meth <- merge(x = subSamples, y = meth, by = "Sequencing_ID", all = FALSE)

# DMR Methylation by Diagnosis ####
diagnosis <- DMRmethLmCell(DMRmeth$DMRid, catVars = "Diagnosis_Alg", contVars = NULL, sampleData = meth,
                           file = "Tables/Differential Methylation by Diagnosis Males Discovery DMRs.txt", 
                           adj = FALSE)
table(abs(diagnosis$Estimate) > 0.05) # 628 / 635 DMRs
table(diagnosis$pvalue < 0.05) # 600 / 635 DMRs

# DMR Methylation by Diagnosis and Cell Type ####
diagnosisCell <- DMRmethLmCell(DMRmeth$DMRid, catVars = "Diagnosis_Alg", contVars = NULL, sampleData = meth,
                               file = "Tables/Differential Methylation by Diagnosis and Cell Proportion Males Discovery DMRs.txt", 
                               adj = TRUE)
diagnosisCellASD <- subset(diagnosisCell, Term == "ASD")
table(abs(diagnosisCellASD$Estimate) > 0.05) # 593 / 635 DMRs
table(diagnosisCellASD$pvalue < 0.05) # 538 / 635 DMRs

# Conclusion ####
table(abs(diagnosisCellASD$Estimate[abs(diagnosis$Estimate) > 0.05]) > 0.05) # 592 / 628 DMRs (94.3%) with estimate > 5%
table(diagnosisCellASD$pvalue[diagnosis$pvalue < 0.05] < 0.05) # 534 / 600 DMRs (89.0%) with p < 0.05

sigDiagnosisDMRs <- diagnosis$Region[abs(diagnosis$Estimate) > 0.05 & diagnosis$pvalue < 0.05] %>% as.character() # 596 / 635 DMRs
sigDiagnosisCellDMRs <- diagnosisCellASD$Region[abs(diagnosisCellASD$Estimate) > 0.05 & diagnosisCellASD$pvalue < 0.05] %>% as.character() # 522 / 635 DMRs
table(sigDiagnosisDMRs %in% sigDiagnosisCellDMRs) # 518 / 596 DMRs (86.9%) with estimate > 5% and p < 0.05

pdf("Figures/DMR Methylation and Cell Proportion Adjustment Plot Males Discovery.pdf", width = 8, height = 4)
par(mfcol = c(1,2), mar = c(5,5,1,1), cex = 0.8)
plot(diagnosis$Estimate, diagnosisCellASD$Estimate) 
plot(diagnosis$log_pvalue, diagnosisCellASD$log_pvalue)
dev.off()

#' With adjustment for all cell types, 87% of DMRs still have methylation difference > 5% and p < 0.05

# Males Replication Analysis ------------------------------------------------
# Load Methylation Data ####
subSamples <- subset(samples, Sex == "M" & Platform == "HiSeq4000", 
                     select = c("Sequencing_ID", "Diagnosis_Alg", "Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC"))
DMRmeth <- loadRegions("DMRs/Replication/Diagnosis Males 50/DMR_raw_methylation_Dx_Replication50_males.txt", 
                       chroms = maleChroms, sort = TRUE, DMRid = TRUE)
meth <- DMRmeth[,grepl("JLCM", colnames(DMRmeth))] %>% t() %>% as.data.frame()
colnames(meth) <- DMRmeth$DMRid
meth$Sequencing_ID <- rownames(meth)
meth <- merge(x = subSamples, y = meth, by = "Sequencing_ID", all = FALSE)

# DMR Methylation by Diagnosis ####
diagnosis <- DMRmethLmCell(DMRmeth$DMRid, catVars = "Diagnosis_Alg", contVars = NULL, sampleData = meth,
                           file = "Tables/Differential Methylation by Diagnosis Males Replication DMRs.txt", 
                           adj = FALSE)
table(abs(diagnosis$Estimate) > 0.05) # 4589 / 4650 DMRs
table(diagnosis$pvalue < 0.05) # 4311 / 4650 DMRs

# DMR Methylation by Diagnosis and Cell Type ####
diagnosisCell <- DMRmethLmCell(DMRmeth$DMRid, catVars = "Diagnosis_Alg", contVars = NULL, sampleData = meth,
                               file = "Tables/Differential Methylation by Diagnosis and Cell Proportion Males Replication DMRs.txt", 
                               adj = TRUE)
diagnosisCellASD <- subset(diagnosisCell, Term == "ASD")
table(abs(diagnosisCellASD$Estimate) > 0.05) # 4410 / 4650 DMRs
table(diagnosisCellASD$pvalue < 0.05) # 3616 / 4650 DMRs

# Conclusion ####
table(abs(diagnosisCellASD$Estimate[abs(diagnosis$Estimate) > 0.05]) > 0.05) # 4399 / 4589 DMRs (95.9%) with estimate > 5%
table(diagnosisCellASD$pvalue[diagnosis$pvalue < 0.05] < 0.05) # 3562 / 4311 DMRs (82.6%) with p < 0.05

sigDiagnosisDMRs <- diagnosis$Region[abs(diagnosis$Estimate) > 0.05 & diagnosis$pvalue < 0.05] %>% as.character() # 4284 / 4650 DMRs
sigDiagnosisCellDMRs <- diagnosisCellASD$Region[abs(diagnosisCellASD$Estimate) > 0.05 & diagnosisCellASD$pvalue < 0.05] %>% as.character() # 3569 / 4650 DMRs
table(sigDiagnosisDMRs %in% sigDiagnosisCellDMRs) # 3512 / 4284 DMRs (82.0%) with estimate > 5% and p < 0.05

pdf("Figures/DMR Methylation and Cell Proportion Adjustment Plot Males Replication.pdf", width = 8, height = 4)
par(mfcol = c(1,2), mar = c(5,5,1,1), cex = 0.8)
plot(diagnosis$Estimate, diagnosisCellASD$Estimate) 
plot(diagnosis$log_pvalue, diagnosisCellASD$log_pvalue)
dev.off()

#' With adjustment for all cell types, 82% of DMRs still have methylation difference > 5% and p < 0.05

# Females Discovery Analysis ------------------------------------------------
# Load Methylation Data ####
subSamples <- subset(samples, Sex == "F" & Platform == "HiSeqX10", 
                     select = c("Sequencing_ID", "Diagnosis_Alg", "Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC"))
femaleChroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrM")
DMRmeth <- loadRegions("DMRs/Discovery/Diagnosis Females 50/DMR_raw_methylation_Dx_Discovery50_females.txt", 
                       chroms = femaleChroms, sort = TRUE, DMRid = TRUE)
meth <- DMRmeth[,grepl("JLCM", colnames(DMRmeth))] %>% t() %>% as.data.frame()
colnames(meth) <- DMRmeth$DMRid
meth$Sequencing_ID <- rownames(meth)
meth <- merge(x = subSamples, y = meth, by = "Sequencing_ID", all = FALSE)

# DMR Methylation by Diagnosis ####
diagnosis <- DMRmethLmCell(DMRmeth$DMRid, catVars = "Diagnosis_Alg", contVars = NULL, sampleData = meth,
                           file = "Tables/Differential Methylation by Diagnosis Females Discovery DMRs.txt", 
                           adj = FALSE)
table(abs(diagnosis$Estimate) > 0.05) # 1926 / 1952 DMRs
table(diagnosis$pvalue < 0.05) # 1715 / 1952 DMRs

# DMR Methylation by Diagnosis and Cell Type ####
diagnosisCell <- DMRmethLmCell(DMRmeth$DMRid, catVars = "Diagnosis_Alg", contVars = NULL, sampleData = meth,
                               file = "Tables/Differential Methylation by Diagnosis and Cell Proportion Females Discovery DMRs.txt", 
                               adj = TRUE)
diagnosisCellASD <- subset(diagnosisCell, Term == "ASD")
table(abs(diagnosisCellASD$Estimate) > 0.05) # 1920 / 1952 DMRs
table(diagnosisCellASD$pvalue < 0.05) # 1651 / 1952 DMRs

# Conclusion ####
table(abs(diagnosisCellASD$Estimate[abs(diagnosis$Estimate) > 0.05]) > 0.05) # 1907 / 1926 DMRs (99.0%) with estimate > 5%
table(diagnosisCellASD$pvalue[diagnosis$pvalue < 0.05] < 0.05) # 1561 / 1715 DMRs (91.0%) with p < 0.05

sigDiagnosisDMRs <- diagnosis$Region[abs(diagnosis$Estimate) > 0.05 & diagnosis$pvalue < 0.05] %>% as.character() # 1705 / 1952 DMRs
sigDiagnosisCellDMRs <- diagnosisCellASD$Region[abs(diagnosisCellASD$Estimate) > 0.05 & diagnosisCellASD$pvalue < 0.05] %>% as.character() # 1647 / 1952 DMRs
table(sigDiagnosisDMRs %in% sigDiagnosisCellDMRs) # 1553 / 1705 DMRs (91.1%) with estimate > 5% and p < 0.05

pdf("Figures/DMR Methylation and Cell Proportion Adjustment Plot Females Discovery.pdf", width = 8, height = 4)
par(mfcol = c(1,2), mar = c(5,5,1,1), cex = 0.8)
plot(diagnosis$Estimate, diagnosisCellASD$Estimate) 
plot(diagnosis$log_pvalue, diagnosisCellASD$log_pvalue)
dev.off()

#' With adjustment for all cell types, 91% of DMRs still have methylation difference > 5% and p < 0.05

# Females Replication Analysis ------------------------------------------------
# Load Methylation Data ####
subSamples <- subset(samples, Sex == "F" & Platform == "HiSeq4000", 
                     select = c("Sequencing_ID", "Diagnosis_Alg", "Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC"))
DMRmeth <- loadRegions("DMRs/Replication/Diagnosis Females 100/DMR_raw_methylation_Dx_Replication100_females.txt", 
                       chroms = femaleChroms, sort = TRUE, DMRid = TRUE)
meth <- DMRmeth[,grepl("JLCM", colnames(DMRmeth))] %>% t() %>% as.data.frame()
colnames(meth) <- DMRmeth$DMRid
meth$Sequencing_ID <- rownames(meth)
meth <- merge(x = subSamples, y = meth, by = "Sequencing_ID", all = FALSE)

# DMR Methylation by Diagnosis ####
diagnosis <- DMRmethLm(DMRmeth$DMRid, catVars = "Diagnosis_Alg", contVars = NULL, sampleData = meth,
                       file = "Tables/Differential Methylation by Diagnosis Females Replication DMRs.txt", 
                       adj = NULL)
table(abs(diagnosis$Estimate) > 0.05) # 8704 / 8728 DMRs
table(diagnosis$pvalue < 0.05) # 6311 / 8728 DMRs
sigDiagnosisDMRs <- diagnosis$Region[abs(diagnosis$Estimate) > 0.05 & diagnosis$pvalue < 0.05] %>% as.character() # 6308 / 8728 DMRs

# DMR Methylation by Diagnosis and Cell Type ####
# Not enough samples to adjust for all cell types simultaneously, so adjust individually

# Adjust for nRBCs only
diagnosis_nRBC <- DMRmethLm(DMRmeth$DMRid, catVars = "Diagnosis_Alg", contVars = NULL, sampleData = meth,
                           file = "Tables/Differential Methylation by Diagnosis and nRBC Females Replication DMRs.txt", 
                           adj = "nRBC")
table(abs(diagnosis_nRBC$Estimate) > 0.05) # 8696 / 8728 DMRs
table(diagnosis_nRBC$pvalue < 0.05) # 5798 / 8728 DMRs

table(abs(diagnosis_nRBC$Estimate[abs(diagnosis$Estimate) > 0.05]) > 0.05) # 8693 / 8704 DMRs (99.9%) with estimate > 5%
table(diagnosis_nRBC$pvalue[diagnosis$pvalue < 0.05] < 0.05) # 5481 / 6311 DMRs (86.8%) with p < 0.05

sigDiagnosis_nRBC_DMRs <- diagnosis_nRBC$Region[abs(diagnosis_nRBC$Estimate) > 0.05 & diagnosis_nRBC$pvalue < 0.05] %>% as.character() # 5794 / 8728 DMRs
table(sigDiagnosisDMRs %in% sigDiagnosis_nRBC_DMRs) # 5477 / 6308 DMRs (86.8%) with estimate > 5% and p < 0.05

pdf("Figures/DMR Methylation and nRBC Adjustment Plot Females Replication.pdf", width = 8, height = 4)
par(mfcol = c(1,2), mar = c(5,5,1,1), cex = 0.8)
plot(diagnosis$Estimate, diagnosis_nRBC$Estimate) 
plot(diagnosis$log_pvalue, diagnosis_nRBC$log_pvalue, xlim = c(0,12), ylim = c(0,12))
dev.off()

# Adjust for Bcells only
diagnosis_Bcell <- DMRmethLm(DMRmeth$DMRid, catVars = "Diagnosis_Alg", contVars = NULL, sampleData = meth,
                            file = "Tables/Differential Methylation by Diagnosis and Bcell Females Replication DMRs.txt", 
                            adj = "Bcell")
table(abs(diagnosis_Bcell$Estimate) > 0.05) # 8697 / 8728 DMRs
table(diagnosis_Bcell$pvalue < 0.05) # 5798 / 8728 DMRs

table(abs(diagnosis_Bcell$Estimate[abs(diagnosis$Estimate) > 0.05]) > 0.05) # 8688 / 8704 DMRs (99.8%) with estimate > 5%
table(diagnosis_Bcell$pvalue[diagnosis$pvalue < 0.05] < 0.05) # 5549 / 6311 DMRs (87.9%) with p < 0.05

sigDiagnosis_Bcell_DMRs <- diagnosis_Bcell$Region[abs(diagnosis_Bcell$Estimate) > 0.05 & diagnosis_Bcell$pvalue < 0.05] %>% as.character() # 5793 / 8728 DMRs
table(sigDiagnosisDMRs %in% sigDiagnosis_Bcell_DMRs) # 5544 / 6308 DMRs (87.9%) with estimate > 5% and p < 0.05

pdf("Figures/DMR Methylation and Bcell Adjustment Plot Females Replication.pdf", width = 8, height = 4)
par(mfcol = c(1,2), mar = c(5,5,1,1), cex = 0.8)
plot(diagnosis$Estimate, diagnosis_Bcell$Estimate) 
plot(diagnosis$log_pvalue, diagnosis_Bcell$log_pvalue, xlim = c(0,12), ylim = c(0,12))
dev.off()

# Adjust for CD4Ts only
diagnosis_CD4T <- DMRmethLm(DMRmeth$DMRid, catVars = "Diagnosis_Alg", contVars = NULL, sampleData = meth,
                             file = "Tables/Differential Methylation by Diagnosis and CD4T Females Replication DMRs.txt", 
                             adj = "CD4T")
table(abs(diagnosis_CD4T$Estimate) > 0.05) # 8684 / 8728 DMRs
table(diagnosis_CD4T$pvalue < 0.05) # 5518 / 8728 DMRs

table(abs(diagnosis_CD4T$Estimate[abs(diagnosis$Estimate) > 0.05]) > 0.05) # 8680 / 8704 DMRs (99.7%) with estimate > 5%
table(diagnosis_CD4T$pvalue[diagnosis$pvalue < 0.05] < 0.05) # 5312 / 6311 DMRs (84.2%) with p < 0.05

sigDiagnosis_CD4T_DMRs <- diagnosis_CD4T$Region[abs(diagnosis_CD4T$Estimate) > 0.05 & diagnosis_CD4T$pvalue < 0.05] %>% as.character() # 5514 / 8728 DMRs
table(sigDiagnosisDMRs %in% sigDiagnosis_CD4T_DMRs) # 5307 / 6308 DMRs (84.1%) with estimate > 5% and p < 0.05

pdf("Figures/DMR Methylation and CD4T Adjustment Plot Females Replication.pdf", width = 8, height = 4)
par(mfcol = c(1,2), mar = c(5,5,1,1), cex = 0.8)
plot(diagnosis$Estimate, diagnosis_CD4T$Estimate) 
plot(diagnosis$log_pvalue, diagnosis_CD4T$log_pvalue, xlim = c(0,12), ylim = c(0,12))
dev.off()

# Adjust for CD8Ts only
diagnosis_CD8T <- DMRmethLm(DMRmeth$DMRid, catVars = "Diagnosis_Alg", contVars = NULL, sampleData = meth,
                            file = "Tables/Differential Methylation by Diagnosis and CD8T Females Replication DMRs.txt", 
                            adj = "CD8T")
table(abs(diagnosis_CD8T$Estimate) > 0.05) # 8690 / 8728 DMRs
table(diagnosis_CD8T$pvalue < 0.05) # 5635 / 8728 DMRs

table(abs(diagnosis_CD8T$Estimate[abs(diagnosis$Estimate) > 0.05]) > 0.05) # 8686 / 8704 DMRs (99.8%) with estimate > 5%
table(diagnosis_CD8T$pvalue[diagnosis$pvalue < 0.05] < 0.05) # 5439 / 6311 DMRs (86.2%) with p < 0.05

sigDiagnosis_CD8T_DMRs <- diagnosis_CD8T$Region[abs(diagnosis_CD8T$Estimate) > 0.05 & diagnosis_CD8T$pvalue < 0.05] %>% as.character() # 5631 / 8728 DMRs
table(sigDiagnosisDMRs %in% sigDiagnosis_CD8T_DMRs) # 5436 / 6308 DMRs (86.2%) with estimate > 5% and p < 0.05

pdf("Figures/DMR Methylation and CD8T Adjustment Plot Females Replication.pdf", width = 8, height = 4)
par(mfcol = c(1,2), mar = c(5,5,1,1), cex = 0.8)
plot(diagnosis$Estimate, diagnosis_CD8T$Estimate) 
plot(diagnosis$log_pvalue, diagnosis_CD8T$log_pvalue, xlim = c(0,12), ylim = c(0,12))
dev.off()

# Adjust for Grans only
diagnosis_Gran <- DMRmethLm(DMRmeth$DMRid, catVars = "Diagnosis_Alg", contVars = NULL, sampleData = meth,
                            file = "Tables/Differential Methylation by Diagnosis and Gran Females Replication DMRs.txt", 
                            adj = "Gran")
table(abs(diagnosis_Gran$Estimate) > 0.05) # 8692 / 8728 DMRs
table(diagnosis_Gran$pvalue < 0.05) # 5673 / 8728 DMRs

table(abs(diagnosis_Gran$Estimate[abs(diagnosis$Estimate) > 0.05]) > 0.05) # 8684 / 8704 DMRs (99.8%) with estimate > 5%
table(diagnosis_Gran$pvalue[diagnosis$pvalue < 0.05] < 0.05) # 5439 / 6311 DMRs (86.2%) with p < 0.05

sigDiagnosis_Gran_DMRs <- diagnosis_Gran$Region[abs(diagnosis_Gran$Estimate) > 0.05 & diagnosis_Gran$pvalue < 0.05] %>% as.character() # 5668 / 8728 DMRs
table(sigDiagnosisDMRs %in% sigDiagnosis_Gran_DMRs) # 5434 / 6308 DMRs (86.1%) with estimate > 5% and p < 0.05

pdf("Figures/DMR Methylation and Gran Adjustment Plot Females Replication.pdf", width = 8, height = 4)
par(mfcol = c(1,2), mar = c(5,5,1,1), cex = 0.8)
plot(diagnosis$Estimate, diagnosis_Gran$Estimate) 
plot(diagnosis$log_pvalue, diagnosis_Gran$log_pvalue, xlim = c(0,12), ylim = c(0,12))
dev.off()

# Adjust for Monos only
diagnosis_Mono <- DMRmethLm(DMRmeth$DMRid, catVars = "Diagnosis_Alg", contVars = NULL, sampleData = meth,
                            file = "Tables/Differential Methylation by Diagnosis and Mono Females Replication DMRs.txt", 
                            adj = "Mono")
table(abs(diagnosis_Mono$Estimate) > 0.05) # 8703 / 8728 DMRs
table(diagnosis_Mono$pvalue < 0.05) # 5938 / 8728 DMRs

table(abs(diagnosis_Mono$Estimate[abs(diagnosis$Estimate) > 0.05]) > 0.05) # 8699 / 8704 DMRs (99.9%) with estimate > 5%
table(diagnosis_Mono$pvalue[diagnosis$pvalue < 0.05] < 0.05) # 5678 / 6311 DMRs (90.0%) with p < 0.05

sigDiagnosis_Mono_DMRs <- diagnosis_Mono$Region[abs(diagnosis_Mono$Estimate) > 0.05 & diagnosis_Mono$pvalue < 0.05] %>% as.character() # 5935 / 8728 DMRs
table(sigDiagnosisDMRs %in% sigDiagnosis_Mono_DMRs) # 5675 / 6308 DMRs (90.0%) with estimate > 5% and p < 0.05

pdf("Figures/DMR Methylation and Mono Adjustment Plot Females Replication.pdf", width = 8, height = 4)
par(mfcol = c(1,2), mar = c(5,5,1,1), cex = 0.8)
plot(diagnosis$Estimate, diagnosis_Mono$Estimate) 
plot(diagnosis$log_pvalue, diagnosis_Mono$log_pvalue, xlim = c(0,12), ylim = c(0,12))
dev.off()

# Adjust for NKs only
diagnosis_NK <- DMRmethLm(DMRmeth$DMRid, catVars = "Diagnosis_Alg", contVars = NULL, sampleData = meth,
                            file = "Tables/Differential Methylation by Diagnosis and NK Females Replication DMRs.txt", 
                            adj = "NK")
table(abs(diagnosis_NK$Estimate) > 0.05) # 8617 / 8728 DMRs
table(diagnosis_NK$pvalue < 0.05) # 4764 / 8728 DMRs

table(abs(diagnosis_NK$Estimate[abs(diagnosis$Estimate) > 0.05]) > 0.05) # 8608 / 8704 DMRs (98.9%) with estimate > 5%
table(diagnosis_NK$pvalue[diagnosis$pvalue < 0.05] < 0.05) # 4592 / 6311 DMRs (72.8%) with p < 0.05

sigDiagnosis_NK_DMRs <- diagnosis_NK$Region[abs(diagnosis_NK$Estimate) > 0.05 & diagnosis_NK$pvalue < 0.05] %>% as.character() # 4756 / 8728 DMRs
table(sigDiagnosisDMRs %in% sigDiagnosis_NK_DMRs) # 4585 / 6308 DMRs (72.7%) with estimate > 5% and p < 0.05

pdf("Figures/DMR Methylation and NK Adjustment Plot Females Replication.pdf", width = 8, height = 4)
par(mfcol = c(1,2), mar = c(5,5,1,1), cex = 0.8)
plot(diagnosis$Estimate, diagnosis_NK$Estimate) 
plot(diagnosis$log_pvalue, diagnosis_NK$log_pvalue, xlim = c(0,12), ylim = c(0,12))
dev.off()

# Conclusions ####
sigAllDMRs <- intersect(sigDiagnosisDMRs, sigDiagnosis_Bcell_DMRs) %>% intersect(sigDiagnosis_CD4T_DMRs) %>%
        intersect(sigDiagnosis_CD8T_DMRs) %>% intersect(sigDiagnosis_Gran_DMRs) %>% 
        intersect(sigDiagnosis_Mono_DMRs) %>% intersect(sigDiagnosis_NK_DMRs) %>% 
        intersect(sigDiagnosis_nRBC_DMRs) # 4287 / 6308 DMRs (68.0%)
#' With individual adjustment for any of these 7 cell types, 
#' 68% of DMRs still have methylation difference > 5% and p < 0.05

sigAllDMRs <- intersect(sigDiagnosisDMRs, sigDiagnosis_Bcell_DMRs) %>% intersect(sigDiagnosis_CD4T_DMRs) %>%
        intersect(sigDiagnosis_CD8T_DMRs) %>% intersect(sigDiagnosis_Gran_DMRs) %>% 
        intersect(sigDiagnosis_Mono_DMRs) %>% intersect(sigDiagnosis_nRBC_DMRs) # 4807 / 6308 DMRs 76%
#' With individual adjustment for any of these 7 cell types except NK cells, 
#' 76% of DMRs still have methylation difference > 5% and p < 0.05

#' DMRs still with methylation difference > 5% and p < 0.05 with adjustment for
#' Bcells, 88% of DMRs
#' CD4Ts, 84% of DMRs
#' CD8Ts, 86% of DMRs
#' Grans, 86% of DMRs
#' Monos, 90% of DMRs
#' NKs, 73% of DMRs
#' nRBCs, 87% of DMRs

