# Global Methylation and Covariate Analysis ####
# Autism Cord Blood Methylation Project
# Charles Mordaunt
# 5/31/19
# Excluded JLCM032B and JLCM050B

# Packages ####
sapply(c("reshape2", "scales", "tidyverse"), require, character.only = TRUE)

# Functions ####
catFisher <- function(variables, sampleData, diagnosis, numHypotheses = length(variables)){
        # Analyzes table of categorical covariates for differences by ASD vs TD diagnosis using Fisher's exact test
        stats <- NULL
        for(i in 1:length(variables)){
                counts <- table(sampleData[,variables[i]], sampleData[,diagnosis])
                p <- fisher.test(sampleData[,variables[i]], sampleData[,diagnosis], workspace=2e7)$p.value
                temp <- cbind(variables[i], rownames(counts), counts, p)
                stats <- rbind(stats, temp)
        }
        rownames(stats) <- 1:nrow(stats)
        colnames(stats) <- c("Variable", "Value", "TD", "ASD", "pvalue")
        stats <- as.data.frame(stats, stringsAsFactors = FALSE)
        stats[,c("Variable", "Value")] <- lapply(stats[,c("Variable", "Value")], as.factor)
        stats[,c("TD", "ASD")] <- lapply(stats[,c("TD", "ASD")], as.integer)
        stats$pvalue <- as.numeric(stats$pvalue)
        varSums <- aggregate(as.matrix(stats[,c("TD", "ASD")]) ~ Variable, data=stats, FUN=sum)
        stats$per_TD <- stats$TD / varSums$TD[match(stats$Variable, varSums$Variable)]
        stats$per_ASD <- stats$ASD / varSums$ASD[match(stats$Variable, varSums$Variable)]
        stats$nper_TD <- paste(stats$TD, " (", round(stats$per_TD*100,0), ")", sep = "")
        stats$nper_ASD <- paste(stats$ASD, " (", round(stats$per_ASD*100,0), ")", sep = "")
        stats_padj <- unique(stats[,c("Variable", "pvalue")])
        stats_padj$qvalue <- p.adjust(stats_padj$pvalue, method = "fdr", n = numHypotheses)
        stats$qvalue <- stats_padj$qvalue[match(stats$pvalue, stats_padj$pvalue)]
        stats <- stats[,c("Variable", "Value", "TD", "ASD", "per_TD", "per_ASD", "nper_TD", "nper_ASD", "pvalue", "qvalue")]
        return(stats)
}

contANOVA <- function(variables, sampleData, diagnosis, numHypotheses = length(variables)){
        # Compares continuous variables for differences by diagnosis using one-way ANOVA
        stats <- NULL
        for(i in 1:length(variables)){
                means <- aggregate(sampleData[,variables[i]] ~ sampleData[,diagnosis], FUN=mean)[,2]
                sds <- aggregate(sampleData[,variables[i]] ~ sampleData[,diagnosis], FUN=sd)[,2]
                p <- summary(aov(sampleData[,variables[i]] ~ sampleData[,diagnosis]))[[1]][1, "Pr(>F)"]
                temp <- cbind(rep(variables[i], length(means)), levels(sampleData[,diagnosis]), means, sds, rep(p, length(means)))
                stats <- rbind(stats, temp)
        }
        colnames(stats) <- c("Variable", "Diagnosis", "Mean", "SD", "pvalue")
        stats <- as.data.frame(stats, stringsAsFactors = FALSE)
        stats[,c("Variable", "Diagnosis")] <- lapply(stats[,c("Variable", "Diagnosis")], as.factor)
        stats[,c("Mean", "SD", "pvalue")] <- lapply(stats[,c("Mean", "SD", "pvalue")], as.numeric)
        stats$MeanSD <- paste(round(stats$Mean, 2), " (", round(stats$SD, 2), ")", sep = "")
        stats_padj <- unique(stats[,c("Variable", "pvalue")])
        stats_padj$qvalue <- p.adjust(stats_padj$pvalue, method = "fdr", n = numHypotheses)
        stats$qvalue <- stats_padj$qvalue[match(stats$pvalue, stats_padj$pvalue)]
        stats <- stats[,c("Variable", "Diagnosis", "Mean", "SD", "MeanSD", "pvalue", "qvalue")]
        return(stats)
}

ggBoxPlot <- function(data, x, y, fill, xlab, ylab, breaks = c("TD", "ASD"), 
                      values = c("TD"="#3366CC", "ASD"="#FF3366"), file, axis.title.x = element_blank()) {
        # Plots a boxplot comparing a continous variable by diagnosis
        g <- ggplot(data = data)
        g + 
                geom_boxplot(aes(x=x, y=y, fill = fill), size = 1.1) +
                theme_bw(base_size = 25) +
                theme(legend.direction = 'horizontal', legend.position = c(0.87,1.03), panel.grid.major = element_blank(), 
                      panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25),
                      legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank(),
                      axis.text = element_text(color = "black"), legend.background = element_blank(), 
                      plot.margin = unit(c(2,1,1,1), "lines"), axis.title.x = axis.title.x) +
                scale_y_continuous(breaks=pretty_breaks(n=4)) +
                xlab(xlab) +
                ylab(ylab) +
                scale_fill_manual(breaks = breaks, values = values)
        ggsave(file = file, dpi = 600, width = 7, height = 7, units = "in")
}

methLm <- function(catVars, contVars, sampleData, globalMeth){
        # Analyzes categorical and continuous variables for association with global methylation, using linear regression
        stats <- NULL
        methDataCol <- as.numeric(sampleData[,globalMeth])
        for(i in 1:length(catVars)){
                sampleDataCol <- sampleData[,catVars[i]]
                temp <- summary(lm(methDataCol ~ sampleDataCol))$coefficients[-1,]
                if(length(temp) == 4){
                        temp <- c(catVars[i], levels(sampleDataCol)[2], temp)
                } else {
                        temp <- cbind(rep(catVars[i], nrow(temp)), rownames(temp), temp)
                }
                stats <- rbind(stats, temp)
        }
        for(i in 1:length(contVars)){
                sampleDataCol <- as.numeric(sampleData[,contVars[i]])
                sampleDataCol <- sampleDataCol/sd(sampleDataCol, na.rm = TRUE)
                temp <- summary(lm(methDataCol ~ sampleDataCol))$coefficients[-1,]
                temp <- c(contVars[i], contVars[i], temp)
                stats <- rbind(stats, temp)
        }
        rownames(stats) <- 1:nrow(stats)
        colnames(stats) <- c("Variable", "Term", "Estimate", "StdError", "tvalue", "pvalue")
        stats <- as.data.frame(stats, stringsAsFactors = FALSE)
        stats[,c("Variable")] <- as.factor(stats[,"Variable"])
        stats[,c("Estimate", "StdError", "tvalue", "pvalue")] <- lapply(stats[,c("Estimate", "StdError", "tvalue", "pvalue")], as.numeric)
        stats$qvalue <- p.adjust(stats$pvalue, method = "fdr")
        return(stats)
}

methLmPlatform <- function(catVars, contVars, sampleData, globalMeth, platform){
        # Analyzes categorical and continuous variables for association with global methylation, adjusting for platform using linear regression
        stats <- NULL
        methDataCol <- as.numeric(sampleData[,globalMeth])
        for(i in 1:length(catVars)){
                sampleDataCol <- sampleData[,catVars[i]]
                temp <- summary(lm(methDataCol ~ sampleDataCol + sampleData[,platform]))$coefficients[-1,]
                temp <- cbind(rep(catVars[i], nrow(temp)), rownames(temp), temp)
                stats <- rbind(stats, temp)
        }
        for(i in 1:length(contVars)){
                sampleDataCol <- as.numeric(sampleData[,contVars[i]])
                sampleDataCol <- sampleDataCol/sd(sampleDataCol, na.rm = TRUE)
                temp <- summary(lm(methDataCol ~ sampleDataCol + sampleData[,platform]))$coefficients[-1,]
                temp <- cbind(rep(contVars[i], nrow(temp)), rownames(temp), temp)
                stats <- rbind(stats, temp)
        }
        rownames(stats) <- 1:nrow(stats)
        colnames(stats) <- c("Variable", "Term", "Estimate", "StdError", "tvalue", "pvalue")
        stats <- as.data.frame(stats, stringsAsFactors = FALSE)
        stats[,c("Variable")] <- as.factor(stats[,"Variable"])
        stats[,c("Estimate", "StdError", "tvalue", "pvalue")] <- lapply(stats[,c("Estimate", "StdError", "tvalue", "pvalue")], as.numeric)
        stats$qvalue <- p.adjust(stats$pvalue, method = "fdr")
        return(stats)
}

methLmAdj <- function(catVars, contVars, sampleData, globalMeth, adjVar1, adjVar2 = NULL, adjVar3 = NULL){
        # Analyzes categorical and continuous variables for association with global methylation, adjusting for 1 or more covariates using linear regression
        stats <- NULL
        methDataCol <- as.numeric(sampleData[,globalMeth])
        if(is.null(adjVar2) & is.null(adjVar3)){
                for(i in 1:length(catVars)){
                        sampleDataCol <- sampleData[,catVars[i]]
                        temp <- summary(lm(methDataCol ~ sampleDataCol + sampleData[,adjVar1]))$coefficients[-1,]
                        temp <- cbind(rep(catVars[i], nrow(temp)), rownames(temp), temp)
                        stats <- rbind(stats, temp)
                }
                for(i in 1:length(contVars)){
                        sampleDataCol <- as.numeric(sampleData[,contVars[i]])
                        sampleDataCol <- sampleDataCol/sd(sampleDataCol, na.rm = TRUE)
                        temp <- summary(lm(methDataCol ~ sampleDataCol + sampleData[,adjVar1]))$coefficients[-1,]
                        temp <- cbind(rep(contVars[i], nrow(temp)), rownames(temp), temp)
                        stats <- rbind(stats, temp)
                }
        } else {
                if(is.null(adjVar3)){
                        for(i in 1:length(catVars)){
                                sampleDataCol <- sampleData[,catVars[i]]
                                temp <- summary(lm(methDataCol ~ sampleDataCol + sampleData[,adjVar1] + sampleData[,adjVar2]))$coefficients[-1,]
                                temp <- cbind(rep(catVars[i], nrow(temp)), rownames(temp), temp)
                                stats <- rbind(stats, temp)
                        }
                        for(i in 1:length(contVars)){
                                sampleDataCol <- as.numeric(sampleData[,contVars[i]])
                                sampleDataCol <- sampleDataCol/sd(sampleDataCol, na.rm = TRUE)
                                temp <- summary(lm(methDataCol ~ sampleDataCol + sampleData[,adjVar1] + sampleData[,adjVar2]))$coefficients[-1,]
                                temp <- cbind(rep(contVars[i], nrow(temp)), rownames(temp), temp)
                                stats <- rbind(stats, temp)
                        }
                } else {
                        for(i in 1:length(catVars)){
                                sampleDataCol <- sampleData[,catVars[i]]
                                temp <- summary(lm(methDataCol ~ sampleDataCol + sampleData[,adjVar1] + sampleData[,adjVar2] + sampleData[,adjVar3]))$coefficients[-1,]
                                temp <- cbind(rep(catVars[i], nrow(temp)), rownames(temp), temp)
                                stats <- rbind(stats, temp)
                        }
                        for(i in 1:length(contVars)){
                                sampleDataCol <- as.numeric(sampleData[,contVars[i]])
                                sampleDataCol <- sampleDataCol/sd(sampleDataCol, na.rm = TRUE)
                                temp <- summary(lm(methDataCol ~ sampleDataCol + sampleData[,adjVar1] + sampleData[,adjVar2] + sampleData[,adjVar3]))$coefficients[-1,]
                                temp <- cbind(rep(contVars[i], nrow(temp)), rownames(temp), temp)
                                stats <- rbind(stats, temp)
                        }
                }
                
        }
        rownames(stats) <- 1:nrow(stats)
        colnames(stats) <- c("Variable", "Term", "Estimate", "StdError", "tvalue", "pvalue")
        stats <- as.data.frame(stats, stringsAsFactors = FALSE)
        stats[,c("Variable")] <- as.factor(stats[,"Variable"])
        stats[,c("Estimate", "StdError", "tvalue", "pvalue")] <- lapply(stats[,c("Estimate", "StdError", "tvalue", "pvalue")], as.numeric)
        stats <- subset(stats, !grepl("adjVar", x = stats$Term, fixed = TRUE))
        stats$qvalue <- p.adjust(stats$pvalue, method = "fdr")
        return(stats)
}

ggScatterPlot <- function(x, y, groupVar, fileName, xlab, ylab, xlim = NULL, ylim = NULL, legendPos = c(0.87,1.03)){
        # Plots 2 continuous variables colored by a grouping variable and writes the file
        g <- ggplot()
        g + 
                geom_smooth(aes(x=x, y=y), method="lm") +
                geom_point(aes(x=x, y=y, color=groupVar), size=3) +
                theme_bw(base_size = 25) +
                theme(legend.direction = 'horizontal', legend.position = legendPos, panel.grid.major = element_blank(), 
                      panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
                      legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank(),
                      axis.text = element_text(color = "black"), legend.background = element_blank(), 
                      plot.margin = unit(c(2,1,1,1), "lines")) +
                coord_cartesian(xlim = xlim, ylim = ylim) +
                scale_x_continuous(breaks=pretty_breaks(n=5)) +
                scale_y_continuous(breaks=pretty_breaks(n=5)) +
                xlab(xlab) +
                ylab(ylab) +
                scale_color_manual(breaks = c(levels(groupVar)[1], levels(groupVar)[2]), values = c("#3366CC", "#FF3366"))
        ggsave(fileName, dpi = 600, width = 8, height = 7, units = "in")
}

# Data ####
samples <- read.delim(file = "Merged Database/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo.txt", sep = "\t", header = TRUE, 
                      stringsAsFactors = FALSE)
samples <- subset(samples, !(Cord_Blood_IBC == 101470 & Platform %in% c("HiSeq2500", "HiSeq4000"))) # Remove 101470 Duplicates
samples <- subset(samples, !Sequencing_ID %in% c("JLCM032B", "JLCM050B")) # Remove mislabeled females

# Add in global CpG methylation from bsseq, see 'Chromosome Methylation Analysis.R'
meth_disc <- read.csv("Tables/Chromosome_Meth_CpG_Counts_Discovery.csv", header = TRUE, stringsAsFactors = FALSE)
cov_disc <- read.csv("Tables/Chromosome_Total_CpG_Counts_Discovery.csv", header = TRUE, stringsAsFactors = FALSE)
meth_rep <- read.csv("Tables/Chromosome_Meth_CpG_Counts_Replication.csv", header = TRUE, stringsAsFactors = FALSE)
cov_rep <- read.csv("Tables/Chromosome_Total_CpG_Counts_Replication.csv", header = TRUE, stringsAsFactors = FALSE)
globalMeth <- data.frame(Sample = c(meth_disc$sampleNames.bs., meth_rep$sampleNames.bs.),
                         meth = c(rowSums(meth_disc[,grepl("chr", colnames(meth_disc), fixed = TRUE)]),
                                  rowSums(meth_rep[,grepl("chr", colnames(meth_rep), fixed = TRUE)])),
                         cov = c(rowSums(cov_disc[,grepl("chr", colnames(cov_disc), fixed = TRUE)]),
                                 rowSums(cov_rep[,grepl("chr", colnames(cov_rep), fixed = TRUE)])))
globalMeth <- subset(globalMeth, !Sample %in% c("JLCM032B", "JLCM050B")) # Remove mislabeled females
globalMeth$percent_cpg_meth_bsseq <- globalMeth$meth * 100 / globalMeth$cov
samples <- merge(x = samples, y = globalMeth[,c("Sample", "percent_cpg_meth_bsseq")], by.x = "Sequencing_ID", by.y = "Sample", 
                 all.x = TRUE, all.y = FALSE, sort = FALSE)
# bsseq permeth is on average 0.27% (min = 0.11%, max = 0.35%) higher than bismark permeth
rm(cov_disc, cov_rep, globalMeth, meth_disc, meth_rep)

# Covariates by Diagnosis ####
samples <- samples[,!colnames(samples) == "percent_cpg_meth"] # Remove bismark permeth
samples$Diagnosis_Alg <- factor(samples$Diagnosis_Alg, levels = c("TD", "ASD"))
catVars <- c("Study", "Platform", "Sex", "Site", "MomEdu_detail", "DM1or2", "GDM", "PE", "home_ownership", "marital_status", 
             "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3", "AllEQ_PV_YN_Mo_2", "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", "AllEQ_PV_YN_Mo2",
             "AllEQ_PV_YN_Mo3", "AllEQ_PV_YN_Mo4", "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", 
             "AllEQ_PV_YN_Mo9")
contVars <- colnames(samples)[!colnames(samples) %in% catVars]
(contVars <- contVars[!contVars %in% c("COI_ID", "Sequencing_ID", "Cord_Blood_IBC", "Diagnosis_Alg")])
# [1] "ADOScs"                    "MSLelcStandard36"          "MSLelTscore36"             "MSLfmTscore36"            
# [5] "MSLrlTscore36"             "MSLvrTscore36"             "percent_trimmed"           "percent_aligned"          
# [9] "percent_duplicate"         "dedup_reads_M"             "C_coverage"                "CG_coverage"              
# [13] "percent_chg_meth"          "percent_chh_meth"          "ga_w"                      "bw_g"                     
# [17] "MomAgeYr"                  "Mat_Height_cm"             "Mat_Weight_kg_PrePreg"     "Mat_BMI_PrePreg"          
# [21] "parity"                    "dad_age"                   "cotinine_urine_ngml"       "final_creatinine_mgdl"    
# [25] "AllEQ_tot_All_FA_mcg_Mo_3" "AllEQ_tot_All_FA_mcg_Mo_2" "AllEQ_tot_All_FA_mcg_Mo_1" "AllEQ_tot_All_FA_mcg_Mo1" 
# [29] "AllEQ_tot_All_FA_mcg_Mo2"  "AllEQ_tot_All_FA_mcg_Mo3"  "AllEQ_tot_All_FA_mcg_Mo4"  "AllEQ_tot_All_FA_mcg_Mo5" 
# [33] "AllEQ_tot_All_FA_mcg_Mo6"  "AllEQ_tot_All_FA_mcg_Mo7"  "AllEQ_tot_All_FA_mcg_Mo8"  "AllEQ_tot_All_FA_mcg_Mo9" 
# [37] "percent_cpg_meth_bsseq"
allVars <- c(catVars, contVars)

# HiSeqX10 Only (Discovery)
catVarsX10 <- catVars[!catVars %in% c("Platform")]
contVarsX10 <- contVars
allVarsX10 <- c(catVarsX10, contVarsX10)
x10 <- subset(samples, Platform == "HiSeqX10")

catStats <- catFisher(variables = catVarsX10, sampleData = x10, diagnosis = "Diagnosis_Alg", numHypotheses = length(allVarsX10))
catStats$Variable[catStats$pvalue < 0.05] %>% as.character %>% unique # "home_ownership"
write.table(catStats, "Tables/Categorical Covariate Stats by Diagnosis HiSeqX10 only.txt", sep="\t", quote = FALSE, row.names = FALSE)

contStats <- contANOVA(variables = contVarsX10, sampleData = x10, diagnosis = "Diagnosis_Alg", numHypotheses = length(allVarsX10))
contStats$Variable[contStats$pvalue < 0.05] %>% as.character %>% unique
# [1] "ADOScs"                "MSLelcStandard36"      "MSLelTscore36"         "MSLfmTscore36"         "MSLrlTscore36"         "MSLvrTscore36"        
# [7] "cotinine_urine_ngml"
write.table(contStats, "Tables/Continuous Covariate Stats by Diagnosis HiSeqX10 only.txt", sep="\t", quote = FALSE, row.names = FALSE)

# HiSeq4000 Only (Replication)
catVars4000 <- catVars[!catVars %in% c("Platform", "Site", "Study")]
contVars4000 <- contVars
allVars4000 <- c(catVars4000, contVars4000)
hs4000 <- subset(samples, Platform == "HiSeq4000")

catStats <- catFisher(variables = catVars4000, sampleData = hs4000, diagnosis = "Diagnosis_Alg", numHypotheses = length(allVars4000))
catStats$Variable[catStats$pvalue < 0.05] %>% as.character %>% unique # "GDM"
write.table(catStats, "Tables/Categorical Covariate Stats by Diagnosis HiSeq4000 only.txt", sep="\t", quote = FALSE, row.names = FALSE)

contStats <- contANOVA(variables = contVars4000, sampleData = hs4000, diagnosis = "Diagnosis_Alg", numHypotheses = length(allVars4000))
contStats$Variable[contStats$pvalue < 0.05] %>% as.character %>% unique
# [1] "ADOScs"                 "MSLelcStandard36"       "MSLelTscore36"          "MSLfmTscore36"          "MSLrlTscore36"         
# [6] "MSLvrTscore36"          "ga_w"                   "percent_cpg_meth_bsseq"
write.table(contStats, "Tables/Continuous Covariate Stats by Diagnosis HiSeq4000 only.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Pooled
catStats <- catFisher(variables = catVars, sampleData = samples, diagnosis = "Diagnosis_Alg", numHypotheses = length(allVars))
catStats$Variable[catStats$pvalue < 0.05] %>% as.character %>% unique # "MomEdu_detail" "home_ownership"    "SmokeYN_Pregnancy"
write.table(catStats, "Tables/Categorical Covariate Stats by Diagnosis Pooled.txt", sep="\t", quote = FALSE, row.names = FALSE)

contStats <- contANOVA(variables = contVars, sampleData = samples, diagnosis = "Diagnosis_Alg", numHypotheses = length(allVars))
contStats$Variable[contStats$pvalue < 0.05] %>% as.character %>% unique
# [1] "ADOScs"                 "MSLelcStandard36"       "MSLelTscore36"          "MSLfmTscore36"          "MSLrlTscore36"         
# [6] "MSLvrTscore36"          "Mat_BMI_PrePreg"        "cotinine_urine_ngml"    "final_creatinine_mgdl"  "percent_cpg_meth_bsseq"
write.table(contStats, "Tables/Continuous Covariate Stats by Diagnosis Pooled.txt", sep="\t", quote = FALSE, row.names = FALSE)

rm(catStats, contStats, hs4000, x10, allVars4000, allVarsX10, catVars4000, catVarsX10, contVars4000, contVarsX10)

# Significant Covariate by Diagnosis Plots ####
# Maternal pre-pregnancy BMI ~ Diagnosis, pooled

ggBoxPlot(data = samples, x = samples$Diagnosis_Alg, y = samples$Mat_BMI_PrePreg, fill = samples$Diagnosis_Alg, 
          xlab = NULL, ylab = "Maternal Pre-pregnancy BMI", file = "Figures/Maternal Prepreg BMI by Diagnosis Pooled Boxplot.png")

# Maternal Cotinine ~ Diagnosis, pooled
ggBoxPlot(data = samples, x = samples$Diagnosis_Alg, y = log10(samples$cotinine_urine_ngml), fill = samples$Diagnosis_Alg, 
          xlab = NULL, ylab = "log(Maternal Cotinine)", file = "Figures/Maternal Cotinine by Diagnosis Pooled Boxplot.png")

# Maternal Cotinine ~ Diagnosis + Prenatal Smoking, pooled
cot_smoke <- samples[,c("COI_ID", "Diagnosis_Alg", "SmokeYN_Pregnancy", "cotinine_urine_ngml")]
cot_smoke <- subset(cot_smoke, !is.na(cot_smoke$SmokeYN_Pregnancy) & !is.na(cot_smoke$cotinine_urine_ngml))
cot_smoke$SmokeYN_Pregnancy2 <- NA
cot_smoke$SmokeYN_Pregnancy2[cot_smoke$SmokeYN_Pregnancy == 0] <- "No"
cot_smoke$SmokeYN_Pregnancy2[cot_smoke$SmokeYN_Pregnancy == 1] <- "Yes"
cot_smoke$SmokeYN_Pregnancy2 <- factor(cot_smoke$SmokeYN_Pregnancy2, levels = c("No", "Yes"))

ggBoxPlot(data = cot_smoke, x = cot_smoke$SmokeYN_Pregnancy2, y = log10(cot_smoke$cotinine_urine_ngml), fill = cot_smoke$Diagnosis_Alg, 
          xlab = "Maternal Prenatal Smoking", ylab = "log(Maternal Cotinine)", 
          file = "Figures/Maternal Cotinine by Diagnosis and Prenatal Smoking Pooled Boxplot.png", axis.title.x = element_text())
rm(cot_smoke)

# Maternal Prenatal Smoking ~ Diagnosis, pooled
mat_smoke <- data.frame("Diagnosis_Alg" = factor(rep(c("TD", "ASD"), each = 2), levels = c("TD", "ASD")),
                        "SmokeYN_Pregnancy" = factor(rep(c("No", "Yes"), 2), levels = c("No", "Yes")),
                        "Count" = c(table(samples$SmokeYN_Pregnancy[samples$Diagnosis_Alg == "TD"]),
                                    table(samples$SmokeYN_Pregnancy[samples$Diagnosis_Alg == "ASD"])))
mat_smoke$percent <- mat_smoke$Count * 100 / c(rep(sum(mat_smoke$Count[mat_smoke$Diagnosis_Alg == "TD"]), 2),
                                               rep(sum(mat_smoke$Count[mat_smoke$Diagnosis_Alg == "ASD"]), 2))
g <- ggplot(data = mat_smoke)
g + 
        geom_col(aes(x = SmokeYN_Pregnancy, y = percent, fill = Diagnosis_Alg), size = 1.1, position = "dodge") +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.87,1.03), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25),
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              plot.margin = unit(c(2,1,1,1), "lines")) +
        scale_y_continuous(breaks=pretty_breaks(n=4), limits = c(0,100), expand = c(0,0.3)) +
        xlab("Maternal Prenatal Smoking") +
        ylab("Subjects (%)") +
        scale_fill_manual(breaks = c("TD", "ASD"), values = c("TD"="#3366CC", "ASD"="#FF3366"))
ggsave(file = "Figures/Maternal Prenatal Smoking by Diagnosis Pooled Barplot.png", dpi = 600, width = 7, height = 7, units = "in")
rm(mat_smoke)

# Home Ownership ~ Diagnosis, pooled
home <- data.frame("Diagnosis_Alg" = factor(rep(c("TD", "ASD"), each = 2), levels = c("TD", "ASD")),
                        "home_ownership" = factor(rep(c("No", "Yes"), 2), levels = c("No", "Yes")),
                        "Count" = c(table(samples$home_ownership[samples$Diagnosis_Alg == "TD"])[1:2],
                                    table(samples$home_ownership[samples$Diagnosis_Alg == "ASD"])[1:2]))
home$percent <- home$Count * 100 / c(rep(sum(home$Count[home$Diagnosis_Alg == "TD"]), 2),
                                     rep(sum(home$Count[home$Diagnosis_Alg == "ASD"]), 2))
g <- ggplot(data = home)
g + 
        geom_col(aes(x = home_ownership, y = percent, fill = Diagnosis_Alg), size = 1.1, position = "dodge") +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.87,1.03), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25),
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              plot.margin = unit(c(2,1,1,1), "lines")) +
        scale_y_continuous(breaks=pretty_breaks(n=4), limits = c(0,100), expand = c(0,0.3)) +
        xlab("Home Ownership") +
        ylab("Subjects (%)") +
        scale_fill_manual(breaks = c("TD", "ASD"), values = c("TD"="#3366CC", "ASD"="#FF3366"))
ggsave(file = "Figures/Home Ownership by Diagnosis Pooled Barplot.png", dpi = 600, width = 7, height = 7, units = "in")
rm(home)

# Maternal Education ~ Diagnosis, pooled
mat_edu <- data.frame("Diagnosis_Alg" = factor(rep(c("TD", "ASD"), each = 8), levels = c("TD", "ASD")),
                   "MomEdu_detail" = factor(rep(1:8, 2), levels = 1:8),
                   "Count" = c(table(samples$MomEdu_detail[samples$Diagnosis_Alg == "TD"])[1],0,
                               table(samples$MomEdu_detail[samples$Diagnosis_Alg == "TD"])[2:7],
                               table(samples$MomEdu_detail[samples$Diagnosis_Alg == "ASD"])))
mat_edu$percent <- mat_edu$Count * 100 / c(rep(sum(mat_edu$Count[mat_edu$Diagnosis_Alg == "TD"]), 8),
                                     rep(sum(mat_edu$Count[mat_edu$Diagnosis_Alg == "ASD"]), 8))
g <- ggplot(data = mat_edu)
g + 
        geom_col(aes(x = MomEdu_detail, y = percent, fill = Diagnosis_Alg), size = 1.1, position = "dodge") +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.87,1.03), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25),
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              plot.margin = unit(c(2,1,1,1), "lines")) +
        scale_y_continuous(breaks=pretty_breaks(n=4), limits = c(0,50), expand = c(0,0.15)) +
        xlab("Maternal Education") +
        ylab("Subjects (%)") +
        scale_fill_manual(breaks = c("TD", "ASD"), values = c("TD"="#3366CC", "ASD"="#FF3366"))
ggsave(file = "Figures/Maternal Education by Diagnosis Pooled Barplot.png", dpi = 600, width = 7, height = 7, units = "in")
rm(mat_edu, g)

# Global Methylation by Covariates ####
samples$Study <- factor(samples$Study, levels = c("MARBLES", "EARLI"))
samples$Platform <- factor(samples$Platform, levels = c("HiSeqX10", "HiSeq4000"))
samples$Sex <- factor(samples$Sex, levels = c("M", "F"))
samples$Site <- factor(samples$Site, levels = c("UC Davis", "Kaiser Permanente", "Drexel", "Johns Hopkins University"))
samples$Diagnosis_Alg <- factor(samples$Diagnosis_Alg, levels = c("TD", "ASD"))
samples$MomEdu_detail <- factor(samples$MomEdu_detail, levels = c(6, 1:5, 7,8))
samples$home_ownership[samples$home_ownership == 99] <- NA
samples$marital_status[samples$marital_status == 99] <- NA
factorCols <- c("DM1or2", "GDM", "PE", "marital_status", "home_ownership", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3", "AllEQ_PV_YN_Mo_2", "AllEQ_PV_YN_Mo_1",
                "AllEQ_PV_YN_Mo1", "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", "AllEQ_PV_YN_Mo4", "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6",
                "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", "AllEQ_PV_YN_Mo9")
samples[,factorCols] <- lapply(samples[,factorCols], as.factor)

# Pooled
catVarsPooled <- c("Diagnosis_Alg", catVars)
contVarsPooled <- contVars[!contVars == "percent_cpg_meth_bsseq"]
statsPooled <- methLm(catVars = catVarsPooled, contVars = contVarsPooled, sampleData = samples, globalMeth = "percent_cpg_meth_bsseq")
write.table(statsPooled, "Tables/Covariate Stats by Global Meth Pooled.txt", sep="\t", quote = FALSE, row.names = FALSE)

# X10 Only, Discovery
catVarsX10 <- catVarsPooled[!catVarsPooled == "Platform"]
contVarsX10 <- contVarsPooled
samplesX10 <- subset(samples, Platform == "HiSeqX10")
statsX10 <- methLm(catVars = catVarsX10, contVars = contVarsX10, sampleData = samplesX10, globalMeth = "percent_cpg_meth_bsseq")
write.table(statsX10, "Tables/Covariate Stats by Global Meth X10.txt", sep="\t", quote = FALSE, row.names = FALSE)

# 4000 Only, Replication
catVars4000 <- catVarsPooled[!catVarsPooled %in% c("Study", "Platform", "Site")]
contVars4000 <- contVarsPooled
samples4000 <- subset(samples, Platform == "HiSeq4000")
stats4000 <- methLm(catVars = catVars4000, contVars = contVars4000, sampleData = samples4000, globalMeth = "percent_cpg_meth_bsseq")
write.table(stats4000, "Tables/Covariate Stats by Global Meth 4000.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Pooled Adjusting for Platform
catVarsPlatform <- catVarsPooled[!catVarsPooled == "Platform"]
contVarsPlatform <- contVarsPooled
samplesPlatform <- samples
statsPlatform <- methLmPlatform(catVars = catVarsPlatform, contVars = contVarsPlatform, sampleData = samplesPlatform, 
                                globalMeth = "percent_cpg_meth_bsseq", platform = "Platform")
write.table(statsPlatform, "Tables/Covariate Stats by Global Meth Pooled Platform Adjusted.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Global Methylation by Covariates and PCR Duplicates ####

# Discovery
contVarsX10Dups <- contVarsX10[!contVarsX10 == "percent_duplicate"]
statsX10Dups <- methLmAdj(catVars = catVarsX10, contVars = contVarsX10Dups, sampleData = samplesX10, 
                          globalMeth = "percent_cpg_meth_bsseq", adjVar1 = "percent_duplicate")
write.table(statsX10Dups, "Tables/Covariate Stats by Global Meth X10 PCR Duplicate Adjusted.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Replication
contVars4000Dups <- contVars4000[!contVars4000 == "percent_duplicate"]
stats4000Dups <- methLmAdj(catVars = catVars4000, contVars = contVars4000Dups, sampleData = samples4000, 
                           globalMeth = "percent_cpg_meth_bsseq", adjVar1 = "percent_duplicate")
write.table(stats4000Dups, "Tables/Covariate Stats by Global Meth 4000 PCR Duplicate Adjusted.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Pooled + Platform
contVarsPlatformDups <- contVarsPlatform[!contVarsPlatform == "percent_duplicate"]
platformStatsDups <- methLmAdj(catVars = catVarsPlatform, contVars = contVarsPlatformDups, sampleData = samplesPlatform, 
                               globalMeth = "percent_cpg_meth_bsseq", adjVar1 = "percent_duplicate", adjVar2 = "Platform")
write.table(platformStatsDups, "Tables/Covariate Stats by Global Meth Pooled, Platform and PCR Duplicate Adjusted.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Global Methylation by Covariates, Males Only ####
# Discovery
samplesX10M <- subset(samplesX10, Sex == "M")
catVarsX10M <- catVarsX10[!catVarsX10 %in% c("Sex", "PE")]
statsX10M <- methLm(catVars = catVarsX10M, contVars = contVarsX10, sampleData = samplesX10M, globalMeth = "percent_cpg_meth_bsseq")
write.table(statsX10M, "Tables/Covariate Stats by Global Meth X10 Males Only.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Replication
samples4000M <- subset(samples4000, Sex == "M")
catVars4000M <- catVars4000[!catVars4000 %in% c("Sex", "DM1or2")]
stats4000M <- methLm(catVars = catVars4000M, contVars = contVars4000, sampleData = samples4000M, globalMeth = "percent_cpg_meth_bsseq")
write.table(stats4000M, "Tables/Covariate Stats by Global Meth 4000 Males Only.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Pooled + Platform
samplesM <- subset(samples, Sex == "M")
catVarsPlatformM <- catVarsPlatform[!catVarsPlatform == "Sex"]
platformStatsM <- methLmPlatform(catVars = catVarsPlatformM, contVars = contVarsPlatform, sampleData = samplesM, 
                                 globalMeth = "percent_cpg_meth_bsseq", platform = "Platform")
write.table(platformStatsM, "Tables/Covariate Stats by Global Meth Pooled Platform Adjusted Males Only.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Global Methylation by Covariates and PCR Duplicates, Males Only ####

# Discovery, Males only
statsX10DupsM <- methLmAdj(catVars = catVarsX10M, contVars = contVarsX10Dups, sampleData = samplesX10M, 
                           globalMeth = "percent_cpg_meth_bsseq", adjVar1 = "percent_duplicate")
write.table(statsX10DupsM, "Tables/Covariate Stats by Global Meth X10 PCR Duplicate Adjusted, Males Only.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Replication, Males only
stats4000DupsM <- methLmAdj(catVars = catVars4000M, contVars = contVars4000Dups, sampleData = samples4000M, 
                            globalMeth = "percent_cpg_meth_bsseq", adjVar1 = "percent_duplicate")
write.table(stats4000DupsM, "Tables/Covariate Stats by Global Meth 4000 PCR Duplicate Adjusted, Males Only.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Pooled + Platform, Males only
platformStatsDupsM <- methLmAdj(catVars = catVarsPlatformM, contVars = contVarsPlatformDups, sampleData = samplesM, 
                                globalMeth = "percent_cpg_meth_bsseq", adjVar1 = "percent_duplicate", adjVar2 = "Platform")
write.table(platformStatsDupsM, "Tables/Covariate Stats by Global Meth Platform PCR Duplicate Adjusted, Males Only.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Global Methylation by Covariates, Females Only ####
# Discovery
samplesX10F <- subset(samplesX10, Sex == "F")
catVarsX10F <- catVarsX10[!catVarsX10 %in% c("Sex")]
statsX10F <- methLm(catVars = catVarsX10F, contVars = contVarsX10, sampleData = samplesX10F, globalMeth = "percent_cpg_meth_bsseq")
write.table(statsX10F, "Tables/Covariate Stats by Global Meth X10 Females Only.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Replication
samples4000F <- subset(samples4000, Sex == "F")
catVars4000F <- catVars4000[!catVars4000 %in% c("Sex", "PE", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3", "AllEQ_PV_YN_Mo_2",
                                                "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1")]
stats4000F <- methLm(catVars = catVars4000F, contVars = contVars4000, sampleData = samples4000F, globalMeth = "percent_cpg_meth_bsseq")
write.table(stats4000F, "Tables/Covariate Stats by Global Meth 4000 Females Only.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Pooled + Platform
samplesF <- subset(samples, Sex == "F")
catVarsPlatformF <- catVarsPlatform[!catVarsPlatform == "Sex"]
platformStatsF <- methLmPlatform(catVars = catVarsPlatformF, contVars = contVarsPlatform, sampleData = samplesF, 
                                 globalMeth = "percent_cpg_meth_bsseq", platform = "Platform")
write.table(platformStatsF, "Tables/Covariate Stats by Global Meth Pooled Platform Adjusted Females Only.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Global Methylation by Covariates and PCR Duplicates, Females Only ####

# Discovery, Females only
statsX10DupsF <- methLmAdj(catVars = catVarsX10F, contVars = contVarsX10Dups, sampleData = samplesX10F, 
                           globalMeth = "percent_cpg_meth_bsseq", adjVar1 = "percent_duplicate")
write.table(statsX10DupsF, "Tables/Covariate Stats by Global Meth X10 PCR Duplicate Adjusted, Females Only.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Replication, Females only
stats4000DupsF <- methLmAdj(catVars = catVars4000F, contVars = contVars4000Dups, sampleData = samples4000F, 
                            globalMeth = "percent_cpg_meth_bsseq", adjVar1 = "percent_duplicate")
write.table(stats4000DupsF, "Tables/Covariate Stats by Global Meth 4000 PCR Duplicate Adjusted, Females Only.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Pooled + Platform, Females only
platformStatsDupsF <- methLmAdj(catVars = catVarsPlatformF, contVars = contVarsPlatformDups, sampleData = samplesF, 
                                globalMeth = "percent_cpg_meth_bsseq", adjVar1 = "percent_duplicate", adjVar2 = "Platform")
write.table(platformStatsDupsF, "Tables/Covariate Stats by Global Meth Platform PCR Duplicate Adjusted, Females Only.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Cleanup
rm(platformStatsF, platformStatsM, samples4000F, 
   samplesF, samplesM, samplesPlatform, samplesX10F, stats4000,
   stats4000F, stats4000M, statsPlatform, statsPooled, statsX10, statsX10F, statsX10M,
   allVars, catVars, catVars4000, catVars4000F, catVars4000M, catVarsPlatform, catVarsPlatformF, catVarsPlatformM, catVarsPooled,
   catVarsX10, catVarsX10F, catVarsX10M, contVars, contVars4000, contVars4000Dups, contVarsPlatform, contVarsPlatformDups, contVarsPooled,
   contVarsX10, contVarsX10Dups, factorCols)

# Global Methylation Box Plots ####
# mCpG ~ Diagnosis + Platform
g <- ggplot(data = samples)
g + 
        geom_boxplot(aes(x=Platform, y=percent_cpg_meth_bsseq, fill = Diagnosis_Alg), size = 1.1) +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.87,1.03), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              plot.margin = unit(c(2,1,1,1), "lines"), axis.title.x = element_blank()) +
        scale_y_continuous(breaks=pretty_breaks(n=4)) +
        ylab("Global CpG Methylation (%)") +
        scale_fill_manual(breaks = c("TD", "ASD"), values = c("TD"="#3366CC", "ASD"="#FF3366"))
ggsave("Figures/Global Methylation by Diagnosis and Platform Boxplot.png", dpi = 600, width = 7, height = 7, units = "in")

# mCpG ~ Sex + Platform
g <- ggplot(data = samples)
g + 
        geom_boxplot(aes(x=Platform, y=percent_cpg_meth_bsseq, fill = Sex), size = 1.1) +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.92,1.03), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              plot.margin = unit(c(2,1,1,1), "lines"), axis.title.x = element_blank()) +
        scale_y_continuous(breaks=pretty_breaks(n=4)) +
        ylab("Global CpG Methylation (%)") +
        scale_fill_manual(breaks = c("M", "F"), values = c("M"="#3366CC", "F"="#FF3366"))
ggsave("Figures/Global Methylation by Sex and Platform Boxplot.png", dpi = 600, width = 7, height = 7, units = "in")

# mCpG ~ Diagnosis + Sex + Platform
g <- ggplot(data = samples)
g + 
        geom_boxplot(aes(x=Sex, y=percent_cpg_meth_bsseq, fill = Diagnosis_Alg), size = 1.1) +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.89,1.12), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              plot.margin = unit(c(2,1,1,1), "lines"), axis.title.x = element_blank(), strip.background = element_blank()) +
        scale_y_continuous(breaks=pretty_breaks(n=4)) +
        ylab("Global CpG Methylation (%)") +
        facet_wrap(vars(Platform)) +
        scale_fill_manual(breaks = c("TD", "ASD"), values = c("TD"="#3366CC", "ASD"="#FF3366"))
ggsave("Figures/Global Methylation by Diagnosis, Sex, and Platform Boxplot.png", dpi = 600, width = 8, height = 7, units = "in")

# Global Methylation Dot Plots ####
# mCpG ~ Diagnosis + Sex + Platform
samplesDot <- samples
samplesDot$Sex <- ifelse(samplesDot$Sex == "M", yes = "Males", no = "Females") %>% factor(levels = c("Males", "Females"))
samplesDot$Set <- ifelse(samplesDot$Platform == "HiSeqX10", yes = "Discovery", no = "Replication") %>% 
        factor(levels = c("Discovery", "Replication"))

g <- ggplot(data = samplesDot)
g + 
        stat_summary(aes(x = Sex, y = percent_cpg_meth_bsseq, group = Diagnosis_Alg), fun.data = "mean_cl_boot", 
                     geom = "crossbar", color = "black", size = 0.7, position = "dodge2") +
        geom_dotplot(aes(x = Sex, y = percent_cpg_meth_bsseq, fill = Diagnosis_Alg, color = Diagnosis_Alg), binwidth = 0.18, 
                     binaxis = "y", stackdir = "center", position = "dodge", stackratio = 1.12, dotsize = 0.75) + 
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.895, 0.96), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              strip.text = element_text(size = 22), plot.margin = unit(c(0.3, 1, 1, 1), "lines"), 
              axis.title.x = element_blank(), strip.background = element_blank()) +
        scale_y_continuous(breaks = pretty_breaks(n = 4)) +
        ylab("Global Methylation (%)") +
        facet_wrap(vars(Set)) +
        scale_fill_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366")) +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366"))
ggsave("Figures/Global Methylation by Diagnosis, Sex, and Platform Dotplot with mean CL.png", dpi = 600, width = 9, height = 7, units = "in")

# Global Methylation ~ PCR Duplicates Scatterplots ####

# Global Methylation ~ Diagnosis + PCR Duplicates, Discovery Samples
ggScatterPlot(x = samplesX10$percent_duplicate, y = samplesX10$percent_cpg_meth_bsseq, groupVar = samplesX10$Diagnosis_Alg, 
              fileName = "Figures/Global Methylation by Diagnosis and PCR Duplicates Discovery Scatterplot.png",
              xlab = "PCR Duplicates (%)", ylab = "Global CpG Methylation (%)", legendPos = c(0.89, 1.035))

# Global Methylation ~ Diagnosis + PCR Duplicates, Replication Samples
ggScatterPlot(x = samples4000$percent_duplicate, y = samples4000$percent_cpg_meth_bsseq, groupVar = samples4000$Diagnosis_Alg, 
              fileName = "Figures/Global Methylation by Diagnosis and PCR Duplicates Replication Scatterplot.png",
              xlab = "PCR Duplicates (%)", ylab = "Global CpG Methylation (%)", legendPos = c(0.89, 1.035))

# Global Methylation ~ Diagnosis + Platform + PCR Duplicates, All Samples
ggScatterPlot(x = samples$percent_duplicate, y = samples$percent_cpg_meth_bsseq, groupVar = samples$Platform, 
              fileName = "Figures/Global Methylation by Platform and PCR Duplicates Scatterplot.png",
              xlab = "PCR Duplicates (%)", ylab = "Global CpG Methylation (%)", legendPos = c(0.78, 1.035))

ggScatterPlot(x = samples$percent_duplicate, y = samples$percent_cpg_meth_bsseq, groupVar = samples$Diagnosis_Alg, 
              fileName = "Figures/Global Methylation by Diagnosis and PCR Duplicates All Samples Scatterplot.png",
              xlab = "PCR Duplicates (%)", ylab = "Global CpG Methylation (%)", legendPos = c(0.89, 1.035))

# Global Methylation ~ Diagnosis + PCR Duplicates, Discovery Samples, Males
ggScatterPlot(x = samplesX10M$percent_duplicate, y = samplesX10M$percent_cpg_meth_bsseq, groupVar = samplesX10M$Diagnosis_Alg, 
              fileName = "Figures/Global Methylation by Diagnosis and PCR Duplicates Discovery Males Scatterplot.png",
              xlab = "PCR Duplicates (%)", ylab = "Global CpG Methylation (%)", legendPos = c(0.89, 1.035))

# Global Methylation ~ Diagnosis + PCR Duplicates, Replication Samples, Males
ggScatterPlot(x = samples4000M$percent_duplicate, y = samples4000M$percent_cpg_meth_bsseq, groupVar = samples4000M$Diagnosis_Alg, 
              fileName = "Figures/Global Methylation by Diagnosis and PCR Duplicates Replication Males Scatterplot.png",
              xlab = "PCR Duplicates (%)", ylab = "Global CpG Methylation (%)", legendPos = c(0.89, 1.035))

# Global Methylation ~ MullenELC Scatterplots ####
# Global Methylation ~ MullenELC + PCR Duplicates, Discovery Samples
ggScatterPlot(x = samplesX10$MSLelcStandard36, y = samplesX10$percent_cpg_meth_bsseq, groupVar = samplesX10$Diagnosis_Alg, 
              fileName = "Figures/Global Methylation by Mullen ELC Discovery Scatterplot.png",
              xlab = "Mullen Early Learning Composite", ylab = "Global CpG Methylation (%)", legendPos = c(0.89, 1.035))

# Global Methylation ~ MullenELC + PCR Duplicates, Replication Samples
ggScatterPlot(x = samples4000$MSLelcStandard36, y = samples4000$percent_cpg_meth_bsseq, groupVar = samples4000$Diagnosis_Alg, 
              fileName = "Figures/Global Methylation by Mullen ELC Replication Scatterplot.png",
              xlab = "Mullen Early Learning Composite", ylab = "Global CpG Methylation (%)", legendPos = c(0.89, 1.035))

# Global Methylation ~ MullenELC + Platform + PCR Duplicates, All Samples
ggScatterPlot(x = samples$MSLelcStandard36, y = samples$percent_cpg_meth_bsseq, groupVar = samples$Diagnosis_Alg, 
              fileName = "Figures/Global Methylation by Mullen ELC and Diagnosis All Samples Scatterplot.png",
              xlab = "Mullen Early Learning Composite", ylab = "Global CpG Methylation (%)", legendPos = c(0.89, 1.035))

ggScatterPlot(x = samples$MSLelcStandard36, y = samples$percent_cpg_meth_bsseq, groupVar = samples$Platform, 
              fileName = "Figures/Global Methylation by Mullen ELC and Platform All Samples Scatterplot.png",
              xlab = "Mullen Early Learning Composite", ylab = "Global CpG Methylation (%)", legendPos = c(0.78, 1.035))

# Global Methylation ~ MullenELC + PCR Duplicates, Discovery Samples, Males
ggScatterPlot(x = samplesX10M$MSLelcStandard36, y = samplesX10M$percent_cpg_meth_bsseq, groupVar = samplesX10M$Diagnosis_Alg, 
              fileName = "Figures/Global Methylation by Mullen ELC Discovery Males Scatterplot.png",
              xlab = "Mullen Early Learning Composite", ylab = "Global CpG Methylation (%)", legendPos = c(0.89, 1.035))

# Global Methylation ~ MullenELC + PCR Duplicates, Replication Samples, Males
ggScatterPlot(x = samples4000M$MSLelcStandard36, y = samples4000M$percent_cpg_meth_bsseq, groupVar = samples4000M$Diagnosis_Alg, 
              fileName = "Figures/Global Methylation by Mullen ELC Replication Males Scatterplot.png",
              xlab = "Mullen Early Learning Composite", ylab = "Global CpG Methylation (%)", legendPos = c(0.89, 1.035))

# Global Methylation ~ MullenELC by Sex and Platform
samples$Set <- ifelse(samples$Platform == "HiSeqX10", yes = "Discovery", no = "Replication") %>% 
        factor(levels = c("Discovery", "Replication"))
samples$Sex <- ifelse(samples$Sex == "M", yes = "Males", no = "Females") %>% factor(levels = c("Males", "Females"))
g <- ggplot(data = samples)
g + 
        geom_smooth(aes(x = MSLelcStandard36, y = percent_cpg_meth_bsseq), method = "lm") +
        geom_point(aes(x = MSLelcStandard36, y = percent_cpg_meth_bsseq, color = Diagnosis_Alg), size = 3) +
        facet_grid(rows = vars(Set), cols = vars(Sex), scales = "free") +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.89, 1.08), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(size = 16, color = "black"), axis.title = element_text(size = 20),
              legend.background = element_blank(), strip.background = element_blank(), 
              plot.margin = unit(c(1,0.25,0.5,0.5), "lines"), strip.text = element_text(size = 20)) +
        scale_x_continuous(breaks = pretty_breaks(n = 5)) +
        scale_y_continuous(breaks = pretty_breaks(n = 5)) +
        xlab("Mullen Early Learning Composite") +
        ylab("Global CpG Methylation (%)") +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("#3366CC", "#FF3366"))
ggsave("Figures/Global Methylation by Mullen ELC split by Sex and Platform.png", dpi = 600, width = 8, height = 8, units = "in")

# Global Methylation ~ ADOS Scatterplots ####
# Global Methylation ~ ADOS + PCR Duplicates, Discovery Samples
ggScatterPlot(x = samplesX10$ADOScs, y = samplesX10$percent_cpg_meth_bsseq, groupVar = samplesX10$Diagnosis_Alg, 
              fileName = "Figures/Global Methylation by ADOS Discovery Scatterplot.png",
              xlab = "ADOS Comparison Score", ylab = "Global CpG Methylation (%)", legendPos = c(0.89, 1.035))

# Global Methylation ~ ADOS + PCR Duplicates, Replication Samples
ggScatterPlot(x = samples4000$ADOScs, y = samples4000$percent_cpg_meth_bsseq, groupVar = samples4000$Diagnosis_Alg, 
              fileName = "Figures/Global Methylation by ADOS Replication Scatterplot.png",
              xlab = "ADOS Comparison Score", ylab = "Global CpG Methylation (%)", legendPos = c(0.89, 1.035))

# Global Methylation ~ ADOS + PCR Duplicates, Discovery Samples, Males
ggScatterPlot(x = samplesX10M$ADOScs, y = samplesX10M$percent_cpg_meth_bsseq, groupVar = samplesX10M$Diagnosis_Alg, 
              fileName = "Figures/Global Methylation by ADOS Discovery Males Scatterplot.png",
              xlab = "ADOS Comparison Score", ylab = "Global CpG Methylation (%)", legendPos = c(0.89, 1.035))

# Global Methylation ~ ADOS + PCR Duplicates, Replication Samples, Males
ggScatterPlot(x = samples4000M$ADOScs, y = samples4000M$percent_cpg_meth_bsseq, groupVar = samples4000M$Diagnosis_Alg, 
              fileName = "Figures/Global Methylation by ADOS Replication Males Scatterplot.png",
              xlab = "ADOS Comparison Score", ylab = "Global CpG Methylation (%)", legendPos = c(0.89, 1.035))

# Global Methylation ~ ADOS by Sex and Platform
g <- ggplot(data = samples)
g + 
        geom_smooth(aes(x = ADOScs, y = percent_cpg_meth_bsseq), method = "lm") +
        geom_point(aes(x = ADOScs, y = percent_cpg_meth_bsseq, color = Diagnosis_Alg), size = 3) +
        facet_grid(rows = vars(Set), cols = vars(Sex), scales = "free") +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.89, 1.08), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(size = 16, color = "black"), axis.title = element_text(size = 20),
              legend.background = element_blank(), strip.background = element_blank(), 
              plot.margin = unit(c(1,0.25,0.5,0.5), "lines"), strip.text = element_text(size = 20)) +
        scale_x_continuous(breaks = pretty_breaks(n = 5)) +
        scale_y_continuous(breaks = pretty_breaks(n = 5)) +
        xlab("ADOS Comparison Score") +
        ylab("Global CpG Methylation (%)") +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("#3366CC", "#FF3366"))
ggsave("Figures/Global Methylation by ADOS split by Sex and Platform.png", dpi = 600, width = 8, height = 8, units = "in")

# Global Methylation ~ Covariates Heatmap ####
# Plot effect size for discovery, replication, males, and females
heatStats <- rbind(statsX10DupsM, stats4000DupsM, platformStatsDupsM, statsX10DupsF, stats4000DupsF, platformStatsDupsF)
heatStats$Set <- c(rep("Discovery", nrow(statsX10DupsM)), 
                   rep("Replication", nrow(stats4000DupsM)), 
                   rep("Pooled", nrow(platformStatsDupsM)),
                   rep("Discovery", nrow(statsX10DupsF)),
                   rep("Replication", nrow(stats4000DupsF)),
                   rep("Pooled", nrow(platformStatsDupsF))) %>%
        factor(levels = c("Discovery", "Replication", "Pooled"))
heatStats$Sex <- c(rep("Males", (nrow(statsX10DupsM) + nrow(stats4000DupsM) + nrow(platformStatsDupsM))),
                   rep("Females", (nrow(statsX10DupsF) + nrow(stats4000DupsF) + nrow(platformStatsDupsF)))) %>%
        factor(levels = c("Males", "Females"))
vars <- intersect(statsX10DupsM$Variable, statsX10DupsF$Variable) %>% 
        intersect(stats4000DupsM$Variable) %>% intersect(stats4000DupsF$Variable) %>% 
        intersect(platformStatsDupsM$Variable) %>% intersect(platformStatsDupsF$Variable) %>% unique
vars <- vars[!grepl("AllEQ", vars, fixed = TRUE) & !vars %in% c("MomEdu_detail", "final_creatinine_mgdl")] # Remove folate variables
heatStats <- subset(heatStats, Variable %in% vars)
varNames <- as.character(heatStats$Variable) %>% 
        str_replace_all(c("Diagnosis_Alg" = "Diagnosis", "GDM" = "Gestational Diabetes", "home_ownership" = "Own Home",
                          "marital_status" = "Married", "ADOScs" = "ADOS", "MSLelcStandard36" = "Mullen Composite",
                          "MSLelTscore36" = "Mullen Expressive Language", "MSLfmTscore36" = "Mullen Fine Motor",
                          "MSLrlTscore36" = "Mullen Receptive Language", "MSLvrTscore36" = "Mullen Visual Reception",
                          "percent_trimmed" = "Bases Trimmed", "percent_aligned" = "Aligned Reads", 
                          "dedup_reads_M" = "Unique Reads", "C_coverage" = "C Coverage", "CG_coverage" = "CpG Coverage",
                          "percent_chg_meth" = "CHG Methylation", "percent_chh_meth" = "CHH Methylation", 
                          "ga_w" = "Gestational Age", "bw_g" = "Birthweight", "MomAgeYr" = "Maternal Age",
                          "Mat_Height_cm" = "Maternal Height", "Mat_Weight_kg_PrePreg" = "Maternal Weight",
                          "Mat_BMI_PrePreg" = "Maternal BMI", "parity" = "Parity", "dad_age" = "Paternal Age",
                          "cotinine_urine_ngml" = "Urine Cotinine"))
heatStats$Variable <- factor(varNames, 
                             levels = rev(c("Diagnosis", "ADOS", "Mullen Composite", "Mullen Expressive Language",
                                            "Mullen Fine Motor", "Mullen Receptive Language", "Mullen Visual Reception",
                                            "Gestational Age", "Birthweight", "Paternal Age", "Maternal Age",
                                            "Maternal Height", "Maternal Weight", "Maternal BMI", 
                                            "Gestational Diabetes", "Parity", "Urine Cotinine", "Married", "Own Home",
                                            "Bases Trimmed", "Aligned Reads", "Unique Reads", "C Coverage", 
                                            "CpG Coverage", "CHG Methylation", "CHH Methylation")))
heatStats$Significant <- (heatStats$qvalue < 0.05) %>% factor(levels = c("TRUE", "FALSE"))

g <- ggplot(data = heatStats)
g + 
        geom_tile(aes(x = Set, y = Variable, fill = Estimate, color = Estimate)) + 
        geom_text(aes(x = Set, y = Variable, alpha = Significant, label = "*"), color = "white", size = 15, nudge_y = -0.42) +
        facet_grid(cols = vars(Sex)) +
        scale_fill_gradientn("Estimate", colors = c("#3366CC", "Black", "#FF0000"), values = c(0, 1), na.value = "#FF0000", 
                             limits = c(-1.21, 1.21), breaks = pretty_breaks(n = 3)) +
        scale_color_gradientn("Estimate", colors = c("#3366CC", "Black", "#FF0000"), values = c(0, 1), na.value = "#FF0000", 
                             limits = c(-1.21, 1.21), breaks = pretty_breaks(n = 3)) +
        scale_alpha_manual(breaks = c("TRUE", "FALSE"), values = c(1, 0), guide = FALSE) +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'vertical', legend.position = c(1.13, 0.9), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "Black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_text(size = 24),
              axis.text.x = element_text(color = "Black", angle = 30, hjust = 1), 
              axis.text.y = element_text(color = "Black", size = 20),
              legend.background = element_blank(), 
              strip.text = element_text(size = 26), plot.margin = unit(c(0.5, 8, 1, 1), "lines"), 
              axis.title = element_blank(), strip.background = element_blank())
ggsave("Figures/Global Methylation by Covariates Heatmap.png", dpi = 600, width = 12, height = 9.75, units = "in")

