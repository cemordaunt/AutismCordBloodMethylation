# Global Methylation and Covariate Analysis ####
# Autism Cord Blood Methylation Project
# Charles Mordaunt
# 11/18/18

# Packages ####
library(reshape2)
library(dplyr)
library(scales)
library(ggplot2)
library(ggbiplot)

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

methLmAdj <- function(catVars, contVars, sampleData, globalMeth, adjVar1, adjVar2 = NULL){
        # Analyzes categorical and continuous variables for association with global methylation, adjusting for 1 or more covariates using linear regression
        stats <- NULL
        methDataCol <- as.numeric(sampleData[,globalMeth])
        if(is.null(adjVar2)){
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
        }
        rownames(stats) <- 1:nrow(stats)
        colnames(stats) <- c("Variable", "Term", "Estimate", "StdError", "tvalue", "pvalue")
        stats <- as.data.frame(stats, stringsAsFactors = FALSE)
        stats[,c("Variable")] <- as.factor(stats[,"Variable"])
        stats[,c("Estimate", "StdError", "tvalue", "pvalue")] <- lapply(stats[,c("Estimate", "StdError", "tvalue", "pvalue")], as.numeric)
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

# Questions ####
# Diagnosis
# Stats of Diagnosis with all categorical and continuous variables
#       Is diagnosis confounded by other variables?
#       Platforms pooled and Stratified by platform
#       Variable ~ Diagnosis
#       Categorical: Fisher’s exact test
#       Continuous: One-way ANOVA
# 
# Global Methylation
# Stats with all categorical and continuous variables
#       What variables are associated with global methylation?
#       Platforms pooled, Stratified by platform
#       Global mCpG ~ variable
#       Categorical: One-way ANOVA
#       Continuous: Linear regression
# Global methylation and diagnosis stats and boxplots
#       Is global methylation different by diagnosis?
#       What happens when we adjust for platform and or sex?
#       Global mCpG ~ Diagnosis, Global mCpG ~ Platform, Global mCpG ~ Sex
# 	Global mCpG ~ Diagnosis + Platform, Global mCpG ~ Diagnosis + Sex, Global mCpG ~ Diagnosis + Platform + Sex
# 	Linear regression, set contrasts

# Covariates by Diagnosis ####
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
# [13] "percent_cpg_meth"          "percent_chg_meth"          "percent_chh_meth"          "ga_w"                     
# [17] "bw_g"                      "MomAgeYr"                  "Mat_Height_cm"             "Mat_Weight_kg_PrePreg"    
# [21] "Mat_BMI_PrePreg"           "parity"                    "dad_age"                   "cotinine_urine_ngml"      
# [25] "final_creatinine_mgdl"     "AllEQ_tot_All_FA_mcg_Mo_3" "AllEQ_tot_All_FA_mcg_Mo_2" "AllEQ_tot_All_FA_mcg_Mo_1"
# [29] "AllEQ_tot_All_FA_mcg_Mo1"  "AllEQ_tot_All_FA_mcg_Mo2"  "AllEQ_tot_All_FA_mcg_Mo3"  "AllEQ_tot_All_FA_mcg_Mo4" 
# [33] "AllEQ_tot_All_FA_mcg_Mo5"  "AllEQ_tot_All_FA_mcg_Mo6"  "AllEQ_tot_All_FA_mcg_Mo7"  "AllEQ_tot_All_FA_mcg_Mo8" 
# [37] "AllEQ_tot_All_FA_mcg_Mo9" 
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
# [7] "cotinine_urine_ngml"   "final_creatinine_mgdl"
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
# [1] "ADOScs"           "MSLelcStandard36" "MSLelTscore36"    "MSLfmTscore36"    "MSLrlTscore36"    "MSLvrTscore36"    "percent_cpg_meth"
# [8] "ga_w" 
write.table(contStats, "Tables/Continuous Covariate Stats by Diagnosis HiSeq4000 only.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Pooled
catStats <- catFisher(variables = catVars, sampleData = samples, diagnosis = "Diagnosis_Alg", numHypotheses = length(allVars))
catStats$Variable[catStats$pvalue < 0.05] %>% as.character %>% unique # "MomEdu_detail" "home_ownership"    "SmokeYN_Pregnancy"
write.table(catStats, "Tables/Categorical Covariate Stats by Diagnosis Pooled.txt", sep="\t", quote = FALSE, row.names = FALSE)

contStats <- contANOVA(variables = contVars, sampleData = samples, diagnosis = "Diagnosis_Alg", numHypotheses = length(allVars))
contStats$Variable[contStats$pvalue < 0.05] %>% as.character %>% unique
# [1] "ADOScs"                "MSLelcStandard36"      "MSLelTscore36"         "MSLfmTscore36"         "MSLrlTscore36"         "MSLvrTscore36"        
# [7] "percent_cpg_meth"      "Mat_BMI_PrePreg"       "cotinine_urine_ngml"   "final_creatinine_mgdl"
write.table(contStats, "Tables/Continuous Covariate Stats by Diagnosis Pooled.txt", sep="\t", quote = FALSE, row.names = FALSE)

rm(catStats, contStats, hs4000, x10, allVars4000, allVarsX10, catVars4000, catVarsX10, contVars4000, contVarsX10)

# Significant Covariate by Diagnosis Plots ####
# Maternal pre-pregnancy BMI ~ Diagnosis, pooled

ggBoxPlot(data = samples, x = samples$Diagnosis_Alg, y = samples$Mat_BMI_PrePreg, fill = samples$Diagnosis_Alg, 
          ylab = "Maternal Pre-pregnancy BMI", file = "Figures/Maternal Prepreg BMI by Diagnosis Pooled Boxplot.png")

# Maternal Cotinine ~ Diagnosis, pooled
ggBoxPlot(data = samples, x = samples$Diagnosis_Alg, y = log10(samples$cotinine_urine_ngml), fill = samples$Diagnosis_Alg, 
          ylab = "log(Maternal Cotinine)", file = "Figures/Maternal Cotinine by Diagnosis Pooled Boxplot.png")

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
rm(mat_edu)


# Global Methylation by Covariates ####
samples$Diagnosis_Alg <- factor(samples$Diagnosis_Alg, levels = c("TD", "ASD"))
samples$Platform <- factor(samples$Platform, levels = c("HiSeqX10", "HiSeq4000"))
samples$Study <- factor(samples$Study, levels = c("MARBLES", "EARLI"))
samples$Sex <- factor(samples$Sex, levels = c("M", "F"))
samples$Site <- factor(samples$Site, levels = c("UC Davis", "Kaiser Permanente", "Drexel", "Johns Hopkins University"))
samples$GDM_EQ <- factor(samples$GDM_EQ, levels = c(0, 1))
samples$PE_EQ <- factor(samples$PE_EQ, levels = c(0, 1))
samples$SmokeYN_Pregnancy <- factor(samples$SmokeYN_Pregnancy, levels = c(0, 1))
samples$Supp_mv_mo_1 <- factor(samples$Supp_mv_mo_1, levels = c(0, 1))
samples$Supp_mv_mo1 <- factor(samples$Supp_mv_mo1, levels = c(0, 1))

# Pooled
catVarsPooled <- c("Diagnosis_Alg", catVars)
contVarsPooled <- contVars[!contVars == "percent_cpg_meth"]
pooledStats <- methLm(catVars = catVarsPooled, contVars = contVarsPooled, sampleData = samples, globalMeth = "percent_cpg_meth")
write.table(pooledStats, "Tables/Covariate Stats by Global Meth Pooled.txt", sep="\t", quote = FALSE, row.names = FALSE)

# X10 Only, Discovery
catVarsX10 <- catVarsPooled[!catVarsPooled == "Platform"]
contVarsX10 <- contVarsPooled
samplesX10 <- subset(samples, Platform == "HiSeqX10")
x10stats <- methLm(catVars = catVarsX10, contVars = contVarsX10, sampleData = samplesX10, globalMeth = "percent_cpg_meth")
write.table(x10stats, "Tables/Covariate Stats by Global Meth X10.txt", sep="\t", quote = FALSE, row.names = FALSE)

# 4000 Only, Replication
catVars4000 <- catVarsPooled[!catVarsPooled %in% c("Study", "Platform", "Site")]
contVars4000 <- contVarsPooled
samples4000 <- subset(samples, Platform == "HiSeq4000")
stats4000 <- methLm(catVars = catVars4000, contVars = contVars4000, sampleData = samples4000, globalMeth = "percent_cpg_meth")
write.table(stats4000, "Tables/Covariate Stats by Global Meth 4000.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Pooled Adjusting for Platform
catVarsPlatform <- catVarsPooled[!catVarsPooled == "Platform"]
contVarsPlatform <- contVarsPooled
samplesPlatform <- samples
platformStats <- methLmPlatform(catVars = catVarsPlatform, contVars = contVarsPlatform, sampleData = samplesPlatform, 
                                globalMeth = "percent_cpg_meth", platform = "Platform")
write.table(platformStats, "Tables/Covariate Stats by Global Meth Pooled Platform Adjusted.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Global Methylation Box Plots ####
# Global Methylation ~ Diagnosis, Discovery Samples
summary(lm(percent_cpg_meth ~ Diagnosis_Alg, data = samplesX10))$coefficients
#                    Estimate Std. Error    t value      Pr(>|t|)
# Diagnosis_AlgASD -0.2166847  0.1831888  -1.182849  2.394663e-01

g <- ggplot(data = samplesX10)
g + 
        geom_boxplot(aes(x=Diagnosis_Alg, y=percent_cpg_meth, fill = Diagnosis_Alg), size = 1.1) +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.87,1.03), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              plot.margin = unit(c(2,1,1,1), "lines"), axis.title.x = element_blank()) +
        scale_y_continuous(breaks=pretty_breaks(n=4)) +
        ylab("Global CpG Methylation (%)") +
        scale_fill_manual(breaks = c("TD", "ASD"), values = c("TD"="#3366CC", "ASD"="#FF3366"))
ggsave("Figures/Global Methylation by Diagnosis, Discovery Boxplot.png", dpi = 600, width = 7, height = 7, units = "in")

# Global Methylation ~ Diagnosis, Replication Samples
summary(lm(percent_cpg_meth ~ Diagnosis_Alg, data = samples4000))$coefficients
#                    Estimate Std. Error    t value      Pr(>|t|)
# Diagnosis_AlgASD -0.7387492  0.3316281  -2.227644 3.106585e-02

g <- ggplot(data = samples4000)
g + 
        geom_boxplot(aes(x=Diagnosis_Alg, y=percent_cpg_meth, fill = Diagnosis_Alg), size = 1.1) +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.87,1.03), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              plot.margin = unit(c(2,1,1,1), "lines"), axis.title.x = element_blank()) +
        scale_y_continuous(breaks=pretty_breaks(n=4)) +
        ylab("Global CpG Methylation (%)") +
        scale_fill_manual(breaks = c("TD", "ASD"), values = c("TD"="#3366CC", "ASD"="#FF3366"))
ggsave("Figures/Global Methylation by Diagnosis, Replication Boxplot.png", dpi = 600, width = 7, height = 7, units = "in")

# Global Methylation ~ Diagnosis + Platform, All Samples
summary(lm(percent_cpg_meth ~ Diagnosis_Alg + Platform, data = samples))$coefficients
#                     Estimate Std. Error    t value      Pr(>|t|)
# Diagnosis_AlgASD  -0.3688064  0.1623863  -2.271167  2.453165e-02
# PlatformHiSeq4000 -1.6966585  0.1780024  -9.531659  3.406750e-17

g <- ggplot(data = samples)
g + 
        geom_boxplot(aes(x=Platform, y=percent_cpg_meth, fill = Diagnosis_Alg), size = 1.1) +
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

# Global Methylation Scatterplots ####

# Global Methylation ~ Diagnosis + PCR Duplicates, Discovery Samples
summary(lm(percent_cpg_meth ~ Diagnosis_Alg + percent_duplicate, data = samplesX10))$coefficients
#                     Estimate  Std. Error    t value      Pr(>|t|)
# Diagnosis_AlgASD  -0.2332074 0.181109791  -1.287657  2.006433e-01
# percent_duplicate -0.0131507 0.006799923  -1.933948  5.576102e-02

ggScatterPlot(x = samplesX10$percent_duplicate, y = samplesX10$percent_cpg_meth, groupVar = samplesX10$Diagnosis_Alg, 
              fileName = "Figures/Global Methylation by Diagnosis and PCR Duplicates Discovery Scatterplot.png",
              xlab = "PCR Duplicates (%)", ylab = "Global CpG Methylation (%)", legendPos = c(0.89, 1.035))

# Global Methylation ~ Diagnosis + PCR Duplicates, Replication Samples
summary(lm(percent_cpg_meth ~ Diagnosis_Alg + percent_duplicate, data = samples4000))$coefficients
#                     Estimate  Std. Error    t value      Pr(>|t|)
# Diagnosis_AlgASD  -0.70403673 0.29785097  -2.363721 2.268113e-02
# percent_duplicate -0.06634262 0.01947115  -3.407226 1.434686e-03

ggScatterPlot(x = samples4000$percent_duplicate, y = samples4000$percent_cpg_meth, groupVar = samples4000$Diagnosis_Alg, 
              fileName = "Figures/Global Methylation by Diagnosis and PCR Duplicates Replication Scatterplot.png",
              xlab = "PCR Duplicates (%)", ylab = "Global CpG Methylation (%)", legendPos = c(0.89, 1.035))

# Global Methylation ~ Diagnosis + Platform + PCR Duplicates, All Samples
summary(lm(percent_cpg_meth ~ Diagnosis_Alg + Platform + percent_duplicate, data = samples))$coefficients
#                     Estimate  Std. Error    t value      Pr(>|t|)
# Diagnosis_AlgASD  -0.38343343 0.158405214  -2.420586  1.667457e-02
# PlatformHiSeq4000 -1.55984607 0.179485718  -8.690642  5.368313e-15
# percent_duplicate -0.01982391 0.006629877  -2.990087  3.254209e-03

ggScatterPlot(x = samples$percent_duplicate, y = samples$percent_cpg_meth, groupVar = samples$Platform, 
              fileName = "Figures/Global Methylation by Platform and PCR Duplicates Scatterplot.png",
              xlab = "PCR Duplicates (%)", ylab = "Global CpG Methylation (%)", legendPos = c(0.78, 1.035))

ggScatterPlot(x = samples$percent_duplicate, y = samples$percent_cpg_meth, groupVar = samples$Diagnosis_Alg, 
              fileName = "Figures/Global Methylation by Diagnosis and PCR Duplicates All Samples Scatterplot.png",
              xlab = "PCR Duplicates (%)", ylab = "Global CpG Methylation (%)", legendPos = c(0.89, 1.035))

# Global Methylation ~ MullenELC + PCR Duplicates, Discovery Samples
summary(lm(percent_cpg_meth ~ MSLelcStandard36 + percent_duplicate, data = samplesX10))$coefficients
#                     Estimate  Std. Error    t value      Pr(>|t|)
# MSLelcStandard36   0.004681534 0.003849800   1.216046  2.267741e-01
# percent_duplicate -0.013071049 0.006476512  -2.018224  4.619326e-02

ggScatterPlot(x = samplesX10$MSLelcStandard36, y = samplesX10$percent_cpg_meth, groupVar = samplesX10$Diagnosis_Alg, 
              fileName = "Figures/Global Methylation by Mullen ELC Discovery Scatterplot.png",
              xlab = "Mullen Early Learning Composite", ylab = "Global CpG Methylation (%)", legendPos = c(0.89, 1.035))

# Global Methylation ~ MullenELC + PCR Duplicates, Replication Samples
summary(lm(percent_cpg_meth ~ MSLelcStandard36 + percent_duplicate, data = samples4000))$coefficients
#                     Estimate  Std. Error    t value      Pr(>|t|)
# MSLelcStandard36   0.02302187 0.006090009  3.780268 4.887965e-04
# percent_duplicate -0.07267041 0.022457760 -3.235871 2.367658e-03

ggScatterPlot(x = samples4000$MSLelcStandard36, y = samples4000$percent_cpg_meth, groupVar = samples4000$Diagnosis_Alg, 
              fileName = "Figures/Global Methylation by Mullen ELC Replication Scatterplot.png",
              xlab = "Mullen Early Learning Composite", ylab = "Global CpG Methylation (%)", legendPos = c(0.89, 1.035))

# Global Methylation ~ MullenELC + Platform + PCR Duplicates, All Samples
summary(lm(percent_cpg_meth ~ MSLelcStandard36 + Platform + percent_duplicate, data = samples))$coefficients
#                     Estimate  Std. Error    t value      Pr(>|t|)
# MSLelcStandard36   0.009995785 0.003375352   2.961405  3.575668e-03
# PlatformHiSeq4000 -1.565851597 0.172691790  -9.067319  7.483491e-16
# percent_duplicate -0.017718184 0.006500379  -2.725716  7.201553e-03

ggScatterPlot(x = samples$MSLelcStandard36, y = samples$percent_cpg_meth, groupVar = samples$Diagnosis_Alg, 
              fileName = "Figures/Global Methylation by Mullen ELC and Diagnosis All Samples Scatterplot.png",
              xlab = "Mullen Early Learning Composite", ylab = "Global CpG Methylation (%)", legendPos = c(0.89, 1.035))

ggScatterPlot(x = samples$MSLelcStandard36, y = samples$percent_cpg_meth, groupVar = samples$Platform, 
              fileName = "Figures/Global Methylation by Mullen ELC and Platform All Samples Scatterplot.png",
              xlab = "Mullen Early Learning Composite", ylab = "Global CpG Methylation (%)", legendPos = c(0.78, 1.035))

# Global Methylation by Covariates, Males Only ####
# Discovery
samplesX10M <- subset(samplesX10, Sex == "M")
catVarsX10M <- catVarsX10[!catVarsX10 %in% c("Sex", "PE_EQ")]
x10statsM <- methLm(catVars = catVarsX10M, contVars = contVarsX10, sampleData = samplesX10M, globalMeth = "percent_cpg_meth")
write.table(x10statsM, "Tables/Covariate Stats by Global Meth X10 Males Only.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Replication
samples4000M <- subset(samples4000, Sex == "M")
catVars4000M <- catVars4000[!catVars4000 %in% c("Sex")]
stats4000M <- methLm(catVars = catVars4000M, contVars = contVars4000, sampleData = samples4000M, globalMeth = "percent_cpg_meth")
write.table(stats4000M, "Tables/Covariate Stats by Global Meth 4000 Males Only.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Pooled
samplesM <- subset(samples, Sex == "M")
catVarsPooledM <- catVarsPooled[!catVarsPooled %in% c("Sex")]
pooledStatsM <- methLm(catVars = catVarsPooledM, contVars = contVarsPooled, sampleData = samplesM, globalMeth = "percent_cpg_meth")
write.table(pooledStatsM, "Tables/Covariate Stats by Global Meth Pooled Males Only.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Pooled + Platform
catVarsPlatformM <- catVarsPlatform[!catVarsPlatform == "Sex"]
platformStatsM <- methLmPlatform(catVars = catVarsPlatformM, contVars = contVarsPlatform, sampleData = samplesM, 
                                globalMeth = "percent_cpg_meth", platform = "Platform")
write.table(platformStatsM, "Tables/Covariate Stats by Global Meth Pooled Platform Adjusted Males Only.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Global Methylation by Covariates and PCR Duplicates ####

# Discovery
contVarsX10Dups <- contVarsX10[!contVarsX10 == "percent_duplicate"]
x10statsDups <- methLmAdj(catVars = catVarsX10, contVars = contVarsX10Dups, sampleData = samplesX10, 
                          globalMeth = "percent_cpg_meth", adjVar1 = "percent_duplicate")
write.table(x10statsDups, "Tables/Covariate Stats by Global Meth X10 PCR Duplicate Adjusted.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Discovery, Males only
x10statsDupsM <- methLmAdj(catVars = catVarsX10M, contVars = contVarsX10Dups, sampleData = samplesX10M, 
                          globalMeth = "percent_cpg_meth", adjVar1 = "percent_duplicate")
write.table(x10statsDupsM, "Tables/Covariate Stats by Global Meth X10 PCR Duplicate Adjusted, Males Only.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Replication
contVars4000Dups <- contVars4000[!contVars4000 == "percent_duplicate"]
stats4000Dups <- methLmAdj(catVars = catVars4000, contVars = contVars4000Dups, sampleData = samples4000, 
                          globalMeth = "percent_cpg_meth", adjVar1 = "percent_duplicate")
write.table(stats4000Dups, "Tables/Covariate Stats by Global Meth 4000 PCR Duplicate Adjusted.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Replication, Males only
stats4000DupsM <- methLmAdj(catVars = catVars4000M, contVars = contVars4000Dups, sampleData = samples4000M, 
                           globalMeth = "percent_cpg_meth", adjVar1 = "percent_duplicate")
write.table(stats4000DupsM, "Tables/Covariate Stats by Global Meth 4000 PCR Duplicate Adjusted, Males Only.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Pooled + Platform
contVarsPlatformDups <- contVarsPlatform[!contVarsPlatform == "percent_duplicate"]
platformStatsDups <- methLmAdj(catVars = catVarsPlatform, contVars = contVarsPlatformDups, sampleData = samplesPlatform, 
                           globalMeth = "percent_cpg_meth", adjVar1 = "percent_duplicate", adjVar2 = "Platform")
write.table(platformStatsDups, "Tables/Covariate Stats by Global Meth Pooled, Platform and PCR Duplicate Adjusted.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Pooled + Platform, Males only
platformStatsDupsM <- methLmAdj(catVars = catVarsPlatformM, contVars = contVarsPlatformDups, sampleData = samplesM, 
                            globalMeth = "percent_cpg_meth", adjVar1 = "percent_duplicate", adjVar2 = "Platform")
write.table(platformStatsDupsM, "Tables/Covariate Stats by Global Meth Platform PCR Duplicate Adjusted, Males Only.txt", sep="\t", quote = FALSE, row.names = FALSE)

