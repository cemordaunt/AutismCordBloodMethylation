# Global Methylation and Covariate Analysis ####
# Autism Cord Blood Methylation Project
# Charles Mordaunt
# 10/3/18

# Packages ####
library(reshape2)
library(dplyr)
library(scales)
library(ggplot2)
library(ggbiplot)

# Functions ####
catFisher <- function(variables, sampleData, diagnosis){
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
        return(stats)
}

contANOVA <- function(variables, sampleData, diagnosis){
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
        return(stats)
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

PCAcatANOVA <- function(variables, sampleData, PCdata){
        # Analyzes table of principal components for differences by categorical covariates using one-way ANOVA
        stats <- NULL
        PCtemp <- NULL
        PCs <- colnames(PCdata)[grep("PC", colnames(PCdata))]
        for(i in 1:length(PCs)){
                varTemp <- NULL
                for(j in 1:length(variables)){
                        means <- aggregate(PCdata[,PCs[i]] ~ sampleData[,variables[j]], FUN=mean)[,2]
                        sds <- aggregate(PCdata[,PCs[i]] ~ sampleData[,variables[j]], FUN=sd)[,2]
                        p <- summary(aov(PCdata[,PCs[i]] ~ sampleData[,variables[j]]))[[1]][1, "Pr(>F)"]
                        temp <- cbind(rep(PCs[i], length(means)), rep(variables[j], length(means)), levels(sampleData[,variables[j]]), means, sds, rep(p, length(means)))
                        varTemp <- rbind(varTemp, temp)
                }
                PCtemp <- rbind(PCtemp, varTemp)
        }
        stats <- rbind(stats, PCtemp)
        colnames(stats) <- c("PC", "Variable", "Value", "Mean", "SD", "pvalue")
        stats <- as.data.frame(stats, stringsAsFactors = FALSE)
        stats[,c("PC", "Variable", "Value")] <- lapply(stats[,c("PC", "Variable", "Value")], as.factor)
        stats[,c("Mean", "SD", "pvalue")] <- lapply(stats[,c("Mean", "SD", "pvalue")], as.numeric)
        return(stats)
}

PCAcontLm <- function(variables, sampleData, PCdata){
        # Analyzes table of principal components for differences by continuous covariates using linear regression
        stats <- NULL
        PCtemp <- NULL
        PCs <- colnames(PCdata)[grep("PC", colnames(PCdata))]
        for(i in 1:length(PCs)){
                PCdataCol <- as.numeric(PCdata[,PCs[i]])
                varTemp <- NULL
                for(j in 1:length(variables)){
                        sampleDataCol <- as.numeric(sampleData[,variables[j]])
                        temp <- c(PCs[i], variables[j], summary(lm(PCdataCol ~ sampleDataCol))$coefficients["sampleDataCol",], use.names = FALSE)
                        varTemp <- rbind(varTemp, temp)
                }
                PCtemp <- rbind(PCtemp, varTemp)
        }
        stats <- rbind(stats, PCtemp)
        rownames(stats) <- 1:nrow(stats)
        colnames(stats) <- c("PC", "Variable", "Effect", "SD", "tstat", "pvalue")
        stats <- as.data.frame(stats, stringsAsFactors = FALSE)
        stats[,c("PC", "Variable")] <- lapply(stats[,c("PC", "Variable")], as.factor)
        stats[,c("Effect", "SD", "tstat", "pvalue")] <- lapply(stats[,c("Effect", "SD", "tstat", "pvalue")], as.numeric)
        return(stats)
}

PCAggbiplot <- function(data.pca, groups, pc = 1:2, fileName, xlim = NULL, ylim = NULL, breaks = c("TD", "ASD"), values = c("TD"="#3366CC", "ASD"="#FF3366")){
        # Plots principal components colored by grouping variable and writes file
        pc1 <- summary(data.pca)$importance["Proportion of Variance", paste("PC", pc[1], sep = "")]*100
        pc2 <- summary(data.pca)$importance["Proportion of Variance", paste("PC", pc[2], sep = "")]*100
        g <- ggbiplot(data.pca, obs.scale = 1, var.scale = 1, groups = groups, ellipse = TRUE, circle = FALSE, 
                      var.axes = FALSE, varname.abbrev = FALSE, choices = pc, ellipse.prob = 0.95)
        g + 
                theme_bw(base_size = 25) +
                theme(legend.direction = 'vertical', legend.position = c(0.85, 0.92), panel.grid.major = element_blank(), 
                      panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
                      legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
                      axis.text = element_text(color = "black"), legend.background = element_blank(), legend.spacing.x = unit(0.5,"lines")) +
                coord_cartesian(xlim = xlim, ylim = ylim) +
                xlab(paste("PC", pc[1], " (", round(pc1,1), "% of Variance)", sep = "")) +
                ylab(paste("PC", pc[2], " (", round(pc2,1), "% of Variance)", sep = "")) +
                scale_color_manual(breaks = breaks , values = values) +
                scale_x_continuous(breaks=pretty_breaks(6)) +
                scale_y_continuous(breaks=pretty_breaks(6)) +
                geom_point(aes(color = groups), size=3)
        ggsave(fileName, dpi = 600, width = 8, height = 8, units = "in")
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

windowANOVA <- function(methData, sampleData, formula, diagnosis){
        # Analyzes windowed methylation data for differences by diagnosis with covariates
        sampleData$meth <- methData
        means <- aggregate(sampleData$meth ~ sampleData[,diagnosis], FUN=mean)[,2]
        sds <- aggregate(sampleData$meth ~ sampleData[,diagnosis], FUN=sd)[,2]
        test <- summary(lm(formula, data = sampleData))$coefficients
        temp <- c(means, sds, test["Diagnosis_AlgASD",], test["PlatformHiSeqX10",], test["SexF",], 
                  use.names = FALSE)
        return(temp)
}

# Data ####
samples <- read.delim(file = "Samples/Merged Cord Blood WGBS Database.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
permeth <- read.delim(file = "Tables/windows_10kb_methylation_ASD_CordBlood2.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# 156 samples x 263618 windows

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
# 
# 10kb Window Methylation (PCA)
# PC1-3 stats with all categorical and continuous variables
# 	What variables are associated with methylation overall?
# 	What variables are associated with specific methylation PCs?
# 	PC ~ variable
# 	Linear regression, set contrasts for categorical variables
# Significant variable PC scatterplots
# 	How do samples cluster by PCs and variables?
# 	Categorical: Diagnosis, Platform, Sex
# 	Continuous: Mullen Composite, Global mCpG, gestational age?
# 	Other significant variables?

# Covariates by Diagnosis ####
samples$Diagnosis_Alg <- factor(samples$Diagnosis_Alg, levels = c("TD", "ASD"))
catVars <- c("Study", "Platform", "Sex", "Site", "GDM_EQ", "PE_EQ", "SmokeYN_Pregnancy", "Supp_mv_mo_1", "Supp_mv_mo1")
contVars <- c(colnames(samples)[14:ncol(samples)])
# [1] "MSLvrTscore36"         "MSLfmTscore36"         "MSLrlTscore36"         "MSLelTscore36"         "MSLelcStandard36"     
# [6] "ADOScs"                "percent_trimmed"       "percent_aligned"       "percent_duplicate"     "dedup_reads"          
# [11] "C_coverage"            "CG_coverage"           "percent_chg_meth"      "percent_chh_meth"      "percent_cpg_meth"     
# [16] "MomAgeYr"              "Mat_Height_cm"         "Mat_Weight_kg_PrePreg" "Mat_BMI_PrePreg"   
allVars <- c(catVars, contVars)

# Pooled
catStats <- catFisher(variables = catVars, sampleData = samples, diagnosis = "Diagnosis_Alg")
contStats <- contANOVA(variables = contVars, sampleData = samples, diagnosis = "Diagnosis_Alg")
catStats$qvalue <- p.adjust(catStats$pvalue, method = "fdr", n = length(allVars))
cont_padj <- unique(contStats[,c("Variable", "pvalue")])
cont_padj$qvalue <- p.adjust(cont_padj$pvalue, method = "fdr", n = length(allVars))
contStats$qvalue <- cont_padj$qvalue[match(contStats$pvalue, cont_padj$pvalue)]
catStats$Variable[catStats$pvalue < 0.05] %>% as.character %>% unique #  "SmokeYN_Pregnancy"
contStats$Variable[contStats$pvalue < 0.05] %>% as.character %>% unique
# "MSLvrTscore36"    "MSLfmTscore36"    "MSLrlTscore36"    "MSLelTscore36"    "MSLelcStandard36" 
# "ADOScs"           "percent_cpg_meth"
write.table(catStats, "Tables/Categorical Covariate Stats by Diagnosis Pooled.txt", sep="\t", quote = FALSE, row.names = FALSE)
write.table(contStats, "Tables/Continuous Covariate Stats by Diagnosis Pooled.txt", sep="\t", quote = FALSE, row.names = FALSE)

# HiSeqX10 Only 
catVarsX10 <- catVars[!catVars %in% c("Platform")]
x10 <- subset(samples, Platform == "HiSeqX10")
catStats <- catFisher(variables = catVarsX10, sampleData = x10, diagnosis = "Diagnosis_Alg")
contStats <- contANOVA(variables = contVars, sampleData = x10, diagnosis = "Diagnosis_Alg")
catStats$qvalue <- p.adjust(catStats$pvalue, method = "fdr", n = length(allVars))
cont_padj <- unique(contStats[,c("Variable", "pvalue")])
cont_padj$qvalue <- p.adjust(cont_padj$pvalue, method = "fdr", n = length(allVars))
contStats$qvalue <- cont_padj$qvalue[match(contStats$pvalue, cont_padj$pvalue)]
catStats$Variable[catStats$pvalue < 0.05] %>% as.character %>% unique # None
contStats$Variable[contStats$pvalue < 0.05] %>% as.character %>% unique
# "MSLvrTscore36"    "MSLfmTscore36"    "MSLrlTscore36"    "MSLelTscore36"    "MSLelcStandard36" "ADOScs"
write.table(catStats, "Tables/Categorical Covariate Stats by Diagnosis HiSeqX10 only.txt", sep="\t", quote = FALSE, row.names = FALSE)
write.table(contStats, "Tables/Continuous Covariate Stats by Diagnosis HiSeqX10 only.txt", sep="\t", quote = FALSE, row.names = FALSE)

# HiSeq4000 Only 
catVars4000 <- catVars[!catVars %in% c("Platform", "Study", "Site")]
hs4000 <- subset(samples, Platform == "HiSeq4000")
catStats <- catFisher(variables = catVars4000, sampleData = hs4000, diagnosis = "Diagnosis_Alg")
contStats <- contANOVA(variables = contVars, sampleData = hs4000, diagnosis = "Diagnosis_Alg")
catStats$qvalue <- p.adjust(catStats$pvalue, method = "fdr", n = length(allVars))
cont_padj <- unique(contStats[,c("Variable", "pvalue")])
cont_padj$qvalue <- p.adjust(cont_padj$pvalue, method = "fdr", n = length(allVars))
contStats$qvalue <- cont_padj$qvalue[match(contStats$pvalue, cont_padj$pvalue)]
catStats$Variable[catStats$pvalue < 0.05] %>% as.character %>% unique # None
contStats$Variable[contStats$pvalue < 0.05] %>% as.character %>% unique
# "MSLvrTscore36"    "MSLfmTscore36"    "MSLrlTscore36"    "MSLelTscore36"    "MSLelcStandard36" 
# "ADOScs"      "percent_cpg_meth"        
write.table(catStats, "Tables/Categorical Covariate Stats by Diagnosis HiSeq4000 only.txt", sep="\t", quote = FALSE, row.names = FALSE)
write.table(contStats, "Tables/Continuous Covariate Stats by Diagnosis HiSeq4000 only.txt", sep="\t", quote = FALSE, row.names = FALSE)
rm(catStats, cont_padj, contStats)

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

# X10 Only
catVarsX10 <- catVarsPooled[!catVarsPooled == "Platform"]
contVarsX10 <- contVarsPooled
samplesX10 <- subset(samples, Platform == "HiSeqX10")
x10stats <- methLm(catVars = catVarsX10, contVars = contVarsX10, sampleData = samplesX10, globalMeth = "percent_cpg_meth")
write.table(x10stats, "Tables/Covariate Stats by Global Meth X10.txt", sep="\t", quote = FALSE, row.names = FALSE)

# 4000 Only
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
methLmAdj <- function(catVars, contVars, sampleData, globalMeth, adjVar1, adjVar2 = NULL){
        # Analyzes categorical and continuous variables for association with global methylation, adjusting for platform using linear regression
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

# Start here

# Window Methylation by Diagnosis ####
# PCA
table(samples$Sequencing_ID == colnames(permeth)[4:ncol(permeth)]) # All TRUE
table(is.na(permeth)) # All FALSE
data <- t(as.matrix(permeth[,4:ncol(permeth)]))
data.pca <- prcomp(data, center = TRUE, scale. = TRUE) 
ggbiplotPCA(data.pca = data.pca, groups = samples$Diagnosis_Alg, pc = c(1,2),
            fileName = "Figures/ASD Cord Window Methylation by Diagnosis PC1 PC2 Plot.png", xlim = c(-750,550), ylim = c(-650,650))
ggbiplotPCA(data.pca = data.pca, groups = samples$Diagnosis_Alg, pc = c(1,3),
            fileName = "Figures/ASD Cord Window Methylation by Diagnosis PC1 PC3 Plot.png", xlim = c(-750,550), ylim = c(-650,650))
ggbiplotPCA(data.pca = data.pca, groups = samples$Diagnosis_Alg, pc = c(2,3),
            fileName = "Figures/ASD Cord Window Methylation by Diagnosis PC2 PC3 Plot.png", xlim = c(-650,650), ylim = c(-650,650))

# Get PCs 1-3
PCs <- data.pca$x[,c("PC1", "PC2", "PC3")] %>% as.data.frame
PCs$Sequencing_ID <- rownames(PCs)
rownames(PCs) <- 1:nrow(PCs)
PCs <- PCs[,c("Sequencing_ID", "PC1", "PC2", "PC3")]

# Methylation PCs by Covariate Stats ####
# Categorical Covariate Stats with 1-way ANOVA
catVars <- c("Diagnosis_Alg", "Study", "Platform", "Sex", "Site")
samples[,c("Study", "Platform", "Sex", "Site")] <- lapply(samples[,c("Study", "Platform", "Sex", "Site")], as.factor)
PCAcatANOVAstats <- PCAcatANOVA(variables = catVars, sampleData = samples, PCdata = PCs)

# Continuous Covariate Stats with linear regression
contVars <- c(colnames(samples)[9:ncol(samples)])
PCAcontLmStats <- PCAcontLm(variables = contVars, sampleData = samples, PCdata = PCs)

# Adjust p-values
cat_padj <- unique(PCAcatANOVAstats[,c("PC", "Variable", "pvalue")])
cont_padj <- unique(PCAcontLmStats[,c("PC", "Variable", "pvalue")])

cat_padj$qvalue <- p.adjust(cat_padj$pvalue, method = "fdr", n = nrow(cat_padj) + nrow(cont_padj))
cont_padj$qvalue <- p.adjust(cont_padj$pvalue, method = "fdr", n = nrow(cat_padj) + nrow(cont_padj))

PCAcatANOVAstats$qvalue <- cat_padj$qvalue[match(PCAcatANOVAstats$pvalue, cat_padj$pvalue)]
PCAcontLmStats$qvalue <- cont_padj$qvalue[match(PCAcontLmStats$pvalue, cont_padj$pvalue)]

write.table(PCAcatANOVAstats, "Tables/Methylation PCs and Categorical Covariate Stats by ANOVA.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(PCAcontLmStats, "Tables/Methylation PCs and Continuous Covariate Stats by lm.txt", sep = "\t", quote = FALSE, row.names = FALSE)

unique(PCAcatANOVAstats[PCAcatANOVAstats$qvalue < 0.05 ,c("PC", "Variable")])
# PC2: Study, Platform, Sex, Site
# PC3: Study, Sex

unique(PCAcontLmStats[PCAcontLmStats$qvalue < 0.05 ,c("PC", "Variable")])
# PC1: MSLvrTscore36, MSLfmTscore36, MSLrlTscore36, MSLelTscore36, MSLelcStandard36, percent_cpg_meth, 
# PC2: percent_trimmed, percent_aligned, percent_duplicate, C_coverage, CG_coverage, percent_cpg_meth
# PC3: percent_trimmed, percent_aligned, percent_duplicate, dedup_reads, C_coverage, CG_coverage, percent_chg_meth, percent_cpg_meth

# Methylation PCs by Covariate Plots ####
# PCA Colored by Platform
pc1 <- summary(data.pca)$importance["Proportion of Variance", "PC1"]*100
pc2 <- summary(data.pca)$importance["Proportion of Variance","PC2"]*100
ggbiplotPCA(data.pca = data.pca, groups = samples$Platform, pc = c(1,2), 
            fileName = "Figures/ASD Cord Window Methylation by Platform PC1 PC2 Plot.png", xlim = c(-750,550), ylim = c(-650,650), 
            breaks = c("HiSeqX10", "HiSeq4000", "HiSeq2500"), 
            values = c("HiSeqX10"="#3366CC", "HiSeq4000"="#FF3366", "HiSeq2500"="#009933"))
ggbiplotPCA(data.pca = data.pca, groups = samples$Platform, pc = c(1,3), 
            fileName = "Figures/ASD Cord Window Methylation by Platform PC1 PC3 Plot.png", xlim = c(-750,550), ylim = c(-650,650), 
            breaks = c("HiSeqX10", "HiSeq4000", "HiSeq2500"), 
            values = c("HiSeqX10"="#3366CC", "HiSeq4000"="#FF3366", "HiSeq2500"="#009933"))
ggbiplotPCA(data.pca = data.pca, groups = samples$Platform, pc = c(2,3), 
            fileName = "Figures/ASD Cord Window Methylation by Platform PC2 PC3 Plot.png", xlim = c(-650,650), ylim = c(-650,650), 
            breaks = c("HiSeqX10", "HiSeq4000", "HiSeq2500"), 
            values = c("HiSeqX10"="#3366CC", "HiSeq4000"="#FF3366", "HiSeq2500"="#009933"))

# PCA Colored by Sex
ggbiplotPCA(data.pca = data.pca, groups = samples$Sex, pc = c(1,2), 
            fileName = "Figures/ASD Cord Window Methylation by Sex PC1 PC2 Plot.png", xlim = c(-750,550), ylim = c(-650,650), 
            breaks = c("M", "F"), values = c("M"="#3366CC", "F"="#FF3366"))
ggbiplotPCA(data.pca = data.pca, groups = samples$Sex, pc = c(1,3), 
            fileName = "Figures/ASD Cord Window Methylation by Sex PC1 PC3 Plot.png", xlim = c(-750,550), ylim = c(-650,650), 
            breaks = c("M", "F"), values = c("M"="#3366CC", "F"="#FF3366"))
ggbiplotPCA(data.pca = data.pca, groups = samples$Sex, pc = c(2,3), 
            fileName = "Figures/ASD Cord Window Methylation by Sex PC2 PC3 Plot.png", xlim = c(-650,650), ylim = c(-650,650), 
            breaks = c("M", "F"), values = c("M"="#3366CC", "F"="#FF3366"))

# PCA Colored by Study
ggbiplotPCA(data.pca = data.pca, groups = samples$Study, pc = c(1,2), 
            fileName = "Figures/ASD Cord Window Methylation by Study PC1 PC2 Plot.png", xlim = c(-750,550), ylim = c(-650,650), 
            breaks = c("MARBLES", "EARLI"), values = c("MARBLES"="#3366CC", "EARLI"="#FF3366"))
ggbiplotPCA(data.pca = data.pca, groups = samples$Study, pc = c(1,3), 
            fileName = "Figures/ASD Cord Window Methylation by Study PC1 PC3 Plot.png", xlim = c(-750,550), ylim = c(-650,650), 
            breaks = c("MARBLES", "EARLI"), values = c("MARBLES"="#3366CC", "EARLI"="#FF3366"))
ggbiplotPCA(data.pca = data.pca, groups = samples$Study, pc = c(2,3), 
            fileName = "Figures/ASD Cord Window Methylation by Study PC2 PC3 Plot.png", xlim = c(-650,650), ylim = c(-650,650), 
            breaks = c("MARBLES", "EARLI"), values = c("MARBLES"="#3366CC", "EARLI"="#FF3366"))

# PC1 vs Mullen Composite
ggScatterPlot(x = samples$MSLelcStandard36, y = PCs$PC1, diagnosis = samples$Diagnosis_Alg, ylim = c(-750,550),
              fileName = "Figures/ASD Cord Window Methylation by Mullen Composite PC1 Plot.png",
              xlab = "Mullen Composite Score", ylab = "Window Methylation PC1")
ggScatterPlot(x = samples$MSLelcStandard36[samples$Diagnosis_Alg == "ASD"], y = PCs$PC1[samples$Diagnosis_Alg == "ASD"], 
              diagnosis = samples$Diagnosis_Alg[samples$Diagnosis_Alg == "ASD"], ylim = c(-750,550),
              fileName = "Figures/ASD Cord Window Methylation by Mullen Composite PC1 Plot ASD only.png",
              xlab = "Mullen Composite Score", ylab = "Window Methylation PC1")
ggScatterPlot(x = samples$MSLelcStandard36[samples$Diagnosis_Alg == "TD"], y = PCs$PC1[samples$Diagnosis_Alg == "TD"], 
              diagnosis = samples$Diagnosis_Alg[samples$Diagnosis_Alg == "TD"], ylim = c(-750,550),
              fileName = "Figures/ASD Cord Window Methylation by Mullen Composite PC1 Plot TD only.png",
              xlab = "Mullen Composite Score", ylab = "Window Methylation PC1")

# PCs vs Global mCpG
ggScatterPlot(x = samples$percent_cpg_meth, y = PCs$PC1, diagnosis = samples$Diagnosis_Alg,
              fileName = "Figures/ASD Cord Window Methylation by mCpG PC1 Plot.png",
              xlab = "Global mCpG (%)", ylab = "Window Methylation PC1")
ggScatterPlot(x = samples$percent_cpg_meth, y = PCs$PC2, diagnosis = samples$Diagnosis_Alg,
              fileName = "Figures/ASD Cord Window Methylation by mCpG PC2 Plot.png",
              xlab = "Global mCpG (%)", ylab = "Window Methylation PC2")
ggScatterPlot(x = samples$percent_cpg_meth, y = PCs$PC3, diagnosis = samples$Diagnosis_Alg,
              fileName = "Figures/ASD Cord Window Methylation by mCpG PC3 Plot.png",
              xlab = "Global mCpG (%)", ylab = "Window Methylation PC3")




