# Window Methylation Analysis ####
# Autism Cord Blood Methylation Project
# Charles Mordaunt
# 9/24/18

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
        colnames(stats) <- c("Variable", "Value", "ASD", "TD", "pvalue")
        stats <- as.data.frame(stats, stringsAsFactors = FALSE)
        stats[,c("Variable", "Value")] <- lapply(stats[,c("Variable", "Value")], as.factor)
        stats[,c("ASD", "TD")] <- lapply(stats[,c("ASD", "TD")], as.integer)
        stats$pvalue <- as.numeric(stats$pvalue)
        varSums <- aggregate(as.matrix(stats[,c("ASD", "TD")]) ~ Variable, data=stats, FUN=sum)
        stats$per_ASD <- stats$ASD / varSums$ASD[match(stats$Variable, varSums$Variable)]
        stats$per_TD <- stats$TD / varSums$TD[match(stats$Variable, varSums$Variable)]
        return(stats)
}

contANOVA <- function(variables, sampleData, diagnosis){
        # Analyzes table of continuous covariates for differences by ASD vs TD diagnosis using one-way ANOVA
        stats <- NULL
        for(i in 1:length(variables)){
                means <- aggregate(sampleData[,variables[i]] ~ sampleData[,diagnosis], FUN=mean)[,2]
                sds <- aggregate(sampleData[,variables[i]] ~ sampleData[,diagnosis], FUN=sd)[,2]
                p <- summary(aov(sampleData[,variables[i]] ~ sampleData[,diagnosis]))[[1]][1, "Pr(>F)"]
                temp <- cbind(rep(variables[i], 2), levels(sampleData[,diagnosis]), means, sds, rep(p, 2))
                stats <- rbind(stats, temp)
        }
        rownames(stats) <- 1:nrow(stats)
        colnames(stats) <- c("Variable", "Diagnosis", "Mean", "SD", "pvalue")
        stats <- as.data.frame(stats, stringsAsFactors = FALSE)
        stats[,c("Variable", "Diagnosis")] <- lapply(stats[,c("Variable", "Diagnosis")], as.factor)
        stats[,c("Mean", "SD", "pvalue")] <- lapply(stats[,c("Mean", "SD", "pvalue")], as.numeric)
        return(stats)
}

ggbiplotPCA <- function(data.pca, groups, pc = 1:2, fileName, xlim = NULL, ylim = NULL, breaks = c("TD", "ASD"), values = c("TD"="#3366CC", "ASD"="#FF3366")){
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

meth_ttest <- function(methData, sampleData, group){
        # Analyzes windowed methylation data for differences by a grouping variable with 2 levels
        test <- t.test(methData ~ sampleData[,group])
        temp <- c(test$estimate[1], test$estimate[2], test$estimate[1] - test$estimate[2], test$conf.int[1], test$conf.int[2], 
                  test$p.value, use.names = FALSE)
        return(temp)
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

ggScatterPlot <- function(x, y, diagnosis, fileName, xlab, ylab, xlim = NULL, ylim = NULL){
        # Plots 2 continuous variables colored by ASD vs TD diagnosis and writes the file
        g <- ggplot()
        g + 
                geom_smooth(aes(x=x, y=y), method="lm") +
                geom_point(aes(x=x, y=y, color=diagnosis), size=3) +
                theme_bw(base_size = 25) +
                theme(legend.direction = 'horizontal', legend.position = c(0.89, 1.035), panel.grid.major = element_blank(), 
                      panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
                      legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank(),
                      axis.text = element_text(color = "black"), legend.background = element_blank(), 
                      plot.margin = unit(c(2,1,1,1), "lines")) +
                coord_cartesian(xlim = xlim, ylim = ylim) +
                scale_x_continuous(breaks=pretty_breaks(n=5)) +
                scale_y_continuous(breaks=pretty_breaks(n=5)) +
                xlab(xlab) +
                ylab(ylab) +
                scale_color_manual(breaks = c("TD", "ASD"), values = c("TD"="#3366CC", "ASD"="#FF3366"))
        ggsave(fileName, dpi = 600, width = 8, height = 7, units = "in")
}

meth_anova <- function(methData, sampleData, formula, diagnosis){
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

# Questions ####
# How do covariates differ by diagnosis?
#       Covariate ~ ASD vs TD
#       Categorical: Fisher's exact test
#       Continuous: One-way ANOVA
# How does methylation differ by diagnosis?
#       Methylation ~ ASD vs TD
#       One-way ANOVA
#       PCA, color by diagnosis
# How does methylation differ by covariates?
#       PCA ~ covariate
#       Categorical: One-way ANOVA
#       Continuous: Linear regression
#       PCA, color by significant covariates

# Covariates by Diagnosis ####
samples$Diagnosis_Alg <- factor(samples$Diagnosis_Alg, levels = c("ASD", "TD"), ordered = TRUE)

# Categorical
catVars <- c("Study", "Platform", "Sex", "Site")
catStats <- catFisher(variables = catVars, sampleData = samples, diagnosis = "Diagnosis_Alg")

# Continuous
contVars <- c(colnames(samples)[9:ncol(samples)])
# [1] "MSLvrTscore36"     "MSLfmTscore36"     "MSLrlTscore36"     "MSLelTscore36"     "MSLelcStandard36"  "ADOScs"           
# [7] "percent_trimmed"   "percent_aligned"   "percent_duplicate" "dedup_reads"       "C_coverage"        "CG_coverage"      
# [13] "percent_chg_meth"  "percent_chh_meth"  "percent_cpg_meth" 
contStats <- contANOVA(variables = contVars, sampleData = samples, diagnosis = "Diagnosis_Alg")

# Adjust pvalues
allVars <- c(catVars, contVars)
catStats$qvalue <- p.adjust(catStats$pvalue, method = "fdr", n = length(allVars))
cont_padj <- unique(contStats[,c("Variable", "pvalue")])
cont_padj$qvalue <- p.adjust(cont_padj$pvalue, method = "fdr", n = length(allVars))
contStats$qvalue <- cont_padj$qvalue[match(contStats$pvalue, cont_padj$pvalue)]

# Write results
#write.table(catStats, "Tables/Categorical Covariate Stats by Diagnosis.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
#write.table(contStats, "Tables/Continuous Covariate Stats by Diagnosis.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

catStats$Variable[catStats$pvalue < 0.05] %>% as.character %>% unique # None
contStats$Variable[contStats$pvalue < 0.05] %>% as.character %>% unique
# "MSLvrTscore36"    "MSLfmTscore36"    "MSLrlTscore36"    "MSLelTscore36"    "MSLelcStandard36" 
# "ADOScs"           "percent_cpg_meth"
contStats$Variable[contStats$qvalue < 0.05] %>% as.character %>% unique
# "MSLvrTscore36"    "MSLfmTscore36"    "MSLrlTscore36"    "MSLelTscore36"    "MSLelcStandard36" "ADOScs"          
subset(contStats, Variable == "percent_cpg_meth")
# Diagnosis     Mean       SD    pvalue     qvalue
#       ASD 76.39468 1.415202 0.0217335 0.05899092
#        TD 76.86798 1.116904 0.0217335 0.05899092
# Only behavioral outcomes and global CpG methylation differs by diagnosis, not other covariates
rm(cont_padj, allVars, catVars, contVars)

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

# t-test
methData <- permeth[,4:ncol(permeth)] %>% as.matrix
table(colnames(methData) == samples$Sequencing_ID) # All TRUE
methStats <- apply(methData, 1, meth_ttest, samples, "Diagnosis_Alg")
methStats <- as.data.frame(t(methStats))
colnames(methStats) <- c("meanASD", "meanTD", "meanDiff", "confIntLow", "confIntHigh", "pValue")
methStats <- cbind(permeth[,c("chr", "start", "end")], methStats)
methStats$qValue <- p.adjust(methStats$pValue, method = "fdr")
methStats$Differential <- abs(methStats$meanDiff) > 1 & methStats$pValue < 0.05
table(methStats$Differential) # 15634 5.9%
table(methStats$qValue < 0.05) # None
write.table(methStats, "Tables/Window Methylation Stats by Diagnosis.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

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

# Global Methylation by Platform ####
# Global mCpG ~ Platform
mCpG_platform_data <- subset(samples, Platform %in% c("HiSeq4000", "HiSeqX10"), select = c(Sequencing_ID, Platform, percent_cpg_meth))
mCpG_platform_data$Platform <- mCpG_platform_data$Platform %>% as.character %>% as.factor
mCpG_platform <- aov(mCpG_platform_data$percent_cpg_meth ~ mCpG_platform_data$Platform)
summary(mCpG_platform)[[1]][1,"Pr(>F)"] # 7.789803e-17
TukeyHSD(mCpG_platform)
#                      diff      lwr      upr p adj
# HiSeqX10-HiSeq4000 1.6955 1.338953 2.052047     0
aggregate(mCpG_platform_data$percent_cpg_meth, by = list(mCpG_platform_data$Platform), FUN = mean)
# HiSeq4000 75.45824
#  HiSeqX10 77.15374

g <- ggplot(data = mCpG_platform_data)
g + 
        geom_boxplot(aes(x=Platform, y=percent_cpg_meth, fill = Platform), size = 1.1) +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = "none", panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              plot.margin = unit(c(2,1,1,1), "lines"), axis.title.x = element_blank()) +
        scale_y_continuous(breaks=pretty_breaks(n=4)) +
        ylab("Global CpG Methylation (%)") +
        scale_fill_manual(breaks = c("HiSeqX10", "HiSeq4000"), values = c("HiSeqX10"="#3366CC", "HiSeq4000"="#FF3366"))
ggsave("Figures/ASD Cord Global Methylation by Platform Boxplot.png", dpi = 600, width = 7, height = 7, units = "in")

# Global mCpG ~ Platform + Diagnosis
mCpG_platformDiag_data <- subset(samples, Platform %in% c("HiSeq4000", "HiSeqX10"), select = c(Sequencing_ID, Platform, Diagnosis_Alg, percent_cpg_meth))
mCpG_platformDiag_data$Platform <- mCpG_platformDiag_data$Platform %>% as.character %>% as.factor
mCpG_platformDiag <- aov(percent_cpg_meth ~ Platform * Diagnosis_Alg, data = mCpG_platformDiag_data)
summary(mCpG_platformDiag)
#                         Df Sum Sq Mean Sq F value Pr(>F)    
# Platform                 1  94.14   94.14  91.829 <2e-16 ***
# Diagnosis_Alg            1   5.93    5.93   5.783 0.0174 *  
# Platform:Diagnosis_Alg   1   2.47    2.47   2.407 0.1229  

TukeyHSD(mCpG_platformDiag)
#                                  diff          lwr        upr     p adj
# HiSeqX10:ASD-HiSeq4000:ASD  1.9256250  1.293825963  2.5574240 0.0000000
# HiSeq4000:TD-HiSeq4000:ASD  0.7773216  0.005578767  1.5490645 0.0476425
# HiSeqX10:TD-HiSeq4000:ASD   2.1516317  1.527399423  2.7758640 0.0000000
# HiSeq4000:TD-HiSeqX10:ASD  -1.1483034 -1.828399223 -0.4682076 0.0001257
# HiSeqX10:TD-HiSeqX10:ASD    0.2260067 -0.280559455  0.7325729 0.6535297
# HiSeqX10:TD-HiSeq4000:TD    1.3743101  0.701237838  2.0473824 0.0000024
aggregate(percent_cpg_meth ~ Platform + Diagnosis_Alg, data = mCpG_platformDiag_data, FUN = mean)
#  Platform Diagnosis_Alg percent_cpg_meth
# HiSeq4000           ASD         75.11093
#  HiSeqX10           ASD         77.03656
# HiSeq4000            TD         75.88825
#  HiSeqX10            TD         77.26256

fisher.test(mCpG_platformDiag_data$Platform, mCpG_platformDiag_data$Diagnosis_Alg, workspace=2e7)
# odds ratio    1.330844 
# 95% CI        0.6340182 2.8184452
# p-value       0.4854

g <- ggplot(data = mCpG_platformDiag_data)
g + 
        geom_boxplot(aes(x=Platform, y=percent_cpg_meth, fill = Diagnosis_Alg), size = 1.1) +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.88, y=1.03), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              plot.margin = unit(c(2,1,1,1), "lines"), axis.title.x = element_blank()) +
        scale_y_continuous(breaks=pretty_breaks(n=4)) +
        ylab("Global CpG Methylation (%)") +
        scale_fill_manual(breaks = c("TD", "ASD"), values = c("TD"="#3366CC", "ASD"="#FF3366"))
ggsave("Figures/ASD Cord Global Methylation by Platform and Diagnosis Boxplot.png", dpi = 600, width = 7, height = 7, units = "in")

# Global mCG by Diagnosis
mCpG_Diag <- aov(percent_cpg_meth ~ Diagnosis_Alg, data = mCpG_platformDiag_data)
summary(mCpG_Diag)
#                Df Sum Sq Mean Sq F value Pr(>F)  
# Diagnosis_Alg   1   9.42    9.42   5.814 0.0171 *
TukeyHSD(mCpG_Diag)
#             diff        lwr       upr     p adj
# TD-ASD 0.4930699 0.08906898 0.8970708 0.0170879
aggregate(percent_cpg_meth ~ Diagnosis_Alg, data = mCpG_platformDiag_data, FUN = mean)
# Diagnosis_Alg percent_cpg_meth
#           ASD         76.39468
#            TD         76.88775

g <- ggplot(data = mCpG_platformDiag_data)
g + 
        geom_boxplot(aes(x=Diagnosis_Alg, y=percent_cpg_meth, fill = Diagnosis_Alg), size = 1.1) +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.88, y=1.03), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              plot.margin = unit(c(2,1,1,1), "lines"), axis.title.x = element_blank()) +
        scale_y_continuous(breaks=pretty_breaks(n=4)) +
        ylab("Global CpG Methylation (%)") +
        scale_fill_manual(breaks = c("TD", "ASD"), values = c("TD"="#3366CC", "ASD"="#FF3366"))
ggsave("Figures/ASD Cord Global Methylation by Diagnosis Boxplot.png", dpi = 600, width = 7, height = 7, units = "in")

# Window Methylation by Platform ####
# HiSeq 4000 vs X10 t-test
table(colnames(methData) == samples$Sequencing_ID) # All TRUE
methData <- methData[,!colnames(methData) %in% samples$Sequencing_ID[samples$Platform == "HiSeq2500"]] # Remove 2500 sample
samplesPlatform <- subset(samples, !Platform == "HiSeq2500")
table(colnames(methData) == samplesPlatform$Sequencing_ID) # All TRUE
methStatsPlatform <- apply(methData, 1, meth_ttest, samplesPlatform, "Platform")
rm(methData)
methStatsPlatform <- as.data.frame(t(methStatsPlatform))
colnames(methStatsPlatform) <- c("mean4000", "meanX10", "meanDiff", "confIntLow", "confIntHigh", "pValue")
methStatsPlatform <- cbind(permeth[,c("chr", "start", "end")], methStatsPlatform)
methStatsPlatform$qValue <- p.adjust(methStatsPlatform$pValue, method = "fdr")
methStatsPlatform$Differential <- abs(methStatsPlatform$meanDiff) > 1 & methStatsPlatform$pValue < 0.05
table(methStatsPlatform$Differential) # TRUE 66862, 25%
table(methStatsPlatform$Differential & methStatsPlatform$qValue < 0.05) # TRUE 45912, 19%
table(methStatsPlatform$qValue < 0.05 & methStatsPlatform$meanDiff > 1) # TRUE 9023 4000 > X10
table(methStatsPlatform$qValue < 0.05 & methStatsPlatform$meanDiff < -1) # TRUE 36889 X10 > 4000
write.table(methStatsPlatform, "Tables/Window Methylation Stats by Platform.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# Percent differential windows at q < 0.05 by chromosome
table(methStatsPlatform$Differential & methStatsPlatform$qValue < 0.05, methStatsPlatform$chr)["TRUE",]/table(methStatsPlatform$chr)
# chr1     chr10     chr11     chr12     chr13     chr14     chr15     chr16     chr17     chr18     chr19      chr2     chr20 
# 0.2019841 0.1670230 0.1839813 0.1906263 0.1413128 0.1808132 0.1921542 0.2100000 0.2636732 0.1408018 0.3780812 0.1587645 0.2043747 
# chr21     chr22      chr3      chr4      chr5      chr6      chr7      chr8      chr9      chrX      chrY 
# 0.1536578 0.2691729 0.1562030 0.1358117 0.1500543 0.1596256 0.1723808 0.1473843 0.1820283 0.1211864 0.5238095 

# Window Methylation by Sex ####
# M vs F t-test
methData <- permeth[,4:ncol(permeth)] %>% as.matrix
table(colnames(methData) == samples$Sequencing_ID) # All TRUE
methStatsSex <- apply(methData, 1, meth_ttest, samples, "Sex")
rm(methData)
methStatsSex <- as.data.frame(t(methStatsSex))
colnames(methStatsSex) <- c("meanF", "meanM", "meanDiff", "confIntLow", "confIntHigh", "pValue")
methStatsSex <- cbind(permeth[,c("chr", "start", "end")], methStatsSex)
methStatsSex$qValue <- p.adjust(methStatsSex$pValue, method = "fdr")
methStatsSex$Differential <- abs(methStatsSex$meanDiff) > 1 & methStatsSex$pValue < 0.05
table(methStatsSex$Differential) # TRUE 19548, 7%
table(methStatsSex$Differential & methStatsSex$qValue < 0.05) # TRUE 11346, 4%, 10547 on chrX
table(methStatsSex$qValue < 0.05 & methStatsSex$meanDiff > 1) # TRUE 1224 F > M, 640 on chrX
table(methStatsSex$qValue < 0.05 & methStatsSex$meanDiff < -1) # TRUE 10122 M > F, 9907 on chrX
write.table(methStatsSex, "Tables/Window Methylation Stats by Sex.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# Percent differential windows at q < 0.05 by chromosome
table(methStatsSex$Differential & methStatsSex$qValue < 0.05, methStatsSex$chr)["TRUE",]/table(methStatsSex$chr)
#        chr1       chr10       chr11       chr12       chr13       chr14       chr15       chr16       chr17       chr18       chr19 
# 0.003573276 0.002296666 0.003710276 0.002424242 0.002006018 0.002702371 0.002064694 0.003424658 0.005315429 0.003352424 0.008153204 
#        chr2       chr20       chr21       chr22        chr3        chr4        chr5        chr6        chr7        chr8        chr9 
# 0.002827979 0.011521122 0.002754821 0.008120301 0.001717107 0.001841515 0.002835767 0.002209801 0.002805721 0.001656003 0.003882777 
#        chrX        chrY 
# 0.896320218 0.666666667 

# Window Methylation by Diagnosis, Adjusting for Platform and Sex ####
methData <- permeth[,4:ncol(permeth)] %>% as.matrix
table(colnames(methData) == samples$Sequencing_ID) # All TRUE
methData <- methData[,!colnames(methData) %in% samples$Sequencing_ID[samples$Platform == "HiSeq2500"]] # Remove 2500 sample
samplesPlatform <- subset(samples, !Platform == "HiSeq2500")
table(colnames(methData) == samplesPlatform$Sequencing_ID) # All TRUE

samplesPlatform$Diagnosis_Alg <- samplesPlatform$Diagnosis_Alg %>% as.character %>% factor(levels = c("TD", "ASD"))
samplesPlatform$Platform <- samplesPlatform$Platform %>% as.character %>% factor(levels = c("HiSeq4000", "HiSeqX10"))
samplesPlatform$Sex <- samplesPlatform$Sex %>% as.character %>% factor(levels = c("M", "F"))

formula <- meth ~ Diagnosis_Alg + Platform + Sex
methANOVAstats <- apply(methData, 1, meth_anova, sampleData = samplesPlatform, formula = formula, diagnosis = "Diagnosis_Alg")
methANOVAstats <- as.data.frame(t(methANOVAstats))
colnames(methANOVAstats) <- c("meanTD", "meanASD", "sdTD", "sdASD", "DiagnosisEstimate", "DiagnosisSE", "DiagnosisTvalue", "DiagnosisPvalue",
                     "PlatformEstimate", "PlatformSE", "PlatformTvalue", "PlatformPvalue", "SexEstimate", "SexSE", "SexTvalue", "SexPvalue")
methANOVAstats <- cbind(permeth[,c("chr", "start", "end")], methANOVAstats)
methANOVAstats$DiagnosisQvalue <- p.adjust(methANOVAstats$DiagnosisPvalue, method = "fdr")
methANOVAstats$PlatformQvalue <- p.adjust(methANOVAstats$PlatformPvalue, method = "fdr")
methANOVAstats$SexQvalue <- p.adjust(methANOVAstats$SexPvalue, method = "fdr")
methANOVAstats$DiagnosisDifferential <- abs(methANOVAstats$DiagnosisEstimate) > 1 & methANOVAstats$DiagnosisPvalue < 0.05
methANOVAstats$PlatformDifferential <- abs(methANOVAstats$PlatformEstimate) > 1 & methANOVAstats$PlatformPvalue < 0.05
methANOVAstats$SexDifferential <- abs(methANOVAstats$SexEstimate) > 1 & methANOVAstats$SexPvalue < 0.05

table(methANOVAstats$DiagnosisDifferential, methANOVAstats$DiagnosisQvalue < 0.05)
#        FALSE
# FALSE 247696
# TRUE   15922
table(methANOVAstats$PlatformDifferential, methANOVAstats$PlatformQvalue < 0.05)
#        FALSE   TRUE
# FALSE 188356   4581
# TRUE   17947  52734
table(methANOVAstats$SexDifferential, methANOVAstats$SexQvalue < 0.05)
#        FALSE   TRUE
# FALSE 244162     45
# TRUE    8132  11279

write.table(methANOVAstats, "Tables/Window Methylation ANOVA Stats with model.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# Next: how does adding in covariates change diagnosis differential windows?
# Can we adjust for platform and/or sex? or should we stratify?




