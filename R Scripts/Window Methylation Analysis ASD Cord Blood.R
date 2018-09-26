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

ggbiplotPCA <- function(data.pca, diagnosis, pc = 1:2, fileName, xlim = NULL, ylim = NULL){
        # Plots principal components colored by diagnosis and writes file
        pc1 <- summary(data.pca)$importance["Proportion of Variance", paste("PC", pc[1], sep = "")]*100
        cat(paste("[ggbiplotPCA] PC", pc[1], " is ", pc1, "% of variance\n", sep=""))
        pc2 <- summary(data.pca)$importance["Proportion of Variance", paste("PC", pc[2], sep = "")]*100
        cat(paste("[ggbiplotPCA] PC", pc[2], " is ", pc2, "% of variance\n", sep=""))
        g <- ggbiplot(data.pca, obs.scale = 1, var.scale = 1, groups = diagnosis, ellipse = TRUE, circle = FALSE, 
                      var.axes = FALSE, varname.abbrev = FALSE, choices = pc, ellipse.prob = 0.95)
        g + 
                theme_bw(base_size = 25) +
                theme(legend.direction = 'horizontal', legend.position = c(0.83, 0.95), panel.grid.major = element_blank(), 
                      panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
                      legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
                      axis.text = element_text(color = "black"), legend.background = element_blank(), legend.spacing.x = unit(0.5,"lines")) +
                coord_cartesian(xlim = xlim, ylim = ylim) +
                xlab(paste("PC", pc[1], " (", round(pc1,0), "% of Variance)", sep = "")) +
                ylab(paste("PC", pc[2], " (", round(pc2,0), "% of Variance)", sep = "")) +
                scale_color_manual(breaks = c("TD", "ASD"), values = c("TD"="#3366CC", "ASD"="#FF3366")) +
                scale_x_continuous(breaks=pretty_breaks(5)) +
                scale_y_continuous(breaks=pretty_breaks(5)) +
                geom_point(aes(color = diagnosis), size=3)
        ggsave(fileName, dpi = 600, width = 8, height = 8, units = "in")
}

meth_ttest <- function(methData, sampleData, diagnosis){
        # Analyzes windowed methylation data for differences by diagnosis
        test <- t.test(methData ~ sampleData[,diagnosis])
        temp <- c(test$estimate[1], test$estimate[2], test$estimate[1] - test$estimate[2], test$conf.int[1], test$conf.int[2], 
                  test$p.value, use.names = FALSE)
        return(temp)
}

# Data ####
samples <- read.delim(file = "Samples/Merged Cord Blood WGBS Database.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
permeth <- read.delim(file = "Tables/windows_10kb_methylation_ASD_CordBlood.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

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
write.table(catStats, "Tables/Categorical Covariate Stats by Diagnosis.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(contStats, "Tables/Continuous Covariate Stats by Diagnosis.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

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

# Methylation by Diagnosis ####
# PCA
table(samples$Sequencing_ID == colnames(permeth)[4:ncol(permeth)]) # All TRUE
permeth <- permeth[1:(nrow(permeth)-1),] # remove last row, had na values
table(is.na(permeth)) # All FALSE
data <- t(as.matrix(permeth[,4:ncol(permeth)]))
data.pca <- prcomp(data, center = TRUE, scale. = TRUE) 
ggbiplotPCA(data.pca = data.pca, diagnosis = samples$Diagnosis_Alg, pc = c(1,2),
            fileName = "Figures/ASD Cord Window Methylation PC1 PC2 Plot by Diagnosis.png", xlim = c(-600,400), ylim = c(-500,500))
ggbiplotPCA(data.pca = data.pca, diagnosis = samples$Diagnosis_Alg, pc = c(1,3),
            fileName = "Figures/ASD Cord Window Methylation PC1 PC3 Plot by Diagnosis.png", xlim = c(-600,400), ylim = c(-500,500))
ggbiplotPCA(data.pca = data.pca, diagnosis = samples$Diagnosis_Alg, pc = c(2,3),
            fileName = "Figures/ASD Cord Window Methylation PC2 PC3 Plot by Diagnosis.png", xlim = c(-500,500), ylim = c(-500,500))

# Get PC 1-3
pc123 <- data.pca$x[,c("PC1", "PC2", "PC3")] %>% as.data.frame
pc123$Sequencing_ID <- rownames(pc123)
rownames(pc123) <- 1:nrow(pc123)
pc123 <- pc123[,c("Sequencing_ID", "PC1", "PC2", "PC3")]
rm(data, data.pca)

# t-test
methData <- permeth[,4:ncol(permeth)] %>% as.matrix
table(colnames(methData) == samples$Sequencing_ID) # All TRUE
methStats <- apply(methData, 1, meth_ttest, samples, "Diagnosis_Alg")
rm(methData)
methStats <- as.data.frame(t(methStats))
colnames(methStats) <- c("meanASD", "meanTD", "meanDiff", "confIntLow", "confIntHigh", "pValue")
methStats <- cbind(permeth[,c("chr", "start", "end")], methStats)
methStats$qValue <- p.adjust(methStats$pValue, method = "fdr")
methStats$Differential <- abs(methStats$meanDiff) > 1 & methStats$pValue < 0.05
table(methStats$Differential) # 9023 5.8%
table(methStats$qValue < 0.05) #1
methStats[methStats$qValue < 0.05,]
#   chr    start      end  meanASD   meanTD  meanDiff confIntLow confIntHigh       pValue     qValue Differential
# chr16 65520000 65530000 64.54697 67.64191 -3.094936  -4.227428   -1.962444 2.494654e-07 0.03850349         TRUE
write.table(methStats, "Tables/Window Methylation Stats by Diagnosis.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# Methylation by Covariates ####


