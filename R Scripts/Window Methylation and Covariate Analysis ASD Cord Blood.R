# Window Methylation and Covariate Analysis ####
# Autism Cord Blood Methylation
# Charles Mordaunt
# 10/2/18

# Packages ####
library(reshape2)
library(dplyr)
library(scales)
library(ggplot2)

# Functions ####
meth_ttest <- function(methData, sampleData, group){
        # Analyzes windowed methylation data for differences by a grouping variable with 2 levels
        test <- t.test(methData ~ sampleData[,group])
        temp <- c(test$estimate[1], test$estimate[2], test$estimate[1] - test$estimate[2], test$conf.int[1], test$conf.int[2], 
                  test$p.value, use.names = FALSE)
        return(temp)
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

# Window Methylation by Diagnosis ####
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

table(samplesPlatform$Platform, samplesPlatform$Sex, samplesPlatform$Diagnosis_Alg)
# TD (77)
#            M  F     M     F
# HiSeq4000 18  3 0.234 0.039
# HiSeqX10  39 17 0.506 0.221

# ASD (78)
#            M  F     M     F
# HiSeq4000 21  5 0.269 0.064
# HiSeqX10  37 15 0.474 0.192
