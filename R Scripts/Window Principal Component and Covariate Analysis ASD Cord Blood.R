# Window Principal Component and Covariate Analysis ####
# Autism Cord Blood Methylation Project
# Charles Mordaunt
# 1/30/19

# Packages ####
library(reshape2)
library(dplyr)
library(scales)
library(ggplot2)
library(ggbiplot)
library(sm)

# Functions ####
ggbiplotPCA <- function(data.pca, groups, pc = 1:2, fileName, xlim = NULL, ylim = NULL, breaks = c("TD", "ASD"), 
                        values = c("TD"="#3366CC", "ASD"="#FF3366"), legend.position = c(0.89, 0.94)){
         # Plots principal components colored by grouping variable and writes file
         pc1 <- summary(data.pca)$importance["Proportion of Variance", paste("PC", pc[1], sep = "")]*100
         pc2 <- summary(data.pca)$importance["Proportion of Variance", paste("PC", pc[2], sep = "")]*100
         g <- ggbiplot(data.pca, obs.scale = 1, var.scale = 1, groups = groups, ellipse = TRUE, circle = FALSE, 
                       var.axes = FALSE, varname.abbrev = FALSE, choices = pc, ellipse.prob = 0.95)
         g + 
                 theme_bw(base_size = 25) +
                 theme(legend.direction = 'vertical', legend.position = legend.position, panel.grid.major = element_blank(), 
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

ggScatterPlot <- function(x, y, groupVar, fileName, xlab, ylab, xlim = NULL, ylim = NULL, breaks = c("TD", "ASD"), 
                          values = c("TD"="#3366CC", "ASD"="#FF3366"), legend.position = c(0.87,1.03)){
        # Plots 2 continuous variables colored by a grouping variable and writes the file
        g <- ggplot()
        g + 
                geom_smooth(aes(x=x, y=y), method="lm") +
                geom_point(aes(x=x, y=y, color=groupVar), size=3) +
                theme_bw(base_size = 25) +
                theme(legend.direction = 'horizontal', legend.position = legend.position, panel.grid.major = element_blank(), 
                      panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
                      legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank(),
                      axis.text = element_text(color = "black"), legend.background = element_blank(), 
                      plot.margin = unit(c(2,1,1,1), "lines")) +
                coord_cartesian(xlim = xlim, ylim = ylim) +
                scale_x_continuous(breaks=pretty_breaks(n=5)) +
                scale_y_continuous(breaks=pretty_breaks(n=5)) +
                xlab(xlab) +
                ylab(ylab) +
                scale_color_manual(breaks = breaks, values = values)
        ggsave(fileName, dpi = 600, width = 8, height = 7, units = "in")
}

# Data ####
samples <- read.delim(file = "Merged Database/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo.txt", sep = "\t", header = TRUE, 
                      stringsAsFactors = FALSE)
samples101470 <- subset(samples, Cord_Blood_IBC == 101470)
samples <- subset(samples, !(Cord_Blood_IBC == 101470 & Platform %in% c("HiSeq2500", "HiSeq4000"))) # Remove 101470 Duplicates

permeth <- read.delim(file = "Tables/windows_10kb_methylation_ASD_CordBlood2.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
permeth101470 <- cbind(permeth[,c("chr", "start", "end")], permeth[,colnames(permeth) %in% samples101470$Sequencing_ID])
permeth <- permeth[,!colnames(permeth) %in% c("JLCM006C", "JLDS061B")] # Remove 101470 Duplicates

# 156 samples x 263618 windows

# Questions ####
# 10kb Window Methylation (PCA)
# PC1-3 stats with all categorical and continuous variables
# 	What variables are associated with specific methylation PCs?
# 	PC ~ variable
# 	Linear regression
# Significant variable PC scatterplots
# 	How do samples cluster by PCs and variables?
# 	Categorical: Diagnosis, Platform, Sex, Study
# 	Continuous: Mullen Composite, Global mCpG, gestational age
# 	Other significant variables?

# Window Methylation in 101470 Across Platforms ####
# Correlation Tests
cor.test(permeth101470$JLCM006C, permeth101470$JLCM069B) # HiSeq 4000 vs X10
# t = 1799.9, df = 263620, p-value < 2.2e-16
# 95 percent confidence interval:
# 0.9613498 0.9619243
# cor:
# 0.9616381 

cor.test(permeth101470$JLCM006C, permeth101470$JLDS061B) # HiSeq 4000 vs 2500
# t = 1330.7, df = 263620, p-value < 2.2e-16
# 95 percent confidence interval:
# 0.9324642 0.9334535
# cor:
# 0.9329606 

cor.test(permeth101470$JLCM069B, permeth101470$JLDS061B) # HiSeq X10 vs 2500
# t = 1381.5, df = 263620, p-value < 2.2e-16
# 95 percent confidence interval:
# 0.9368883 0.9378149
# cor:
# 0.9373532 

# Plots
permeth101470$Density <- sm.density(permeth101470[,c("JLCM006C", "JLCM069B")], 
                                    eval.points=permeth101470[,c("JLCM006C", "JLCM069B")], display = "none", 
                                    eval.grid = FALSE, panel = FALSE, verbose = FALSE)$estimate
gg <- ggplot()
gg + 
        geom_point(data = permeth101470, aes(JLCM006C, JLCM069B, color = Density), size = 1.5) +
        geom_abline(slope = 1, size = 1.25, linetype = "longdash", color = "white") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.25), panel.grid.minor = element_blank(),
              legend.position = "none", axis.text = element_text(color = "Black")) +
        xlab("HiSeq 4000 Methylation (%)") +
        ylab("HiSeq X10 Methylation (%)") +
        coord_cartesian(xlim = c(0, 100), ylim = c(0,100))
ggsave("Figures/Sample 101470 HiSeq 4000 vs X10 10kb Window Methylation.png", dpi = 600, height = 6, width = 7)

permeth101470$Density <- sm.density(permeth101470[,c("JLCM006C", "JLDS061B")], 
                                    eval.points=permeth101470[,c("JLCM006C", "JLDS061B")], display = "none", 
                                    eval.grid = FALSE, panel = FALSE, verbose = FALSE)$estimate
gg <- ggplot()
gg + 
        geom_point(data = permeth101470, aes(JLCM006C, JLDS061B, color = Density), size = 1.5) +
        geom_abline(slope = 1, size = 1.25, linetype = "longdash", color = "white") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.25), panel.grid.minor = element_blank(),
              legend.position = "none", axis.text = element_text(color = "Black")) +
        xlab("HiSeq 4000 Methylation (%)") +
        ylab("HiSeq 2500 Methylation (%)") +
        coord_cartesian(xlim = c(0, 100), ylim = c(0,100))
ggsave("Figures/Sample 101470 HiSeq 4000 vs 2500 10kb Window Methylation.png", dpi = 600, height = 6, width = 7)

permeth101470$Density <- sm.density(permeth101470[,c("JLCM069B", "JLDS061B")], 
                                    eval.points=permeth101470[,c("JLCM069B", "JLDS061B")], display = "none", 
                                    eval.grid = FALSE, panel = FALSE, verbose = FALSE)$estimate
gg <- ggplot()
gg + 
        geom_point(data = permeth101470, aes(JLCM069B, JLDS061B, color = Density), size = 1.5) +
        geom_abline(slope = 1, size = 1.25, linetype = "longdash", color = "white") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.25), panel.grid.minor = element_blank(),
              legend.position = "none", axis.text = element_text(color = "Black")) +
        xlab("HiSeq X10 Methylation (%)") +
        ylab("HiSeq 2500 Methylation (%)") +
        coord_cartesian(xlim = c(0, 100), ylim = c(0,100))
ggsave("Figures/Sample 101470 HiSeq X10 vs 2500 10kb Window Methylation.png", dpi = 600, height = 6, width = 7)
rm(gg, permeth101470, samples101470)

# Window Methylation PCA ####
samples <- samples[match(colnames(permeth)[4:ncol(permeth)], samples$Sequencing_ID),] # Make Same Order
table(samples$Sequencing_ID == colnames(permeth)[4:ncol(permeth)]) # All TRUE
table(is.na(permeth)) # All FALSE
data <- permeth[,4:ncol(permeth)] %>% as.matrix %>% t
data.pca <- prcomp(data, center = TRUE, scale. = TRUE) 

# Get PCs 1-3
PCs <- data.pca$x[,c("PC1", "PC2", "PC3")] %>% as.data.frame
PCs$Sequencing_ID <- rownames(PCs)
rownames(PCs) <- 1:nrow(PCs)
PCs <- PCs[,c("Sequencing_ID", "PC1", "PC2", "PC3")]

# Methylation PCs by Covariate Stats ####
# Prep Sample Data
samples <- merge(x = samples, y = PCs, by = "Sequencing_ID", all = FALSE, sort = FALSE)
catVars <- c("Study", "Platform", "Sex", "Site", "Diagnosis_Alg", "MomEdu_detail", "DM1or2", "GDM", "PE", 
             "home_ownership", "marital_status", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3", "AllEQ_PV_YN_Mo_2", 
             "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", "AllEQ_PV_YN_Mo4", 
             "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", "AllEQ_PV_YN_Mo9")
samples$Study <- factor(samples$Study, levels = c("MARBLES", "EARLI"))
samples$Platform <- factor(samples$Platform, levels = c("HiSeqX10", "HiSeq4000"))
samples$Sex <- factor(samples$Sex, levels = c("M", "F"))
samples$Site <- factor(samples$Site, levels = c("UC Davis", "Kaiser Permanente", "Drexel", "Johns Hopkins University"))
samples$Diagnosis_Alg <- factor(samples$Diagnosis_Alg, levels = c("TD", "ASD"))
samples$MomEdu_detail <- factor(samples$MomEdu_detail, levels = c(6, 1:5, 7,8))
samples$home_ownership[samples$home_ownership == 99] <- NA
samples$marital_status[samples$marital_status == 99] <- NA
factorCols <- c("DM1or2", "GDM", "PE", "marital_status", "home_ownership", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3", 
                "AllEQ_PV_YN_Mo_2", "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", 
                "AllEQ_PV_YN_Mo4", "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", 
                "AllEQ_PV_YN_Mo9")
samples[,factorCols] <- lapply(samples[,factorCols], as.factor)

contVars <- colnames(samples)[!colnames(samples) %in% catVars]
contVars <- contVars[!contVars %in% c("COI_ID", "Sequencing_ID", "Cord_Blood_IBC", "PC1", "PC2", "PC3")]
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

# Get Stats
statsPC1 <- methLm(catVars = catVars, contVars = contVars, sampleData = samples, globalMeth = "PC1")
write.table(statsPC1, "Tables/Window Methylation PC1 by Covariate Stats.txt", sep="\t", quote = FALSE, row.names = FALSE)
unique(statsPC1$Variable[statsPC1$qvalue < 0.05])
# percent_cpg_meth

statsPC2 <- methLm(catVars = catVars, contVars = contVars, sampleData = samples, globalMeth = "PC2")
write.table(statsPC2, "Tables/Window Methylation PC2 by Covariate Stats.txt", sep="\t", quote = FALSE, row.names = FALSE)
unique(statsPC2$Variable[statsPC2$qvalue < 0.05])
# [1] Study             Platform          Sex               Site              percent_trimmed   percent_aligned  
# [7] percent_duplicate C_coverage        CG_coverage       percent_cpg_meth 

statsPC3 <- methLm(catVars = catVars, contVars = contVars, sampleData = samples, globalMeth = "PC3")
write.table(statsPC3, "Tables/Window Methylation PC3 by Covariate Stats.txt", sep="\t", quote = FALSE, row.names = FALSE)
unique(statsPC3$Variable[statsPC3$qvalue < 0.05])
# [1] Study             Platform          Sex               percent_trimmed   percent_duplicate dedup_reads_M    
# [7] C_coverage        CG_coverage       percent_cpg_meth  percent_chg_meth 

# Methylation PCs by Covariate Plots ####
# PCA Colored by Diagnosis
ggbiplotPCA(data.pca = data.pca, groups = samples$Diagnosis_Alg, pc = c(1,2),
            fileName = "Figures/ASD Cord Window Methylation by Diagnosis PC1 PC2 Plot.png",
            xlim = c(-750,550), ylim = c(-650,650))
ggbiplotPCA(data.pca = data.pca, groups = samples$Diagnosis_Alg, pc = c(1,3),
            fileName = "Figures/ASD Cord Window Methylation by Diagnosis PC1 PC3 Plot.png", 
            xlim = c(-750,550), ylim = c(-650,650))
ggbiplotPCA(data.pca = data.pca, groups = samples$Diagnosis_Alg, pc = c(2,3),
            fileName = "Figures/ASD Cord Window Methylation by Diagnosis PC2 PC3 Plot.png", 
            xlim = c(-650,650), ylim = c(-650,650))

# PCA Colored by Platform
ggbiplotPCA(data.pca = data.pca, groups = samples$Platform, pc = c(1,2), 
            fileName = "Figures/ASD Cord Window Methylation by Platform PC1 PC2 Plot.png", 
            xlim = c(-750,550), ylim = c(-650,650), legend.position = c(0.85, 0.94),
            breaks = c("HiSeqX10", "HiSeq4000"), values = c("HiSeqX10"="#3366CC", "HiSeq4000"="#FF3366"))
ggbiplotPCA(data.pca = data.pca, groups = samples$Platform, pc = c(1,3), 
            fileName = "Figures/ASD Cord Window Methylation by Platform PC1 PC3 Plot.png", 
            xlim = c(-750,550), ylim = c(-650,650), legend.position = c(0.85, 0.94),
            breaks = c("HiSeqX10", "HiSeq4000"), values = c("HiSeqX10"="#3366CC", "HiSeq4000"="#FF3366"))
ggbiplotPCA(data.pca = data.pca, groups = samples$Platform, pc = c(2,3), 
            fileName = "Figures/ASD Cord Window Methylation by Platform PC2 PC3 Plot.png", 
            xlim = c(-650,650), ylim = c(-650,650),legend.position = c(0.85, 0.94),
            breaks = c("HiSeqX10", "HiSeq4000"), values = c("HiSeqX10"="#3366CC", "HiSeq4000"="#FF3366"))

# PCA Colored by Sex
ggbiplotPCA(data.pca = data.pca, groups = samples$Sex, pc = c(1,2), 
            fileName = "Figures/ASD Cord Window Methylation by Sex PC1 PC2 Plot.png", 
            xlim = c(-750,550), ylim = c(-650,650), legend.position = c(0.93, 0.94),
            breaks = c("M", "F"), values = c("M"="#3366CC", "F"="#FF3366"))
ggbiplotPCA(data.pca = data.pca, groups = samples$Sex, pc = c(1,3), 
            fileName = "Figures/ASD Cord Window Methylation by Sex PC1 PC3 Plot.png", 
            xlim = c(-750,550), ylim = c(-650,650), legend.position = c(0.93, 0.94),
            breaks = c("M", "F"), values = c("M"="#3366CC", "F"="#FF3366"))
ggbiplotPCA(data.pca = data.pca, groups = samples$Sex, pc = c(2,3), 
            fileName = "Figures/ASD Cord Window Methylation by Sex PC2 PC3 Plot.png", 
            xlim = c(-650,650), ylim = c(-650,650), legend.position = c(0.93, 0.94),
            breaks = c("M", "F"), values = c("M"="#3366CC", "F"="#FF3366"))

# PCA Colored by Study
ggbiplotPCA(data.pca = data.pca, groups = samples$Study, pc = c(1,2), 
            fileName = "Figures/ASD Cord Window Methylation by Study PC1 PC2 Plot.png", 
            xlim = c(-750,550), ylim = c(-650,650), legend.position = c(0.85, 0.94),
            breaks = c("MARBLES", "EARLI"), values = c("MARBLES"="#3366CC", "EARLI"="#FF3366"))
ggbiplotPCA(data.pca = data.pca, groups = samples$Study, pc = c(1,3), 
            fileName = "Figures/ASD Cord Window Methylation by Study PC1 PC3 Plot.png", 
            xlim = c(-750,550), ylim = c(-650,650), legend.position = c(0.85, 0.94),
            breaks = c("MARBLES", "EARLI"), values = c("MARBLES"="#3366CC", "EARLI"="#FF3366"))
ggbiplotPCA(data.pca = data.pca, groups = samples$Study, pc = c(2,3), 
            fileName = "Figures/ASD Cord Window Methylation by Study PC2 PC3 Plot.png", 
            xlim = c(-650,650), ylim = c(-650,650), legend.position = c(0.85, 0.94),
            breaks = c("MARBLES", "EARLI"), values = c("MARBLES"="#3366CC", "EARLI"="#FF3366"))

# PC1 vs Mullen Composite
ggScatterPlot(x = samples$MSLelcStandard36, y = PCs$PC1, groupVar = samples$Diagnosis_Alg, ylim = c(-750,550),
              fileName = "Figures/ASD Cord Window Methylation by Mullen Composite PC1 Plot.png",
              xlab = "Mullen Composite Score", ylab = "Window Methylation PC1")
ggScatterPlot(x = samples$MSLelcStandard36[samples$Diagnosis_Alg == "ASD"], y = PCs$PC1[samples$Diagnosis_Alg == "ASD"], 
              groupVar = samples$Diagnosis_Alg[samples$Diagnosis_Alg == "ASD"], ylim = c(-750,550),
              fileName = "Figures/ASD Cord Window Methylation by Mullen Composite PC1 Plot ASD only.png",
              xlab = "Mullen Composite Score", ylab = "Window Methylation PC1")
ggScatterPlot(x = samples$MSLelcStandard36[samples$Diagnosis_Alg == "TD"], y = PCs$PC1[samples$Diagnosis_Alg == "TD"], 
              groupVar = samples$Diagnosis_Alg[samples$Diagnosis_Alg == "TD"], ylim = c(-750,550),
              fileName = "Figures/ASD Cord Window Methylation by Mullen Composite PC1 Plot TD only.png",
              xlab = "Mullen Composite Score", ylab = "Window Methylation PC1")

# PCs vs Global mCpG
ggScatterPlot(x = samples$percent_cpg_meth, y = PCs$PC1, groupVar = samples$Diagnosis_Alg,
              fileName = "Figures/ASD Cord Window Methylation by mCpG PC1 Plot.png",
              xlab = "Global mCpG (%)", ylab = "Window Methylation PC1")
ggScatterPlot(x = samples$percent_cpg_meth, y = PCs$PC2, groupVar = samples$Platform,
              fileName = "Figures/ASD Cord Window Methylation by mCpG PC2 Plot.png",
              xlab = "Global mCpG (%)", ylab = "Window Methylation PC2", legend.position = c(0.75,1.03),
              breaks = c("HiSeqX10", "HiSeq4000"), values = c("HiSeqX10"="#3366CC", "HiSeq4000"="#FF3366"))
ggScatterPlot(x = samples$percent_cpg_meth, y = PCs$PC3, groupVar = samples$Sex,
              fileName = "Figures/ASD Cord Window Methylation by mCpG PC3 Plot.png",
              xlab = "Global mCpG (%)", ylab = "Window Methylation PC3", legend.position = c(0.92,1.03),
              breaks = c("M", "F"), values = c("M"="#3366CC", "F"="#FF3366"))

# PCs vs PCR Duplicates
ggScatterPlot(x = samples$percent_duplicate, y = PCs$PC1, groupVar = samples$Diagnosis_Alg,
              fileName = "Figures/ASD Cord Window Methylation by PCR Dups PC1 Plot.png",
              xlab = "PCR Duplicates (%)", ylab = "Window Methylation PC1")
ggScatterPlot(x = samples$percent_duplicate, y = PCs$PC2, groupVar = samples$Platform,
              fileName = "Figures/ASD Cord Window Methylation by PCR Dups PC2 Plot.png",
              xlab = "PCR Duplicates (%)", ylab = "Window Methylation PC2", legend.position = c(0.75,1.03),
              breaks = c("HiSeqX10", "HiSeq4000"), values = c("HiSeqX10"="#3366CC", "HiSeq4000"="#FF3366"))
ggScatterPlot(x = samples$percent_duplicate, y = PCs$PC3, groupVar = samples$Sex,
              fileName = "Figures/ASD Cord Window Methylation by PCR Dups PC3 Plot.png",
              xlab = "PCR Duplicates (%)", ylab = "Window Methylation PC3", legend.position = c(0.92,1.03),
              breaks = c("M", "F"), values = c("M"="#3366CC", "F"="#FF3366"))

# PCs vs Gestational Age
ggScatterPlot(x = samples$ga_w, y = PCs$PC1, groupVar = samples$Diagnosis_Alg,
              fileName = "Figures/ASD Cord Window Methylation by Gestational Age PC1 Plot.png",
              xlab = "Gestational Age (Weeks)", ylab = "Window Methylation PC1")
ggScatterPlot(x = samples$ga_w, y = PCs$PC2, groupVar = samples$Platform,
              fileName = "Figures/ASD Cord Window Methylation by Gestational Age PC2 Plot.png",
              xlab = "Gestational Age (Weeks)", ylab = "Window Methylation PC2", legend.position = c(0.75,1.03),
              breaks = c("HiSeqX10", "HiSeq4000"), values = c("HiSeqX10"="#3366CC", "HiSeq4000"="#FF3366"))
ggScatterPlot(x = samples$ga_w, y = PCs$PC3, groupVar = samples$Sex,
              fileName = "Figures/ASD Cord Window Methylation by Gestataional Age PC3 Plot.png",
              xlab = "Gestational Age (Weeks)", ylab = "Window Methylation PC3", legend.position = c(0.92,1.03),
              breaks = c("M", "F"), values = c("M"="#3366CC", "F"="#FF3366"))


