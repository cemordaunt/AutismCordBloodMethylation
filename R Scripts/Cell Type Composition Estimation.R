# Cell Type Composition Estimation ----------------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 3/5/19

# Packages ####
packages <- c("Biostrings", "AnnotationDbi", "openssl", "FlowSorted.CordBlood.450k", "minfi", 
              "IlluminaHumanMethylation450kmanifest", "IlluminaHumanMethylation450kanno.ilmn12.hg19", "data.table", "tidyverse",
              "reshape2", "scales", "wesanderson", "RColorBrewer", "rlist")
sapply(packages, require, character.only = TRUE)
rm(packages)

# Functions ####
source("R Scripts/DMR Analysis Functions.R")

# Estimate Cell Composition --------------------------------------------------------------------
# Get 450K Probe Methylation in WGBS Samples ####
samples <- read.delim(file = "Merged Database/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo.txt", 
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE)
samples <- subset(samples, !(Cord_Blood_IBC == 101470 & Platform %in% c("HiSeq2500", "HiSeq4000")))
probes <- read.delim(file = "UCSC Tracks/HM450_hg38.bed", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
probes_dup <- probes[duplicated(probes[,1:3]),] # 16 locations map 2 probes
probes <- probes[!duplicated(probes[,1:3]),]
write.table(probes, "UCSC Tracks/HM450_hg38_unique.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
rm(probes)

# Get Cell-Type Specific 450K probes ####
data(FlowSorted.CordBlood.450k.ModelPars)
ModelPars <- as.data.frame(FlowSorted.CordBlood.450k.ModelPars)
rm(FlowSorted.CordBlood.450k.ModelPars)

# Get Cell-Type Specific 450K Probe Methylation in WGBS Samples ####
methCounts <- read.delim(file = "Tables/probe_450k_methylation_ASD_CordBlood.txt.2col", sep = "\t", header = TRUE, 
                         stringsAsFactors = FALSE) # 447978 probes
methCounts <- subset(methCounts, Name %in% rownames(ModelPars)) # 668 / 700 probes covered in >= 1 sample
meth <- methCounts[,grepl("meth", colnames(methCounts), fixed = TRUE)]
cov <- methCounts[,grepl("total", colnames(methCounts), fixed = TRUE)]
naCount <- apply(cov, 1, function(x) sum(is.na(x)))
table(naCount < nrow(samples) / 4) # 636 probes with na in < 25% of samples
table(naCount < nrow(samples) / 10) # 560 probes with na in < 10% of samples
table(naCount == 0) # 127 probes covered in all samples

cov <- subset(cov, naCount < nrow(samples) / 10) %>% as.matrix # Subset for probes covered in 90% of samples
meth <- subset(meth, naCount < nrow(samples) / 10) %>% as.matrix
permeth <- meth / cov
table(is.na(permeth)) # 2062
for(i in 1:nrow(permeth)){
        temp <- permeth[i,]
        temp[is.na(temp)] <- mean(temp, na.rm = TRUE) # Replace missing values with probe mean value
        permeth[i,] <- temp
}
table(is.na(permeth)) # 0
rownames(permeth) <- methCounts$Name[naCount < nrow(samples) / 10] # probe IDs
colnames(permeth) <- gsub("_meth", "", colnames(permeth), fixed = TRUE) # Sample SeqIDs
methCounts <- subset(methCounts, Name %in% rownames(permeth))
colnames(probes_dup) <- c("Chromosome", "Start", "End")
table(duplicated(rbind(methCounts[,1:3], probes_dup[,1:3]))) # probes used for cell estimation assigned to the same location
table(methCounts$Chromosome) # No probes on chrX or chrY
quantile(as.matrix(methCounts[,grepl("total", colnames(methCounts), fixed = TRUE)]), probes = c(0,0.25, 0.5, 0.75, 1), na.rm = TRUE)
# 0%  25%  50%  75% 100% 
# 1    4    6    9   32 probe coverage by sample
cellprobes <- data.frame("probe" = rownames(ModelPars), 
                         "cellType" = rep(c("Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC"), each = 100))
cellprobes <- subset(cellprobes, probe %in% methCounts$Name)
table(cellprobes$cellType)

rm(meth, cov, naCount, temp, i, probes_dup)

# Get Estimated Cell Counts ####
ModelPars <- subset(ModelPars, rownames(ModelPars) %in% rownames(permeth)) %>% as.matrix
cellCounts <- minfi:::projectCellType(permeth[rownames(ModelPars),], ModelPars) %>% as.data.frame
cellCounts$Sum <- rowSums2(as.matrix(cellCounts))
cellCounts$Sample <- rownames(cellCounts)
cellCounts <- merge(x = cellCounts, y = samples[,c("Sequencing_ID", "Diagnosis_Alg", "Sex", "Platform")], 
                    by.x = "Sample", by.y = "Sequencing_ID", all = FALSE)
cellCounts_scaled <- (as.matrix(cellCounts[,c("Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC")]) * 100 / 
        cellCounts$Sum) %>% as.data.frame
colnames(cellCounts_scaled) <- paste(colnames(cellCounts_scaled), "scaled", sep = "_")
cellCounts_scaled$Sum_scaled <- rowSums2(as.matrix(cellCounts_scaled))
cellCounts <- cbind(cellCounts[,c("Sample", "Diagnosis_Alg", "Sex", "Platform", "Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC", "Sum")], 
                    cellCounts_scaled)
rm(cellCounts_scaled)

# Analyze Cell Composition ------------------------------------------------
# Plots ####
# Composition by Sample Stacked Barplot
cellCounts$Diagnosis_Alg <- factor(cellCounts$Diagnosis_Alg, levels = c("TD", "ASD"))
cellCounts$Sex <- factor(cellCounts$Sex, levels = c("M", "F"))
cellCounts$Platform <- factor(cellCounts$Platform, levels = c("HiSeqX10", "HiSeq4000"))
sapply(cellCounts[,grepl("scaled", colnames(cellCounts), fixed = TRUE)], mean)
# Bcell_scaled  CD4T_scaled  CD8T_scaled  Gran_scaled  Mono_scaled    NK_scaled  nRBC_scaled   Sum_scaled 
#     6.322674    20.521061     7.258283    47.597714    10.039982     1.180118     7.080168   100.000000 
# Mean Proportion Gran > CD4T > Mono > CD8T > nRBC > Bcell > NK

cellCounts_m <- cellCounts[,c("Sample", "Diagnosis_Alg", "Sex", "Platform", "Bcell_scaled", "CD4T_scaled", "CD8T_scaled",
                              "Gran_scaled", "Mono_scaled", "NK_scaled", "nRBC_scaled")]
colnames(cellCounts_m) <- c("Sample", "Diagnosis", "Sex", "Platform", "B cells", "CD4 T cells", "CD8 T cells", "Granulocytes", 
                            "Monocytes", "NK cells", "nRBCs")
cellCounts_m <- melt(cellCounts_m, id.vars = c("Sample", "Diagnosis", "Sex", "Platform"))
colnames(cellCounts_m) <- c("Sample", "Diagnosis", "Sex", "Platform", "CellType", "Percent")

gg <- ggplot(cellCounts_m, aes(x = Sample, y = Percent, fill = CellType, color = CellType))
gg + 
        geom_bar(stat="identity") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks.x = element_blank(), legend.key = element_blank(), panel.grid.minor = element_blank(),
              legend.position = c(1.15, 0.82), legend.background = element_blank(), axis.text.x = element_blank(), 
              legend.key.size = unit(0.8, "cm"), strip.text.x = element_text(size = 19), 
              axis.ticks.y = element_line(size = 1.25, color = "black"), 
              legend.text = element_text(size = 16, margin = unit(c(0, 0, 0, 0.5), "lines")),
              strip.background = element_blank(), legend.direction = "vertical", panel.spacing.y = unit(0, "lines"), 
              plot.margin = unit(c(0.5, 10, 1, 1), "lines"), axis.title.x = element_blank(), 
              axis.text.y = element_text(size = 17, color = "black"), legend.title = element_blank()) +
        ylab("Cell Populations (%)") +
        scale_fill_manual(values = wes_palette(n = 7, name = "FantasticFox1", type = "continuous")) +
        scale_color_manual(values = wes_palette(n = 7, name = "FantasticFox1", type = "continuous")) +
        coord_cartesian(ylim = c(0, 100)) +
        scale_y_continuous(expand = c(0.003, 0)) +
        facet_grid(cols = vars(Diagnosis), scales = "free")
ggsave("Figures/Cell composition by sample stacked barplot.png", dpi = 600, width = 10, height = 7, units = "in")

# Composition by Diagnosis and Sex Stacked Bar Plot
cellCounts_agg <- aggregate(Percent ~ CellType + Diagnosis + Sex, data = cellCounts_m, FUN = mean)
gg <- ggplot(cellCounts_agg, aes(x = Sex, y = Percent, fill = CellType, color = CellType))
gg + 
        geom_bar(stat="identity") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), panel.grid.minor = element_blank(),
              legend.position = c(1.21, 0.81), legend.background = element_blank(), 
              legend.key.size = unit(0.8, "cm"), strip.text.x = element_text(size = 19), 
              axis.ticks = element_line(size = 1.25, color = "black"), 
              legend.text = element_text(size = 16, margin = unit(c(0, 0, 0, 0.5), "lines")),
              strip.background = element_blank(), legend.direction = "vertical", 
              panel.spacing.y = unit(0, "lines"), 
              plot.margin = unit(c(0.5, 10, 1, 1), "lines"), axis.title.x = element_blank(), 
              axis.text = element_text(size = 17, color = "black"), legend.title = element_blank()) +
        ylab("Cell Populations (%)") +
        scale_fill_manual(values = wes_palette(n = 7, name = "FantasticFox1", type = "continuous")) +
        scale_color_manual(values = wes_palette(n = 7, name = "FantasticFox1", type = "continuous")) +
        coord_cartesian(ylim = c(0, 100)) +
        scale_y_continuous(expand = c(0.003, 0)) +
        facet_grid(cols = vars(Diagnosis), scales = "free")
ggsave("Figures/Cell composition by diagnosis and sex stacked barplot.png", dpi = 600, width = 8, height = 7, units = "in")

# Cell Type by Diagnosis Box Plot
ggBoxPlot(data = cellCounts_m, x = cellCounts_m$Diagnosis, y = cellCounts_m$Percent, 
          fill = cellCounts_m$Diagnosis, ylab = "Cell Populations (%)", legend.name = "Diagnosis", 
          facet = vars(CellType), file = "Figures/Cell composition by diagnosis boxplot.png")

# Cell Type by Sex Box Plot
ggBoxPlot(data = cellCounts_m, x = cellCounts_m$Sex, y = cellCounts_m$Percent, 
          fill = cellCounts_m$Sex, ylab = "Cell Populations (%)", legend.name = "Sex", 
          facet = vars(CellType), file = "Figures/Cell composition by sex boxplot.png",
          legend.position = c(0.81, 0.38))

# Cell Type by Diagnosis and Sex Box Plot
ggBoxPlot(data = cellCounts_m, x = cellCounts_m$Sex, y = cellCounts_m$Percent, 
          fill = cellCounts_m$Diagnosis, ylab = "Cell Populations (%)", legend.name = "Diagnosis", 
          facet = vars(CellType), file = "Figures/Cell composition by diagnosis and Sex boxplot.png",
          axis.ticks.x = element_line(size = 1.25), axis.text.x = element_text(size = 16, color = "black"),
          legend.position = c(0.84, 0.37))

# Cell Type by Platform Box Plot
ggBoxPlot(data = cellCounts_m, x = cellCounts_m$Platform, y = cellCounts_m$Percent, 
          fill = cellCounts_m$Platform, ylab = "Cell Populations (%)", legend.name = "Platform", 
          facet = vars(CellType), file = "Figures/Cell composition by platform boxplot.png",
          legend.position = c(0.865, 0.38))

rm(cellCounts_agg, gg, ModelPars, permeth, cellCounts_m)

# Covariate Association ####
# Prep Data
cellCov <- cellCounts[,c("Sample", "Bcell_scaled", "CD4T_scaled", "CD8T_scaled", "Gran_scaled", "Mono_scaled", "NK_scaled",
                         "nRBC_scaled")]
colnames(cellCov) <- c("Sample", "Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC")
samples_cov <- merge(x = samples, y = cellCov, by.x = "Sequencing_ID", by.y = "Sample", all = FALSE, sort = FALSE)
catVars <- c("Study", "Platform", "Sex", "Site", "Diagnosis_Alg", "MomEdu_detail", "DM1or2", "GDM", "PE", 
             "home_ownership", "marital_status", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3", "AllEQ_PV_YN_Mo_2", 
             "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", "AllEQ_PV_YN_Mo4", 
             "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", "AllEQ_PV_YN_Mo9")
catVars <- catVars[!catVars == "Platform"]
samples_cov$Study <- factor(samples_cov$Study, levels = c("MARBLES", "EARLI"))
samples_cov$Platform <- factor(samples_cov$Platform, levels = c("HiSeqX10", "HiSeq4000"))
samples_cov$Sex <- factor(samples_cov$Sex, levels = c("M", "F"))
samples_cov$Site <- factor(samples_cov$Site, levels = c("UC Davis", "Kaiser Permanente", "Drexel", "Johns Hopkins University"))
samples_cov$Diagnosis_Alg <- factor(samples_cov$Diagnosis_Alg, levels = c("TD", "ASD"))
samples_cov$MomEdu_detail <- factor(samples_cov$MomEdu_detail, levels = c(6, 1:5, 7,8))
samples_cov$home_ownership[samples_cov$home_ownership == 99] <- NA
samples_cov$marital_status[samples_cov$marital_status == 99] <- NA
factorCols <- c("DM1or2", "GDM", "PE", "marital_status", "home_ownership", "SmokeYN_Pregnancy", "AllEQ_PV_YN_Mo_3", 
                "AllEQ_PV_YN_Mo_2", "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", 
                "AllEQ_PV_YN_Mo4", "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", 
                "AllEQ_PV_YN_Mo9")
samples_cov[,factorCols] <- lapply(samples_cov[,factorCols], as.factor)
contVars <- colnames(samples_cov)[!colnames(samples_cov) %in% catVars & !colnames(samples_cov) %in% colnames(cellCov) &
                                          !colnames(samples_cov) %in% c("COI_ID", "Sequencing_ID", "Cord_Blood_IBC",
                                                                        "Platform")]

# Get Stats
covStats <- DMRmethLm(DMRs = c("Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC"), catVars = catVars, 
                      contVars = contVars, sampleData = samples_cov, adj = "Platform",
                      file = "Tables/Estimated Cell Populations by Covariate Stats.txt")
covSum <- DMRmethLmSum(covStats, file = "Tables/Estimated Cell Populations by Covariate Summary.txt")

# Plot Heatmap
variables<-c("Diagnosis_Alg", "Sex", "Study", "Site_Drexel", "Site_Johns_Hopkins_University", "Site_Kaiser_Permanente", 
             "ADOScs", "MSLelcStandard36", "MSLelTscore36", "MSLfmTscore36", "MSLrlTscore36", "MSLvrTscore36", "ga_w", 
             "bw_g", "percent_trimmed", "percent_aligned", "percent_duplicate", "dedup_reads_M", "C_coverage", 
             "CG_coverage", "percent_cpg_meth", "percent_chg_meth", "percent_chh_meth", "MomEdu_detail_1", 
             "MomEdu_detail_2", "MomEdu_detail_3", "MomEdu_detail_4", "MomEdu_detail_5", "MomEdu_detail_7", 
             "MomEdu_detail_8", "home_ownership", "marital_status", "MomAgeYr", "Mat_Height_cm", "Mat_Weight_kg_PrePreg",
             "Mat_BMI_PrePreg", "DM1or2", "GDM", "PE", "parity", "dad_age", "SmokeYN_Pregnancy","cotinine_urine_ngml", 
             "final_creatinine_mgdl", "AllEQ_PV_YN_Mo_3", "AllEQ_PV_YN_Mo_2", "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", 
             "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", "AllEQ_PV_YN_Mo4", "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", 
             "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", "AllEQ_PV_YN_Mo9", "AllEQ_tot_All_FA_mcg_Mo_3", 
             "AllEQ_tot_All_FA_mcg_Mo_2", "AllEQ_tot_All_FA_mcg_Mo_1", "AllEQ_tot_All_FA_mcg_Mo1", 
             "AllEQ_tot_All_FA_mcg_Mo2", "AllEQ_tot_All_FA_mcg_Mo3", "AllEQ_tot_All_FA_mcg_Mo4", 
             "AllEQ_tot_All_FA_mcg_Mo5", "AllEQ_tot_All_FA_mcg_Mo6", "AllEQ_tot_All_FA_mcg_Mo7", 
             "AllEQ_tot_All_FA_mcg_Mo8","AllEQ_tot_All_FA_mcg_Mo9")
covHeatmap(covStats, variableOrdering = "manual", regionOrdering = "variable", variables = variables,
           sortVariable = "Diagnosis_Alg", probs = 0.997, axis.ticks.y = element_line(size = 1.25),
           axis.text.y = element_text(size = 11, color = "Black"), width = 11, height = 6, 
           legend.position = c(1.06, 0.835),
           file = "Figures/Estimated Cell Populations by Covariate Heatmap Sorted by Diagnosis.png")
covHeatmap(covStats, variableOrdering = "hierarchical", regionOrdering = "hierarchical", variables = variables,
           sortVariable = "Diagnosis_Alg", probs = 0.997, axis.ticks.y = element_line(size = 1.25),
           axis.text.y = element_text(size = 11, color = "Black"), width = 11, height = 6, 
           legend.position = c(1.06, 0.835),
           file = "Figures/Estimated Cell Populations by Covariate Heatmap Clustered.png")

# Plot nRBCs by Global mCG
ggScatterPlot(x = samples_cov$nRBC, y = samples_cov$percent_cpg_meth, groupVar = samples_cov$Diagnosis_Alg,
              fileName = "Figures/Estimated nRBCs by Global mCpG.png", xlab = "Estimated nRBCs (%)",
              ylab = "Global CpG Methylation (%)", legendPos = c(0.9, 1.03))

# Plot nRBCs by Mullen ELC
ggScatterPlot(x = samples_cov$nRBC, y = samples_cov$MSLelcStandard36, groupVar = samples_cov$Diagnosis_Alg,
              fileName = "Figures/Estimated nRBCs by Mullen Composite.png", xlab = "Estimated nRBCs (%)",
              ylab = "Mullen Composite Score", legendPos = c(0.9, 1.03), xlim = c(0, 23), ylim = c(40, 145))

# Plot nRBCs by Diagnosis and Sex
ggBoxPlot(data = samples_cov, x = samples_cov$Diagnosis_Alg, y = samples_cov$nRBC, fill = samples_cov$Diagnosis_Alg,
          ylab = "Estimated nRBCs (%)", legend.name = NULL, legend.position = "none", facet = vars(Sex), 
          file = "Figures/Estimated nRBCs by Diagnosis and Sex.png", width = 8, height = 7, nrow = NULL, ncol = 2,
          axis.ticks.x = element_line(size = 1.25), axis.text.x = element_text(size = 16, color = "black"),
          ylim = c(0, 40), outlier.size = 1.5)

summary(lm(nRBC ~ Diagnosis_Alg + Platform, data = samples_cov))$coefficients
#                   Estimate Std. Error  t value     Pr(>|t|)
# Diagnosis_AlgASD  1.675600  0.7866926 2.129930 3.479822e-02
# PlatformHiSeq4000 1.199557  0.8593466 1.395894 1.647955e-01

summary(lm(nRBC ~ Diagnosis_Alg + Platform, data = subset(samples_cov, Sex == "M")))$coefficients
#                   Estimate Std. Error  t value     Pr(>|t|)
# Diagnosis_AlgASD  2.429736  0.8800973 2.760758 6.749086e-03
# PlatformHiSeq4000 1.029075  0.9333405 1.102572 2.725981e-01

summary(lm(nRBC ~ Diagnosis_Alg + Platform, data = subset(samples_cov, Sex == "F")))$coefficients
#                     Estimate Std. Error    t value     Pr(>|t|)
# Diagnosis_AlgASD  -0.5075624   1.709197 -0.2969596 7.681582e-01
# PlatformHiSeq4000  1.6908438   2.136496  0.7914097 4.337493e-01

# Plot Gestational Age by Cell Populations
gaCells <- samples_cov[,c("Sequencing_ID", "Diagnosis_Alg", "ga_w", "Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC")]
colnames(gaCells) <- c("Sample", "Diagnosis", "Gestational Age", "B cells", "CD4 T cells", "CD8 T cells",
                       "Granulocytes", "Monocytes", "NK cells", "nRBCs")
gaCells <- melt(gaCells, id.vars = c("Sample", "Diagnosis", "Gestational Age"))
colnames(gaCells) <- c("Sample", "Diagnosis", "GestationalAge", "CellType", "Percent")

g <- ggplot(gaCells, aes(x = GestationalAge, y = Percent))
g + 
        geom_smooth(method="lm") +
        geom_point(aes(color = Diagnosis), size = 2) +
        theme_bw(base_size = 25) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), panel.grid.minor = element_blank(),
              legend.position = c(0.83, 0.34), legend.background = element_blank(),
              legend.key.size = unit(0.8, "cm"), strip.text.x = element_text(size = 19), 
              axis.ticks = element_line(size = 1.25), legend.title = element_text(size = 22),
              strip.background = element_blank(), legend.direction = "vertical", panel.spacing.y = unit(0, "lines"), 
              plot.margin = unit(c(0,1,1,0.4), "lines"), axis.title = element_text(size = 22, color = "black"),
              axis.text = element_text(size = 16, color = "black")) +
        scale_x_continuous(breaks = pretty_breaks(n = 5)) +
        scale_y_continuous(breaks = pretty_breaks(n = 5)) +
        xlab("Gestational Age (weeks)") +
        ylab("Cell Populations (%)") +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("#3366CC", "#FF3366")) +
        facet_wrap(vars(CellType), nrow = 2, scales = "free")
ggsave("Figures/Estimated Cell Populations by Gestational Age Scatterplots.png", dpi = 600, width = 10, height = 7, units = "in")

