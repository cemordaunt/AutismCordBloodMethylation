# Cell Type Composition Estimation ----------------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 6/1/19
# Excluded JLCM032B and JLCM050B
# Removed folate and added in bsseq global meth

# Packages ####
sapply(c("Biostrings", "AnnotationDbi", "openssl", "FlowSorted.CordBlood.450k", "minfi",
         "IlluminaHumanMethylation450kmanifest", "IlluminaHumanMethylation450kanno.ilmn12.hg19", "data.table", "tidyverse",
         "reshape2", "scales", "wesanderson", "RColorBrewer", "rlist"), require, character.only = TRUE)

# Functions ####
source("R Scripts/DMR Analysis Functions.R")

# Estimate Cell Composition --------------------------------------------------------------------
# Get 450K Probe Methylation in WGBS Samples ####
samples <- read.delim(file = "Merged Database/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo.txt", 
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE)
samples <- subset(samples, !(Cord_Blood_IBC == 101470 & Platform %in% c("HiSeq2500", "HiSeq4000")))
samples <- subset(samples, !Sequencing_ID %in% c("JLCM032B", "JLCM050B")) # Remove mislabeled females
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
methCounts <- methCounts[,!grepl("JLCM032B", colnames(methCounts), fixed = TRUE) & 
                                 !grepl("JLCM050B", colnames(methCounts), fixed = TRUE)] # Remove mislabeled females
meth <- methCounts[,grepl("meth", colnames(methCounts), fixed = TRUE)]
cov <- methCounts[,grepl("total", colnames(methCounts), fixed = TRUE)]
naCount <- apply(cov, 1, function(x) sum(is.na(x)))
table(naCount < nrow(samples) / 4) # 635 probes with na in < 25% of samples
table(naCount < nrow(samples) / 10) # 562 probes with na in < 10% of samples
table(naCount == 0) # 128 probes covered in all samples

cov <- subset(cov, naCount < nrow(samples) / 10) %>% as.matrix # Subset for probes covered in 90% of samples
meth <- subset(meth, naCount < nrow(samples) / 10) %>% as.matrix
permeth <- meth / cov
table(is.na(permeth)) # 2086
for(i in 1:nrow(permeth)){
        temp <- permeth[i,]
        temp[is.na(temp)] <- mean(temp, na.rm = TRUE) # Replace missing values with probe mean value
        permeth[i,] <- temp
}
table(is.na(permeth)) # 0
rownames(permeth) <- methCounts$Name[naCount < nrow(samples) / 10] # probe IDs
colnames(permeth) <- gsub("_meth", "", colnames(permeth), fixed = TRUE) # Sample SeqIDs
methCounts <- subset(methCounts, Name %in% rownames(permeth))
write.table(methCounts[,c("Chromosome", "Start", "End", "Name")], file = "Tables/Cord Cell Type Probes for Minfi.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
colnames(probes_dup) <- c("Chromosome", "Start", "End")
table(duplicated(rbind(methCounts[,1:3], probes_dup[,1:3]))) # probes used for cell estimation assigned to the same location
table(methCounts$Chromosome) # No probes on chrX or chrY
quantile(as.matrix(methCounts[,grepl("total", colnames(methCounts), fixed = TRUE)]), probes = c(0,0.25, 0.5, 0.75, 1), na.rm = TRUE)
# 0%  25%  50%  75% 100% 
# 1    3    6    9   32 probe coverage by sample
cellprobes <- data.frame("probe" = rownames(ModelPars), 
                         "cellType" = rep(c("Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC"), each = 100))
cellprobes <- subset(cellprobes, probe %in% methCounts$Name)
table(cellprobes$cellType)
# Bcell  CD4T  CD8T  Gran  Mono    NK  nRBC 
#    84    80    76    78    80    79    85
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

# Get Estimated Cell Counts in Reference ####
FlowSorted.CordBlood.450k <- updateObject(FlowSorted.CordBlood.450k)
reference <- preprocessIllumina(FlowSorted.CordBlood.450k) %>% mapToGenome(mergeManifest = FALSE) # hg19
refMeth <- reference[rownames(ModelPars),] %>% getBeta() 
refCellCounts <- minfi:::projectCellType(refMeth, ModelPars)
refCellCounts <- (refCellCounts / rowSums(refCellCounts)) %>% as.data.frame() # Scale to 100%
refpData <- pData(reference) %>% as.data.frame()
table(refpData$X == rownames(refCellCounts)) # TRUE
refCellCounts <- cbind(refpData, refCellCounts)
write.table(refCellCounts, file = "Tables/Cord Cell Type Reference Counts minfi.txt", sep = "\t", quote = FALSE, 
            row.names = FALSE)

cellTypes <- unique(refCellCounts$CellType) %>% as.character() %>% .[!. == "WholeBlood"]

# Proportion
sapply(cellTypes, function(x) mean(refCellCounts[refCellCounts$CellType == x,x]))
#     Bcell     CD4T     CD8T     Gran     Mono       NK     nRBC 
# 0.8981972 0.9004195 0.7309882 0.9505342 0.9398066 0.9170121 0.9257731 

# Whole Blood Proportion
sapply(cellTypes, function(x) mean(refCellCounts[refCellCounts$CellType == "WholeBlood",x]))
#      Bcell       CD4T       CD8T       Gran       Mono         NK       nRBC 
# 0.07624771 0.14065901 0.14484327 0.43910732 0.07744317 0.01164557 0.11005395

# Error Detail
expected <- c(0.08,0.15,0.1,0.45,0.1,0.02,0.1)
mapply(function(x,y) mean(abs(refCellCounts[refCellCounts$CellType == "WholeBlood",x] - y)), x = cellTypes, y = expected)
#      Bcell       CD4T       CD8T       Gran       Mono         NK       nRBC 
# 0.02195389 0.03851564 0.04779539 0.07607101 0.02788162 0.02364557 0.06178837

# Average Error Detail
sapply(cellTypes, function(x) mean(refCellCounts[refCellCounts$CellType == "WholeBlood",x])) - expected
#        Bcell         CD4T         CD8T         Gran         Mono           NK         nRBC 
# -0.003752287 -0.009340993  0.044843267 -0.010892678 -0.022556827 -0.008354428  0.010053947 

# Analyze Cell Composition ------------------------------------------------
# Plots ####
# Composition by Sample Stacked Barplot
cellCounts$Diagnosis_Alg <- factor(cellCounts$Diagnosis_Alg, levels = c("TD", "ASD"))
cellCounts$Sex <- factor(cellCounts$Sex, levels = c("M", "F"))
cellCounts$Platform <- factor(cellCounts$Platform, levels = c("HiSeqX10", "HiSeq4000"))
sapply(cellCounts[,grepl("scaled", colnames(cellCounts), fixed = TRUE)], mean)
# Bcell_scaled  CD4T_scaled  CD8T_scaled  Gran_scaled  Mono_scaled    NK_scaled  nRBC_scaled   Sum_scaled 
#     6.343706    20.605626     7.236536    47.518511    10.028071     1.171758     7.095792   100.000000 
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

# Composition by Diagnosis, Sex, and Platform Stacked Bar Plot
cellCounts_agg <- aggregate(Percent ~ CellType + Diagnosis + Sex + Platform, data = cellCounts_m, FUN = mean)
cellCounts_agg$SampleSet <- ifelse(cellCounts_agg$Platform == "HiSeqX10", yes = "Discovery", no = "Replication") %>%
        factor(levels = c("Discovery", "Replication"))
gg <- ggplot(cellCounts_agg, aes(x = Diagnosis, y = Percent, fill = CellType, color = CellType))
gg + 
        geom_bar(stat="identity") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), panel.grid.minor = element_blank(),
              legend.position = c(1.175, 0.748), legend.background = element_blank(), 
              legend.key.size = unit(0.8, "cm"), strip.text.x = element_text(size = 19), 
              axis.ticks = element_line(size = 1.25, color = "black"), 
              legend.text = element_text(size = 16, margin = unit(c(0, 0, 0, 0.5), "lines")),
              strip.background = element_blank(), legend.direction = "vertical", 
              panel.spacing.y = unit(0, "lines"), axis.title.y = element_text(size = 21),
              plot.margin = unit(c(0.5, 10, 1, 1), "lines"), axis.title.x = element_blank(), 
              axis.text = element_text(size = 17, color = "black"), legend.title = element_blank()) +
        ylab("Mean Cell Populations (%)") +
        scale_fill_manual(values = wes_palette(n = 7, name = "FantasticFox1", type = "continuous")) +
        scale_color_manual(values = wes_palette(n = 7, name = "FantasticFox1", type = "continuous")) +
        coord_cartesian(ylim = c(0, 100)) +
        scale_y_continuous(expand = c(0.004, 0)) +
        facet_grid(cols = vars(SampleSet, Sex), scales = "free")
ggsave("Figures/Cell composition by diagnosis, sex, and platform stacked barplot.png", dpi = 600, width = 9, height = 6, units = "in")

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

# Covariate Association (Update this to exclude folate variables and include bsseq global meth) ####

# Prep Data
cellCov <- cellCounts[,c("Sample", "Bcell_scaled", "CD4T_scaled", "CD8T_scaled", "Gran_scaled", "Mono_scaled", "NK_scaled",
                         "nRBC_scaled")]
colnames(cellCov) <- c("Sample", "Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC")
write.csv(cellCov, file = "Tables/Estimated Cell Proportions by Sample Discovery and Replication.csv", quote = FALSE, 
          row.names = FALSE)
samples_cov <- read.csv("Tables/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo and Cell Type.csv", header = TRUE,
                        stringsAsFactors = FALSE) 
# Covariates merged with cell type in "Global Methylation and Covariate Analysis ASD Cord Blood.R"
catVars <- c("Study", "Platform", "Sex", "Site", "Diagnosis_Alg", "MomEdu_detail", "DM1or2", "GDM", "PE", "home_ownership", 
             "marital_status", "SmokeYN_Pregnancy")
catVars <- catVars[!catVars == "Platform"]
samples_cov$Study <- factor(samples_cov$Study, levels = c("MARBLES", "EARLI"))
samples_cov$Platform <- factor(samples_cov$Platform, levels = c("HiSeqX10", "HiSeq4000"))
samples_cov$Sex <- factor(samples_cov$Sex, levels = c("M", "F"))
samples_cov$Site <- factor(samples_cov$Site, levels = c("UC Davis", "Kaiser Permanente", "Drexel", "Johns Hopkins University"))
samples_cov$Diagnosis_Alg <- factor(samples_cov$Diagnosis_Alg, levels = c("TD", "ASD"))
samples_cov$MomEdu_detail <- factor(samples_cov$MomEdu_detail, levels = c(6, 1:5, 7,8))
samples_cov$home_ownership[samples_cov$home_ownership == 99] <- NA
samples_cov$marital_status[samples_cov$marital_status == 99] <- NA
factorCols <- c("DM1or2", "GDM", "PE", "marital_status", "home_ownership", "SmokeYN_Pregnancy")
samples_cov[,factorCols] <- lapply(samples_cov[,factorCols], as.factor)
celltypes <- c("Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC")
contVars <- colnames(samples_cov)[!colnames(samples_cov) %in% catVars & 
                                          !colnames(samples_cov) %in% celltypes &
                                          !colnames(samples_cov) %in% c("COI_ID", "Sequencing_ID", "Cord_Blood_IBC", "Platform")]

# Get Stats
covStats <- DMRmethLm(DMRs = celltypes, catVars = catVars, contVars = contVars, sampleData = samples_cov, adj = "Platform",
                      file = "Tables/Estimated Cell Populations by Covariate Stats.txt")

# Plot Heatmap
variables<-c("Diagnosis_Alg", "Sex", "Study", "Site_Drexel", "Site_Johns_Hopkins_University", "Site_Kaiser_Permanente", 
             "ADOScs", "MSLelcStandard36", "MSLelTscore36", "MSLfmTscore36", "MSLrlTscore36", "MSLvrTscore36", "ga_w", 
             "bw_g", "percent_trimmed", "percent_aligned", "percent_duplicate", "dedup_reads_M", "C_coverage", 
             "CG_coverage", "percent_cpg_meth_bsseq", "percent_chg_meth", "percent_chh_meth", "MomEdu_detail_1", 
             "MomEdu_detail_2", "MomEdu_detail_3", "MomEdu_detail_4", "MomEdu_detail_5", "MomEdu_detail_7", 
             "MomEdu_detail_8", "home_ownership", "marital_status", "MomAgeYr", "Mat_Height_cm", "Mat_Weight_kg_PrePreg",
             "Mat_BMI_PrePreg", "DM1or2", "GDM", "PE", "parity", "dad_age", "SmokeYN_Pregnancy","cotinine_urine_ngml", 
             "final_creatinine_mgdl")
covHeatmap(covStats, variableOrdering = "manual", regionOrdering = "variable", variables = variables,
           sortVariable = "Diagnosis_Alg", probs = 0.997, axis.ticks = element_line(size = 0.9),
           axis.text.y = element_text(size = 11, color = "Black"), width = 11, height = 6, 
           legend.position = c(1.06, 0.835),
           file = "Figures/Estimated Cell Populations by Covariate Heatmap Sorted by Diagnosis.png")

# nRBCs Plots ####
# Plot nRBCs by Diagnosis and Sex
ggBoxPlot(data = samples_cov, x = samples_cov$Diagnosis_Alg, y = samples_cov$nRBC, fill = samples_cov$Diagnosis_Alg,
          ylab = "Estimated nRBCs (%)", legend.name = NULL, legend.position = "none", facet = vars(Sex), 
          file = "Figures/Estimated nRBCs by Diagnosis and Sex.png", width = 8, height = 7, nrow = NULL, ncol = 2,
          axis.ticks.x = element_line(size = 1.25), axis.text.x = element_text(size = 16, color = "black"),
          ylim = c(0, 40), outlier.size = 1.5)

# nRBC Dot Plot by Diagnosis, Sex, and Platform
samples_cov$SampleSet <- ifelse(samples_cov$Platform == "HiSeqX10", yes = "Discovery", no = "Replication") %>%
        factor(levels = c("Discovery", "Replication"))
samples_cov_dot <- samples_cov
samples_cov_dot$Sex <- ifelse(samples_cov_dot$Sex == "M", yes = "Males", no = "Females") %>% factor(levels = c("Males", "Females"))
g <- ggplot(data = samples_cov_dot)
g + 
        stat_summary(aes(x = Sex, y = nRBC, group = Diagnosis_Alg), fun.data = "mean_cl_boot", 
                     geom = "crossbar", color = "black", size = 0.7, position = "dodge2") +
        geom_dotplot(aes(x = Sex, y = nRBC, fill = Diagnosis_Alg, color = Diagnosis_Alg), binwidth = 0.6,
                     binaxis = "y", stackdir = "center", position = "dodge", stackratio = 1.15, dotsize = 0.8) + 
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.895, 0.96), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              strip.text = element_text(size = 22), plot.margin = unit(c(0.3, 1, 1, 1), "lines"), 
              axis.title.x = element_blank(), strip.background = element_blank()) +
        scale_y_continuous(breaks = pretty_breaks(n = 5)) +
        ylab("Estimated nRBCs (%)") +
        facet_wrap(vars(SampleSet)) +
        scale_fill_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366")) +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366"))
ggsave("Figures/Estimated nRBCs by Diagnosis, Sex, and Platform Dotplot Mean CL.png", dpi = 600, width = 9, 
       height = 7, units = "in")
rm(samples_cov_dot)

# Diagnosis, Global Methylation, and nRBCs Stats ####
samples_cov_pooled <- list(All = samples_cov, 
                           Males = subset(samples_cov, Sex == "M"), 
                           Females = subset(samples_cov, Sex == "F"))
samples_cov_disc <- list(All = subset(samples_cov, Platform == "HiSeqX10"), 
                         Males = subset(samples_cov, Sex == "M" & Platform == "HiSeqX10"), 
                         Females = subset(samples_cov, Sex == "F" & Platform == "HiSeqX10"))
samples_cov_rep <- list(All = subset(samples_cov, Platform == "HiSeq4000"), 
                        Males = subset(samples_cov, Sex == "M" & Platform == "HiSeq4000"), 
                        Females = subset(samples_cov, Sex == "F" & Platform == "HiSeq4000"))

# nRBCs ~ Diagnosis
nRBC_Dx <- list(
        Disc = sapply(samples_cov_disc, function(x){
                summary(lm(nRBC ~ Diagnosis_Alg, 
                           data = x))$coefficients["Diagnosis_AlgASD", c("Estimate", "Pr(>|t|)")]}),
        Rep = sapply(samples_cov_rep, function(x){
                summary(lm(nRBC ~ Diagnosis_Alg, 
                           data = x))$coefficients["Diagnosis_AlgASD", c("Estimate", "Pr(>|t|)")]}),
        Pool = sapply(samples_cov_pooled, function(x){
                summary(lm(nRBC ~ Diagnosis_Alg + Platform, 
                           data = x))$coefficients["Diagnosis_AlgASD", c("Estimate", "Pr(>|t|)")]}))

# nRBCs ~ Global Methylation
nRBC_Meth <- list(
        Disc = sapply(samples_cov_disc, function(x){
                summary(lm(nRBC ~ percent_cpg_meth_bsseq,
                           data = x))$coefficients["percent_cpg_meth_bsseq", c("Estimate", "Pr(>|t|)")]}),
        Rep = sapply(samples_cov_rep, function(x){
                summary(lm(nRBC ~ percent_cpg_meth_bsseq,
                           data = x))$coefficients["percent_cpg_meth_bsseq", c("Estimate", "Pr(>|t|)")]}),
        Pool = sapply(samples_cov_pooled, function(x){
                summary(lm(nRBC ~ percent_cpg_meth_bsseq + Platform, 
                           data = x))$coefficients["percent_cpg_meth_bsseq", c("Estimate", "Pr(>|t|)")]}))

# Global Methylation ~ Diagnosis + PCR Duplicates Stats ####
Meth_DxDup = list(
        Disc = sapply(samples_cov_disc, function(x){
                summary(lm(percent_cpg_meth_bsseq ~ Diagnosis_Alg + percent_duplicate,
                           data = x))$coefficients["Diagnosis_AlgASD", c("Estimate", "Pr(>|t|)")]}),
        Rep = sapply(samples_cov_rep, function(x){
                summary(lm(percent_cpg_meth_bsseq ~ Diagnosis_Alg + percent_duplicate,
                           data = x))$coefficients["Diagnosis_AlgASD", c("Estimate", "Pr(>|t|)")]}),
        Pool = sapply(samples_cov_pooled, function(x){
                summary(lm(percent_cpg_meth_bsseq ~ Diagnosis_Alg + percent_duplicate + Platform,
                           data = x))$coefficients["Diagnosis_AlgASD", c("Estimate", "Pr(>|t|)")]}))

# Global Methylation ~ Diagnosis + PCR Duplicates + nRBCs Stats ####
Meth_DxDup_nRBC = list(
        Disc = sapply(samples_cov_disc, function(x){
                summary(lm(percent_cpg_meth_bsseq ~ Diagnosis_Alg + percent_duplicate + nRBC,
                           data = x))$coefficients["Diagnosis_AlgASD", c("Estimate", "Pr(>|t|)")]}),
        Rep = sapply(samples_cov_rep, function(x){
                summary(lm(percent_cpg_meth_bsseq ~ Diagnosis_Alg + percent_duplicate + nRBC,
                           data = x))$coefficients["Diagnosis_AlgASD", c("Estimate", "Pr(>|t|)")]}),
        Pool = sapply(samples_cov_pooled, function(x){
                summary(lm(percent_cpg_meth_bsseq ~ Diagnosis_Alg + percent_duplicate + Platform + nRBC,
                           data = x))$coefficients["Diagnosis_AlgASD", c("Estimate", "Pr(>|t|)")]}))

# Make Table
estimates <- sapply(list(Meth_DxDup = Meth_DxDup, Meth_DxDup_nRBC = Meth_DxDup_nRBC, nRBC_Dx = nRBC_Dx, nRBC_Meth = nRBC_Meth), 
                    function(x){
                            sapply(x, function(y){
                                    y["Estimate",]})}) %>% t
colnames(estimates) <- paste(rep(c("Disc", "Rep", "Pool"), each = 3), rep(c("All", "Males", "Females"), times = 3), 
                             "Estimate", sep = "_")
estimates <- rbind(estimates[,grepl("Disc", colnames(estimates), fixed = TRUE)],
                   estimates[,grepl("Rep", colnames(estimates), fixed = TRUE)],
                   estimates[,grepl("Pool", colnames(estimates), fixed = TRUE)])

pvalues <- sapply(list(Meth_DxDup = Meth_DxDup, Meth_DxDup_nRBC = Meth_DxDup_nRBC, nRBC_Dx = nRBC_Dx, nRBC_Meth = nRBC_Meth), 
                    function(x){
                            sapply(x, function(y){
                                    y["Pr(>|t|)",]})}) %>% t
colnames(pvalues) <- paste(rep(c("Disc", "Rep", "Pool"), each = 3), rep(c("All", "Males", "Females"), times = 3), 
                           "Pvalue", sep = "_")
pvalues <- rbind(pvalues[,grepl("Disc", colnames(pvalues), fixed = TRUE)],
                   pvalues[,grepl("Rep", colnames(pvalues), fixed = TRUE)],
                   pvalues[,grepl("Pool", colnames(pvalues), fixed = TRUE)])

Dx_Meth_nRBC_stats <- cbind("Model" = rownames(estimates), "Set" = rep(c("Disc", "Rep", "Pool"), each = 4), 
                            estimates, pvalues) %>% as.data.frame
Dx_Meth_nRBC_stats <- Dx_Meth_nRBC_stats[order(Dx_Meth_nRBC_stats$Model),]
write.csv(Dx_Meth_nRBC_stats, file = "Tables/Diagnosis Global Meth and nRBCs Stats.csv", quote = FALSE, row.names = FALSE)

# nRBC ~ Mullen Composite Stats ####
# Discovery
sapply(samples_cov_disc, function(x){
        summary(lm(nRBC ~ MSLelcStandard36, data = x))$coefficients["MSLelcStandard36", c("Estimate", "Pr(>|t|)")]})
#                  All       Males      Females
# Estimate -0.03810434 -0.04735098 -0.004838788
# Pr(>|t|)  0.06479838  0.03362920  0.921925417

# Replication
sapply(samples_cov_rep, function(x){
        summary(lm(nRBC ~ MSLelcStandard36, data = x))$coefficients["MSLelcStandard36", c("Estimate", "Pr(>|t|)")]})
#                  All       Males      Females
# Estimate -0.03810434 -0.04735098 -0.004838788
# Pr(>|t|)  0.06479838  0.03362920  0.921925417

# Pooled
sapply(samples_cov_pooled, function(x){
        summary(lm(nRBC ~ MSLelcStandard36 + Platform, data = x))$coefficients["MSLelcStandard36", c("Estimate", "Pr(>|t|)")]})
#                   All        Males     Females
# Estimate -0.044744041 -0.050752655 -0.01577715
# Pr(>|t|)  0.005435093  0.003061058  0.71084293

# Test for reduced males ####
table(samples_cov$Sex)
#   M   F 
# 112  40 
samples_cov_males <- subset(samples_cov, Sex == "M")

reducedMales <- NULL
for(i in 1:1000){
        temp <- summary(lm(nRBC ~ Diagnosis_Alg + Platform, 
                           data = samples_cov_males[sample(1:nrow(samples_cov_males), 40),])
                        )$coefficients["Diagnosis_AlgASD",]
        reducedMales <- rbind(reducedMales, temp)
}
reducedMales <- as.data.frame(reducedMales)
rownames(reducedMales) <- 1:nrow(reducedMales)
colnames(reducedMales) <- c("Estimate", "StdError", "tvalue", "pvalue")

p <- 2 * ecdf(as.numeric(unlist(reducedMales$Estimate)))(-0.4697397) # 0.004
p <- 2 * (1-ecdf(as.numeric(unlist(reducedMales$pvalue)))(0.7847955)) # 0.034

# Plot Reduced Males Estimate
gg <- ggplot()
gg +
        geom_histogram(data = reducedMales, aes(x = Estimate), fill = "#3366CC", bins = 50) +
        geom_vline(xintercept = c(2.65, -0.47), color = "#FF3366", size = 1.25) +
        annotate("text", x = c(2.7, -0.41), y = 58, label = c("True Males", "True Females"), size = 6, hjust = 0) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), legend.key.size = unit(1, "lines"),
              panel.grid.minor = element_blank(), legend.position = c(0.9, 0.87), 
              legend.background = element_blank(), legend.text = element_text(size = 16, color = "Black"),
              plot.margin = unit(c(1,1,1,1), "lines"), 
              axis.text = element_text(color = "black", size = 16),
              axis.ticks = element_line(color = "black", size = 1.05),
              axis.title = element_text(color = "black", size = 18), 
              legend.title = element_text(color = "black", size = 18)) +
        xlab("ASD - TD nRBC Estimate (%)") +
        ylab("Random Male Samples") +
        scale_x_continuous(breaks = pretty_breaks(n = 4)) +
        scale_y_continuous(breaks = pretty_breaks(n = 4), expand = c(0.002, 0)) +
        coord_cartesian(ylim = c(0, 60))
ggsave("Figures/Reduced Males nRBC Difference Estimate Histogram.png", dpi = 600, width = 8, height = 8, units = "in")

# Plot Reduced Males pvalue
gg <- ggplot()
gg +
        geom_histogram(data = reducedMales, aes(x = pvalue), fill = "#3366CC", bins = 50) +
        geom_vline(xintercept = c(0.003, 0.785), color = "#FF3366", size = 1.25) +
        annotate("text", x = c(0.04, 0.8), y = 125, label = c("True Males", "True Females"), size = 6, hjust = 0) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), legend.key.size = unit(1, "lines"),
              panel.grid.minor = element_blank(), legend.position = c(0.9, 0.87), 
              legend.background = element_blank(), legend.text = element_text(size = 16, color = "Black"),
              plot.margin = unit(c(1,1,1,1), "lines"), 
              axis.text = element_text(color = "black", size = 16),
              axis.ticks = element_line(color = "black", size = 1.05),
              axis.title = element_text(color = "black", size = 18), 
              legend.title = element_text(color = "black", size = 18)) +
        xlab("ASD - TD nRBC p-value") +
        ylab("Random Male Samples") +
        scale_x_continuous(breaks = pretty_breaks(n = 4)) +
        scale_y_continuous(breaks = pretty_breaks(n = 6), expand = c(0.002, 0)) +
        coord_cartesian(xlim = c(0,1), ylim = c(0, 130))
ggsave("Figures/Reduced Males nRBC Difference pvalue Histogram.png", dpi = 600, width = 8, height = 8, units = "in")

# Covariate Plots ####
# Plot nRBCs by Global mCG
ggScatterPlot(x = samples_cov$nRBC, y = samples_cov$percent_cpg_meth_bsseq, groupVar = samples_cov$Diagnosis_Alg,
              fileName = "Figures/Estimated nRBCs by Global mCpG.png", xlab = "Estimated nRBCs (%)",
              ylab = "Global CpG Methylation (%)", legendPos = c(0.9, 1.03))

# Plot nRBCs by Global mCG and Sex
g <- ggplot(samples_cov, aes(x = nRBC, y = percent_cpg_meth_bsseq))
g + 
        geom_smooth(method = "lm") +
        geom_point(aes(color = Diagnosis_Alg), size = 3) +
        facet_grid(cols = vars(Sex)) +
        theme_bw(base_size = 25) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), panel.grid.minor = element_blank(),
              legend.position = c(0.92, 1.05), legend.background = element_blank(),
              legend.key.size = unit(0.8, "cm"), strip.text.x = element_text(size = 22), 
              axis.ticks = element_line(size = 1.25), legend.title = element_blank(),
              strip.background = element_blank(), legend.direction = "horizontal", panel.spacing.y = unit(0, "lines"), 
              plot.margin = unit(c(0,1,1,0.4), "lines"), axis.title = element_text(size = 22, color = "black"),
              axis.text = element_text(size = 17, color = "black")) +
        scale_x_continuous(breaks = pretty_breaks(n = 5)) +
        scale_y_continuous(breaks = pretty_breaks(n = 5)) +
        xlab("Estimated nRBCs (%)") +
        ylab("Global CpG Methylation (%)") +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("#3366CC", "#FF3366"))
ggsave("Figures/Estimated nRBCs by Global mCpG and Sex.png", dpi = 600, width = 10, height = 6, units = "in")

# Plot nRBCs by Global mCG and Platform
g <- ggplot(samples_cov, aes(x = nRBC, y = percent_cpg_meth_bsseq))
g + 
        geom_smooth(method = "lm") +
        geom_point(aes(color = Diagnosis_Alg), size = 3) +
        facet_grid(cols = vars(SampleSet)) +
        theme_bw(base_size = 25) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), panel.grid.minor = element_blank(),
              legend.position = c(0.92, 1.05), legend.background = element_blank(),
              legend.key.size = unit(0.8, "cm"), strip.text.x = element_text(size = 22), 
              axis.ticks = element_line(size = 1.25), legend.title = element_blank(),
              strip.background = element_blank(), legend.direction = "horizontal", panel.spacing.y = unit(0, "lines"), 
              plot.margin = unit(c(0,1,1,0.4), "lines"), axis.title = element_text(size = 22, color = "black"),
              axis.text = element_text(size = 17, color = "black")) +
        scale_x_continuous(breaks = pretty_breaks(n = 5)) +
        scale_y_continuous(breaks = pretty_breaks(n = 5)) +
        xlab("Estimated nRBCs (%)") +
        ylab("Global CpG Methylation (%)") +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("#3366CC", "#FF3366"))
ggsave("Figures/Estimated nRBCs by Global mCpG and Platform.png", dpi = 600, width = 10, height = 6, units = "in")

# Plot nRBCs by Global mCG and Platform, Males Only
g <- ggplot(subset(samples_cov, Sex == "M"), aes(x = nRBC, y = percent_cpg_meth_bsseq))
g + 
        geom_smooth(method = "lm") +
        geom_point(aes(color = Diagnosis_Alg), size = 3) +
        facet_grid(cols = vars(SampleSet)) +
        theme_bw(base_size = 25) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), panel.grid.minor = element_blank(),
              legend.position = c(0.895, 1.12), legend.background = element_blank(),
              legend.key.size = unit(0.8, "cm"), strip.text.x = element_text(size = 22), 
              axis.ticks = element_line(size = 1.25), legend.title = element_blank(),
              strip.background = element_blank(), legend.direction = "horizontal", panel.spacing.y = unit(0, "lines"), 
              plot.margin = unit(c(1,1,1,0.4), "lines"), axis.title = element_text(size = 22, color = "black"),
              axis.text = element_text(size = 17, color = "black")) +
        scale_x_continuous(breaks = pretty_breaks(n = 5)) +
        scale_y_continuous(breaks = pretty_breaks(n = 5)) +
        xlab("Estimated nRBCs (%)") +
        ylab("Global CpG Methylation (%)") +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("#3366CC", "#FF3366"))
ggsave("Figures/Estimated nRBCs by Global mCpG and Platform, Males.png", dpi = 600, width = 8.5, height = 6.07, units = "in")

# Plot nRBCs by Global mCG, Sex, and Platform
samples_cov_plot <- samples_cov
samples_cov_plot$Sex <- ifelse(samples_cov_plot$Sex == "M", yes = "Males", no = "Females") %>% 
        factor(levels = c("Males", "Females"))
g <- ggplot(samples_cov_plot, aes(x = nRBC, y = percent_cpg_meth_bsseq))
g + 
        geom_smooth(method = "lm") +
        geom_point(aes(color = Diagnosis_Alg), size = 3) +
        facet_grid(rows = vars(SampleSet), cols = vars(Sex), scales = "free") +
        theme_bw(base_size = 25) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), panel.grid.minor = element_blank(),
              legend.position = c(0.89, 1.08), legend.background = element_blank(),
              legend.key.size = unit(0.8, "cm"), strip.text.x = element_text(size = 20), 
              axis.ticks = element_line(size = 1.25), legend.title = element_blank(),
              strip.background = element_blank(), legend.direction = "horizontal", panel.spacing.y = unit(1, "lines"), 
              plot.margin = unit(c(1,0.25,0.5,0.5), "lines"), axis.title = element_text(size = 20, color = "black"),
              axis.text = element_text(size = 16, color = "black")) +
        scale_x_continuous(breaks = pretty_breaks(n = 5)) +
        scale_y_continuous(breaks = pretty_breaks(n = 5), expand = c(0.08, 0)) +
        xlab("Estimated nRBCs (%)") +
        ylab("Global CpG Methylation (%)") +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("#3366CC", "#FF3366"))
ggsave("Figures/Estimated nRBCs by Global mCpG, Sex, and Platform.png", dpi = 600, width = 8, height = 8, units = "in")

# Plot nRBCs by Mullen ELC
ggScatterPlot(x = samples_cov$nRBC, y = samples_cov$MSLelcStandard36, groupVar = samples_cov$Diagnosis_Alg,
              fileName = "Figures/Estimated nRBCs by Mullen Composite.png", xlab = "Estimated nRBCs (%)",
              ylab = "Mullen Composite Score", legendPos = c(0.9, 1.03), xlim = c(0, 23), ylim = c(40, 145))

# Plot nRBCs by Mullen ELC and Sex
g <- ggplot(samples_cov, aes(x = nRBC, y = MSLelcStandard36))
g + 
        geom_smooth(method = "lm") +
        geom_point(aes(color = Diagnosis_Alg), size = 3) +
        facet_grid(cols = vars(Sex)) +
        theme_bw(base_size = 25) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), panel.grid.minor = element_blank(),
              legend.position = c(0.92, 1.05), legend.background = element_blank(),
              legend.key.size = unit(0.8, "cm"), strip.text.x = element_text(size = 22), 
              axis.ticks = element_line(size = 1.25), legend.title = element_blank(),
              strip.background = element_blank(), legend.direction = "horizontal", panel.spacing.y = unit(0, "lines"), 
              plot.margin = unit(c(0,1,1,0.4), "lines"), axis.title = element_text(size = 22, color = "black"),
              axis.text = element_text(size = 17, color = "black")) +
        scale_x_continuous(breaks = pretty_breaks(n = 5)) +
        scale_y_continuous(breaks = pretty_breaks(n = 5)) +
        coord_cartesian(xlim = c(0, 23), ylim = c(40, 145)) +
        xlab("Estimated nRBCs (%)") +
        ylab("Mullen Composite Score") +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("#3366CC", "#FF3366"))
ggsave("Figures/Estimated nRBCs by Mullen Composite and Sex.png", dpi = 600, width = 10, height = 6, units = "in")

# Plot nRBCs by Mullen ELC and Platform, Males Only
g <- ggplot(subset(samples_cov, Sex == "M"), aes(x = nRBC, y = MSLelcStandard36))
g + 
        geom_smooth(method = "lm") +
        geom_point(aes(color = Diagnosis_Alg), size = 3) +
        facet_grid(cols = vars(SampleSet)) +
        theme_bw(base_size = 25) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), panel.grid.minor = element_blank(),
              legend.position = c(0.92, 1.05), legend.background = element_blank(),
              legend.key.size = unit(0.8, "cm"), strip.text.x = element_text(size = 22), 
              axis.ticks = element_line(size = 1.25), legend.title = element_blank(),
              strip.background = element_blank(), legend.direction = "horizontal", panel.spacing.y = unit(0, "lines"), 
              plot.margin = unit(c(0,1,1,0.4), "lines"), axis.title = element_text(size = 22, color = "black"),
              axis.text = element_text(size = 17, color = "black")) +
        scale_x_continuous(breaks = pretty_breaks(n = 5)) +
        scale_y_continuous(breaks = pretty_breaks(n = 5)) +
        xlab("Estimated nRBCs (%)") +
        ylab("Mullen Composite Score") +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("#3366CC", "#FF3366"))
ggsave("Figures/Estimated nRBCs by Mullen Composite and Platform, Males.png", dpi = 600, width = 10, height = 6, units = "in")

# Plot nRBCs by Mullen ELC, Sex, and Platform
g <- ggplot(samples_cov_plot, aes(x = nRBC, y = MSLelcStandard36))
g + 
        geom_smooth(method = "lm") +
        geom_point(aes(color = Diagnosis_Alg), size = 3) +
        facet_grid(rows = vars(SampleSet), cols = vars(Sex), scales = "free") +
        theme_bw(base_size = 25) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), panel.grid.minor = element_blank(),
              legend.position = c(0.89, 1.08), legend.background = element_blank(),
              legend.key.size = unit(0.8, "cm"), strip.text.x = element_text(size = 20), 
              axis.ticks = element_line(size = 1.25), legend.title = element_blank(),
              strip.background = element_blank(), legend.direction = "horizontal", panel.spacing.y = unit(1, "lines"), 
              plot.margin = unit(c(1,0.25,0.5,0.5), "lines"), axis.title = element_text(size = 20, color = "black"),
              axis.text = element_text(size = 16, color = "black")) +
        scale_x_continuous(breaks = pretty_breaks(n = 5)) +
        scale_y_continuous(breaks = pretty_breaks(n = 5), expand = c(0.08, 0)) +
        xlab("Estimated nRBCs (%)") +
        ylab("Mullen Early Learning Composite") +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("#3366CC", "#FF3366"))
ggsave("Figures/Estimated nRBCs by Mullen Composite, Sex, and Platform.png", dpi = 600, width = 8, height = 8, units = "in")

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

# Covariate Association, Males Only ####
# Data
samples_cov_males <- subset(samples_cov, Sex == "M")
catVars <- catVars[!catVars == "Sex"]
variables <- variables[!variables == "Sex"]

# Get Stats
covStats_males <- DMRmethLm(DMRs = celltypes, catVars = catVars, contVars = contVars, sampleData = samples_cov_males, 
                            adj = "Platform", file = "Tables/Estimated Cell Populations by Covariate Stats Males Only.txt")

# Get Stats for nRBCs only
covStats_males_nRBCs <- DMRmethLm(DMRs = "nRBC", catVars = catVars, contVars = contVars, sampleData = samples_cov_males, 
                            adj = "Platform", file = "Tables/Estimated nRBCs by Covariate Stats Males Only.txt")

# Covariate Association, Females Only ####
# Data
samples_cov_females <- subset(samples_cov, Sex == "F")

# Get Stats
covStats_females <- DMRmethLm(DMRs = celltypes, catVars = catVars, 
                            contVars = contVars, sampleData = samples_cov_females, adj = "Platform",
                            file = "Tables/Estimated Cell Populations by Covariate Stats Females Only.txt")

# Get Stats for nRBCs only
covStats_females_nRBCs <- DMRmethLm(DMRs = "nRBC", catVars = catVars, 
                              contVars = contVars, sampleData = samples_cov_females, adj = "Platform",
                              file = "Tables/Estimated nRBCs by Covariate Stats Females Only.txt")

rm(gaCells, samples_cov_disc, samples_cov_rep, samples_cov_pooled)

# Covariate Association, Discovery Samples --------------------------------
# All Discovery Samples ####
# Data
samples_cov_disc <- subset(samples_cov, Platform == "HiSeqX10")
catVars <- c("Study", "Platform", "Sex", "Site", "Diagnosis_Alg", "MomEdu_detail", "DM1or2", "GDM", "PE", 
             "home_ownership", "marital_status", "SmokeYN_Pregnancy")
catVars <- catVars[!catVars == "Platform"]

# Get Stats
covStats_disc <- DMRmethLm(DMRs = celltypes, catVars = catVars, contVars = contVars, 
                           sampleData = samples_cov_disc, adj = NULL,
                           file = "Tables/Estimated Cell Populations by Covariate Stats Discovery.txt")

# Discovery Males Only ####
# Data
samples_cov_disc_males <- subset(samples_cov_disc, Sex == "M")
catVars <- catVars[!catVars == "Sex"]

# Get Stats
covStats_disc_males <- DMRmethLm(DMRs = celltypes, catVars = catVars, contVars = contVars, 
                                 sampleData = samples_cov_disc_males, adj = NULL,
                                 file = "Tables/Estimated Cell Populations by Covariate Stats Discovery Males Only.txt")
# Error with Bcell and PE
# Error with CD4T and PE
# Error with CD8T and PE
# Error with Gran and PE
# Error with Mono and PE
# Error with NK and PE
# Error with nRBC and PE

# Discovery Females Only ####
# Data
samples_cov_disc_females <- subset(samples_cov_disc, Sex == "F")

# Get Stats
covStats_disc_females <- DMRmethLm(DMRs = celltypes, catVars = catVars, contVars = contVars, 
                                   sampleData = samples_cov_disc_females, adj = NULL,
                                   file = "Tables/Estimated Cell Populations by Covariate Stats Discovery Females Only.txt")

# Plots ####
# Plot nRBCs by Diagnosis and Sex
ggBoxPlot(data = samples_cov_disc, x = samples_cov_disc$Diagnosis_Alg, y = samples_cov_disc$nRBC, 
          fill = samples_cov_disc$Diagnosis_Alg, ylab = "Estimated nRBCs (%)", legend.name = NULL, 
          legend.position = "none", facet = vars(Sex), file = "Figures/Estimated nRBCs by Diagnosis and Sex Discovery.png", 
          width = 8, height = 7, nrow = NULL, ncol = 2, axis.ticks.x = element_line(size = 1.25), 
          axis.text.x = element_text(size = 16, color = "black"), ylim = c(0, 40), outlier.size = 1.5)

# Plot nRBCs by Global mCG
ggScatterPlot(x = samples_cov_disc$nRBC, y = samples_cov_disc$percent_cpg_meth_bsseq, groupVar = samples_cov_disc$Diagnosis_Alg,
              fileName = "Figures/Estimated nRBCs by Global mCpG Discovery.png", xlab = "Estimated nRBCs (%)",
              ylab = "Global CpG Methylation (%)", legendPos = c(0.9, 1.03))

# Plot nRBCs by Global mCG and Sex
g <- ggplot(samples_cov_disc, aes(x = nRBC, y = percent_cpg_meth_bsseq))
g + 
        geom_smooth(method = "lm") +
        geom_point(aes(color = Diagnosis_Alg), size = 3) +
        facet_grid(cols = vars(Sex)) +
        theme_bw(base_size = 25) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), panel.grid.minor = element_blank(),
              legend.position = c(0.92, 1.05), legend.background = element_blank(),
              legend.key.size = unit(0.8, "cm"), strip.text.x = element_text(size = 22), 
              axis.ticks = element_line(size = 1.25), legend.title = element_blank(),
              strip.background = element_blank(), legend.direction = "horizontal", panel.spacing.y = unit(0, "lines"), 
              plot.margin = unit(c(0,1,1,0.4), "lines"), axis.title = element_text(size = 22, color = "black"),
              axis.text = element_text(size = 17, color = "black")) +
        scale_x_continuous(breaks = pretty_breaks(n = 5)) +
        scale_y_continuous(breaks = pretty_breaks(n = 5)) +
        xlab("Estimated nRBCs (%)") +
        ylab("Global CpG Methylation (%)") +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("#3366CC", "#FF3366"))
ggsave("Figures/Estimated nRBCs by Global mCpG and Sex Discovery.png", dpi = 600, width = 10, height = 6, units = "in")

# Plot Gestational Age by Cell Populations
gaCells <- samples_cov_disc[,c("Sequencing_ID", "Diagnosis_Alg", "ga_w", "Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC")]
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
ggsave("Figures/Estimated Cell Populations by Gestational Age Scatterplots Discovery.png", dpi = 600, width = 10, height = 7, units = "in")

# Covariate Association, Replication Samples --------------------------------
# All Replication Samples ####
# Data
samples_cov_rep <- subset(samples_cov, Platform == "HiSeq4000")
catVars <- c("Study", "Platform", "Sex", "Site", "Diagnosis_Alg", "MomEdu_detail", "DM1or2", "GDM", "PE", 
             "home_ownership", "marital_status", "SmokeYN_Pregnancy")
catVars <- catVars[!catVars %in% c("Platform", "Study", "Site")]

# Get Stats
covStats_rep <- DMRmethLm(DMRs = celltypes, catVars = catVars, contVars = contVars, 
                          sampleData = samples_cov_rep, adj = NULL,
                          file = "Tables/Estimated Cell Populations by Covariate Stats Replication.txt")

# Replication Males Only ####
# Data
samples_cov_rep_males <- subset(samples_cov_rep, Sex == "M")
catVars <- catVars[!catVars == "Sex"]

# Get Stats
covStats_rep_males <- DMRmethLm(DMRs = celltypes, catVars = catVars, contVars = contVars, 
                                sampleData = samples_cov_rep_males, adj = NULL,
                                file = "Tables/Estimated Cell Populations by Covariate Stats Replication Males Only.txt")
# Error with DM1or2 for all

# Replication Females Only ####
# Data
samples_cov_rep_females <- subset(samples_cov_rep, Sex == "F")

# Get Stats
covStats_rep_females <- DMRmethLm(DMRs = celltypes, catVars = catVars, contVars = contVars, 
                                  sampleData = samples_cov_rep_females, adj = NULL,
                                  file = "Tables/Estimated Cell Populations by Covariate Stats Replication Females Only.txt")
# Error with PE, SmokeYN_Pregnancy for all

# Plots ####
# Plot nRBCs by Diagnosis and Sex
ggBoxPlot(data = samples_cov_rep, x = samples_cov_rep$Diagnosis_Alg, y = samples_cov_rep$nRBC, 
          fill = samples_cov_rep$Diagnosis_Alg, ylab = "Estimated nRBCs (%)", legend.name = NULL, 
          legend.position = "none", facet = vars(Sex), file = "Figures/Estimated nRBCs by Diagnosis and Sex Replication.png", 
          width = 8, height = 7, nrow = NULL, ncol = 2, axis.ticks.x = element_line(size = 1.25), 
          axis.text.x = element_text(size = 16, color = "black"), ylim = c(0, 40), outlier.size = 1.5)

# Plot nRBCs by Global mCG
ggScatterPlot(x = samples_cov_rep$nRBC, y = samples_cov_rep$percent_cpg_meth_bsseq, groupVar = samples_cov_rep$Diagnosis_Alg,
              fileName = "Figures/Estimated nRBCs by Global mCpG Replication.png", xlab = "Estimated nRBCs (%)",
              ylab = "Global CpG Methylation (%)", legendPos = c(0.9, 1.03))

# Plot nRBCs by Global mCG and Sex
g <- ggplot(samples_cov_rep, aes(x = nRBC, y = percent_cpg_meth_bsseq))
g + 
        geom_smooth(method = "lm") +
        geom_point(aes(color = Diagnosis_Alg), size = 3) +
        facet_grid(cols = vars(Sex)) +
        theme_bw(base_size = 25) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), panel.grid.minor = element_blank(),
              legend.position = c(0.92, 1.05), legend.background = element_blank(),
              legend.key.size = unit(0.8, "cm"), strip.text.x = element_text(size = 22), 
              axis.ticks = element_line(size = 1.25), legend.title = element_blank(),
              strip.background = element_blank(), legend.direction = "horizontal", panel.spacing.y = unit(0, "lines"), 
              plot.margin = unit(c(0,1,1,0.4), "lines"), axis.title = element_text(size = 22, color = "black"),
              axis.text = element_text(size = 17, color = "black")) +
        scale_x_continuous(breaks = pretty_breaks(n = 5)) +
        scale_y_continuous(breaks = pretty_breaks(n = 5)) +
        xlab("Estimated nRBCs (%)") +
        ylab("Global CpG Methylation (%)") +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("#3366CC", "#FF3366"))
ggsave("Figures/Estimated nRBCs by Global mCpG and Sex Replication.png", dpi = 600, width = 10, height = 6, units = "in")

# Plot Gestational Age by Cell Populations
gaCells <- samples_cov_rep[,c("Sequencing_ID", "Diagnosis_Alg", "ga_w", "Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC")]
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
        scale_x_continuous(breaks = pretty_breaks(n = 3)) +
        scale_y_continuous(breaks = pretty_breaks(n = 5)) +
        xlab("Gestational Age (weeks)") +
        ylab("Cell Populations (%)") +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("#3366CC", "#FF3366")) +
        facet_wrap(vars(CellType), nrow = 2, scales = "free")
ggsave("Figures/Estimated Cell Populations by Gestational Age Scatterplots Replication.png", dpi = 600, width = 10, height = 7,
       units = "in")
