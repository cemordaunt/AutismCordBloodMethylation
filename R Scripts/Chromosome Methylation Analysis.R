# Chromosome Methylation Analysis ---------------------------------------------------
# ASD Cord Blood Methylation
# Charles Mordaunt
# 3/1/19

# Packages ####
.libPaths("/share/lasallelab/programs/DMRichR/R_3.5")
sapply(c("tidyverse", "openxlsx", "bsseq", "dmrseq", "DMRichR", "rlist", "reshape2", "scales"), require, character.only = TRUE)

# Functions ####
source("R Scripts/DMR Analysis Functions.R")

# Get Methylation Data ----------------------------------------------------------
# Discovery ####
# Load Methylation Data
meta <- read.xlsx("sample_info.xlsx", colNames = TRUE) %>% mutate_if(is.character,as.factor)
files <- list.files(path = getwd(), pattern = "*.txt.gz")
files <- files[pmatch(meta$Name, files)]
bs <- read.bismark(files = files, rmZeroCov = FALSE, strandCollapse = TRUE, verbose = TRUE, 
                   BPPARAM = MulticoreParam(workers = 30, progressbar = TRUE), nThread = 1)
sampleNames(bs) <- gsub("_.*$", "", sampleNames(bs))
meta <- meta[order(match(meta[, 1], sampleNames(bs))), ]
pData(bs) <- cbind(pData(bs), meta[2:length(meta)])
pData(bs)[["Diagnosis"]] <- as.factor(pData(bs)[["Diagnosis"]])
bs <- GenomeInfoDb::keepStandardChromosomes(bs, pruning.mode = "coarse")

# Get Methylated and Total CpG Counts by Chromosome
bs_split <- split(bs, f = seqnames(bs))
bs_cov <- sapply(bs_split, function(x){
        DelayedMatrixStats::colSums2(getCoverage(BSseq = x, type = "Cov", what = "perBase"))
        })
bs_cov <- as.data.frame(cbind(sampleNames(bs), pData(bs), bs_cov))
write.csv(bs_cov, file = "Chromosome_Total_CpG_Counts_Discovery.csv", quote = FALSE, row.names = FALSE)

bs_meth <- sapply(bs_split, function(x){
        DelayedMatrixStats::colSums2(getCoverage(BSseq = x, type = "M", what = "perBase"))
})
bs_meth <- as.data.frame(cbind(sampleNames(bs), pData(bs), bs_meth))
write.csv(bs_meth, file = "Chromosome_Meth_CpG_Counts_Discovery.csv", quote = FALSE, row.names = FALSE)

# Replication ####
# Load Methylation Data
meta <- read.xlsx("sample_info.xlsx", colNames = TRUE) %>% mutate_if(is.character,as.factor)
files <- list.files(path = getwd(), pattern = "*.txt.gz")
files <- files[pmatch(meta$Name, files)]
bs <- read.bismark(files = files, rmZeroCov = FALSE, strandCollapse = TRUE, verbose = TRUE, 
                   BPPARAM = MulticoreParam(workers = 30, progressbar = TRUE), nThread = 1)
sampleNames(bs) <- gsub("_.*$", "", sampleNames(bs))
meta <- meta[order(match(meta[, 1], sampleNames(bs))), ]
pData(bs) <- cbind(pData(bs), meta[2:length(meta)])
pData(bs)[["Diagnosis"]] <- as.factor(pData(bs)[["Diagnosis"]])
bs <- GenomeInfoDb::keepStandardChromosomes(bs, pruning.mode = "coarse")

# Get Methylated and Total CpG Counts by Chromosome
bs_split <- split(bs, f = seqnames(bs))
bs_cov <- sapply(bs_split, function(x){
        DelayedMatrixStats::colSums2(getCoverage(BSseq = x, type = "Cov", what = "perBase"))
})
bs_cov <- as.data.frame(cbind(sampleNames(bs), pData(bs), bs_cov))
write.csv(bs_cov, file = "Chromosome_Total_CpG_Counts_Replication.csv", quote = FALSE, row.names = FALSE)

bs_meth <- sapply(bs_split, function(x){
        DelayedMatrixStats::colSums2(getCoverage(BSseq = x, type = "M", what = "perBase"))
})
bs_meth <- as.data.frame(cbind(sampleNames(bs), pData(bs), bs_meth))
write.csv(bs_meth, file = "Chromosome_Meth_CpG_Counts_Replication.csv", quote = FALSE, row.names = FALSE)

# Methylation by Diagnosis and Sex ----------------------------------------
# Discovery ####
# Data
meth <- read.csv("Tables/Chromosome_Meth_CpG_Counts_Discovery.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(meth)[1] <- "Sample"
meth$Diagnosis[meth$Diagnosis == "Ctrl_TD"] <- "TD"
meth$Diagnosis[meth$Diagnosis == "Exp_ASD"] <- "ASD"
meth$Diagnosis <- factor(meth$Diagnosis, levels = c("TD", "ASD"))
meth$Sex <- factor(meth$Sex, levels = c("M", "F"))

cov <- read.csv("Tables/Chromosome_Total_CpG_Counts_Discovery.csv", header = TRUE, stringsAsFactors = FALSE)
permeth <- cbind(meth[,c("Sample", "Diagnosis", "Sex", "percent_duplicate")],
                 meth[,grepl("chr", colnames(meth), fixed = TRUE)] * 100 / cov[,grepl("chr", colnames(cov), fixed = TRUE)])

global <- cbind(meth[,c("Sample", "Diagnosis", "Sex", "percent_duplicate")],
                "meth" = rowSums(meth[,grepl("chr", colnames(meth), fixed = TRUE)]), 
                "cov" = rowSums(cov[,grepl("chr", colnames(cov), fixed = TRUE)]))
global$permeth <- global$meth * 100 / global$cov
rm(cov, meth)

# Global Stats
stats_global <- rbind(c("DxAll", summary(lm(permeth ~ Diagnosis + percent_duplicate, 
                                            data = global))$coefficients["DiagnosisASD",]),
                      c("DxAllSex", summary(lm(permeth ~ Diagnosis + Sex + percent_duplicate, 
                                               data = global))$coefficients["DiagnosisASD",]),
                      c("DxMales", summary(lm(permeth ~ Diagnosis + percent_duplicate, 
                                              data = subset(global, Sex == "M")))$coefficients["DiagnosisASD",]),
                      c("DxFemales", summary(lm(permeth ~ Diagnosis + percent_duplicate, 
                                                data = subset(global, Sex == "F")))$coefficients["DiagnosisASD",])) %>%
        as.data.frame(stringsAsFactors = FALSE)
colnames(stats_global) <- c("Comparison", "Estimate", "StdError", "tvalue", "pvalue")
stats_global[,c("Estimate", "StdError", "tvalue", "pvalue")] <- lapply(stats_global[,c("Estimate", "StdError", "tvalue", "pvalue")], as.numeric)

# Chromosomal Stats
stats_chrom <- DMRmethLm(DMRs = colnames(permeth)[grepl("chr", colnames(permeth), fixed = TRUE)],
                         catVars = c("Diagnosis", "Sex"), contVars = NULL, sampleData = permeth, 
                         file = "Tables/Chromosome Methylation by Diagnosis All Samples Discovery Stats.txt",
                         adj = "percent_duplicate")
stats_chrom_males <- DMRmethLm(DMRs = colnames(permeth)[grepl("chr", colnames(permeth), fixed = TRUE)],
                               catVars = c("Diagnosis"), contVars = NULL, sampleData = subset(permeth, Sex == "M"), 
                               file = "Tables/Chromosome Methylation by Diagnosis Males Discovery Stats.txt",
                               adj = "percent_duplicate")
stats_chrom_females <- DMRmethLm(DMRs = colnames(permeth)[grepl("chr", colnames(permeth), fixed = TRUE)],
                                 catVars = c("Diagnosis"), contVars = NULL, sampleData = subset(permeth, Sex == "F"), 
                                 file = "Tables/Chromosome Methylation by Diagnosis Females Discovery Stats.txt",
                                 adj = "percent_duplicate")

# Chromosome Plot
permeth_m <- melt(permeth, id.vars = c("Sample", "Diagnosis", "Sex", "percent_duplicate"))
colnames(permeth_m) <- c("Sample", "Diagnosis", "Sex", "percent_duplicate", "Chromosome", "Methylation")
permeth_m <- subset(permeth_m, !Chromosome %in% c("chrY", "chrM"))

gg <- ggplot(data = permeth_m)
gg +
        geom_boxplot(aes(x = Chromosome, y = Methylation, fill = Diagnosis), outlier.size = 1.1) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              axis.ticks = element_line(size = 1.25), legend.direction = "horizontal",
              legend.key = element_blank(), legend.position = c(0.92, 1.025), 
              legend.background = element_blank(), legend.title = element_blank(), 
              plot.margin = unit(c(2, 0.5, 0.5, 0.5), "lines"),  
              axis.text.x = element_text(size = 15, color = "black", angle = 90, hjust = 1, vjust = 0.5),
              axis.text.y = element_text(size = 15, color = "black"),
              axis.title = element_text(size = 16, color = "black"),
              strip.background = element_blank()) +
        scale_fill_manual(breaks = c("TD", "ASD"), values  = c("TD" = "#3366CC", "ASD" = "#FF3366")) +
        coord_cartesian(ylim = c(70, 82)) +
        scale_y_continuous(breaks = pretty_breaks(6)) +
        ylab("Methylation (%)") +
        facet_grid(rows = vars(Sex))
ggsave("Figures/Chromosome Methylation Males, Females Discovery Boxplot.png", dpi = 600, width = 10, height = 9, units = "in")

# Replication ####
# Data
meth <- read.csv("Tables/Chromosome_Meth_CpG_Counts_Replication.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(meth)[1] <- "Sample"
meth$Diagnosis[meth$Diagnosis == "Ctrl_TD"] <- "TD"
meth$Diagnosis[meth$Diagnosis == "Exp_ASD"] <- "ASD"
meth$Diagnosis <- factor(meth$Diagnosis, levels = c("TD", "ASD"))
meth$Sex <- factor(meth$Sex, levels = c("M", "F"))

cov <- read.csv("Tables/Chromosome_Total_CpG_Counts_Replication.csv", header = TRUE, stringsAsFactors = FALSE)
permeth_rep <- cbind(meth[,c("Sample", "Diagnosis", "Sex", "percent_duplicate")],
                     meth[,grepl("chr", colnames(meth), fixed = TRUE)] * 100 / cov[,grepl("chr", colnames(cov), fixed = TRUE)])

global_rep <- cbind(meth[,c("Sample", "Diagnosis", "Sex", "percent_duplicate")],
                    "meth" = rowSums(meth[,grepl("chr", colnames(meth), fixed = TRUE)]), 
                    "cov" = rowSums(cov[,grepl("chr", colnames(cov), fixed = TRUE)]))
global_rep$permeth <- global_rep$meth * 100 / global_rep$cov
rm(cov, meth)

# Global Stats
stats_global_rep <- rbind(c("DxAll", summary(lm(permeth ~ Diagnosis + percent_duplicate, 
                                            data = global_rep))$coefficients["DiagnosisASD",]),
                      c("DxAllSex", summary(lm(permeth ~ Diagnosis + Sex + percent_duplicate, 
                                               data = global_rep))$coefficients["DiagnosisASD",]),
                      c("DxMales", summary(lm(permeth ~ Diagnosis + percent_duplicate, 
                                              data = subset(global_rep, Sex == "M")))$coefficients["DiagnosisASD",]),
                      c("DxFemales", summary(lm(permeth ~ Diagnosis + percent_duplicate, 
                                                data = subset(global_rep, Sex == "F")))$coefficients["DiagnosisASD",])) %>%
        as.data.frame(stringsAsFactors = FALSE)
colnames(stats_global_rep) <- c("Comparison", "Estimate", "StdError", "tvalue", "pvalue")
stats_global_rep[,c("Estimate", "StdError", "tvalue", "pvalue")] <- lapply(stats_global_rep[,c("Estimate", "StdError", "tvalue", "pvalue")], as.numeric)

# Chromosomal Stats
stats_chrom_rep <- DMRmethLm(DMRs = colnames(permeth_rep)[grepl("chr", colnames(permeth_rep), fixed = TRUE)],
                         catVars = c("Diagnosis", "Sex"), contVars = NULL, sampleData = permeth_rep, 
                         file = "Tables/Chromosome Methylation by Diagnosis All Samples Replication Stats.txt",
                         adj = "percent_duplicate")
stats_chrom_males_rep <- DMRmethLm(DMRs = colnames(permeth_rep)[grepl("chr", colnames(permeth_rep), fixed = TRUE)],
                               catVars = c("Diagnosis"), contVars = NULL, sampleData = subset(permeth_rep, Sex == "M"), 
                               file = "Tables/Chromosome Methylation by Diagnosis Males Replication Stats.txt",
                               adj = "percent_duplicate")
stats_chrom_females_rep <- DMRmethLm(DMRs = colnames(permeth_rep)[grepl("chr", colnames(permeth_rep), fixed = TRUE)],
                                 catVars = c("Diagnosis"), contVars = NULL, sampleData = subset(permeth_rep, Sex == "F"), 
                                 file = "Tables/Chromosome Methylation by Diagnosis Females Replication Stats.txt",
                                 adj = "percent_duplicate")

# Chromosome Plot
permeth_rep_m <- melt(permeth_rep, id.vars = c("Sample", "Diagnosis", "Sex", "percent_duplicate"))
colnames(permeth_rep_m) <- c("Sample", "Diagnosis", "Sex", "percent_duplicate", "Chromosome", "Methylation")
permeth_rep_m <- subset(permeth_rep_m, !Chromosome %in% c("chrY", "chrM"))

gg <- ggplot(data = permeth_rep_m)
gg +
        geom_boxplot(aes(x = Chromosome, y = Methylation, fill = Diagnosis), outlier.size = 1.1) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              axis.ticks = element_line(size = 1.25), legend.direction = "horizontal",
              legend.key = element_blank(), legend.position = c(0.92, 1.025), 
              legend.background = element_blank(), legend.title = element_blank(), 
              plot.margin = unit(c(2, 0.5, 0.5, 0.5), "lines"),  
              axis.text.x = element_text(size = 15, color = "black", angle = 90, hjust = 1, vjust = 0.5),
              axis.text.y = element_text(size = 15, color = "black"),
              axis.title = element_text(size = 16, color = "black"),
              strip.background = element_blank()) +
        scale_fill_manual(breaks = c("TD", "ASD"), values  = c("TD" = "#3366CC", "ASD" = "#FF3366")) +
        coord_cartesian(ylim = c(68, 80)) +
        scale_y_continuous(breaks = pretty_breaks(6)) +
        ylab("Methylation (%)") +
        facet_grid(rows = vars(Sex))
ggsave("Figures/Chromosome Methylation Males, Females Replication Boxplot.png", dpi = 600, width = 10, height = 9, units = "in")

# Discovery vs Replication Stats and Plots ####
# Global
global_pooled <- rbind(global, global_rep)
global_pooled$Set <- factor(c(rep("Discovery", nrow(global)), rep("Replication", nrow(global_rep))), 
                            levels = c("Discovery", "Replication"))

stats_global_pooled <- rbind(c("DxAll", summary(lm(permeth ~ Diagnosis + percent_duplicate + Set, 
                                                data = global_pooled))$coefficients["DiagnosisASD",]),
                          c("DxAllSex", summary(lm(permeth ~ Diagnosis + Sex + percent_duplicate + Set, 
                                                   data = global_pooled))$coefficients["DiagnosisASD",]),
                          c("DxMales", summary(lm(permeth ~ Diagnosis + percent_duplicate + Set, 
                                                  data = subset(global_pooled, Sex == "M")))$coefficients["DiagnosisASD",]),
                          c("DxFemales", summary(lm(permeth ~ Diagnosis + percent_duplicate + Set, 
                                                    data = subset(global_pooled, Sex == "F")))$coefficients["DiagnosisASD",])) %>%
        as.data.frame(stringsAsFactors = FALSE)
colnames(stats_global_pooled) <- c("Comparison", "Estimate", "StdError", "tvalue", "pvalue")
stats_global_pooled[,c("Estimate", "StdError", "tvalue", "pvalue")] <- lapply(stats_global_pooled[,c("Estimate", "StdError", "tvalue", "pvalue")], as.numeric)

gg <- ggplot(data = global_pooled)
gg +
        geom_boxplot(aes(x = Sex, y = permeth, fill = Diagnosis), outlier.size = 1.1) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              axis.ticks = element_line(size = 1.25), legend.direction = "horizontal",
              legend.key = element_blank(), legend.position = c(0.83, 0.95), 
              legend.background = element_blank(), legend.title = element_blank(), 
              legend.text = element_text(size = 15), strip.background = element_blank(),
              plot.margin = unit(c(0, 0.5, 0.5, 0.5), "lines"), strip.text = element_text(size = 17), 
              axis.text.x = element_text(size = 15, color = "black"),
              axis.text.y = element_text(size = 15, color = "black"),
              axis.title = element_text(size = 16, color = "black")) +
        scale_fill_manual(breaks = c("TD", "ASD"), values  = c("TD" = "#3366CC", "ASD" = "#FF3366")) +
        coord_cartesian(ylim = c(70, 80)) +
        scale_y_continuous(breaks = pretty_breaks(6)) +
        ylab("Methylation (%)") +
        facet_grid(cols = vars(Set))
ggsave("Figures/Global Methylation from BSseq Discovery and Replication Boxplot.png", dpi = 600, width = 5, height = 5, units = "in")

# Discovery chrX Stats
chrX_permeth <- subset(permeth_m, Chromosome == "chrX")
stats_chrX <- rbind(c("DxAll", summary(lm(Methylation ~ Diagnosis + percent_duplicate, 
                                                 data = chrX_permeth))$coefficients["DiagnosisASD",]),
                           c("DxAllSex", summary(lm(Methylation ~ Diagnosis + Sex + percent_duplicate, 
                                                    data = chrX_permeth))$coefficients["DiagnosisASD",]),
                           c("DxMales", summary(lm(Methylation ~ Diagnosis + percent_duplicate, 
                                                   data = subset(chrX_permeth, Sex == "M")))$coefficients["DiagnosisASD",]),
                           c("DxFemales", summary(lm(Methylation ~ Diagnosis + percent_duplicate, 
                                                     data = subset(chrX_permeth, Sex == "F")))$coefficients["DiagnosisASD",])) %>%
        as.data.frame(stringsAsFactors = FALSE)
colnames(stats_chrX) <- c("Comparison", "Estimate", "StdError", "tvalue", "pvalue")
stats_chrX[,c("Estimate", "StdError", "tvalue", "pvalue")] <- lapply(stats_chrX[,c("Estimate", "StdError", "tvalue", "pvalue")], as.numeric)

# Replication chrX Stats
chrX_permeth_rep <- subset(permeth_rep_m, Chromosome == "chrX")
stats_chrX_rep <- rbind(c("DxAll", summary(lm(Methylation ~ Diagnosis + percent_duplicate, 
                                          data = chrX_permeth_rep))$coefficients["DiagnosisASD",]),
                    c("DxAllSex", summary(lm(Methylation ~ Diagnosis + Sex + percent_duplicate, 
                                             data = chrX_permeth_rep))$coefficients["DiagnosisASD",]),
                    c("DxMales", summary(lm(Methylation ~ Diagnosis + percent_duplicate, 
                                            data = subset(chrX_permeth_rep, Sex == "M")))$coefficients["DiagnosisASD",]),
                    c("DxFemales", summary(lm(Methylation ~ Diagnosis + percent_duplicate, 
                                              data = subset(chrX_permeth_rep, Sex == "F")))$coefficients["DiagnosisASD",])) %>%
        as.data.frame(stringsAsFactors = FALSE)
colnames(stats_chrX_rep) <- c("Comparison", "Estimate", "StdError", "tvalue", "pvalue")
stats_chrX_rep[,c("Estimate", "StdError", "tvalue", "pvalue")] <- lapply(stats_chrX_rep[,c("Estimate", "StdError", "tvalue", "pvalue")], as.numeric)

# Pooled chrX Stats
chrX_pooled <- rbind(subset(permeth_m, Chromosome == "chrX"), subset(permeth_rep_m, Chromosome == "chrX"))
chrX_pooled$Set <- factor(c(rep("Discovery", nrow(permeth)), rep("Replication", nrow(permeth_rep))), 
                            levels = c("Discovery", "Replication"))

stats_chrX_pooled <- rbind(c("DxAll", summary(lm(Methylation ~ Diagnosis + percent_duplicate + Set, 
                                                   data = chrX_pooled))$coefficients["DiagnosisASD",]),
                             c("DxAllSex", summary(lm(Methylation ~ Diagnosis + Sex + percent_duplicate + Set, 
                                                      data = chrX_pooled))$coefficients["DiagnosisASD",]),
                             c("DxMales", summary(lm(Methylation ~ Diagnosis + percent_duplicate + Set, 
                                                     data = subset(chrX_pooled, Sex == "M")))$coefficients["DiagnosisASD",]),
                             c("DxFemales", summary(lm(Methylation ~ Diagnosis + percent_duplicate + Set, 
                                                       data = subset(chrX_pooled, Sex == "F")))$coefficients["DiagnosisASD",])) %>%
        as.data.frame(stringsAsFactors = FALSE)
colnames(stats_chrX_pooled) <- c("Comparison", "Estimate", "StdError", "tvalue", "pvalue")
stats_chrX_pooled[,c("Estimate", "StdError", "tvalue", "pvalue")] <- lapply(stats_chrX_pooled[,c("Estimate", "StdError", "tvalue", "pvalue")], as.numeric)

# chrX Plot
gg <- ggplot(data = chrX_pooled)
gg +
        geom_boxplot(aes(x = Sex, y = Methylation, fill = Diagnosis), outlier.size = 1.1) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              axis.ticks = element_line(size = 1.25), legend.direction = "horizontal",
              legend.key = element_blank(), legend.position = c(0.83, 0.95), 
              legend.background = element_blank(), legend.title = element_blank(), 
              legend.text = element_text(size = 15), strip.background = element_blank(),
              plot.margin = unit(c(0, 0.5, 0.5, 0.5), "lines"), strip.text = element_text(size = 17), 
              axis.text.x = element_text(size = 15, color = "black"),
              axis.text.y = element_text(size = 15, color = "black"),
              axis.title = element_text(size = 16, color = "black")) +
        scale_fill_manual(breaks = c("TD", "ASD"), values  = c("TD" = "#3366CC", "ASD" = "#FF3366")) +
        coord_cartesian(ylim = c(68, 82)) +
        scale_y_continuous(breaks = pretty_breaks(6)) +
        ylab("Methylation (%)") +
        facet_grid(cols = vars(Set))
ggsave("Figures/chrX Methylation from BSseq Discovery and Replication Boxplot.png", dpi = 600, width = 5, height = 5, units = "in")



