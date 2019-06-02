# Chromosome Methylation Analysis ---------------------------------------------------
# ASD Cord Blood Methylation
# Charles Mordaunt
# 5/31/19
# Excluded JLCM032B and JLCM050B

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
permeth <- subset(permeth, !Sample %in% c("JLCM032B", "JLCM050B")) # Exclude mislabeled females
cov_disc <- cov
rm(meth, cov)

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
cov_rep <- cov
rm(meth, cov)

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

# All Chromosomes Pooled Stats
permeth_pooled <- rbind(permeth, permeth_rep)
permeth_pooled$Platform <- c(rep("HiSeqX10", nrow(permeth)), rep("HiSeq4000", nrow(permeth_rep))) %>% 
        factor(levels = c("HiSeqX10", "HiSeq4000"))

stats_chrom_pooled <- sapply(permeth_pooled[,grepl("chr", colnames(permeth_pooled), fixed = TRUE)], function(x){
        summary(lm(x ~ permeth_pooled$Diagnosis + permeth_pooled$percent_duplicate + permeth_pooled$Platform))$coefficients["permeth_pooled$DiagnosisASD",]
}) %>% t %>% as.data.frame
colnames(stats_chrom_pooled) <- c("Estimate", "StdError", "tvalue", "pvalue")
stats_chrom_pooled$Region <- rownames(stats_chrom_pooled)
stats_chrom_pooled$qvalue <- p.adjust(stats_chrom_pooled$pvalue, method = "fdr")
stats_chrom_pooled <- stats_chrom_pooled[, c("Region", "Estimate", "StdError", "tvalue", "pvalue", "qvalue")]
write.table(stats_chrom_pooled, "Tables/Chromosome Methylation by Diagnosis Pooled Stats.txt", sep = "\t", quote = FALSE, 
            row.names = FALSE)

permeth_pooled_males <- subset(permeth_pooled, Sex == "M")
stats_chrom_pooled_males <- sapply(permeth_pooled_males[,grepl("chr", colnames(permeth_pooled_males), fixed = TRUE)], function(x){
        summary(lm(x ~ permeth_pooled_males$Diagnosis + permeth_pooled_males$percent_duplicate + permeth_pooled_males$Platform))$coefficients["permeth_pooled_males$DiagnosisASD",]
}) %>% t %>% as.data.frame
colnames(stats_chrom_pooled_males) <- c("Estimate", "StdError", "tvalue", "pvalue")
stats_chrom_pooled_males$Region <- rownames(stats_chrom_pooled_males)
stats_chrom_pooled_males$qvalue <- p.adjust(stats_chrom_pooled_males$pvalue, method = "fdr")
stats_chrom_pooled_males <- stats_chrom_pooled_males[, c("Region", "Estimate", "StdError", "tvalue", "pvalue", "qvalue")]
write.table(stats_chrom_pooled_males, "Tables/Chromosome Methylation by Diagnosis Pooled Stats Males.txt", sep = "\t", quote = FALSE, 
            row.names = FALSE)

permeth_pooled_females <- subset(permeth_pooled, Sex == "F")
stats_chrom_pooled_females <- sapply(permeth_pooled_females[,grepl("chr", colnames(permeth_pooled_females), fixed = TRUE)], function(x){
        summary(lm(x ~ permeth_pooled_females$Diagnosis + permeth_pooled_females$percent_duplicate + permeth_pooled_females$Platform))$coefficients["permeth_pooled_females$DiagnosisASD",]
}) %>% t %>% as.data.frame
colnames(stats_chrom_pooled_females) <- c("Estimate", "StdError", "tvalue", "pvalue")
stats_chrom_pooled_females$Region <- rownames(stats_chrom_pooled_females)
stats_chrom_pooled_females$qvalue <- p.adjust(stats_chrom_pooled_females$pvalue, method = "fdr")
stats_chrom_pooled_females <- stats_chrom_pooled_females[, c("Region", "Estimate", "StdError", "tvalue", "pvalue", "qvalue")]
write.table(stats_chrom_pooled_females, "Tables/Chromosome Methylation by Diagnosis Pooled Stats Females.txt", sep = "\t", quote = FALSE, 
            row.names = FALSE)

# Chromosome Methylation ~ Diagnosis Heatmap
heatStats <- rbind(stats_chrom_males[,c("Region", "Estimate", "pvalue")], stats_chrom_males_rep[,c("Region", "Estimate", "pvalue")],
                   stats_chrom_pooled_males[,c("Region", "Estimate", "pvalue")], stats_chrom_females[,c("Region", "Estimate", "pvalue")],
                   stats_chrom_females_rep[,c("Region", "Estimate", "pvalue")], stats_chrom_pooled_females[,c("Region", "Estimate", "pvalue")])
heatStats$Set <- rep(c("Discovery", "Replication", "Pooled", "Discovery", "Replication", "Pooled"), each = 25) %>%
        factor(levels = c("Discovery", "Replication", "Pooled"))
heatStats$Sex <- rep(c("Males", "Females"), each = 75) %>% factor(levels = c("Males", "Females"))
heatStats <- subset(heatStats, !Region %in% c("chrY", "chrM"))
heatStats$Significant <- (heatStats$pvalue < 0.05) %>% factor(levels = c("TRUE", "FALSE"))
heatStats$Region <- factor(heatStats$Region, levels = rev(unique(heatStats$Region)))

g <- ggplot(data = heatStats)
g + 
        geom_tile(aes(x = Set, y = Region, fill = Estimate, color = Estimate)) + 
        geom_text(aes(x = Set, y = Region, alpha = Significant, label = "*"), color = "white", size = 15, nudge_y = -0.35) +
        facet_grid(cols = vars(Sex)) +
        scale_fill_gradientn("Estimate", colors = c("#3366CC", "Black", "#FF0000"), values = c(0, 1), na.value = "#FF0000", 
                             limits = c(-1, 1), breaks = pretty_breaks(n = 2)) +
        scale_color_gradientn("Estimate", colors = c("#3366CC", "Black", "#FF0000"), values = c(0, 1), na.value = "#FF0000", 
                              limits = c(-1, 1), breaks = pretty_breaks(n = 2)) +
        scale_alpha_manual(breaks = c("TRUE", "FALSE"), values = c(1, 0), guide = FALSE) +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'vertical', legend.position = c(1.11, 0.902), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "Black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_text(size = 24),
              axis.text.x = element_text(color = "Black", angle = 30, hjust = 1), 
              axis.text.y = element_text(color = "Black", size = 20),
              legend.background = element_blank(), 
              strip.text = element_text(size = 26), plot.margin = unit(c(1, 8, 1, 1), "lines"), 
              axis.title = element_blank(), strip.background = element_blank())
ggsave("Figures/Chromosome Methylation by Diagnosis Heatmap.png", dpi = 600, width = 10, height = 10.25, units = "in")

# All Chromosomes, Males Boxplot
all_chroms_males <- rbind(subset(permeth_m, Sex == "M"), subset(permeth_rep_m, Sex == "M"))
all_chroms_males$Set <- c(rep("Discovery", nrow(subset(permeth_m, Sex == "M"))),
                          rep("Replication", nrow(subset(permeth_rep_m, Sex == "M")))) %>% 
        factor(levels = c("Discovery", "Replication"))
gg <- ggplot(data = all_chroms_males)
gg +
        geom_boxplot(aes(x = Chromosome, y = Methylation, fill = Diagnosis), outlier.size = 1.1) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              axis.ticks = element_line(size = 1.25), legend.direction = "horizontal",
              legend.key = element_blank(), legend.position = c(0.92, 1.025), 
              legend.background = element_blank(), legend.title = element_blank(), 
              plot.margin = unit(c(2, 0.5, 0.5, 0.5), "lines"), strip.text = element_text(size = 20), 
              axis.text.x = element_text(size = 16, color = "black", angle = 90, hjust = 1, vjust = 0.5),
              axis.text.y = element_text(size = 16, color = "black"),
              axis.title = element_text(size = 20, color = "black"),
              strip.background = element_blank()) +
        scale_fill_manual(breaks = c("TD", "ASD"), values  = c("TD" = "#3366CC", "ASD" = "#FF3366")) +
        scale_y_continuous(breaks = pretty_breaks(6)) +
        ylab("Methylation (%)") +
        facet_grid(rows = vars(Set), scales = "free")
ggsave("Figures/Chromosome Methylation Males Discovery and Replication Boxplot.png", dpi = 600, width = 10, height = 9, units = "in")

# All Chromosomes, Males Mean CL
gg <- ggplot(data = all_chroms_males)
gg +
        stat_summary(aes(x = Chromosome, y = Methylation, color = Diagnosis), fun.data = "mean_cl_boot",
                     geom = "crossbar", size = 0.8, position = position_dodge2(width = 0.75, padding = 0.2)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              axis.ticks = element_line(size = 1.25), legend.direction = "horizontal",
              legend.key = element_blank(), legend.position = c(0.92, 1.025), legend.text = element_text(size = 22),
              legend.background = element_blank(), legend.title = element_blank(), 
              plot.margin = unit(c(2, 0.5, 0.5, 0.5), "lines"), strip.text = element_text(size = 24), 
              axis.text.x = element_text(size = 20, color = "black", angle = 90, hjust = 1, vjust = 0.5),
              axis.text.y = element_text(size = 20, color = "black"),
              axis.title = element_text(size = 24, color = "black"),
              strip.background = element_blank()) +
        scale_color_manual(breaks = c("TD", "ASD"), values  = c("TD" = "#3366CC", "ASD" = "#FF3366")) +
        scale_y_continuous(breaks = pretty_breaks(6)) +
        ylab("Methylation (%)") +
        facet_grid(rows = vars(Set), scales = "free")
ggsave("Figures/Chromosome Methylation Males Discovery and Replication Mean CL.png", dpi = 600, width = 10, height = 9, units = "in")

# All Chromosomes, Males Line Plot
gg <- ggplot(data = all_chroms_males)
gg +
        geom_line(aes(x = Chromosome, y = Methylation, color = Diagnosis, group = Sample), size = 0.85, alpha = 0.6) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              axis.ticks = element_line(size = 1.25), legend.direction = "horizontal",
              legend.key = element_blank(), legend.position = c(0.92, 1.025), legend.text = element_text(size = 22),
              legend.background = element_blank(), legend.title = element_blank(), 
              plot.margin = unit(c(2, 0.5, 0.5, 0.5), "lines"), strip.text = element_text(size = 24), 
              axis.text.x = element_text(size = 20, color = "black", angle = 90, hjust = 1, vjust = 0.5),
              axis.text.y = element_text(size = 20, color = "black"),
              axis.title = element_text(size = 24, color = "black"),
              strip.background = element_blank()) +
        scale_color_manual(breaks = c("TD", "ASD"), values  = c("TD" = "#3366CC", "ASD" = "#FF3366")) +
        scale_y_continuous(breaks = pretty_breaks(6)) +
        ylab("Methylation (%)") +
        facet_grid(rows = vars(Set), scales = "free")
ggsave("Figures/Chromosome Methylation Males Discovery and Replication Line Plot.png", dpi = 600, width = 10, height = 9, units = "in")

# All Chromosomes, Females
all_chroms_females <- rbind(subset(permeth_m, Sex == "F"), subset(permeth_rep_m, Sex == "F"))
all_chroms_females$Set <- c(rep("Discovery", nrow(subset(permeth_m, Sex == "F"))),
                          rep("Replication", nrow(subset(permeth_rep_m, Sex == "F")))) %>% 
        factor(levels = c("Discovery", "Replication"))
gg <- ggplot(data = all_chroms_females)
gg +
        geom_boxplot(aes(x = Chromosome, y = Methylation, fill = Diagnosis), outlier.size = 1.1) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              axis.ticks = element_line(size = 1.25), legend.direction = "horizontal",
              legend.key = element_blank(), legend.position = c(0.92, 1.025), 
              legend.background = element_blank(), legend.title = element_blank(), 
              plot.margin = unit(c(2, 0.5, 0.5, 0.5), "lines"), strip.text = element_text(size = 20), 
              axis.text.x = element_text(size = 16, color = "black", angle = 90, hjust = 1, vjust = 0.5),
              axis.text.y = element_text(size = 16, color = "black"),
              axis.title = element_text(size = 20, color = "black"),
              strip.background = element_blank()) +
        scale_fill_manual(breaks = c("TD", "ASD"), values  = c("TD" = "#3366CC", "ASD" = "#FF3366")) +
        scale_y_continuous(breaks = pretty_breaks(6)) +
        ylab("Methylation (%)") +
        facet_grid(rows = vars(Set), scales = "free")
ggsave("Figures/Chromosome Methylation Females Discovery and Replication Boxplot.png", dpi = 600, width = 10, height = 9, units = "in")

# All Chromosomes, Females Mean CL
gg <- ggplot(data = all_chroms_females)
gg +
        stat_summary(aes(x = Chromosome, y = Methylation, color = Diagnosis), fun.data = "mean_cl_boot",
                     geom = "crossbar", size = 0.8, position = position_dodge2(width = 0.75, padding = 0.2)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              axis.ticks = element_line(size = 1.25), legend.direction = "horizontal",
              legend.key = element_blank(), legend.position = c(0.92, 1.025), 
              legend.background = element_blank(), legend.title = element_blank(), 
              plot.margin = unit(c(2, 0.5, 0.5, 0.5), "lines"), strip.text = element_text(size = 20), 
              axis.text.x = element_text(size = 16, color = "black", angle = 90, hjust = 1, vjust = 0.5),
              axis.text.y = element_text(size = 16, color = "black"),
              axis.title = element_text(size = 20, color = "black"),
              strip.background = element_blank()) +
        scale_color_manual(breaks = c("TD", "ASD"), values  = c("TD" = "#3366CC", "ASD" = "#FF3366")) +
        scale_y_continuous(breaks = pretty_breaks(6)) +
        ylab("Methylation (%)") +
        facet_grid(rows = vars(Set), scales = "free")
ggsave("Figures/Chromosome Methylation Females Discovery and Replication Mean CL.png", dpi = 600, width = 10, height = 9, units = "in")

# All Chromosomes, Females Line Plot
gg <- ggplot(data = all_chroms_females)
gg +
        geom_line(aes(x = Chromosome, y = Methylation, color = Diagnosis, group = Sample), size = 0.85, alpha = 0.6) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              axis.ticks = element_line(size = 1.25), legend.direction = "horizontal",
              legend.key = element_blank(), legend.position = c(0.92, 1.025), legend.text = element_text(size = 22),
              legend.background = element_blank(), legend.title = element_blank(), 
              plot.margin = unit(c(2, 0.5, 0.5, 0.5), "lines"), strip.text = element_text(size = 24), 
              axis.text.x = element_text(size = 20, color = "black", angle = 90, hjust = 1, vjust = 0.5),
              axis.text.y = element_text(size = 20, color = "black"),
              axis.title = element_text(size = 24, color = "black"),
              strip.background = element_blank()) +
        scale_color_manual(breaks = c("TD", "ASD"), values  = c("TD" = "#3366CC", "ASD" = "#FF3366")) +
        scale_y_continuous(breaks = pretty_breaks(6)) +
        ylab("Methylation (%)") +
        facet_grid(rows = vars(Set), scales = "free")
ggsave("Figures/Chromosome Methylation Females Discovery and Replication Line Plot.png", dpi = 600, width = 10, height = 9, units = "in")

# ChrX ####
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

# chrX Dot Plot
chrX_pooled$Sex <- ifelse(chrX_pooled$Sex == "M", yes = "Males", no = "Females") %>% factor(levels = c("Males", "Females"))

g <- ggplot(data = chrX_pooled)
g + 
        stat_summary(aes(x = Sex, y = Methylation, group = Diagnosis), fun.data = "mean_cl_boot", 
                     geom = "crossbar", color = "black", size = 0.7, position = "dodge2") +
        geom_dotplot(aes(x = Sex, y = Methylation, fill = Diagnosis, color = Diagnosis), binwidth = 0.21, 
                     binaxis = "y", stackdir = "center", position = "dodge", stackratio = 1.15, dotsize = 0.82) + 
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.895, 0.96), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              strip.text = element_text(size = 22), plot.margin = unit(c(0.3, 1, 1, 1), "lines"), 
              axis.title.x = element_blank(), strip.background = element_blank()) +
        scale_y_continuous(breaks = pretty_breaks(n = 5)) +
        ylab("ChrX Methylation (%)") +
        facet_wrap(vars(Set)) +
        scale_fill_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366")) +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366"))
ggsave("Figures/ChrX Methylation from BSseq Discovery and Replication Dotplot with Mean CL.png", dpi = 600, width = 9, height = 7, units = "in")

# Chromosome X and Y Coverage ####
cov_XY <- rbind(cov_disc[,c("sampleNames.bs.", "Diagnosis", "Sex", "chrX", "chrY")],
                cov_rep[,c("sampleNames.bs.", "Diagnosis", "Sex", "chrX", "chrY")])
colnames(cov_XY)[colnames(cov_XY) == "sampleNames.bs."] <- "Sample"
cov_XY$Set <- c(rep("Discovery", nrow(cov_disc)), rep("Replication", nrow(cov_rep))) %>% 
        factor(levels = c("Discovery", "Replication"))
cov_XY$Diagnosis <- ifelse(cov_XY$Diagnosis == "Ctrl_TD", yes = "TD", no = "ASD") %>% factor(levels = c("TD", "ASD"))
cov_XY$Sex <- ifelse(cov_XY$Sex == "M", yes = "Males", no = "Females") %>% factor(levels = c("Males", "Females"))
cov_XY <- cov_XY[,c("Sample", "Diagnosis", "Sex", "Set", "chrX", "chrY")]
cov_XY$Total <- c(rowSums(cov_disc[,grepl("chr", colnames(cov_disc), fixed = TRUE)]),
                  rowSums(cov_rep[,grepl("chr", colnames(cov_rep), fixed = TRUE)]))
cov_XY$per_chrX <- cov_XY$chrX * 100 / cov_XY$Total
cov_XY$per_chrY <- cov_XY$chrY * 100 / cov_XY$Total

# ChrX Coverage Dot Plot
g <- ggplot(data = cov_XY)
g + 
        stat_summary(aes(x = Sex, y = per_chrX, group = Diagnosis), fun.data = "mean_cl_boot", 
                     geom = "crossbar", color = "black", size = 0.7, position = "dodge2") +
        geom_dotplot(aes(x = Sex, y = per_chrX, fill = Diagnosis, color = Diagnosis), binwidth = 0.024,
                     binaxis = "y", stackdir = "center", position = "dodge", stackratio = 1.15, dotsize = 0.97) + 
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.895, 0.96), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              strip.text = element_text(size = 22), plot.margin = unit(c(0.3, 1, 1, 1), "lines"), 
              axis.title.x = element_blank(), strip.background = element_blank()) +
        scale_y_continuous(breaks = pretty_breaks(n = 5)) +
        ylab("ChrX Coverage (% of Total)") +
        facet_wrap(vars(Set)) +
        scale_fill_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366")) +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366"))
ggsave("Figures/ChrX Coverage from BSseq Discovery and Replication Dotplot with Mean CL.png", dpi = 600, width = 10, 
       height = 7, units = "in")

# Chromosome Y Coverage Dotplot
g <- ggplot(data = cov_XY)
g + 
        stat_summary(aes(x = Sex, y = per_chrY, group = Diagnosis), fun.data = "mean_cl_boot", 
                     geom = "crossbar", color = "black", size = 0.7, position = "dodge2") +
        geom_dotplot(aes(x = Sex, y = per_chrY, fill = Diagnosis, color = Diagnosis), binwidth = 0.003,
                     binaxis = "y", stackdir = "center", position = "dodge", stackratio = 1.15, dotsize = 0.97) + 
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.895, 0.96), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              strip.text = element_text(size = 22), plot.margin = unit(c(0.3, 1, 1, 1), "lines"), 
              axis.title.x = element_blank(), strip.background = element_blank()) +
        scale_y_continuous(breaks = pretty_breaks(n = 5)) +
        ylab("ChrY Coverage (% of Total)") +
        facet_wrap(vars(Set)) +
        scale_fill_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366")) +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366"))
ggsave("Figures/ChrY Coverage from BSseq Discovery and Replication Dotplot with Mean CL.png", dpi = 600, width = 10, 
       height = 7, units = "in")

# Two Discovery ASD male samples are actually females: JLCM050B and JLCM032B
library(sas7bdat)
samples <- read.delim(file = "Merged Database/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo.txt", sep = "\t", header = TRUE, 
                      stringsAsFactors = FALSE)
samples <- subset(samples, !(Cord_Blood_IBC == 101470 & Platform %in% c("HiSeq2500", "HiSeq4000"))) # Remove 101470 Duplicates
covars <- read.sas7bdat("Merged Database/e_m_covars_v2.sas7bdat")
samples$COI_ID[samples$Sequencing_ID == "JLCM032B"] # 9031905
samples$COI_ID[samples$Sequencing_ID == "JLCM050B"] # 534104

covars$coi_gender[covars$coi_id == 9031905] #M
covars$coi_gender[covars$coi_id == 534104] #M
# Both of these subjects are labeled as males in most recent covariate database
# Exclude these subjects from analyses

# ChrX and ChrY Coverage without Mislabeled Females
# ChrX Coverage Dot Plot
cov_XY <- subset(cov_XY, !Sample %in% c("JLCM032B", "JLCM050B"))
g <- ggplot(data = cov_XY)
g + 
        stat_summary(aes(x = Sex, y = per_chrX, group = Diagnosis), fun.data = "mean_cl_boot", 
                     geom = "crossbar", color = "black", size = 0.7, position = "dodge2") +
        geom_dotplot(aes(x = Sex, y = per_chrX, fill = Diagnosis, color = Diagnosis), binwidth = 0.024,
                     binaxis = "y", stackdir = "center", position = "dodge", stackratio = 1.15, dotsize = 0.97) + 
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.895, 0.96), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              strip.text = element_text(size = 22), plot.margin = unit(c(0.3, 1, 1, 1), "lines"), 
              axis.title.x = element_blank(), strip.background = element_blank()) +
        scale_y_continuous(breaks = pretty_breaks(n = 5)) +
        ylab("ChrX Coverage (% of Total)") +
        facet_wrap(vars(Set)) +
        scale_fill_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366")) +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366"))
ggsave("Figures/ChrX Coverage from BSseq Discovery and Replication Dotplot with Mean CL no Mislabeled Females.png", dpi = 600, width = 10, 
       height = 7, units = "in")

# Chromosome Y Coverage Dotplot
g <- ggplot(data = cov_XY)
g + 
        stat_summary(aes(x = Sex, y = per_chrY, group = Diagnosis), fun.data = "mean_cl_boot", 
                     geom = "crossbar", color = "black", size = 0.7, position = "dodge2") +
        geom_dotplot(aes(x = Sex, y = per_chrY, fill = Diagnosis, color = Diagnosis), binwidth = 0.003,
                     binaxis = "y", stackdir = "center", position = "dodge", stackratio = 1.15, dotsize = 0.97) + 
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.895, 0.96), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              strip.text = element_text(size = 22), plot.margin = unit(c(0.3, 1, 1, 1), "lines"), 
              axis.title.x = element_blank(), strip.background = element_blank()) +
        scale_y_continuous(breaks = pretty_breaks(n = 5)) +
        ylab("ChrY Coverage (% of Total)") +
        facet_wrap(vars(Set)) +
        scale_fill_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366")) +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366"))
ggsave("Figures/ChrY Coverage from BSseq Discovery and Replication Dotplot with Mean CL no Mislabeled Females.png", dpi = 600, width = 10, 
       height = 7, units = "in")





