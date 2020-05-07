# DMR Power Analysis with ssize.fdr ---------------------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 4/30/20

# Setup -------------------------------------------------------------------
# Load Packages ####
setwd("/share/lasallelab/Charles/CM_WGBS_ASD_CordBlood/Bismark_Reports")
.libPaths("/share/lasallelab/Charles/R")
sapply(c("tidyverse", "effectsize", "ssize.fdr", "qvalue", "bsseq"), require, character.only = TRUE)

# Load Functions ####
loadRegions <- function(file, chroms = c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM"), sort = TRUE, 
                        DMRid = FALSE, DMBid = FALSE){
        if(grepl("txt", file, fixed = TRUE)){
                regions <- read.delim(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
        }
        else{
                regions <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
        }
        if(nrow(regions) == 0){
                stop(c(paste("There are no regions in", file, sep = " ")))
        }
        if("seqnames" %in% colnames(regions)){
                colnames(regions)[colnames(regions) == "seqnames"] <- "chr"
        }
        regions <- subset(regions, chr %in% chroms)
        regions$chr <- factor(regions$chr, levels = chroms)
        if(sort){
                regions <- regions[order(regions$chr, regions$start),]
        }
        if(DMRid){
                regions$DMRid <- paste("DMR", 1:nrow(regions), sep = "_")
        }
        if(DMBid){
                regions$DMBid <- paste("DMB", 1:nrow(regions), sep = "_")
        }
        return(regions)
}

makeGRange <- function(DMRs, direction = c("all", "hyper", "hypo")){
        direction <- match.arg(direction)
        if(direction == "hyper"){
                DMRs <- subset(DMRs, percentDifference > 0)
        }
        if(direction == "hypo"){
                DMRs <- subset(DMRs, percentDifference < 0)
        }
        if("chr" %in% colnames(DMRs)){
                GR <- GRanges(seqnames = DMRs$chr, ranges = IRanges(start = DMRs$start, end = DMRs$end))
        } else {
                GR <- GRanges(seqnames = DMRs$seqnames, ranges = IRanges(start = DMRs$start, end = DMRs$end))
        }
        return(GR)
}

# Power Analysis for Discovery Males --------------------------------------
# Estimate Parameters ####
# Do t-tests on background regions (cluster)
bs <- readRDS("Discovery/Dx_Males/Filtered_BSseq_Discovery50_males.rds")
pData(bs) # TD first
samples <- read.csv("Discovery/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo and Cell Type.csv", 
                    stringsAsFactors = FALSE)
samples <- samples[match(sampleNames(bs), samples$Sequencing_ID),]
samples$Diagnosis_Alg <- factor(samples$Diagnosis_Alg, levels = c("TD", "ASD"))

chroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")
background <- loadRegions("Discovery/Dx_Males/bsseq_background_Discovery50_males.csv", chroms = chroms)
background_GR <- makeGRange(background)
meth <- getMeth(bs, regions = background_GR, type = "raw", what = "perRegion") %>% "*"(100)
table(samples$Sequencing_ID == colnames(meth)) # All TRUE

stats <- sapply(1:nrow(meth), function(x) summary(lm(meth[x,] ~ samples$Diagnosis_Alg))$coefficients[2,]) %>% t() %>% 
        cbind(background, .)
stats$stdev <- sapply(1:nrow(meth), function(x) sd_pooled(x = meth[x,], y = samples$Diagnosis_Alg))
write.csv(stats, file = "Discovery/Dx_Males/bsseq_background_vs_Diagnosis_Stats_Discovery50_males.csv",
          quote = FALSE, row.names = FALSE)

# pi0, proportion of true negatives, from background (laptop)
stats <- read.csv("Tables/bsseq_background_vs_Diagnosis_Stats_Discovery50_males.csv", stringsAsFactors = FALSE)
colnames(stats) <- str_replace_all(colnames(stats), pattern = c("Estimate" = "estimate", "Std..Error" = "std_error",
                                                                "t.value" = "t", "Pr...t.." = "p"))
nrow(stats) # 195304
table(is.na(stats$t)) # TRUE 1
stats <- subset(stats, !is.na(stats$t))
pi0 <- pi0est(stats$p)$pi0 # 0.9387358 

# 90th percentile of pooled sd, from background
stdev90 <- quantile(stats$stdev, probs = 0.9) # 13.16744 

# delta, methylation difference of DMRs
DMRs <- loadRegions("DMRs/Discovery/Diagnosis Males 50/DMRs_Dx_Discovery50_males.csv", DMRid = TRUE)
meth <- loadRegions("DMRs/Discovery/Diagnosis Males 50/DMR_raw_methylation_Dx_Discovery50_males.txt") %>%
        .[,str_subset(colnames(.), pattern = "JLCM")] %>% as.matrix() %>% "*"(100)
samples <- read.csv("Tables/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo and Cell Type.csv", 
                    stringsAsFactors = FALSE)
TD <- with(samples, Sequencing_ID[Sequencing_ID %in% colnames(meth) & Diagnosis_Alg == "TD"])
ASD <- with(samples, Sequencing_ID[Sequencing_ID %in% colnames(meth) & Diagnosis_Alg == "ASD"])
DMRs$delta <- rowMeans(meth[,ASD], na.rm = TRUE) - rowMeans(meth[,TD], na.rm = TRUE)
summary(abs(DMRs$delta))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4.149   7.676   9.602  10.359  12.215  25.839 

# Power Analysis ####
# Get Sample Size
n <- 2:300
power <- sapply(c(5,10,15,20), function(x) ssize.twoSamp(delta = x, sigma = stdev90, pi0 = pi0, maxN = max(n))$power[,2])

# Minimum samples per group
apply(power, 2, function(x) n[x > 0.8][1])
#  5% 10% 15% 20%
# 207  54  25  16

# Power for n = 37 per group
power[n == 37,] 
#        5%       10%       15%       20%
# 0.0006216 0.5106642 0.9652089 0.9996737

# Plot
pdf("Figures/Discovery Males Power Analysis Plot ssize fdr.pdf", width = 7, height = 7)
par(mar = c(5,5,2,2))
matplot(n, power, type = "b", pch = 16, cex = 0.8, ylim = c(0, 1), ylab = "Power",
        xlab = "Sample Size Per Group")
grid()
legend(x = 250, y = 0.23, legend = c("5%", "10%", "15%", "20%"), col = c(1:4), fill = c(1:4), title = "Delta")
dev.off()

# Power Analysis for Replication Males --------------------------------------
# Estimate Parameters ####
# Do t-tests on background regions (cluster)
bs <- readRDS("Replication/Dx_Males/Filtered_BSseq_Replication50_males.rds")
pData(bs) # TD first
samples <- read.csv("Discovery/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo and Cell Type.csv", 
                    stringsAsFactors = FALSE)
samples <- samples[match(sampleNames(bs), samples$Sequencing_ID),]
samples$Diagnosis_Alg <- factor(samples$Diagnosis_Alg, levels = c("TD", "ASD"))

chroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")
background <- loadRegions("Replication/Dx_Males/bsseq_background_Replication50_males.csv", chroms = chroms)
background_GR <- makeGRange(background)
meth <- getMeth(bs, regions = background_GR, type = "raw", what = "perRegion") %>% "*"(100)
table(samples$Sequencing_ID == colnames(meth)) # All TRUE

stats <- sapply(1:nrow(meth), function(x) summary(lm(meth[x,] ~ samples$Diagnosis_Alg))$coefficients[2,]) %>% t() %>% 
        cbind(background, .)
stats$stdev <- sapply(1:nrow(meth), function(x) sd_pooled(x = meth[x,], y = samples$Diagnosis_Alg))
write.csv(stats, file = "Replication/Dx_Males/bsseq_background_vs_Diagnosis_Stats_Replication50_males.csv",
          quote = FALSE, row.names = FALSE)

# pi0, proportion of true negatives, from background (laptop)
stats <- read.csv("Tables/bsseq_background_vs_Diagnosis_Stats_Replication50_males.csv", stringsAsFactors = FALSE)
colnames(stats) <- str_replace_all(colnames(stats), pattern = c("Estimate" = "estimate", "Std..Error" = "std_error",
                                                                "t.value" = "t", "Pr...t.." = "p"))
nrow(stats) # 270439
table(is.na(stats$t)) # TRUE 2
stats <- subset(stats, !is.na(stats$t))
pi0 <- pi0est(stats$p)$pi0 # 0.8646789

# 90th percentile of pooled sd, from background
stdev90 <- quantile(stats$stdev, probs = 0.9) # 15.71407  

# delta, methylation difference of DMRs
DMRs <- loadRegions("DMRs/Replication/Diagnosis Males 50/DMRs_Dx_Replication50_males.csv", DMRid = TRUE)
meth <- loadRegions("DMRs/Replication/Diagnosis Males 50/DMR_raw_methylation_Dx_Replication50_males.txt") %>%
        .[,str_subset(colnames(.), pattern = "JLCM")] %>% as.matrix() %>% "*"(100)
samples <- read.csv("Tables/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo and Cell Type.csv", 
                    stringsAsFactors = FALSE)
TD <- with(samples, Sequencing_ID[Sequencing_ID %in% colnames(meth) & Diagnosis_Alg == "TD"])
ASD <- with(samples, Sequencing_ID[Sequencing_ID %in% colnames(meth) & Diagnosis_Alg == "ASD"])
DMRs$delta <- rowMeans(meth[,ASD], na.rm = TRUE) - rowMeans(meth[,TD], na.rm = TRUE)
summary(abs(DMRs$delta))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.381   8.309  10.602  11.622  13.830  39.206 

# Power Analysis ####
# Get Sample Size
n <- 2:300
power <- sapply(c(5,10,15, 20), function(x) ssize.twoSamp(delta = x, sigma = stdev90, pi0 = pi0, maxN = max(n))$power[,2])

# Minimum samples per group
apply(power, 2, function(x) n[x > 0.8][1])
#  5% 10% 15% 20%
# 253  65  30  18

# Power for n = 19 per group
power[n == 19,] 
#        5%       10%       15%       20%
# 0.0000000 0.0302486 0.4468698 0.8472361

# Plot
pdf("Figures/Replication Males Power Analysis Plot ssize fdr.pdf", width = 7, height = 7)
par(mar = c(5,5,2,2))
matplot(n, power, type = "b", pch = 16, cex = 0.8, ylim = c(0, 1), ylab = "Power",
        xlab = "Sample Size Per Group")
grid()
legend(x = 250, y = 0.23, legend = c("5%", "10%", "15%", "20%"), col = c(1:4), fill = c(1:4), title = "Delta")
dev.off()

# Power Analysis for Discovery Females --------------------------------------
# Estimate Parameters ####
# Do t-tests on background regions (cluster)
bs <- readRDS("Discovery/Dx_Females/Filtered_BSseq_Discovery50_females.rds")
pData(bs) # TD first
samples <- read.csv("Discovery/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo and Cell Type.csv", 
                    stringsAsFactors = FALSE)
samples <- samples[match(sampleNames(bs), samples$Sequencing_ID),]
samples$Diagnosis_Alg <- factor(samples$Diagnosis_Alg, levels = c("TD", "ASD"))

chroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")
background <- loadRegions("Discovery/Dx_Females/bsseq_background_Discovery50_females.csv", chroms = chroms)
background_GR <- makeGRange(background)
meth <- getMeth(bs, regions = background_GR, type = "raw", what = "perRegion") %>% "*"(100)
table(samples$Sequencing_ID == colnames(meth)) # All TRUE

stats <- sapply(1:nrow(meth), function(x) summary(lm(meth[x,] ~ samples$Diagnosis_Alg))$coefficients[2,]) %>% t() %>% 
        cbind(background, .)
stats$stdev <- sapply(1:nrow(meth), function(x) sd_pooled(x = meth[x,], y = samples$Diagnosis_Alg))
write.csv(stats, file = "Discovery/Dx_Females/bsseq_background_vs_Diagnosis_Stats_Discovery50_females.csv",
          quote = FALSE, row.names = FALSE)

# pi0, proportion of true negatives, from background (laptop)
stats <- read.csv("Tables/bsseq_background_vs_Diagnosis_Stats_Discovery50_females.csv", stringsAsFactors = FALSE)
colnames(stats) <- str_replace_all(colnames(stats), pattern = c("Estimate" = "estimate", "Std..Error" = "std_error",
                                                                "t.value" = "t", "Pr...t.." = "p"))
nrow(stats) # 209099
table(is.na(stats$t)) # TRUE 2
stats <- subset(stats, !is.na(stats$t))
pi0 <- pi0est(stats$p)$pi0 # 1 ? set to 0.95 (default)

# 90th percentile of pooled sd, from background
stdev90 <- quantile(stats$stdev, probs = 0.9) # 13.22445  

# delta, methylation difference of DMRs
DMRs <- loadRegions("DMRs/Discovery/Diagnosis Females 50/DMRs_Dx_Discovery50_females.csv", DMRid = TRUE)
meth <- loadRegions("DMRs/Discovery/Diagnosis Females 50/DMR_raw_methylation_Dx_Discovery50_females.txt") %>%
        .[,str_subset(colnames(.), pattern = "JLCM")] %>% as.matrix() %>% "*"(100)
samples <- read.csv("Tables/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo and Cell Type.csv", 
                    stringsAsFactors = FALSE)
TD <- with(samples, Sequencing_ID[Sequencing_ID %in% colnames(meth) & Diagnosis_Alg == "TD"])
ASD <- with(samples, Sequencing_ID[Sequencing_ID %in% colnames(meth) & Diagnosis_Alg == "ASD"])
DMRs$delta <- rowMeans(meth[,ASD], na.rm = TRUE) - rowMeans(meth[,TD], na.rm = TRUE)
summary(abs(DMRs$delta))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.214   8.657  10.985  11.962  14.076  34.769 

# Power Analysis ####
# Get Sample Size
n <- 2:300
power <- sapply(c(5,10,15,20), function(x) ssize.twoSamp(delta = x, sigma = stdev90, pi0 = 0.95, maxN = max(n))$power[,2])

# Minimum samples per group
apply(power, 2, function(x) n[x > 0.8][1])
#  5% 10% 15% 20%
# 216  56  27  16

# Power for n = 16 per group
power[n == 16,] 
#       5%      10%      15%       20%
# 0.000000 0.008567 0.354333 0.8124937

# Plot
pdf("Figures/Discovery Females Power Analysis Plot ssize fdr.pdf", width = 7, height = 7)
par(mar = c(5,5,2,2))
matplot(n, power, type = "b", pch = 16, cex = 0.8, ylim = c(0, 1), ylab = "Power",
        xlab = "Sample Size Per Group")
grid()
legend(x = 250, y = 0.23, legend = c("5%", "10%", "15%", "20%"), col = c(1:4), fill = c(1:4), title = "Delta")
dev.off()

# Power Analysis for Replication Females --------------------------------------
# Estimate Parameters ####
# Do t-tests on background regions (cluster)
bs <- readRDS("Replication/Dx_Females_100/Filtered_BSseq_Replication100_females.rds")
pData(bs) # TD first
samples <- read.csv("Discovery/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo and Cell Type.csv", 
                    stringsAsFactors = FALSE)
samples <- samples[match(sampleNames(bs), samples$Sequencing_ID),]
samples$Diagnosis_Alg <- factor(samples$Diagnosis_Alg, levels = c("TD", "ASD"))

chroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")
background <- loadRegions("Replication/Dx_Females_100/bsseq_background_Replication100_females.csv", chroms = chroms)
background_GR <- makeGRange(background)
meth <- getMeth(bs, regions = background_GR, type = "raw", what = "perRegion") %>% "*"(100)
table(samples$Sequencing_ID == colnames(meth)) # All TRUE

stats <- sapply(1:nrow(meth), function(x) summary(lm(meth[x,] ~ samples$Diagnosis_Alg))$coefficients[2,]) %>% t() %>% 
        cbind(background, .)
stats$stdev <- sapply(1:nrow(meth), function(x) sd_pooled(x = meth[x,], y = samples$Diagnosis_Alg))
write.csv(stats, file = "Replication/Dx_Females_100/bsseq_background_vs_Diagnosis_Stats_Replication100_females.csv",
          quote = FALSE, row.names = FALSE)

# pi0, proportion of true negatives, from background (laptop)
stats <- read.csv("Tables/bsseq_background_vs_Diagnosis_Stats_Replication100_females.csv", stringsAsFactors = FALSE)
colnames(stats) <- str_replace_all(colnames(stats), pattern = c("Estimate" = "estimate", "Std..Error" = "std_error",
                                                                "t.value" = "t", "Pr...t.." = "p"))
nrow(stats) # 433649
table(is.na(stats$t)) # TRUE 1430
stats <- subset(stats, !is.na(stats$t))
pi0 <- pi0est(stats$p)$pi0 # 0.9250878

# 90th percentile of pooled sd, from background
stdev90 <- quantile(stats$stdev, probs = 0.9) # 13.51373 

# delta, methylation difference of DMRs
DMRs <- loadRegions("DMRs/Replication/Diagnosis Females 100/DMRs_Dx_Replication100_females.csv", DMRid = TRUE)
meth <- loadRegions("DMRs/Replication/Diagnosis Females 100/DMR_raw_methylation_Dx_Replication100_females.txt") %>%
        .[,str_subset(colnames(.), pattern = "JLCM")] %>% as.matrix() %>% "*"(100)
samples <- read.csv("Tables/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo and Cell Type.csv", 
                    stringsAsFactors = FALSE)
TD <- with(samples, Sequencing_ID[Sequencing_ID %in% colnames(meth) & Diagnosis_Alg == "TD"])
ASD <- with(samples, Sequencing_ID[Sequencing_ID %in% colnames(meth) & Diagnosis_Alg == "ASD"])
DMRs$delta <- rowMeans(meth[,ASD], na.rm = TRUE) - rowMeans(meth[,TD], na.rm = TRUE)
summary(abs(DMRs$delta))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.333  11.452  15.243  17.005  20.795  79.881  

# Power Analysis ####
# Get Sample Size
n <- 2:300
power <- sapply(c(5,10,15,20), function(x) ssize.twoSamp(delta = x, sigma = stdev90, pi0 = pi0, maxN = max(n))$power[,2])

# Minimum samples per group
apply(power, 2, function(x) n[x > 0.8][1])
#  5% 10% 15% 20%
# 210  55  26  16

# Power for n = 4 per group
power[n == 4,] 
# 5% 10% 15% 20%
#  0   0   0   0

# Plot
pdf("Figures/Replication Females Power Analysis Plot ssize fdr.pdf", width = 7, height = 7)
par(mar = c(5,5,2,2))
matplot(n, power, type = "b", pch = 16, cex = 0.8, ylim = c(0, 1), ylab = "Power",
        xlab = "Sample Size Per Group")
grid()
legend(x = 250, y = 0.23, legend = c("5%", "10%", "15%", "20%"), col = c(1:4), fill = c(1:4), title = "Delta")
dev.off()


