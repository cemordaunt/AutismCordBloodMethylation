# Methylation Risk Score Analysis -----------------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 9/17/19

# Packages ####
sapply(c("reshape", "scales", "tidyverse", "bsseq", "annotatr", "dmrseq", "rlist"), require, character.only = TRUE)

# Functions ####
source("R Scripts/DMR Analysis Functions.R")

# Data ####
regDomains <- read.delim("Tables/Regulatory domains hg38.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
samples <- read.csv("Tables/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo and Cell Type.csv", header = TRUE,
                    stringsAsFactors = FALSE)
samples$Sequencing_ID <- factor(samples$Sequencing_ID, levels = unique(samples$Sequencing_ID))
samples$Diagnosis_Alg <- factor(samples$Diagnosis_Alg, levels = c("TD", "ASD"))
samples$Sex <- factor(samples$Sex, levels = c("M", "F"))

# Discovery Males Methylation Risk Score Analysis -------------------------
# Data ####
malesDisc <- subset(samples, Sex == "M" & Platform == "HiSeqX10", select = c("Sequencing_ID", "Diagnosis_Alg", "Sex"))
malesDisc$Sequencing_ID <- as.character(malesDisc$Sequencing_ID)
malesDisc$Diagnosis_Alg_01 <- ifelse(malesDisc$Diagnosis_Alg == "TD", yes = 0, no = 1)
maleChroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")
methDisc <- loadRegions("DMRs/Discovery/Diagnosis Males 50/DMR_raw_methylation_Dx_Discovery50_males.txt", 
                        chroms = maleChroms, sort = TRUE, DMRid = TRUE)

# Calculate Discovery MRS in Discovery Samples ####
# Get Scaled Methylation Data
meth <- methDisc[,as.character(malesDisc$Sequencing_ID)]
rownames(meth) <- methDisc$DMRid
table(colnames(meth) == malesDisc$Sequencing_ID) # All TRUE
meth_scaled <- apply(meth, 1, scale) # Subtracts mean and divides by SD
rownames(meth_scaled) <- malesDisc$Sequencing_ID

# Select DMRs and Calculate Estimates
stats <- apply(meth_scaled, 2, function(x){
        summary(glm(malesDisc$Diagnosis_Alg_01 ~ x, family = "binomial", na.action = na.exclude))$coefficient[2,]
}) %>% t() %>% as.data.frame()
estimates <- data.frame(DMRid = rownames(stats)[stats$`Pr(>|z|)` < 0.05],
                        Estimate = stats$Estimate[stats$`Pr(>|z|)` < 0.05])
estimates$DMRid <- as.character(estimates$DMRid)

# Calculate MRS (see "DMR Analysis Functions.R")
MRSdisc <- getMRS(meth = meth, estimates = estimates)

# Calculate Discovery MRS in Replication Samples ####
# Get Methylation Data
malesRep <- subset(samples, Sex == "M" & Platform == "HiSeq4000", select = c("Sequencing_ID", "Diagnosis_Alg", "Sex"))
malesRep$Diagnosis_Alg_01 <- ifelse(malesRep$Diagnosis_Alg == "TD", yes = 0, no = 1)
malesRep$Sequencing_ID <- as.character(malesRep$Sequencing_ID)
methRep <- loadRegions("DMRs/Discovery/Diagnosis Males 50/Dx_Discovery50_males_DMR_raw_methylation_in_Replication.csv", 
                       chroms = maleChroms, sort = TRUE, DMRid = TRUE)
meth <- methRep[,as.character(malesRep$Sequencing_ID)]
rownames(meth) <- methRep$DMRid
meth2 <- meth

# Calculate MRS
MRSrep <- getMRS(meth = meth, estimates = estimates)

# Merge MRS Tables ####
MRSdisc <- merge(x = malesDisc, y = MRSdisc, by = "Sequencing_ID", all = FALSE, sort = FALSE)
MRSrep <- merge(x = malesRep, y = MRSrep, by = "Sequencing_ID", all = FALSE, sort = FALSE)
MRSall <- rbind(MRSdisc, MRSrep)
MRSall$MRS_scaled <- scale(MRSall$MRS)
MRSall$SampleSet <- factor(c(rep("Discovery", nrow(MRSdisc)), rep("Replication", nrow(MRSrep))), 
                           levels = c("Discovery", "Replication"))
MRSall$Diagnosis_Alg_Called <- ifelse(MRSall$MRS > 0, yes = "ASD", no = "TD") %>% factor(levels = c("TD", "ASD"))

table(MRSall$Diagnosis_Alg_Called[MRSall$SampleSet == "Replication"] == MRSall$Diagnosis_Alg[MRSall$SampleSet == "Replication"])
# FALSE  TRUE 
#    16    22 
# 57.89474%

summary(lm(MRS ~ Diagnosis_Alg_01, data = MRSdisc))
#                  Estimate Std. Error t value Pr(>|t|)    
# Diagnosis_Alg_01   385.42      10.69   36.06   <2e-16
        
summary(lm(MRS ~ Diagnosis_Alg_01, data = MRSrep))
#                  Estimate Std. Error t value Pr(>|t|)
# Diagnosis_Alg_01   13.941     11.308   1.233    0.226

g <- ggplot(data = MRSall)
g + 
        geom_hline(yintercept = 0, size = 1.25, lty = 2, color = "black") +
        stat_summary(aes(x = Diagnosis_Alg, y = MRS), fun.data = "mean_cl_boot", geom = "crossbar", color = "black", 
                     fill = "white", size = 0.7, position = "dodge2") +
        geom_dotplot(aes(x = Diagnosis_Alg, y = MRS, fill = Diagnosis_Alg_Called, color = Diagnosis_Alg_Called), 
                     binwidth = 15, binaxis = "y", stackdir = "center", stackratio = 1.15, dotsize = 0.8) + 
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.89, 0.96), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              strip.text = element_text(size = 22), plot.margin = unit(c(0.3, 1, 1, 1), "lines"), 
              axis.title.x = element_blank(), strip.background = element_blank()) +
        scale_y_continuous(breaks = pretty_breaks(n = 6)) +
        ylab("Methylation Risk Score") +
        facet_wrap(vars(SampleSet)) +
        scale_fill_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366")) +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366"))
ggsave("Figures/Discovery Males DMR Methylation Risk Score in Discovery and Replication Dotplot Mean CL.png", dpi = 600, width = 9, 
       height = 7, units = "in")

g <- ggplot(data = MRSall)
g + 
        geom_hline(yintercept = 0, size = 1.25, lty = 2, color = "black") +
        stat_summary(aes(x = Diagnosis_Alg, y = MRS_scaled), fun.data = "mean_cl_boot", geom = "crossbar", color = "black", 
                     fill = "white", size = 0.7, position = "dodge2") +
        geom_dotplot(aes(x = Diagnosis_Alg, y = MRS_scaled, fill = Diagnosis_Alg_Called, color = Diagnosis_Alg_Called), binwidth = 0.09,
                     binaxis = "y", stackdir = "center", stackratio = 1.15, dotsize = 0.8) + 
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.89, 0.96), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              strip.text = element_text(size = 22), plot.margin = unit(c(0.3, 1, 1, 1), "lines"), 
              axis.title.x = element_blank(), strip.background = element_blank()) +
        scale_y_continuous(breaks = pretty_breaks(n = 6)) +
        ylab("Methylation Risk Score") +
        facet_wrap(vars(SampleSet)) +
        scale_fill_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366")) +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366"))
ggsave("Figures/Discovery Males DMR Scaled Methylation Risk Score in Discovery and Replication Dotplot Mean CL.png", dpi = 600, width = 9, 
       height = 7, units = "in")

# MRS with Top Quartile of DMRs ####
methDiscTemp <- methDisc[,as.character(malesDisc$Sequencing_ID)]
rownames(methDiscTemp) <- methDisc$DMRid
methRepTemp <- methRep[,as.character(malesRep$Sequencing_ID)]
rownames(methRepTemp) <- methRep$DMRid

MRStop <- MRSall[,c("Sequencing_ID", "Diagnosis_Alg", "Sex", "SampleSet")]
for(i in seq(0, 0.99, 0.01)){
        estimates_top <- subset(estimates, abs(estimates$Estimate) >= quantile(abs(estimates$Estimate), probs = i))
        tempDisc <- getMRS(meth = methDiscTemp, estimates = estimates_top)
        tempRep <- getMRS(meth = methRepTemp, estimates = estimates_top)
        MRStop$temp <- c(tempDisc$MRS, tempRep$MRS)
        colnames(MRStop) <- c(colnames(MRStop)[!colnames(MRStop) == "temp"], paste("MRS", i, sep = "_"))
}
rm(estimates_top, tempDisc, tempRep)

MRScalls <- sapply(MRStop[,grepl("MRS", colnames(MRStop), fixed = TRUE)], function(x) ifelse(x > 0, yes = "ASD", no = "TD"))
precisionDisc <- apply(MRScalls[MRStop$SampleSet == "Discovery",], 2,
                       function(x) table(x == MRStop$Diagnosis_Alg[MRStop$SampleSet == "Discovery"])["TRUE"] * 100 /
                               length(MRStop$Diagnosis_Alg[MRStop$SampleSet == "Discovery"]))
precisionRep <- apply(MRScalls[MRStop$SampleSet == "Replication",], 2,
                       function(x) table(x == MRStop$Diagnosis_Alg[MRStop$SampleSet == "Replication"])["TRUE"] * 100 /
                               length(MRStop$Diagnosis_Alg[MRStop$SampleSet == "Replication"]))
precision <- data.frame(Quantile = seq(0, 0.99, 0.01), precisionDisc = precisionDisc, precisionRep = precisionRep)
subset(precision, precisionRep == max(precision$precisionRep))
#          Quantile precisionDisc precisionRep
# MRS_0.18     0.18     100.00000     63.15789
# MRS_0.19     0.19     100.00000     63.15789
# MRS_0.95     0.95      98.64865     63.15789
# MRS_0.96     0.96      97.29730     63.15789
# MRS_0.97     0.97      95.94595     63.15789

pdf("Figures/Discovery Males DMR Methylation Risk Score Precision by DMR Quantile Line Plot.pdf", width = 6, height = 6)
par(mar = c(5,5,1,1))
plot(precisionDisc ~ Quantile, data = precision, ylim = c(0, 100), type = "l", lwd = 1.5, ylab = "Precision (%)", col = "blue")
lines(precisionRep ~ Quantile, data = precision, lwd = 1.5, col = "red")
abline(h = 50, lty = 2, lwd = 1.5)
dev.off()

# Best at 0.95, 30 DMRs
estimates_0.95 <- subset(estimates, abs(estimates$Estimate) >= quantile(abs(estimates$Estimate), probs = 0.95))
summary(abs(estimates_0.95$Estimate))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.392   1.440   1.577   1.679   1.689   3.059 

MRSall$MRS_0.95 <- MRStop$MRS_0.95
MRSall$MRS_0.95_scaled <- scale(MRSall$MRS_0.95)
MRSall$Diagnosis_Alg_Called_0.95 <- ifelse(MRSall$MRS_0.95 > 0, yes = "ASD", no = "TD") %>% factor(levels = c("TD", "ASD"))

# Plot
g <- ggplot(data = MRSall)
g + 
        geom_hline(yintercept = 0, size = 1.25, lty = 2, color = "black") +
        stat_summary(aes(x = Diagnosis_Alg, y = MRS_0.95), fun.data = "mean_cl_boot", geom = "crossbar", color = "black", 
                     fill = "white", size = 0.7, position = "dodge2") +
        geom_dotplot(aes(x = Diagnosis_Alg, y = MRS_0.95, fill = Diagnosis_Alg_Called_0.95, color = Diagnosis_Alg_Called_0.95), 
                     binwidth = 2, binaxis = "y", stackdir = "center", stackratio = 1.15, dotsize = 0.8) + 
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.89, 0.96), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              strip.text = element_text(size = 22), plot.margin = unit(c(0.3, 1, 1, 1), "lines"), 
              axis.title.x = element_blank(), strip.background = element_blank()) +
        scale_y_continuous(breaks = pretty_breaks(n = 6)) +
        ylab("Methylation Risk Score") +
        facet_wrap(vars(SampleSet)) +
        scale_fill_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366")) +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366"))
ggsave("Figures/Discovery Males Top DMR Methylation Risk Score in Discovery and Replication Dotplot Mean CL.png", dpi = 600, width = 9, 
       height = 7, units = "in")

g <- ggplot(data = MRSall)
g + 
        geom_hline(yintercept = 0, size = 1.25, lty = 2, color = "black") +
        stat_summary(aes(x = Diagnosis_Alg, y = MRS_0.95_scaled), fun.data = "mean_cl_boot", geom = "crossbar", color = "black", 
                     fill = "white", size = 0.7, position = "dodge2") +
        geom_dotplot(aes(x = Diagnosis_Alg, y = MRS_0.95_scaled, fill = Diagnosis_Alg_Called_0.95, color = Diagnosis_Alg_Called_0.95), binwidth = 0.1,
                     binaxis = "y", stackdir = "center", stackratio = 1.15, dotsize = 0.8) + 
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.89, 0.96), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              strip.text = element_text(size = 22), plot.margin = unit(c(0.3, 1, 1, 1), "lines"), 
              axis.title.x = element_blank(), strip.background = element_blank()) +
        scale_y_continuous(breaks = pretty_breaks(n = 6)) +
        ylab("Methylation Risk Score") +
        facet_wrap(vars(SampleSet)) +
        scale_fill_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366")) +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366"))
ggsave("Figures/Discovery Males Top DMR Scaled Methylation Risk Score in Discovery and Replication Dotplot Mean CL.png", dpi = 600, width = 9, 
       height = 7, units = "in")

# Covariate Correlation ####
# Covariate Data
catVars <- c("Diagnosis_Alg", "Study", "Site", "MomEdu_detail", "DM1or2", "GDM", "PE", "home_ownership", 
             "marital_status", "SmokeYN_Pregnancy")
contVars <- colnames(samples)[!colnames(samples) %in% catVars]
contVars <- contVars[!contVars %in% c("COI_ID", "Sequencing_ID", "Cord_Blood_IBC", "Platform", "Sex")]

samples$Study <- factor(samples$Study, levels = c("MARBLES", "EARLI"))
samples$Site <- factor(samples$Site, levels = c("UC Davis", "Kaiser Permanente", "Drexel", "Johns Hopkins University"))
samples$MomEdu_detail <- factor(samples$MomEdu_detail, levels = c(6, 1:5, 7,8))
samples$home_ownership[samples$home_ownership == 99] <- NA
samples$marital_status[samples$marital_status == 99] <- NA
factorCols <- c("DM1or2", "GDM", "PE", "marital_status", "home_ownership", "SmokeYN_Pregnancy")
samples[,factorCols] <- lapply(samples[,factorCols], as.factor)

malesDisc <- subset(samples, Platform == "HiSeqX10" & Sex == "M") %>%
        merge(y = MRSall[,c("Sequencing_ID", "MRS")], by = "Sequencing_ID", all = FALSE, sort = FALSE)
malesDisc$MRS_scaled <- scale(malesDisc$MRS)

# Run Stats
MRS_stats <- DMRmethLm(DMRs = "MRS_scaled", catVars = catVars, contVars = contVars, sampleData = malesDisc, 
                       file = "Tables/Covariate Stats by MRS Discovery Males.txt")
malesDiscASD <- subset(malesDisc, Diagnosis_Alg == "ASD")
MRS_stats_ASD <- DMRmethLm(DMRs = "MRS_scaled", catVars = catVars, contVars = contVars, sampleData = malesDiscASD,
                           file = "Tables/Covariate Stats by MRS Discovery ASD Males.txt")
malesDiscTD <- subset(malesDisc, Diagnosis_Alg == "TD")
MRS_stats_TD <- DMRmethLm(DMRs = "MRS_scaled", catVars = catVars, contVars = contVars, sampleData = malesDiscTD,
                           file = "Tables/Covariate Stats by MRS Discovery TD Males.txt")
MRS_stats <- merge(x = MRS_stats, y = MRS_stats_TD, by = "Variable", all.x = TRUE, all.y = FALSE, 
                   suffixes = c("", "_TD"), sort = FALSE) %>%
        merge(y = MRS_stats_ASD, by = "Variable", all.x = TRUE, all.y = FALSE, suffixes = c("", "_ASD"), sort = FALSE)

# Plot Results
plotData <- data.frame(Variable = factor(MRS_stats$Variable, levels = rev(levels(MRS_stats$Variable))), 
                       Diagnosis = factor(rep(c("All", "TD", "ASD"), each = nrow(MRS_stats)), 
                                          levels = c("All", "TD", "ASD")),
                       Estimate = c(MRS_stats$Estimate, MRS_stats$Estimate_TD, MRS_stats$Estimate_ASD),
                       qvalue = c(MRS_stats$qvalue, MRS_stats$qvalue_TD, MRS_stats$qvalue_ASD))
plotData$Significant <- (plotData$qvalue < 0.05) %>% factor(levels = c("TRUE", "FALSE"))
plotData <- subset(plotData, !is.na(plotData$Estimate) & !grepl("MomEdu", plotData$Variable, fixed = TRUE) &
                           !plotData$Variable == c("final_creatinine_mgdl"))

varNames <- as.character(plotData$Variable) %>% 
        str_replace_all(c("Diagnosis_Alg" = "Diagnosis", "GDM" = "Gestational Diabetes", "home_ownership" = "Own Home",
                          "marital_status" = "Married", "ADOScs" = "ADOS", "MSLelcStandard36" = "Mullen Composite",
                          "MSLelTscore36" = "Mullen Expressive Language", "MSLfmTscore36" = "Mullen Fine Motor",
                          "MSLrlTscore36" = "Mullen Receptive Language", "MSLvrTscore36" = "Mullen Visual Reception",
                          "percent_trimmed" = "Bases Trimmed", "percent_aligned" = "Aligned Reads", "percent_duplicate" = "PCR Duplicates",
                          "dedup_reads_M" = "Unique Reads", "C_coverage" = "C Coverage", "CG_coverage" = "CpG Coverage",
                          "percent_cpg_meth_bsseq" = "CpG Methylation", "percent_chg_meth" = "CHG Methylation", 
                          "percent_chh_meth" = "CHH Methylation", "ga_w" = "Gestational Age", "bw_g" = "Birthweight", 
                          "MomAgeYr" = "Maternal Age", "Mat_Height_cm" = "Maternal Height", 
                          "Mat_Weight_kg_PrePreg" = "Maternal Weight", "Mat_BMI_PrePreg" = "Maternal BMI", 
                          "parity" = "Parity", "dad_age" = "Paternal Age", "SmokeYN_Pregnancy" = "Maternal Smoking",
                          "cotinine_urine_ngml" = "Urine Cotinine", "DM1or2" = "Maternal Diabetes",
                          "Bcell" = "B Cells", "CD4T" = "CD4 T Cells", "CD8T" = "CD8 T Cells", "Gran" = "Granulocytes", 
                          "Mono" = "Monocytes", "NK" = "NK Cells", "nRBC" = "nRBCs", "Site_" = "", "_" = " ", 
                          " University" = ""))
plotData$Variable <- factor(varNames, 
                             levels = rev(c("Diagnosis", "ADOS", "Mullen Composite", "Mullen Expressive Language",
                                            "Mullen Fine Motor", "Mullen Receptive Language", "Mullen Visual Reception",
                                            "Study", "Kaiser Permanente", "Drexel", "Johns Hopkins",
                                            "Gestational Age", "Birthweight", "Paternal Age", "Maternal Age",
                                            "Maternal Height", "Maternal Weight", "Maternal BMI", "Maternal Diabetes",
                                            "Gestational Diabetes", "Parity", "Maternal Smoking", "Urine Cotinine", "Married", "Own Home",
                                            "Bases Trimmed", "Aligned Reads", "PCR Duplicates", "Unique Reads", "C Coverage", 
                                            "CpG Coverage", "CpG Methylation", "CHG Methylation", "CHH Methylation", "B Cells", "CD4 T Cells",
                                            "CD8 T Cells", "Granulocytes", "Monocytes", "NK Cells", "nRBCs")))
limit <- 0.76
g <- ggplot(data = plotData)
g + 
        geom_tile(aes(x = Diagnosis, y = Variable, fill = Estimate, color = Estimate)) + 
        geom_text(aes(x = Diagnosis, y = Variable, alpha = Significant, label = "*"), color = "white", size = 15, 
                  nudge_y = -0.45) +
        scale_fill_gradientn("Change\nin MRS", colors = c("#3366CC", "Black", "#FF0000"), values = c(0, 1), na.value = "#FF0000", 
                             limits = c(-limit, limit), breaks = pretty_breaks(n = 3)) +
        scale_color_gradientn("Change\nin MRS", colors = c("#3366CC", "Black", "#FF0000"), values = c(0, 1), na.value = "#FF0000", 
                              limits = c(-limit, limit), breaks = pretty_breaks(n = 3)) +
        scale_alpha_manual(breaks = c("TRUE", "FALSE"), values = c(1, 0), guide = FALSE) +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'vertical', legend.position = c(1.18, 0.91), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "Black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_text(size = 24),
              axis.text.x = element_text(color = "Black", size = 24), 
              axis.text.y = element_text(color = "Black", size = 21), legend.text = element_text(size = 22),
              legend.background = element_blank(), panel.background = element_rect(fill = "black"),
              strip.text = element_text(size = 26), plot.margin = unit(c(0.5, 8, 1, 1), "lines"), 
              axis.title = element_blank(), strip.background = element_blank())
ggsave("Figures/MRS Covariate Association Discovery Males Heatmap.png", dpi = 600, width = 10, height = 12, units = "in")

# Discovery Females Methylation Risk Score Analysis -------------------------
# Data ####
femalesDisc <- subset(samples, Sex == "F" & Platform == "HiSeqX10", select = c("Sequencing_ID", "Diagnosis_Alg", "Sex"))
femalesDisc$Sequencing_ID <- as.character(femalesDisc$Sequencing_ID)
femalesDisc$Diagnosis_Alg_01 <- ifelse(femalesDisc$Diagnosis_Alg == "TD", yes = 0, no = 1)
femaleChroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrM")
methDisc <- loadRegions("DMRs/Discovery/Diagnosis Females 50/DMR_raw_methylation_Dx_Discovery50_females.txt", 
                        chroms = femaleChroms, sort = TRUE, DMRid = TRUE)

# Calculate Discovery MRS in Discovery Samples ####
# Get Scaled Methylation Data
meth <- methDisc[,as.character(femalesDisc$Sequencing_ID)]
rownames(meth) <- methDisc$DMRid
table(colnames(meth) == femalesDisc$Sequencing_ID) # All TRUE
meth_scaled <- apply(meth, 1, scale) # Subtracts mean and divides by SD
rownames(meth_scaled) <- femalesDisc$Sequencing_ID

# Select DMRs and Calculate Estimates
stats <- apply(meth_scaled, 2, function(x){
        summary(glm(femalesDisc$Diagnosis_Alg_01 ~ x, family = "binomial", na.action = na.exclude))$coefficient[2,]
}) %>% t() %>% as.data.frame()
estimates <- data.frame(DMRid = rownames(stats)[stats$`Pr(>|z|)` < 0.05],
                        Estimate = stats$Estimate[stats$`Pr(>|z|)` < 0.05])
estimates$DMRid <- as.character(estimates$DMRid)

# Calculate MRS (see "DMR Analysis Functions.R")
MRSdisc <- getMRS(meth = meth, estimates = estimates)

# Calculate Discovery MRS in Replication Samples ####
# Get Methylation Data
femalesRep <- subset(samples, Sex == "F" & Platform == "HiSeq4000", select = c("Sequencing_ID", "Diagnosis_Alg", "Sex"))
femalesRep$Diagnosis_Alg_01 <- ifelse(femalesRep$Diagnosis_Alg == "TD", yes = 0, no = 1)
femalesRep$Sequencing_ID <- as.character(femalesRep$Sequencing_ID)
methRep <- loadRegions("DMRs/Discovery/Diagnosis Females 50/Dx_Discovery50_females_DMR_raw_methylation_in_Replication.csv", 
                       chroms = femaleChroms, sort = TRUE, DMRid = TRUE)
meth <- methRep[,as.character(femalesRep$Sequencing_ID)]
rownames(meth) <- methRep$DMRid

# Calculate MRS
MRSrep <- getMRS(meth = meth, estimates = estimates)

# Merge MRS Tables ####
MRSdisc <- merge(x = femalesDisc, y = MRSdisc, by = "Sequencing_ID", all = FALSE, sort = FALSE)
MRSrep <- merge(x = femalesRep, y = MRSrep, by = "Sequencing_ID", all = FALSE, sort = FALSE)
MRSall <- rbind(MRSdisc, MRSrep)
MRSall$MRS_scaled <- scale(MRSall$MRS)
MRSall$SampleSet <- factor(c(rep("Discovery", nrow(MRSdisc)), rep("Replication", nrow(MRSrep))), 
                           levels = c("Discovery", "Replication"))
MRSall$Diagnosis_Alg_Called <- ifelse(MRSall$MRS > 0, yes = "ASD", no = "TD") %>% factor(levels = c("TD", "ASD"))

table(MRSall$Diagnosis_Alg_Called[MRSall$SampleSet == "Replication"] == MRSall$Diagnosis_Alg[MRSall$SampleSet == "Replication"])
# FALSE  TRUE 
#     5     3 
# 37.5%

summary(lm(MRS ~ Diagnosis_Alg_01, data = MRSdisc))
#                  Estimate Std. Error t value Pr(>|t|)    
# Diagnosis_Alg_01  2511.70      65.67   38.25   <2e-16 ***
        
summary(lm(MRS ~ Diagnosis_Alg_01, data = MRSrep))
#                  Estimate Std. Error t value Pr(>|t|)
# Diagnosis_Alg_01   -46.13      63.85  -0.722    0.497

g <- ggplot(data = MRSall)
g + 
        geom_hline(yintercept = 0, size = 1.25, lty = 2, color = "black") +
        stat_summary(aes(x = Diagnosis_Alg, y = MRS), fun.data = "mean_cl_boot", geom = "crossbar", color = "black", 
                     fill = "white", size = 0.7, position = "dodge2") +
        geom_dotplot(aes(x = Diagnosis_Alg, y = MRS, fill = Diagnosis_Alg_Called, color = Diagnosis_Alg_Called), 
                     binwidth = 75, binaxis = "y", stackdir = "center", stackratio = 1.15, dotsize = 0.8) + 
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.89, 0.96), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              strip.text = element_text(size = 22), plot.margin = unit(c(0.3, 1, 1, 1), "lines"), 
              axis.title.x = element_blank(), strip.background = element_blank()) +
        scale_y_continuous(breaks = pretty_breaks(n = 6)) +
        ylab("Methylation Risk Score") +
        facet_wrap(vars(SampleSet)) +
        scale_fill_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366")) +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366"))
ggsave("Figures/Discovery Females DMR Methylation Risk Score in Discovery and Replication Dotplot Mean CL.png", dpi = 600, width = 9, 
       height = 7, units = "in")

g <- ggplot(data = MRSall)
g + 
        geom_hline(yintercept = 0, size = 1.25, lty = 2, color = "black") +
        stat_summary(aes(x = Diagnosis_Alg, y = MRS_scaled), fun.data = "mean_cl_boot", geom = "crossbar", color = "black", 
                     fill = "white", size = 0.7, position = "dodge2") +
        geom_dotplot(aes(x = Diagnosis_Alg, y = MRS_scaled, fill = Diagnosis_Alg_Called, color = Diagnosis_Alg_Called), binwidth = 0.065,
                     binaxis = "y", stackdir = "center", stackratio = 1.15, dotsize = 0.8) + 
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.89, 0.96), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              strip.text = element_text(size = 22), plot.margin = unit(c(0.3, 1, 1, 1), "lines"), 
              axis.title.x = element_blank(), strip.background = element_blank()) +
        scale_y_continuous(breaks = pretty_breaks(n = 6)) +
        ylab("Methylation Risk Score") +
        facet_wrap(vars(SampleSet)) +
        scale_fill_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366")) +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366"))
ggsave("Figures/Discovery Females DMR Scaled Methylation Risk Score in Discovery and Replication Dotplot Mean CL.png", dpi = 600, width = 9, 
       height = 7, units = "in")

# MRS with Top Quartile of DMRs ####
methDiscTemp <- methDisc[,as.character(femalesDisc$Sequencing_ID)]
rownames(methDiscTemp) <- methDisc$DMRid
methRepTemp <- methRep[,as.character(femalesRep$Sequencing_ID)]
rownames(methRepTemp) <- methRep$DMRid

MRStop <- MRSall[,c("Sequencing_ID", "Diagnosis_Alg", "Sex", "SampleSet")]
for(i in seq(0, 0.99, 0.01)){
        estimates_top <- subset(estimates, abs(estimates$Estimate) >= quantile(abs(estimates$Estimate), probs = i))
        tempDisc <- getMRS(meth = methDiscTemp, estimates = estimates_top)
        tempRep <- getMRS(meth = methRepTemp, estimates = estimates_top)
        MRStop$temp <- c(tempDisc$MRS, tempRep$MRS)
        colnames(MRStop) <- c(colnames(MRStop)[!colnames(MRStop) == "temp"], paste("MRS", i, sep = "_"))
}
rm(estimates_top, tempDisc, tempRep)

MRScalls <- sapply(MRStop[,grepl("MRS", colnames(MRStop), fixed = TRUE)], function(x) ifelse(x > 0, yes = "ASD", no = "TD"))
precisionDisc <- apply(MRScalls[MRStop$SampleSet == "Discovery",], 2,
                       function(x) table(x == MRStop$Diagnosis_Alg[MRStop$SampleSet == "Discovery"])["TRUE"] * 100 /
                               length(MRStop$Diagnosis_Alg[MRStop$SampleSet == "Discovery"]))
precisionRep <- apply(MRScalls[MRStop$SampleSet == "Replication",], 2,
                      function(x) table(x == MRStop$Diagnosis_Alg[MRStop$SampleSet == "Replication"])["TRUE"] * 100 /
                              length(MRStop$Diagnosis_Alg[MRStop$SampleSet == "Replication"]))
precision <- data.frame(Quantile = seq(0, 0.99, 0.01), precisionDisc = precisionDisc, precisionRep = precisionRep)
subset(precision, precisionRep == max(precision$precisionRep))
#          Quantile precisionDisc precisionRep
# MRS_0.73     0.73           100           50
pdf("Figures/Discovery Females DMR Methylation Risk Score Precision by DMR Quantile Line Plot.pdf", width = 6, height = 6)
par(mar = c(5,5,1,1))
plot(precisionDisc ~ Quantile, data = precision, ylim = c(0, 100), type = "l", lwd = 1.5, ylab = "Precision (%)", col = "blue")
lines(precisionRep ~ Quantile, data = precision, lwd = 1.5, col = "red")
abline(h = 50, lty = 2, lwd = 1.5)
dev.off()

# Best at 0.73, 430 DMRs
estimates_0.73 <- subset(estimates, abs(estimates$Estimate) >= quantile(abs(estimates$Estimate), probs = 0.73))
summary(abs(estimates_0.73$Estimate))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.735   1.919   2.098   2.328   2.523   6.782 

MRSall$MRS_0.73 <- MRStop$MRS_0.73
MRSall$MRS_0.73_scaled <- scale(MRSall$MRS_0.73)
MRSall$Diagnosis_Alg_Called_0.73 <- ifelse(MRSall$MRS_0.73 > 0, yes = "ASD", no = "TD") %>% factor(levels = c("TD", "ASD"))

# Plot
g <- ggplot(data = MRSall)
g + 
        geom_hline(yintercept = 0, size = 1.25, lty = 2, color = "black") +
        stat_summary(aes(x = Diagnosis_Alg, y = MRS_0.73), fun.data = "mean_cl_boot", geom = "crossbar", color = "black", 
                     fill = "white", size = 0.7, position = "dodge2") +
        geom_dotplot(aes(x = Diagnosis_Alg, y = MRS_0.73, fill = Diagnosis_Alg_Called_0.73, color = Diagnosis_Alg_Called_0.73), 
                     binwidth = 40, binaxis = "y", stackdir = "center", stackratio = 1.15, dotsize = 0.8) + 
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.89, 0.96), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              strip.text = element_text(size = 22), plot.margin = unit(c(0.3, 1, 1, 1), "lines"), 
              axis.title.x = element_blank(), strip.background = element_blank()) +
        scale_y_continuous(breaks = pretty_breaks(n = 7)) +
        ylab("Methylation Risk Score") +
        facet_wrap(vars(SampleSet)) +
        scale_fill_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366")) +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366"))
ggsave("Figures/Discovery Females Top DMR Methylation Risk Score in Discovery and Replication Dotplot Mean CL.png", dpi = 600, width = 9, 
       height = 7, units = "in")

g <- ggplot(data = MRSall)
g + 
        geom_hline(yintercept = 0, size = 1.25, lty = 2, color = "black") +
        stat_summary(aes(x = Diagnosis_Alg, y = MRS_0.73_scaled), fun.data = "mean_cl_boot", geom = "crossbar", color = "black", 
                     fill = "white", size = 0.7, position = "dodge2") +
        geom_dotplot(aes(x = Diagnosis_Alg, y = MRS_0.73_scaled, fill = Diagnosis_Alg_Called_0.73, color = Diagnosis_Alg_Called_0.73), binwidth = 0.08,
                     binaxis = "y", stackdir = "center", stackratio = 1.15, dotsize = 0.8) + 
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.89, 0.96), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              strip.text = element_text(size = 22), plot.margin = unit(c(0.3, 1, 1, 1), "lines"), 
              axis.title.x = element_blank(), strip.background = element_blank()) +
        scale_y_continuous(breaks = pretty_breaks(n = 6)) +
        ylab("Methylation Risk Score") +
        facet_wrap(vars(SampleSet)) +
        scale_fill_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366")) +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366"))
ggsave("Figures/Discovery Females Top DMR Scaled Methylation Risk Score in Discovery and Replication Dotplot Mean CL.png", dpi = 600, width = 9, 
       height = 7, units = "in")

# Covariate Correlation ####
# Covariate Data
catVars <- c("Diagnosis_Alg", "Study", "Site", "MomEdu_detail", "DM1or2", "GDM", "PE", "home_ownership", 
             "marital_status", "SmokeYN_Pregnancy")
contVars <- colnames(samples)[!colnames(samples) %in% catVars]
contVars <- contVars[!contVars %in% c("COI_ID", "Sequencing_ID", "Cord_Blood_IBC", "Platform", "Sex")]
femalesDisc <- subset(samples, Platform == "HiSeqX10" & Sex == "F") %>%
        merge(y = MRSall[,c("Sequencing_ID", "MRS")], by = "Sequencing_ID", all = FALSE, sort = FALSE)
femalesDisc$MRS_scaled <- scale(femalesDisc$MRS)

# Run Stats
MRS_stats <- DMRmethLm(DMRs = "MRS_scaled", catVars = catVars, contVars = contVars, sampleData = femalesDisc, 
                       file = "Tables/Covariate Stats by MRS Discovery Females.txt")
femalesDiscASD <- subset(femalesDisc, Diagnosis_Alg == "ASD")
MRS_stats_ASD <- DMRmethLm(DMRs = "MRS_scaled", catVars = catVars, contVars = contVars, sampleData = femalesDiscASD,
                           file = "Tables/Covariate Stats by MRS Discovery ASD Females.txt")
femalesDiscTD <- subset(femalesDisc, Diagnosis_Alg == "TD")
MRS_stats_TD <- DMRmethLm(DMRs = "MRS_scaled", catVars = catVars, contVars = contVars, sampleData = femalesDiscTD,
                          file = "Tables/Covariate Stats by MRS Discovery TD Females.txt")
MRS_stats <- merge(x = MRS_stats, y = MRS_stats_TD, by = "Variable", all.x = TRUE, all.y = FALSE, 
                   suffixes = c("", "_TD"), sort = FALSE) %>%
        merge(y = MRS_stats_ASD, by = "Variable", all.x = TRUE, all.y = FALSE, suffixes = c("", "_ASD"), sort = FALSE)

# Plot Results
plotData <- data.frame(Variable = factor(MRS_stats$Variable, levels = rev(levels(MRS_stats$Variable))), 
                       Diagnosis = factor(rep(c("All", "TD", "ASD"), each = nrow(MRS_stats)), 
                                          levels = c("All", "TD", "ASD")),
                       Estimate = c(MRS_stats$Estimate, MRS_stats$Estimate_TD, MRS_stats$Estimate_ASD),
                       qvalue = c(MRS_stats$qvalue, MRS_stats$qvalue_TD, MRS_stats$qvalue_ASD))
plotData$Significant <- (plotData$qvalue < 0.05) %>% factor(levels = c("TRUE", "FALSE"))
plotData <- subset(plotData, !is.na(plotData$Estimate) & !grepl("MomEdu", plotData$Variable, fixed = TRUE) &
                           !plotData$Variable == c("final_creatinine_mgdl"))

varNames <- as.character(plotData$Variable) %>% 
        str_replace_all(c("Diagnosis_Alg" = "Diagnosis", "GDM" = "Gestational Diabetes", "home_ownership" = "Own Home",
                          "marital_status" = "Married", "ADOScs" = "ADOS", "MSLelcStandard36" = "Mullen Composite",
                          "MSLelTscore36" = "Mullen Expressive Language", "MSLfmTscore36" = "Mullen Fine Motor",
                          "MSLrlTscore36" = "Mullen Receptive Language", "MSLvrTscore36" = "Mullen Visual Reception",
                          "percent_trimmed" = "Bases Trimmed", "percent_aligned" = "Aligned Reads", "percent_duplicate" = "PCR Duplicates",
                          "dedup_reads_M" = "Unique Reads", "C_coverage" = "C Coverage", "CG_coverage" = "CpG Coverage",
                          "percent_cpg_meth_bsseq" = "CpG Methylation", "percent_chg_meth" = "CHG Methylation", 
                          "percent_chh_meth" = "CHH Methylation", "ga_w" = "Gestational Age", "bw_g" = "Birthweight", 
                          "MomAgeYr" = "Maternal Age", "Mat_Height_cm" = "Maternal Height", "PE" = "Pre-eclampsia",
                          "Mat_Weight_kg_PrePreg" = "Maternal Weight", "Mat_BMI_PrePreg" = "Maternal BMI", 
                          "parity" = "Parity", "dad_age" = "Paternal Age", "SmokeYN_Pregnancy" = "Maternal Smoking",
                          "cotinine_urine_ngml" = "Urine Cotinine", "DM1or2" = "Maternal Diabetes",
                          "Bcell" = "B Cells", "CD4T" = "CD4 T Cells", "CD8T" = "CD8 T Cells", "Gran" = "Granulocytes", 
                          "Mono" = "Monocytes", "NK" = "NK Cells", "nRBC" = "nRBCs", "Site_" = "", "_" = " ", 
                          " University" = ""))
plotData$Variable <- factor(varNames, 
                            levels = rev(c("Diagnosis", "ADOS", "Mullen Composite", "Mullen Expressive Language",
                                           "Mullen Fine Motor", "Mullen Receptive Language", "Mullen Visual Reception",
                                           "Study", "Kaiser Permanente", "Drexel", "Johns Hopkins",
                                           "Gestational Age", "Birthweight", "Paternal Age", "Maternal Age",
                                           "Maternal Height", "Maternal Weight", "Maternal BMI", "Maternal Diabetes",
                                           "Gestational Diabetes", "Pre-eclampsia", "Parity", "Maternal Smoking", "Urine Cotinine", "Married", "Own Home",
                                           "Bases Trimmed", "Aligned Reads", "PCR Duplicates", "Unique Reads", "C Coverage", 
                                           "CpG Coverage", "CpG Methylation", "CHG Methylation", "CHH Methylation", "B Cells", "CD4 T Cells",
                                           "CD8 T Cells", "Granulocytes", "Monocytes", "NK Cells", "nRBCs")))
limit <- 0.83
g <- ggplot(data = plotData)
g + 
        geom_tile(aes(x = Diagnosis, y = Variable, fill = Estimate, color = Estimate)) + 
        geom_text(aes(x = Diagnosis, y = Variable, alpha = Significant, label = "*"), color = "white", size = 15, 
                  nudge_y = -0.45) +
        scale_fill_gradientn("Change\nin MRS", colors = c("#3366CC", "Black", "#FF0000"), values = c(0, 1), na.value = "#FF0000", 
                             limits = c(-limit, limit), breaks = pretty_breaks(n = 3)) +
        scale_color_gradientn("Change\nin MRS", colors = c("#3366CC", "Black", "#FF0000"), values = c(0, 1), na.value = "#FF0000", 
                              limits = c(-limit, limit), breaks = pretty_breaks(n = 3)) +
        scale_alpha_manual(breaks = c("TRUE", "FALSE"), values = c(1, 0), guide = FALSE) +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'vertical', legend.position = c(1.18, 0.91), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "Black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_text(size = 24),
              axis.text.x = element_text(color = "Black", size = 24), 
              axis.text.y = element_text(color = "Black", size = 21), legend.text = element_text(size = 22),
              legend.background = element_blank(), panel.background = element_rect(fill = "black"),
              strip.text = element_text(size = 26), plot.margin = unit(c(0.5, 8, 1, 1), "lines"), 
              axis.title = element_blank(), strip.background = element_blank())
ggsave("Figures/MRS Covariate Association Discovery Females Heatmap.png", dpi = 600, width = 10, height = 12, units = "in")




