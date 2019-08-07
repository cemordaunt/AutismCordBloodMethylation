# Transcription Factor Motif Enrichment Analysis -----------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 8/5/19

# Packages ####
sapply(c("tidyverse", "GenomicRanges", "LOLA", "reshape2"), require, character.only = TRUE)

# Functions ####
source("R Scripts/DMR Analysis Functions.R")

# Make BED Input Files ----------------------------------------------------
# Load Regions ####
chromsXY <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")
chromsX <- c(paste("chr", 1:22, sep = ""), "chrX", "chrM")
DMRs <- list(MalesDisc = loadRegions("DMRs/Discovery/Diagnosis Males 50/DMRs_Dx_Discovery50_males.csv",
                                     chroms = chromsXY, sort = TRUE, DMRid = TRUE),
             FemalesDisc = loadRegions("DMRs/Discovery/Diagnosis Females 50/DMRs_Dx_Discovery50_females.csv",
                                       chroms = chromsX, sort = TRUE, DMRid = TRUE),
             MalesRep = loadRegions("DMRs/Replication/Diagnosis Males 50/DMRs_Dx_Replication50_males.csv",
                                    chroms = chromsXY, sort = TRUE, DMRid = TRUE),
             FemalesRep = loadRegions("DMRs/Replication/Diagnosis Females 100/DMRs_Dx_Replication100_females.csv",
                                      chroms = chromsX, sort = TRUE, DMRid = TRUE))
Background <- list(MalesDisc = loadRegions("DMRs/Discovery/Diagnosis Males 50/bsseq_background_Discovery50_males.csv",
                                           chroms = chromsXY, sort = TRUE),
                   FemalesDisc = loadRegions("DMRs/Discovery/Diagnosis Females 50/bsseq_background_Discovery50_females.csv",
                                             chroms = chromsX, sort = TRUE),
                   MalesRep = loadRegions("DMRs/Replication/Diagnosis Males 50/bsseq_background_Replication50_males.csv",
                                          chroms = chromsXY, sort = TRUE),
                   FemalesRep = loadRegions("DMRs/Replication/Diagnosis Females 100/bsseq_background_Replication100_females.csv",
                                            chroms = chromsX, sort = TRUE))

# Redefine DMRs in Terms of Background ####
GR_Background <- lapply(Background, makeGRange, direction = "all")

reDMRs <- mapply(function(x, y){
        makeGRange(x, direction = "all") %>% GRangesList() %>% redefineUserSets(userUniverse = y) %>% as.data.frame()
}, x = DMRs, y = GR_Background, SIMPLIFY = FALSE)

reHyperDMRs <- mapply(function(x, y){
        makeGRange(x, direction = "hyper") %>% GRangesList() %>% redefineUserSets(userUniverse = y) %>% as.data.frame()
}, x = DMRs, y = GR_Background, SIMPLIFY = FALSE)

reHypoDMRs <- mapply(function(x, y){
        makeGRange(x, direction = "hypo") %>% GRangesList() %>% redefineUserSets(userUniverse = y) %>% as.data.frame()
}, x = DMRs, y = GR_Background, SIMPLIFY = FALSE)

# Split Autosome and ChrX DMRs ####
Background_auto <- lapply(Background, subset, chr %in% paste("chr", 1:22, sep = ""))
reDMRs_auto <- lapply(reDMRs, subset, seqnames %in% paste("chr", 1:22, sep = ""))
reHyperDMRs_auto <- lapply(reHyperDMRs, subset, seqnames %in% paste("chr", 1:22, sep = ""))
reHypoDMRs_auto <- lapply(reHypoDMRs, subset, seqnames %in% paste("chr", 1:22, sep = ""))

Background_chrX <- lapply(Background, subset, chr == "chrX")
reDMRs_chrX <- lapply(reDMRs, subset, seqnames == "chrX")
reHyperDMRs_chrX <- lapply(reHyperDMRs, subset, seqnames == "chrX")
reHypoDMRs_chrX <- lapply(reHypoDMRs, subset, seqnames == "chrX")

# Write BED Files ####
# All DMRs
mapply(writeBED, regions = reDMRs, 
       file = paste("UCSC Tracks/", c("Discovery_Males", "Discovery_Females", "Replication_Males", "Replication_Females"),
                    "_Diagnosis_DMRs_HOMER.bed", sep = ""))
mapply(writeBED, regions = reHyperDMRs, 
       file = paste("UCSC Tracks/", c("Discovery_Males", "Discovery_Females", "Replication_Males", "Replication_Females"),
                    "_Diagnosis_Hyper_DMRs_HOMER.bed", sep = ""))
mapply(writeBED, regions = reHypoDMRs, 
       file = paste("UCSC Tracks/", c("Discovery_Males", "Discovery_Females", "Replication_Males", "Replication_Females"),
                    "_Diagnosis_Hypo_DMRs_HOMER.bed", sep = ""))
mapply(writeBED, regions = Background, 
       file = paste("UCSC Tracks/", c("Discovery_Males", "Discovery_Females", "Replication_Males", "Replication_Females"),
                    "_background.bed", sep = ""))

# Autosome DMRs
mapply(writeBED, regions = reDMRs_auto, 
       file = paste("UCSC Tracks/", c("Discovery_Males", "Discovery_Females", "Replication_Males", "Replication_Females"),
                    "_Diagnosis_DMRs_auto_HOMER.bed", sep = ""))
mapply(writeBED, regions = reHyperDMRs_auto, 
       file = paste("UCSC Tracks/", c("Discovery_Males", "Discovery_Females", "Replication_Males", "Replication_Females"),
                    "_Diagnosis_Hyper_DMRs_auto_HOMER.bed", sep = ""))
mapply(writeBED, regions = reHypoDMRs_auto, 
       file = paste("UCSC Tracks/", c("Discovery_Males", "Discovery_Females", "Replication_Males", "Replication_Females"),
                    "_Diagnosis_Hypo_DMRs_auto_HOMER.bed", sep = ""))
mapply(writeBED, regions = Background_auto, 
       file = paste("UCSC Tracks/", c("Discovery_Males", "Discovery_Females", "Replication_Males", "Replication_Females"),
                    "_background_auto.bed", sep = ""))

# ChrX DMRs
mapply(writeBED, regions = reDMRs_chrX, 
       file = paste("UCSC Tracks/", c("Discovery_Males", "Discovery_Females", "Replication_Males", "Replication_Females"),
                    "_Diagnosis_DMRs_chrX_HOMER.bed", sep = ""))
mapply(writeBED, regions = reHyperDMRs_chrX, 
       file = paste("UCSC Tracks/", c("Discovery_Males", "Discovery_Females", "Replication_Males", "Replication_Females"),
                    "_Diagnosis_Hyper_DMRs_chrX_HOMER.bed", sep = ""))
mapply(writeBED, regions = reHypoDMRs_chrX, 
       file = paste("UCSC Tracks/", c("Discovery_Males", "Discovery_Females", "Replication_Males", "Replication_Females"),
                    "_Diagnosis_Hypo_DMRs_chrX_HOMER.bed", sep = ""))
mapply(writeBED, regions = Background_chrX, 
       file = paste("UCSC Tracks/", c("Discovery_Males", "Discovery_Females", "Replication_Males", "Replication_Females"),
                    "_background_chrX.bed", sep = ""))

# Get HOMER Enrichments All Chroms ---------------------------------------------------
# Parameters ####
# HOMER Version v4.10
# See Shell Scripts/HOMER_TF_Motif_Enrichment.sh
# Parameters
# genome = hg38
# size = given (size of regions for motifs)
# cpg (Normalized for % CpG content)
# N = 100000 (Randomly pick 100K background regions for comparison, backgrounds are 200K - 400K)
# Queued on barbera 7/11

# Load Results ####
males_disc <- list(All = loadHOMER("DMRs/Discovery/Diagnosis Males 50/HOMER_All/knownResults.txt"), 
                   Hyper = loadHOMER("DMRs/Discovery/Diagnosis Males 50/HOMER_Hyper/knownResults.txt"),
                   Hypo = loadHOMER("DMRs/Discovery/Diagnosis Males 50/HOMER_Hypo/knownResults.txt"))
males_rep <- list(All = loadHOMER("DMRs/Replication/Diagnosis Males 50/HOMER_All/knownResults.txt"), 
                  Hyper = loadHOMER("DMRs/Replication/Diagnosis Males 50/HOMER_Hyper/knownResults.txt"),
                  Hypo = loadHOMER("DMRs/Replication/Diagnosis Males 50/HOMER_Hypo/knownResults.txt"))
females_disc <- list(All = loadHOMER("DMRs/Discovery/Diagnosis Females 50/HOMER_All/knownResults.txt"), 
                     Hyper = loadHOMER("DMRs/Discovery/Diagnosis Females 50/HOMER_Hyper/knownResults.txt"),
                     Hypo = loadHOMER("DMRs/Discovery/Diagnosis Females 50/HOMER_Hypo/knownResults.txt"))
females_rep <- list(All = loadHOMER("DMRs/Replication/Diagnosis Females 100/HOMER_All/knownResults.txt"), 
                    Hyper = loadHOMER("DMRs/Replication/Diagnosis Females 100/HOMER_Hyper/knownResults.txt"),
                    Hypo = loadHOMER("DMRs/Replication/Diagnosis Females 100/HOMER_Hypo/knownResults.txt"))

# Combine ####
id.vars <- colnames(males_disc$All)
homer <- rbind(melt(males_disc, id.vars = id.vars), melt(males_rep, id.vars = id.vars),
               melt(females_disc, id.vars = id.vars), melt(females_rep, id.vars = id.vars))
homer$DMRs <- paste(rep(c("Discovery", "Replication"), each = nrow(males_disc$All) * 3), 
                    rep(c("Males", "Females"), each = nrow(males_disc$All) * 6),
                    rep(c("All", "Hyper", "Hypo"), each = nrow(males_disc$All)))
homer$DMRs <- factor(homer$DMRs, levels = unique(homer$DMRs))
homer <- homer[,c("DMRs", "Transcription_Factor", "Family", "Consensus", "Target_with_Motif", "Target_Sequences",
                  "Percent_Target_with_Motif", "Background_with_Motif", "Background_Sequences", "Percent_Background_with_Motif",
                  "Fold_Enrichment", "pvalue", "log_pvalue", "qvalue", "log_qvalue")]
homer$log_pvalue[is.infinite(homer$log_pvalue)] <- max(homer$log_pvalue[!is.infinite(homer$log_pvalue)])
homer$log_qvalue[is.infinite(homer$log_qvalue)] <- max(homer$log_qvalue[!is.infinite(homer$log_qvalue)])

DMRs <- as.character(homer$DMRs) %>% unique()
homer$Enriched <- NULL
for(i in 1:length(DMRs)){ # Enriched motifs are in the top quartile of fold enrichment and log q-value for that comparison
        Fold_Enrichment <- homer$Fold_Enrichment[homer$DMRs == DMRs[i]]
        log_qvalue <- homer$log_qvalue[homer$DMRs == DMRs[i]]
        qvalue <- homer$qvalue[homer$DMRs == DMRs[i]]
        Enriched <- Fold_Enrichment >= quantile(Fold_Enrichment, 0.75) & 
                log_qvalue >= quantile(log_qvalue, 0.75) &
                qvalue < 0.05
        homer$Enriched[homer$DMRs == DMRs[i]] <- Enriched
}
table(homer$Enriched, homer$DMRs)["TRUE",]
#     Discovery Males All     Discovery Males Hyper      Discovery Males Hypo     Replication Males All 
#                      20                        27                        23                        11 
# Replication Males Hyper    Replication Males Hypo     Discovery Females All   Discovery Females Hyper 
#                      14                        15                        24                        25 
#  Discovery Females Hypo   Replication Females All Replication Females Hyper  Replication Females Hypo 
#                      22                        83                        29                        60 
write.table(homer, file = "Tables/HOMER Motif Enrichments Discovery and Replication DMRs.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)
rm(id.vars, males_disc, males_rep, females_disc, females_rep, i, Enriched, Fold_Enrichment, log_qvalue, qvalue)

# Motif Enrichment Comparison All Chroms ---------------------------------------------
# Hypothesis ####
# ASD DMRs are related to transcription factor-mediated regulation and may be upstream or downstream 
# of altered transcriptional regulation
# What are the top motifs in ASD DMRs?
# Do these differ between hyper and hypo DMRs?
# Do these differ between males and females?
# 413 unique TF motifs, 12 DMR sets

# Enrichment Summary Heatmaps ####
# Fold Enrichment, alphabetical order
homer$Transcription_Factor <- factor(homer$Transcription_Factor, levels = unique(homer$Transcription_Factor) %>% 
                                             sort(decreasing = TRUE))
gg <- ggplot(data = homer)
gg +
        geom_tile(aes(y = Transcription_Factor, x = DMRs, fill = Fold_Enrichment)) +
        scale_fill_gradientn("Fold\nEnrichment", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", 
                             limits = c(0, max(homer$Fold_Enrichment))) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              plot.margin = unit(c(1, 7, 1, 5), "lines"), axis.ticks.x = element_line(size = 1), 
              axis.ticks.y = element_blank(), panel.background = element_rect(fill = "black"),
              axis.text.x = element_text(size = 14, color = "Black", angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_blank(), axis.title = element_blank(), legend.key = element_blank(),  
              legend.position = c(1.125, 0.9), legend.background = element_blank(), 
              legend.key.size = unit(1, "lines"), legend.title = element_text(size = 16), 
              legend.text = element_text(size = 14))
ggsave("Figures/HOMER TF Enrichment Summary Fold Enrichment Heatmap.png", dpi = 600, width = 8, height = 10, units = "in")

# Fold Enrichment, clustered
homer$Transcription_Factor <- as.character(homer$Transcription_Factor)
homer <- homer[order(homer$Transcription_Factor),]
clusterData <- reshape::cast(homer[,c("Transcription_Factor", "DMRs", "Fold_Enrichment")], formula = Transcription_Factor ~ DMRs, 
                             fun.aggregate = mean, value = "Fold_Enrichment", add.missing = TRUE, fill = 0)
clusterOrder <- hclust(dist(clusterData[,2:ncol(clusterData)], method = "euclidean"), method = "ward.D")$order %>% rev()
homer$Transcription_Factor <- factor(homer$Transcription_Factor, levels = unique(homer$Transcription_Factor)[clusterOrder], ordered = TRUE)

gg <- ggplot(data = homer)
gg +
        geom_tile(aes(y = Transcription_Factor, x = DMRs, fill = Fold_Enrichment)) +
        scale_fill_gradientn("Fold\nEnrichment", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", 
                             limits = c(0, max(homer$Fold_Enrichment))) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              plot.margin = unit(c(1, 7, 1, 5), "lines"), axis.ticks.x = element_line(size = 1), 
              axis.ticks.y = element_blank(), panel.background = element_rect(fill = "black"),
              axis.text.x = element_text(size = 14, color = "Black", angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_blank(), axis.title = element_blank(), legend.key = element_blank(),  
              legend.position = c(1.125, 0.9), legend.background = element_blank(), 
              legend.key.size = unit(1, "lines"), legend.title = element_text(size = 16), 
              legend.text = element_text(size = 14))
ggsave("Figures/HOMER TF Enrichment Summary Fold Enrichment Heatmap clustered.png", dpi = 600, width = 8, height = 10, units = "in")

# log q-value, alphabetical order
homer$Transcription_Factor <- factor(homer$Transcription_Factor, levels = unique(homer$Transcription_Factor) %>% 
                                             sort(decreasing = TRUE))
gg <- ggplot(data = homer)
gg +
        geom_tile(aes(y = Transcription_Factor, x = DMRs, fill = log_qvalue)) +
        scale_fill_gradientn("-log(q-value)", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", 
                             limits = c(0, max(homer$log_qvalue))) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              plot.margin = unit(c(1, 7, 1, 5), "lines"), axis.ticks.x = element_line(size = 1), 
              axis.ticks.y = element_blank(), panel.background = element_rect(fill = "black"),
              axis.text.x = element_text(size = 14, color = "Black", angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_blank(), axis.title = element_blank(), legend.key = element_blank(),  
              legend.position = c(1.135, 0.915), legend.background = element_blank(), 
              legend.key.size = unit(1, "lines"), legend.title = element_text(size = 16), 
              legend.text = element_text(size = 14))
ggsave("Figures/HOMER TF Enrichment Summary log qvalue Heatmap.png", dpi = 600, width = 8, height = 10, units = "in")

# Top Enriched Motifs ####
# Enriched if in top quartile of log_qvalue and fold change within that data set
top <- lapply(DMRs, function(x) homer$Transcription_Factor[homer$DMRs == x & homer$Enriched])
names(top) <- DMRs

replicated <- list(Males_Hyper = intersect(top$'Discovery Males Hyper', top$'Replication Males Hyper'),
                   Males_Hypo = intersect(top$'Discovery Males Hypo', top$'Replication Males Hypo'),
                   Females_Hyper = intersect(top$'Discovery Females Hyper', top$'Replication Females Hyper'),
                   Females_Hypo = intersect(top$'Discovery Females Hypo', top$'Replication Females Hypo'))
# Males_Hyper "ETS"
# Males_Hypo "HIF-1A"   "IRF:BATF" "TBX20"    "VDR"      "ZNF675"  
# Females_Hyper "ARE"      "BACH2"    "CRE"      "GATA:SCL" "TCFCP2L1" "ZNF675"  
# Females_Hypo "EBF"      "FOXA1:AR" "HOXA2"    "IRF2"     "JUN-AP1"  "PKNOX1"   "TBX20"    "ZFP281"   "ZFP809"   
#              "ZNF519"   "ZNF675"  

intersect(replicated$Males_Hyper, replicated$Females_Hyper) # None
intersect(replicated$Males_Hypo, replicated$Females_Hypo) # "TBX20"  "ZNF675"
intersect(replicated$Males_Hyper, replicated$Males_Hypo) # None
intersect(replicated$Females_Hyper, replicated$Females_Hypo) # "ZNF675"

# Heatmap of Fold Enrichment, manual order by replicated DMR set
replicated_v <- unlist(replicated, use.names = FALSE) %>% unique()
homer_top <- subset(homer, Transcription_Factor %in% replicated_v & !grepl("All", homer$DMRs, fixed = TRUE))

homer_top$Transcription_Factor <- as.character(homer_top$Transcription_Factor) %>%
        factor(levels = rev(c("ZNF675", "TBX20", "ETS", "VDR", "HIF-1A", "IRF:BATF", "ARE", "GATA:SCL", "TCFCP2L1", "BACH2", "CRE", 
                              "HOXA2", "PKNOX1", "IRF2", "JUN-AP1", "EBF", "ZFP809", "ZNF519", "FOXA1:AR", "ZFP281")), 
               ordered = TRUE)
homer_top$Enriched <- factor(homer_top$Enriched, levels = c("TRUE", "FALSE"))
homer_top$DMRs <- as.character(homer_top$DMRs)
homer_top$Sex <- ifelse(grepl("Males", homer_top$DMRs, fixed = TRUE), yes = "Males", no = "Females") %>%
        factor(levels = c("Males", "Females"))
homer_top$DMRs <- str_replace_all(homer_top$DMRs, pattern = c("Males " = "", "Females " = "")) %>%
        factor(levels = c("Discovery Hyper", "Replication Hyper", "Discovery Hypo", "Replication Hypo"))
gg <- ggplot(data = homer_top)
gg +
        geom_tile(aes(y = Transcription_Factor, x = DMRs, fill = Fold_Enrichment)) +
        geom_point(aes(y = Transcription_Factor, x = DMRs, alpha = Enriched), color = "white", size = 2.5) +
        facet_grid(cols = vars(Sex)) +
        scale_fill_gradientn("Fold\nEnrichment", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", 
                             limits = c(0, max(homer_top$Fold_Enrichment))) +
        scale_alpha_manual(breaks = c("TRUE", "FALSE"), values = c(1, 0), guide = FALSE) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              plot.margin = unit(c(0.25, 8, 0.5, 1), "lines"), axis.ticks = element_line(size = 1), 
              panel.background = element_rect(fill = "black"),
              axis.text.x = element_text(size = 15, color = "black", angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_text(size = 14, color = "black"), axis.title = element_blank(), 
              legend.key = element_blank(), legend.position = c(1.19, 0.86), legend.background = element_blank(), 
              legend.key.size = unit(1, "lines"), legend.title = element_text(size = 18), 
              legend.text = element_text(size = 16), strip.background = element_blank(), 
              strip.text = element_text(size = 20))
ggsave("Figures/HOMER TF Enrichment Top Fold Enrichment Heatmap manual order.png", dpi = 600, width = 7, height = 8, units = "in")

# Get HOMER Enrichments Auto ---------------------------------------------------
# Load Results ####
males_disc <- list(All = loadHOMER("DMRs/Discovery/Diagnosis Males 50/HOMER_Auto_All/knownResults.txt"), 
                   Hyper = loadHOMER("DMRs/Discovery/Diagnosis Males 50/HOMER_Auto_Hyper/knownResults.txt"),
                   Hypo = loadHOMER("DMRs/Discovery/Diagnosis Males 50/HOMER_Auto_Hypo/knownResults.txt"))
males_rep <- list(All = loadHOMER("DMRs/Replication/Diagnosis Males 50/HOMER_Auto_All/knownResults.txt"), 
                  Hyper = loadHOMER("DMRs/Replication/Diagnosis Males 50/HOMER_Auto_Hyper/knownResults.txt"),
                  Hypo = loadHOMER("DMRs/Replication/Diagnosis Males 50/HOMER_Auto_Hypo/knownResults.txt"))
females_disc <- list(All = loadHOMER("DMRs/Discovery/Diagnosis Females 50/HOMER_Auto_All/knownResults.txt"), 
                     Hyper = loadHOMER("DMRs/Discovery/Diagnosis Females 50/HOMER_Auto_Hyper/knownResults.txt"),
                     Hypo = loadHOMER("DMRs/Discovery/Diagnosis Females 50/HOMER_Auto_Hypo/knownResults.txt"))
females_rep <- list(All = loadHOMER("DMRs/Replication/Diagnosis Females 100/HOMER_Auto_All/knownResults.txt"), 
                    Hyper = loadHOMER("DMRs/Replication/Diagnosis Females 100/HOMER_Auto_Hyper/knownResults.txt"),
                    Hypo = loadHOMER("DMRs/Replication/Diagnosis Females 100/HOMER_Auto_Hypo/knownResults.txt"))

# Combine ####
id.vars <- colnames(males_disc$All)
homer <- rbind(melt(males_disc, id.vars = id.vars), melt(males_rep, id.vars = id.vars),
               melt(females_disc, id.vars = id.vars), melt(females_rep, id.vars = id.vars))
homer$DMRs <- paste(rep(c("Discovery", "Replication"), each = nrow(males_disc$All) * 3), 
                    rep(c("Males", "Females"), each = nrow(males_disc$All) * 6),
                    rep(c("All", "Hyper", "Hypo"), each = nrow(males_disc$All)))
homer$DMRs <- factor(homer$DMRs, levels = unique(homer$DMRs))
homer <- homer[,c("DMRs", "Transcription_Factor", "Family", "Consensus", "Target_with_Motif", "Target_Sequences",
                  "Percent_Target_with_Motif", "Background_with_Motif", "Background_Sequences", "Percent_Background_with_Motif",
                  "Fold_Enrichment", "pvalue", "log_pvalue", "qvalue", "log_qvalue")]
homer$log_pvalue[is.infinite(homer$log_pvalue)] <- max(homer$log_pvalue[!is.infinite(homer$log_pvalue)])
homer$log_qvalue[is.infinite(homer$log_qvalue)] <- max(homer$log_qvalue[!is.infinite(homer$log_qvalue)])

DMRs <- as.character(homer$DMRs) %>% unique()
homer$Enriched <- NULL
for(i in 1:length(DMRs)){ # Enriched motifs are in the top quartile of fold enrichment and log q-value for that comparison
        Fold_Enrichment <- homer$Fold_Enrichment[homer$DMRs == DMRs[i]]
        log_qvalue <- homer$log_qvalue[homer$DMRs == DMRs[i]]
        qvalue <- homer$qvalue[homer$DMRs == DMRs[i]]
        Enriched <- Fold_Enrichment >= quantile(Fold_Enrichment, 0.75) & 
                log_qvalue >= quantile(log_qvalue, 0.75) &
                qvalue < 0.05
        homer$Enriched[homer$DMRs == DMRs[i]] <- Enriched
}
table(homer$Enriched, homer$DMRs)["TRUE",]
#     Discovery Males All     Discovery Males Hyper      Discovery Males Hypo     Replication Males All 
#                      26                        34                        32                        15 
# Replication Males Hyper    Replication Males Hypo     Discovery Females All   Discovery Females Hyper 
#                      17                        17                        24                        26 
#  Discovery Females Hypo   Replication Females All Replication Females Hyper  Replication Females Hypo 
#                      27                        82                        33                        56 
write.table(homer, file = "Tables/HOMER Motif Enrichments Discovery and Replication Auto DMRs.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)
rm(id.vars, males_disc, males_rep, females_disc, females_rep, i, Enriched, Fold_Enrichment, log_qvalue, qvalue)

# Motif Enrichment Comparison Auto ---------------------------------------------
# Enrichment Summary Heatmaps ####
# Fold Enrichment, alphabetical order
homer$Transcription_Factor <- factor(homer$Transcription_Factor, levels = unique(homer$Transcription_Factor) %>% 
                                             sort(decreasing = TRUE))
gg <- ggplot(data = homer)
gg +
        geom_tile(aes(y = Transcription_Factor, x = DMRs, fill = Fold_Enrichment)) +
        scale_fill_gradientn("Fold\nEnrichment", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", 
                             limits = c(0, max(homer$Fold_Enrichment))) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              plot.margin = unit(c(1, 7, 1, 5), "lines"), axis.ticks.x = element_line(size = 1), 
              axis.ticks.y = element_blank(), panel.background = element_rect(fill = "black"),
              axis.text.x = element_text(size = 14, color = "Black", angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_blank(), axis.title = element_blank(), legend.key = element_blank(),  
              legend.position = c(1.125, 0.9), legend.background = element_blank(), 
              legend.key.size = unit(1, "lines"), legend.title = element_text(size = 16), 
              legend.text = element_text(size = 14))
ggsave("Figures/HOMER TF Enrichment Summary Fold Enrichment Heatmap Auto.png", dpi = 600, width = 8, height = 10, units = "in")

# Fold Enrichment, clustered
homer$Transcription_Factor <- as.character(homer$Transcription_Factor)
homer <- homer[order(homer$Transcription_Factor),]
clusterData <- reshape::cast(homer[,c("Transcription_Factor", "DMRs", "Fold_Enrichment")], formula = Transcription_Factor ~ DMRs, 
                             fun.aggregate = mean, value = "Fold_Enrichment", add.missing = TRUE, fill = 0)
clusterOrder <- hclust(dist(clusterData[,2:ncol(clusterData)], method = "euclidean"), method = "ward.D")$order %>% rev()
homer$Transcription_Factor <- factor(homer$Transcription_Factor, levels = unique(homer$Transcription_Factor)[clusterOrder], ordered = TRUE)

gg <- ggplot(data = homer)
gg +
        geom_tile(aes(y = Transcription_Factor, x = DMRs, fill = Fold_Enrichment)) +
        scale_fill_gradientn("Fold\nEnrichment", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", 
                             limits = c(0, max(homer$Fold_Enrichment))) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              plot.margin = unit(c(1, 7, 1, 5), "lines"), axis.ticks.x = element_line(size = 1), 
              axis.ticks.y = element_blank(), panel.background = element_rect(fill = "black"),
              axis.text.x = element_text(size = 14, color = "Black", angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_blank(), axis.title = element_blank(), legend.key = element_blank(),  
              legend.position = c(1.125, 0.9), legend.background = element_blank(), 
              legend.key.size = unit(1, "lines"), legend.title = element_text(size = 16), 
              legend.text = element_text(size = 14))
ggsave("Figures/HOMER TF Enrichment Summary Fold Enrichment Heatmap clustered Auto.png", dpi = 600, width = 8, height = 10, units = "in")

# log q-value, alphabetical order
homer$Transcription_Factor <- factor(homer$Transcription_Factor, levels = unique(homer$Transcription_Factor) %>% 
                                             sort(decreasing = TRUE))
gg <- ggplot(data = homer)
gg +
        geom_tile(aes(y = Transcription_Factor, x = DMRs, fill = log_qvalue)) +
        scale_fill_gradientn("-log(q-value)", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", 
                             limits = c(0, max(homer$log_qvalue))) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              plot.margin = unit(c(1, 7, 1, 5), "lines"), axis.ticks.x = element_line(size = 1), 
              axis.ticks.y = element_blank(), panel.background = element_rect(fill = "black"),
              axis.text.x = element_text(size = 14, color = "Black", angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_blank(), axis.title = element_blank(), legend.key = element_blank(),  
              legend.position = c(1.135, 0.915), legend.background = element_blank(), 
              legend.key.size = unit(1, "lines"), legend.title = element_text(size = 16), 
              legend.text = element_text(size = 14))
ggsave("Figures/HOMER TF Enrichment Summary log qvalue Heatmap Auto.png", dpi = 600, width = 8, height = 10, units = "in")

# Top Enriched Motifs ####
# Enriched if in top quartile of log_qvalue and fold change within that data set
top <- lapply(DMRs, function(x) homer$Transcription_Factor[homer$DMRs == x & homer$Enriched])
names(top) <- DMRs

replicated <- list(Males_Hyper = intersect(top$'Discovery Males Hyper', top$'Replication Males Hyper'),
                   Males_Hypo = intersect(top$'Discovery Males Hypo', top$'Replication Males Hypo'),
                   Females_Hyper = intersect(top$'Discovery Females Hyper', top$'Replication Females Hyper'),
                   Females_Hypo = intersect(top$'Discovery Females Hypo', top$'Replication Females Hypo'))
# Males_Hyper "ETS"   "TBX20" "VDR"  
# Males_Hypo "ARE"      "BRN2"     "HIF-1A"   "HOXA2"    "IRF:BATF" "RFX1"     "RORGT"    "TBX20"    "VDR"      "ZNF675"  
# Females_Hyper "ARE"      "CRE"      "GATA:SCL" "REVERB"   "RFX1"     "ZNF519"   "ZNF675"  
# Females_Hypo "BACH2"     "CRE"       "EBF"       "ETS:E-BOX" "HOXA2"     "JUN-AP1"   "PAX5_1"    "PAX6"      "PKNOX1"    
#              "TBX20"    "ZFP809"    "ZNF519"    "ZNF675"   

intersect(replicated$Males_Hyper, replicated$Females_Hyper) # None
intersect(replicated$Males_Hypo, replicated$Females_Hypo) # "HOXA2"  "TBX20"  "ZNF675"
intersect(replicated$Males_Hyper, replicated$Males_Hypo) # "TBX20" "VDR"  
intersect(replicated$Females_Hyper, replicated$Females_Hypo) # "CRE"    "ZNF519" "ZNF675"

# Heatmap of Fold Enrichment, manual order by replicated DMR set
replicated_v <- unlist(replicated, use.names = FALSE) %>% unique()
homer_top <- subset(homer, Transcription_Factor %in% replicated_v & !grepl("All", homer$DMRs, fixed = TRUE))

replicated_v <- str_replace_all(replicated_v, pattern = c("_1" = ""))
homer_top$Transcription_Factor <- as.character(homer_top$Transcription_Factor) %>% 
        str_replace_all(pattern = c("_1" = "")) %>%
        factor(levels = rev(c("TBX20", "ZNF675", "HOXA2", "VDR", "ETS", "ARE", "RFX1", "HIF-1A", "IRF:BATF", "RORGT", "BRN2",    
                              "CRE", "ZNF519", "REVERB", "GATA:SCL", "PKNOX1", "BACH2", "PAX5", "JUN-AP1", "EBF", "ZFP809",
                              "PAX6", "ETS:E-BOX")), ordered = TRUE)
homer_top$Enriched <- factor(homer_top$Enriched, levels = c("TRUE", "FALSE"))
homer_top$DMRs <- as.character(homer_top$DMRs)
homer_top$Sex <- ifelse(grepl("Males", homer_top$DMRs, fixed = TRUE), yes = "Males", no = "Females") %>%
        factor(levels = c("Males", "Females"))
homer_top$DMRs <- str_replace_all(homer_top$DMRs, pattern = c("Males " = "", "Females " = "")) %>%
        factor(levels = c("Discovery Hyper", "Replication Hyper", "Discovery Hypo", "Replication Hypo"))
gg <- ggplot(data = homer_top)
gg +
        geom_tile(aes(y = Transcription_Factor, x = DMRs, fill = Fold_Enrichment)) +
        geom_point(aes(y = Transcription_Factor, x = DMRs, alpha = Enriched), color = "white", size = 2.5) +
        facet_grid(cols = vars(Sex)) +
        scale_fill_gradientn("Fold\nEnrichment", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", 
                             limits = c(0, max(homer_top$Fold_Enrichment))) +
        scale_alpha_manual(breaks = c("TRUE", "FALSE"), values = c(1, 0), guide = FALSE) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              plot.margin = unit(c(0.25, 8, 0.5, 1), "lines"), axis.ticks = element_line(size = 1), 
              panel.background = element_rect(fill = "black"),
              axis.text.x = element_text(size = 15, color = "black", angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_text(size = 14, color = "black"), axis.title = element_blank(), 
              legend.key = element_blank(), legend.position = c(1.19, 0.865), legend.background = element_blank(), 
              legend.key.size = unit(1, "lines"), legend.title = element_text(size = 18), 
              legend.text = element_text(size = 16), strip.background = element_blank(), 
              strip.text = element_text(size = 20))
ggsave("Figures/HOMER TF Enrichment Top Fold Enrichment Heatmap manual order Auto.png", dpi = 600, width = 7, height = 8, units = "in")

# Get HOMER Enrichments ChrX ---------------------------------------------------
# Load Results ####
males_disc <- list(All = loadHOMER("DMRs/Discovery/Diagnosis Males 50/HOMER_ChrX_All/knownResults.txt"), 
                   Hyper = loadHOMER("DMRs/Discovery/Diagnosis Males 50/HOMER_ChrX_Hyper/knownResults.txt"),
                   Hypo = loadHOMER("DMRs/Discovery/Diagnosis Males 50/HOMER_ChrX_Hypo/knownResults.txt"))
males_rep <- list(All = loadHOMER("DMRs/Replication/Diagnosis Males 50/HOMER_ChrX_All/knownResults.txt"), 
                  Hyper = loadHOMER("DMRs/Replication/Diagnosis Males 50/HOMER_ChrX_Hyper/knownResults.txt"),
                  Hypo = loadHOMER("DMRs/Replication/Diagnosis Males 50/HOMER_ChrX_Hypo/knownResults.txt"))
females_disc <- list(All = loadHOMER("DMRs/Discovery/Diagnosis Females 50/HOMER_ChrX_All/knownResults.txt"), 
                     Hyper = loadHOMER("DMRs/Discovery/Diagnosis Females 50/HOMER_ChrX_Hyper/knownResults.txt"),
                     Hypo = loadHOMER("DMRs/Discovery/Diagnosis Females 50/HOMER_ChrX_Hypo/knownResults.txt"))
females_rep <- list(All = loadHOMER("DMRs/Replication/Diagnosis Females 100/HOMER_ChrX_All/knownResults.txt"), 
                    Hyper = loadHOMER("DMRs/Replication/Diagnosis Females 100/HOMER_ChrX_Hyper/knownResults.txt"),
                    Hypo = loadHOMER("DMRs/Replication/Diagnosis Females 100/HOMER_ChrX_Hypo/knownResults.txt"))

# Combine ####
id.vars <- colnames(males_disc$All)
homer <- rbind(melt(males_disc, id.vars = id.vars), melt(males_rep, id.vars = id.vars),
               melt(females_disc, id.vars = id.vars), melt(females_rep, id.vars = id.vars))
homer$DMRs <- paste(rep(c("Discovery", "Replication"), each = nrow(males_disc$All) * 3), 
                    rep(c("Males", "Females"), each = nrow(males_disc$All) * 6),
                    rep(c("All", "Hyper", "Hypo"), each = nrow(males_disc$All)))
homer$DMRs <- factor(homer$DMRs, levels = unique(homer$DMRs))
homer <- homer[,c("DMRs", "Transcription_Factor", "Family", "Consensus", "Target_with_Motif", "Target_Sequences",
                  "Percent_Target_with_Motif", "Background_with_Motif", "Background_Sequences", "Percent_Background_with_Motif",
                  "Fold_Enrichment", "pvalue", "log_pvalue", "qvalue", "log_qvalue")]
homer$log_pvalue[is.infinite(homer$log_pvalue)] <- max(homer$log_pvalue[!is.infinite(homer$log_pvalue)])
homer$log_qvalue[is.infinite(homer$log_qvalue)] <- max(homer$log_qvalue[!is.infinite(homer$log_qvalue)])

DMRs <- as.character(homer$DMRs) %>% unique()
homer$Enriched <- NULL
for(i in 1:length(DMRs)){ # Enriched motifs are in the top quartile of fold enrichment and log q-value for that comparison
        Fold_Enrichment <- homer$Fold_Enrichment[homer$DMRs == DMRs[i]]
        log_qvalue <- homer$log_qvalue[homer$DMRs == DMRs[i]]
        qvalue <- homer$qvalue[homer$DMRs == DMRs[i]]
        Enriched <- Fold_Enrichment >= quantile(Fold_Enrichment, 0.75) & 
                log_qvalue >= quantile(log_qvalue, 0.75) &
                qvalue < 0.05
        homer$Enriched[homer$DMRs == DMRs[i]] <- Enriched
}
table(homer$Enriched, homer$DMRs)["TRUE",]
#     Discovery Males All     Discovery Males Hyper      Discovery Males Hypo     Replication Males All 
#                      41                         0                        56                        14 
# Replication Males Hyper    Replication Males Hypo     Discovery Females All   Discovery Females Hyper 
#                      57                        17                        31                        69 
#  Discovery Females Hypo   Replication Females All Replication Females Hyper  Replication Females Hypo 
#                      34                        23                        25                        27 
write.table(homer, file = "Tables/HOMER Motif Enrichments Discovery and Replication ChrX DMRs.txt", sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)
rm(id.vars, males_disc, males_rep, females_disc, females_rep, i, Enriched, Fold_Enrichment, log_qvalue, qvalue)

# Motif Enrichment Comparison ChrX ---------------------------------------------
# Enrichment Summary Heatmaps ####
# Fold Enrichment, alphabetical order
homer$Transcription_Factor <- factor(homer$Transcription_Factor, levels = unique(homer$Transcription_Factor) %>% 
                                             sort(decreasing = TRUE))
gg <- ggplot(data = homer)
gg +
        geom_tile(aes(y = Transcription_Factor, x = DMRs, fill = Fold_Enrichment)) +
        scale_fill_gradientn("Fold\nEnrichment", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", 
                             limits = c(0, max(homer$Fold_Enrichment))) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              plot.margin = unit(c(1, 7, 1, 5), "lines"), axis.ticks.x = element_line(size = 1), 
              axis.ticks.y = element_blank(), panel.background = element_rect(fill = "black"),
              axis.text.x = element_text(size = 14, color = "Black", angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_blank(), axis.title = element_blank(), legend.key = element_blank(),  
              legend.position = c(1.125, 0.9), legend.background = element_blank(), 
              legend.key.size = unit(1, "lines"), legend.title = element_text(size = 16), 
              legend.text = element_text(size = 14))
ggsave("Figures/HOMER TF Enrichment Summary Fold Enrichment Heatmap ChrX.png", dpi = 600, width = 8, height = 10, units = "in")

# Fold Enrichment, clustered
homer$Transcription_Factor <- as.character(homer$Transcription_Factor)
homer <- homer[order(homer$Transcription_Factor),]
clusterData <- reshape::cast(homer[,c("Transcription_Factor", "DMRs", "Fold_Enrichment")], formula = Transcription_Factor ~ DMRs, 
                             fun.aggregate = mean, value = "Fold_Enrichment", add.missing = TRUE, fill = 0)
clusterOrder <- hclust(dist(clusterData[,2:ncol(clusterData)], method = "euclidean"), method = "ward.D")$order %>% rev()
homer$Transcription_Factor <- factor(homer$Transcription_Factor, levels = unique(homer$Transcription_Factor)[clusterOrder], ordered = TRUE)

gg <- ggplot(data = homer)
gg +
        geom_tile(aes(y = Transcription_Factor, x = DMRs, fill = Fold_Enrichment)) +
        scale_fill_gradientn("Fold\nEnrichment", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", 
                             limits = c(0, max(homer$Fold_Enrichment))) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              plot.margin = unit(c(1, 7, 1, 5), "lines"), axis.ticks.x = element_line(size = 1), 
              axis.ticks.y = element_blank(), panel.background = element_rect(fill = "black"),
              axis.text.x = element_text(size = 14, color = "Black", angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_blank(), axis.title = element_blank(), legend.key = element_blank(),  
              legend.position = c(1.125, 0.9), legend.background = element_blank(), 
              legend.key.size = unit(1, "lines"), legend.title = element_text(size = 16), 
              legend.text = element_text(size = 14))
ggsave("Figures/HOMER TF Enrichment Summary Fold Enrichment Heatmap clustered ChrX.png", dpi = 600, width = 8, height = 10, units = "in")

# log q-value, alphabetical order
homer$Transcription_Factor <- factor(homer$Transcription_Factor, levels = unique(homer$Transcription_Factor) %>% 
                                             sort(decreasing = TRUE))
gg <- ggplot(data = homer)
gg +
        geom_tile(aes(y = Transcription_Factor, x = DMRs, fill = log_qvalue)) +
        scale_fill_gradientn("-log(q-value)", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", 
                             limits = c(0, max(homer$log_qvalue))) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              plot.margin = unit(c(1, 7, 1, 5), "lines"), axis.ticks.x = element_line(size = 1), 
              axis.ticks.y = element_blank(), panel.background = element_rect(fill = "black"),
              axis.text.x = element_text(size = 14, color = "Black", angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_blank(), axis.title = element_blank(), legend.key = element_blank(),  
              legend.position = c(1.135, 0.915), legend.background = element_blank(), 
              legend.key.size = unit(1, "lines"), legend.title = element_text(size = 16), 
              legend.text = element_text(size = 14))
ggsave("Figures/HOMER TF Enrichment Summary log qvalue Heatmap ChrX.png", dpi = 600, width = 8, height = 10, units = "in")

# Top Enriched Motifs ####
# Enriched if in top quartile of log_qvalue and fold change within that data set
top <- lapply(DMRs, function(x) homer$Transcription_Factor[homer$DMRs == x & homer$Enriched])
names(top) <- DMRs

replicated <- list(Males_Hyper = intersect(top$'Discovery Males Hyper', top$'Replication Males Hyper'),
                   Males_Hypo = intersect(top$'Discovery Males Hypo', top$'Replication Males Hypo'),
                   Females_Hyper = intersect(top$'Discovery Females Hyper', top$'Replication Females Hyper'),
                   Females_Hypo = intersect(top$'Discovery Females Hypo', top$'Replication Females Hypo'))
# Males_Hyper None  
# Males_Hypo "GATA3_1" "GRE_2"   "KLF4"    "PAX5_2"  "PGR"    
# Females_Hyper "ARE"      "BACH2"    "CRE"      "E2F4"     "EBF"      "GATA:SCL" "HINFP"    "JUN-AP1"  "MAFK"     
#              "PBX3"     "ZNF322"
# Females_Hypo "E2F1"      "E2F4"      "E2F6"      "PBX3"      "PU.1:IRF8" "RAR:RXR_1" "RORA"      "TBX20"       

intersect(replicated$Males_Hyper, replicated$Females_Hyper) # None
intersect(replicated$Males_Hypo, replicated$Females_Hypo) # None
intersect(replicated$Males_Hyper, replicated$Males_Hypo) # None
intersect(replicated$Females_Hyper, replicated$Females_Hypo) # "E2F4" "PBX3"

# Heatmap of Fold Enrichment, manual order by replicated DMR set
replicated_v <- unlist(replicated, use.names = FALSE) %>% unique()
homer_top <- subset(homer, Transcription_Factor %in% replicated_v & !grepl("All", homer$DMRs, fixed = TRUE))

replicated_v <- str_replace_all(replicated_v, pattern = c("_1" = "", "_2" = ""))
homer_top$Transcription_Factor <- as.character(homer_top$Transcription_Factor) %>% 
        str_replace_all(pattern = c("_1" = "", "_2" = "")) %>%
        factor(levels = rev(c("GATA3", "PAX5", "GRE", "KLF4", "PGR", "PBX3", "E2F4", "JUN-AP1", "ARE", "HINFP", "ZNF322", 
                              "CRE", "EBF", "BACH2", "MAFK", "GATA:SCL", "RORA", "E2F1", "E2F6", "TBX20", "PU.1:IRF8",
                              "RAR:RXR")), ordered = TRUE)
homer_top$Enriched <- factor(homer_top$Enriched, levels = c("TRUE", "FALSE"))
homer_top$DMRs <- as.character(homer_top$DMRs)
homer_top$Sex <- ifelse(grepl("Males", homer_top$DMRs, fixed = TRUE), yes = "Males", no = "Females") %>%
        factor(levels = c("Males", "Females"))
homer_top$DMRs <- str_replace_all(homer_top$DMRs, pattern = c("Males " = "", "Females " = "")) %>%
        factor(levels = c("Discovery Hyper", "Replication Hyper", "Discovery Hypo", "Replication Hypo"))
gg <- ggplot(data = homer_top)
gg +
        geom_tile(aes(y = Transcription_Factor, x = DMRs, fill = Fold_Enrichment)) +
        geom_point(aes(y = Transcription_Factor, x = DMRs, alpha = Enriched), color = "white", size = 2.5) +
        facet_grid(cols = vars(Sex)) +
        scale_fill_gradientn("Fold\nEnrichment", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", 
                             limits = c(0, max(homer_top$Fold_Enrichment))) +
        scale_alpha_manual(breaks = c("TRUE", "FALSE"), values = c(1, 0), guide = FALSE) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              plot.margin = unit(c(0.25, 8, 0.5, 1), "lines"), axis.ticks = element_line(size = 1), 
              panel.background = element_rect(fill = "black"),
              axis.text.x = element_text(size = 15, color = "black", angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_text(size = 14, color = "black"), axis.title = element_blank(), 
              legend.key = element_blank(), legend.position = c(1.19, 0.865), legend.background = element_blank(), 
              legend.key.size = unit(1, "lines"), legend.title = element_text(size = 18), 
              legend.text = element_text(size = 16), strip.background = element_blank(), 
              strip.text = element_text(size = 20))
ggsave("Figures/HOMER TF Enrichment Top Fold Enrichment Heatmap manual order ChrX.png", dpi = 600, width = 7, height = 8, units = "in")
