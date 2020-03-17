# DMR Pyrosequencing Replication ------------------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 3/12/20

# Setup -------------------------------------------------------------------
# Load Packages ####
.libPaths("/share/lasallelab/Charles/R")
sapply(c("tidyverse", "scales", "reshape2", "rlist", "bsseq"), require, character.only = TRUE)

# Load Functions ####
source("R Scripts/DMR Analysis Functions.R")

# Workflow ####
#' What is the rate of technical replication of DMRs identified with WGBS?
#' What proportion of WGBS DMRs can be validated by pyrosequencing?
#' 1. Get locations of DMRs assayed with pyrosequencing
#' 2. Overlap pyrosequenced regions with 4 current DMR comparisons
#' 3. For overlapping pyrosequenced regions, compare pyrosequencing methylation between ASD and TD
#' 4. Extract WGBS methylation just at pyrosequenced regions and compare pyrosequencing methylation between ASD and TD
#' 5. Compare pyrosequencing methylation with WGBS methylation at the same samples

# Overlap Pyrosequencing Regions with DMRs --------------------------------
# Load Regions ####
pyroseq <- read.delim("Tables/DMR Pyrosequencing Regions.txt", stringsAsFactors = FALSE)
pyroseq_GR <- makeGRange(pyroseq)
names(pyroseq_GR) <- pyroseq$gene
maleChroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")
femaleChroms <- c(paste("chr", 1:22, sep = ""), "chrX", "chrM")
DMRs <- list(Males_Discovery = loadRegions("DMRs/Discovery/Diagnosis Males 50/DMRs_Dx_Discovery50_males.csv", 
                                           chroms = maleChroms, DMRid = TRUE),
             Males_Replication = loadRegions("DMRs/Replication/Diagnosis Males 50/DMRs_Dx_Replication50_males.csv", 
                                             chroms = maleChroms, DMRid = TRUE),
             Females_Discovery = loadRegions("DMRs/Discovery/Diagnosis Females 50/DMRs_Dx_Discovery50_females.csv", 
                                             chroms = femaleChroms, DMRid = TRUE),
             Females_Replication = loadRegions("DMRs/Replication/Diagnosis Females 100/DMRs_Dx_Replication100_females.csv", 
                                               chroms = femaleChroms, DMRid = TRUE))
DMRs_GR <- lapply(DMRs, function(x){
        GR <- makeGRange(x)
        names(GR) <- x$DMRid
        return(GR)
})

# Overlap Regions ####
DMRs_sub <- lapply(DMRs_GR, subsetByOverlaps, ranges = pyroseq_GR, maxgap = 500) # Separated by up to 500 b, 3/4 directly overlap
sapply(DMRs_sub, length)
# Males_Discovery   Males_Replication   Females_Discovery Females_Replication 
#               0                   4                   0                   1 
# $Males_Replication
#  DMR_551     chr2 113235266-113236458      
# DMR_2589    chr10     1363157-1363742      
# DMR_4135    chr19   48498509-48498893     
# DMR_4210    chr20   30921928-30922301      
# $Females_Replication
# DMR_7555    chr19   48497534-48499296            

DMRstats_sub <- mapply(FUN = function(x,y){
        subset(x, DMRid %in% names(y), select = c("DMRid", "chr", "start", "end", "width", "L", "percentDifference", "stat", "pval"))
}, x = DMRs, y = DMRs_sub, SIMPLIFY = FALSE) 
# $Males_Replication
#    DMRid   chr     start       end width   L percentDifference     stat        pval
#  DMR_551  chr2 113235266 113236458  1193 102                 8 4.612009 0.002620299
# DMR_2589 chr10   1363157   1363742   586  81                 6 4.780476 0.001728179
# DMR_4135 chr19  48498509  48498893   385  43                 8 3.799209 0.018195845
# DMR_4210 chr20  30921928  30922301   374  41                 6 4.143268 0.008096115
# $Females_Replication
#    DMRid   chr    start      end width   L percentDifference    stat       pval
# DMR_7555 chr19 48497534 48499296  1763 111                11 3.41525 0.02849697

pyroseq_sub <- lapply(DMRs_GR, subsetByOverlaps, x = pyroseq_GR, maxgap = 500)
# $Males_Replication
#  ADARB2    chr10     1363559-1363667 Pyrosequenced to CpG 12, within DMR_2589     
# DEFB115    chr20   30921928-30922031 Pyrosequenced to CpG 10, within DMR_4210     
#   LMTK3    chr19   48499291-48499403 Pyrosequenced to CpG 5, 400 b downstream of DMR_4135      
#    PAX8     chr2 113235460-113235498 Pyrosequenced to CpG 4, within DMR_551     
# $Females_Replication
#   LMTK3    chr19   48499291-48499403 overlaps Females Rep DMR_7555 and continues 100 b downstream 
rm(femaleChroms, maleChroms)

# Compare Pyrosequencing Methylation by Diagnosis -------------------------
# Setup Data ####
samples <- read.csv("Tables/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo and Cell Type.csv", 
                    stringsAsFactors = FALSE)
samples <- samples[,c("Sequencing_ID", "Cord_Blood_IBC", "Platform", "Sex", "Diagnosis_Alg")] %>%
        subset(Platform == "HiSeq4000")
pyroseq_meth <- read.delim("Tables/DMR Pyrosequencing Data.txt", stringsAsFactors = FALSE)
pyroseq_meth <- merge(x = samples, y = pyroseq_meth, by = "Cord_Blood_IBC", all = FALSE)
pyroseq_meth$Diagnosis_Alg <- factor(pyroseq_meth$Diagnosis_Alg, levels = c("TD", "ASD"))
table(pyroseq_meth$Diagnosis_Alg)
# TD ASD 
# 10  12  All replication males

# Compare by Diagnosis ####
pyroseq_DMRs <- c("ADARB2", "DEFB115", "LMTK3", "PAX8")
pyroseq_stats <- sapply(pyroseq_DMRs, function(x){
        summary(lm(pyroseq_meth[,x] ~ pyroseq_meth$Diagnosis_Alg))$coefficients[-1,]
}) %>% t() %>% as.data.frame()
pyroseq_stats$qvalue <- p.adjust(pyroseq_stats$`Pr(>|t|)`, method = "fdr")
#         Estimate Std. Error   t value   Pr(>|t|)     qvalue
# ADARB2  4.628910   4.735684 0.9774534 0.34002071 0.45336094
# DEFB115 5.501964   2.720434 2.0224583 0.05670719 0.11341437
# LMTK3   7.241111   2.783114 2.6018022 0.01706031 0.06824124
# PAX8    5.104241   7.210869 0.7078538 0.48720402 0.48720402
write.table(pyroseq_stats, file = "Tables/Pyrosequencing Methylation by Diagnosis Stats.txt", sep = "\t", quote = FALSE,
            row.names = TRUE, col.names = TRUE)
rm(pyroseq, pyroseq_GR, pyroseq_sub)

# Compare WGBS Methylation by Diagnosis -----------------------------------
# Subset BSseq Object (Cluster) ####
BSmalesRep <- readRDS("Replication/Dx_Males/Filtered_BSseq_Replication50_males.rds")
pyroseq <- read.delim("DMR Pyrosequencing Regions.txt", stringsAsFactors = FALSE)
pyroseq <- subset(pyroseq, gene %in% c("ADARB2", "DEFB115", "LMTK3", "PAX8"))
pyroseq_GR <- with(pyroseq, GRanges(seqnames = chr, ranges = IRanges(start = start, end = end, names = gene)))
BSmalesRepSub <- subsetByOverlaps(BSmalesRep, ranges = pyroseq_GR)
saveRDS(BSmalesRepSub, "Pyrosequencing_Subsetted_BSseq_Replication_Males.rds")

# Get WGBS Methylation ####
options("showHeadLines" = 20)
options("showTailLines" = 20)
BSmalesRepSub <- readRDS("R Objects/Pyrosequencing_Subsetted_BSseq_Replication_Males.rds")
granges(BSmalesRepSub)
pyroseq <- read.delim("Tables/DMR Pyrosequencing Regions 2.txt", stringsAsFactors = FALSE) # New: just overlapping, reduced to only CpGs measured successfully
pyroseq_GR <- with(pyroseq, GRanges(seqnames = chr, ranges = IRanges(start = start, end = end, names = gene)))
wgbs_meth <- getMeth(BSmalesRepSub, regions = pyroseq_GR, type = "raw", what = "perRegion") %>% t() %>%  "*"(100) %>% 
        as.data.frame()
colnames(wgbs_meth) <- pyroseq_DMRs
wgbs_meth$Sequencing_ID <- rownames(wgbs_meth)
wgbs_meth <- merge(x = pyroseq_meth[,c("Cord_Blood_IBC", "Sequencing_ID", "Platform", "Sex", "Diagnosis_Alg")], y = wgbs_meth,
                   by = "Sequencing_ID", all = FALSE, sort = FALSE)
wgbs_meth$Diagnosis_Alg <- factor(wgbs_meth$Diagnosis_Alg, levels = c("TD", "ASD"))

# Compare by Diagnosis ####
wgbs_stats <- sapply(pyroseq_DMRs, function(x){
        summary(lm(wgbs_meth[,x] ~ wgbs_meth$Diagnosis_Alg))$coefficients[-1,]
}) %>% t() %>% as.data.frame()
wgbs_stats$qvalue <- p.adjust(wgbs_stats$`Pr(>|t|)`, method = "fdr")
#          Estimate Std. Error   t value   Pr(>|t|)     qvalue
# ADARB2   4.627891   8.133613 0.5689834 0.57570203 0.57570203
# DEFB115 14.542396   5.869350 2.4776844 0.02226819 0.08907276
# LMTK3    7.894444  10.156888 0.7772503 0.44610788 0.57570203
# PAX8     6.144040  10.434707 0.5888080 0.56258192 0.57570203

write.table(wgbs_stats, file = "Tables/WGBS Methylation at Pyroseq DMRs by Diagnosis Stats.txt", sep = "\t", quote = FALSE,
            row.names = TRUE, col.names = TRUE)

# Compare WGBS and Pyrosequencing -----------------------------------------
# Setup Data ####
pyroseq_vs_wgbs <- merge(x = reshape2::melt(pyroseq_meth[,c("Sequencing_ID", "Diagnosis_Alg", pyroseq_DMRs)], 
                                            id.vars = c("Sequencing_ID", "Diagnosis_Alg")),
                         y = reshape2::melt(wgbs_meth[,c("Sequencing_ID", "Diagnosis_Alg", pyroseq_DMRs)],
                                            id.vars = c("Sequencing_ID", "Diagnosis_Alg")), 
                         by = c("Sequencing_ID", "Diagnosis_Alg", "variable"), all = FALSE)
colnames(pyroseq_vs_wgbs) <- c("Sequencing_ID", "Diagnosis_Alg", "DMR", "Pyrosequencing", "WGBS")

# Correlation Stats ####
# Overall
cor.test(x = pyroseq_vs_wgbs$WGBS, y = pyroseq_vs_wgbs$Pyrosequencing) %>% .[c("estimate", "p.value")] %>% simplify2array()
# estimate.cor      p.value 
#    0.6015974 5.694511e-10 

# By DMR
pyroseq_vs_wgbs_stats <- split(pyroseq_vs_wgbs, f = pyroseq_vs_wgbs$DMR) %>% 
        lapply(function(x) cor.test(x = x$WGBS, y = x$Pyrosequencing) %>% .[c("estimate", "p.value")] %>% simplify2array()) %>%
        list.rbind() %>% as.data.frame()
pyroseq_vs_wgbs_stats$qvalue <- p.adjust(pyroseq_vs_wgbs_stats$p.value, method = "fdr")
#         estimate.cor      p.value       qvalue
# ADARB2     0.6564349 9.066369e-04 0.0018132737
# DEFB115    0.7568527 4.563578e-05 0.0001825431
# LMTK3      0.4259620 4.808243e-02 0.0480824338
# PAX8       0.5674379 5.881672e-03 0.0078422294

# Plot ####
pyroseq_vs_wgbs$DMR <- as.character(pyroseq_vs_wgbs$DMR) %>% 
        str_replace_all(pattern = c("ADARB2" = "DMR_2589", "DEFB115" = "DMR_4210", "LMTK3" = "DMR_4135", "PAX8" = "DMR_551")) %>%
        factor(levels = c("DMR_551", "DMR_2589", "DMR_4135", "DMR_4210"))

gg <- ggplot(data = pyroseq_vs_wgbs, aes(x = WGBS, y = Pyrosequencing)) +
        geom_smooth(method = "lm", formula = y ~ x) +
        geom_point(aes(color = Diagnosis_Alg), size = 3) +
        annotate("text", x = -Inf, y = Inf, label = "*", size = 14, hjust = -0.5, vjust = 1.3) +
        facet_wrap(vars(DMR), nrow = 2, scales = "free") +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'vertical', legend.position = c(1.08,0.98), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), 
              plot.margin = unit(c(1,7,1,1), "lines"), strip.background = element_blank()) +
        scale_x_continuous(breaks = breaks_pretty(n = 4)) +
        scale_y_continuous(breaks = breaks_pretty(n = 4)) +
        xlab("WGBS Methylation (%)") +
        ylab("Pyrosequencing Methylation (%)") +
        scale_color_manual(breaks = c("TD", "ASD"), values = c("TD" = "#3366CC", "ASD" = "#FF3366"))
ggsave("Figures/Pyrosequencing vs WGBS Methylation Scatterplot.pdf", plot = gg, dpi = 600, width = 11, height = 10, 
       units = "in")
