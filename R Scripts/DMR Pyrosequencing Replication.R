# DMR Pyrosequencing Replication ------------------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 3/12/20

# Setup -------------------------------------------------------------------
# Load Packages ####
.libPaths("/share/lasallelab/Charles/R")
sapply(c("tidyverse", "scales", "reshape2", "bsseq"), require, character.only = TRUE)

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
#  ADARB2    chr10     1363559-1363667 Pyrosequenced to CpG 12     
# DEFB115    chr20   30921928-30922031 Pyrosequenced to CpG 10     
#   LMTK3    chr19   48499291-48499403 Pyrosequenced to CpG 5, 400 b downstream of DMR      
#    PAX8     chr2 113235460-113235498 Pyrosequenced to CpG 4     
# $Females_Replication
#   LMTK3    chr19   48499291-48499403 overlaps directly     
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
pyroseq_stats <- sapply(c("ADARB2", "DEFB115", "LMTK3", "PAX8"), function(x){
        summary(lm(pyroseq_meth[,x] ~ pyroseq_meth$Diagnosis_Alg))$coefficients[-1,]
}) %>% t() %>% as.data.frame()
#         Estimate Std. Error   t value   Pr(>|t|)
# ADARB2  4.628910   4.735684 0.9774534 0.34002071
# DEFB115 5.501964   2.720434 2.0224583 0.05670719
# LMTK3   7.241111   2.783114 2.6018022 0.01706031
# PAX8    5.210194   7.190782 0.7245658 0.47711014
