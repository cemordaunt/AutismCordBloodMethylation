# Autism Gene Overlap -----------------------------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 8/14/19

# Packages ####
sapply(c("reshape2", "tidyverse", "GenomicRanges", "annotatr", "GeneOverlap", "rlist", "biomaRt", "rtracklayer",
         "scales"), require, character.only = TRUE)

# Functions ####
source("R Scripts/DMR Analysis Functions.R")
regDomains <- read.delim("Tables/Regulatory domains hg38.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Get DMR EntrezIDs -----------------------------------------------------------
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

# Get Entrez Gene IDs ####
DMRs_entrez <- lapply(DMRs, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_entrezID")
Background_entrez <- lapply(Background, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_entrezID")
chrX_DMRs_entrez <- lapply(DMRs, subset, chr == "chrX") %>% lapply(getDMRgeneList, regDomains = regDomains, direction = "all", 
                                                                   type = "gene_entrezID")

# Write Files
names(DMRs_entrez) <- c("Males_Discovery", "Females_Discovery", "Males_Replication", "Females_Replication")
mapply(write.table, x = DMRs_entrez, file = paste("Entrez ID Lists/EntrezIDs_", names(DMRs_entrez), "_DMRs.txt", sep = ""), 
       MoreArgs = list(sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE))

names(Background_entrez) <- c("Males_Discovery", "Females_Discovery", "Males_Replication", "Females_Replication")
mapply(write.table, x = Background_entrez, file = paste("Entrez ID Lists/EntrezIDs_", names(Background_entrez), "_Background.txt", sep = ""), 
       MoreArgs = list(sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE))

names(chrX_DMRs_entrez) <- c("Males_Discovery", "Females_Discovery", "Males_Replication", "Females_Replication")
mapply(write.table, x = chrX_DMRs_entrez, file = paste("Entrez ID Lists/Cord Blood ChrX DMRs/EntrezIDs_", names(chrX_DMRs_entrez), "_chrX_DMRs.txt", sep = ""), 
       MoreArgs = list(sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE))

# Use entrez IDs instead of ENSGs
# Many DMR genes are lost after converting from entrez to ensembl
# Used entrez IDs for DAVID enrichment
# Used entrez IDs/NCBI refseq for gene assignment

sapply(DMRs_entrez, length)
# MalesDisc FemalesDisc    MalesRep  FemalesRep 
#       949        2912        5429        9496 

sapply(DMRs_ENSGs, nrow)
# MalesDisc FemalesDisc    MalesRep  FemalesRep 
#       591        1984        3735        6890 

# Get Genetic and Expression Study Gene Lists -----------------------------
# Convert ENSGs to Entrez IDs Overlaps from Vogel Ciernia et al 2019 ####
# Read in Files
files <- list.files("ENSG ID Lists")
files <- files[!grepl("CpG", files, fixed = TRUE) & !grepl("CpH", files, fixed = TRUE) & !grepl("H3K", files, fixed = TRUE)]
files <- paste("ENSG ID Lists/", files, sep = "")
ENSG <- lapply(files, read.delim, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
ENSG <- lapply(ENSG, unlist, use.names = FALSE)

# Convert to Entrez IDs
ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
entrezIDs <- lapply(ENSG, getBM, attributes = "entrezgene_id", filters = "ensembl_gene_id", mart = ensembl)
entrezIDs <- lapply(entrezIDs, unlist, use.names = FALSE)

# Write Files
newFiles <- gsub("ENSG ID Lists/ENSGs", replacement = "Entrez ID Lists/EntrezIDs", x = files, fixed = TRUE)
mapply(write.table, x = entrezIDs, file = newFiles, MoreArgs = list(sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE))

# Mordaunt 2019 Cord Expression Gene Lists ####
probes <- list(ASD = read.csv("Tables/Mordaunt 2019 Cord Expression ASD vs TD.csv", header = TRUE, 
                              stringsAsFactors = FALSE),
               NonTD = read.csv("Tables/Mordaunt 2019 Cord Expression NonTD vs TD.csv", header = TRUE,
                                stringsAsFactors = FALSE)) %>%
        lapply(function(x) x$Probe[x$Meta_diff])

ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
entrezIDs <- lapply(probes, getBM, attributes = "entrezgene_id", filters = "affy_hugene_2_0_st_v1", mart = ensembl) %>%
        lapply(unlist, use.names = FALSE)

mapply(write.table, x = entrezIDs, file = c("Entrez ID Lists/EntrezIDs_CordBlood_Mordaunt_ASDvsTD_DEG.txt",
                                            "Entrez ID Lists/EntrezIDs_CordBlood_Mordaunt_NonTDvsTD_DEG.txt"), 
       MoreArgs = list(sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE))

# Get EWAS Gene Lists ----------------------------------------------------------
# Load Regions ####
DMRs_hg18 <- lapply(c("Tables/Nguyen LCL DMRs hg18.csv", "Tables/Wong Blood DMPs hg18.csv"), read.csv, header = TRUE, 
                    stringsAsFactors = FALSE)
names(DMRs_hg18) <- c("LCL_Nguyen", "Blood_Wong")

DMRs_hg19 <- lapply(c("Tables/Nardone BA24 DMPs hg19.csv", "Tables/Nardone BA10 DMPs hg19.csv", 
                      "Tables/Shulha BA10 Neurons DAPs hg19.csv", "Tables/Ladd-Acosta Temporal Cortex DMRs hg19.csv",
                      "Tables/Ladd-Acosta Cerebellum DMRs hg19.csv", "Tables/Berko Buccal DMRs hg19.csv", 
                      "Tables/Feinberg Sperm DMRs hg19.csv", "Tables/Sun Cerebellum DAPs hg19.csv",
                      "Tables/Sun Temporal Cortex DAPs hg19.csv", "Tables/Sun Prefrontal Cortex DAPs hg19.csv",
                      "Tables/Ellis BA19 CpH DMRs hg19.csv", "Tables/Ellis BA19 CpG DMRs hg19.csv",
                      "Tables/Nardone Cortex Neuron DMRs hg19.csv", "Tables/Hannon Bloodspot DMPs hg19.csv",
                      "Tables/Andrews Blood DMPs hg19.csv", "Tables/Wong Cerebellum DMPs hg19.csv",
                      "Tables/Wong Temporal Cortex DMPs hg19.csv", "Tables/Wong BA9 DMPs hg19.csv",
                      "Tables/Wong Cerebellum DMPs Dup15 hg19.csv", "Tables/Wong Temporal Cortex DMPs Dup15 hg19.csv",
                      "Tables/Wong BA9 DMPs Dup15 hg19.csv"),
                    read.csv, header = TRUE, stringsAsFactors = FALSE)
names(DMRs_hg19) <- c("CingulateCortex_Nardone", "PFC_Nardone", "PFCneurons_Shulha_H3K4me3", "TemporalCortex_LaddAcosta",
                      "Cerebellum_LaddAcosta", "Buccal_Berko", "Sperm_Feinberg", "Cerebellum_Sun_H3K27Ac", 
                      "TemporalCortex_Sun_H3K27Ac", "PFC_Sun_H3K27Ac", "VisualCortex_Ellis_CpH", "VisualCortex_Ellis", 
                      "PFCneurons_Nardone", "Bloodspot_Hannon", "Blood_Andrews", "Cerebellum_Wong", "TemporalCortex_Wong", 
                      "PFC_Wong", "Cerebellum_Wong_Dup15", "TemporalCortex_Wong_Dup15", "PFC_Wong_Dup15")
DMRs_hg19$Bloodspot_Hannon$chr <- gsub(" ", replacement = "", DMRs_hg19$Bloodspot_Hannon$chr) # Remove spaces

DMRs_hg38 <- lapply(c("Tables/Zhu Placenta DMRs hg38.csv", "Tables/Vogel Ciernia Brain DMRs hg38.csv",
                      "Tables/Vogel Ciernia Dup15q Brain DMRs hg38.csv", "Tables/Vogel Ciernia RTT Brain DMRs hg38.csv"),
                    read.csv, header = TRUE, stringsAsFactors = FALSE)
names(DMRs_hg38) <- c("Placenta_Zhu", "PFC_VogelCiernia_ASD", "PFC_VogelCiernia_Dup15", "PFC_VogelCiernia_RTT")

# LiftOver to hg38 ####
liftOverDMRs <- function(DMRs, chainFile){
        chain <- import.chain(chainFile)
        DMRs <- subset(DMRs, chr %in% c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")) %>% makeGRange(direction = "all")
        seqlevelsStyle(DMRs) <- "UCSC"
        DMRs <- suppressWarnings(liftOver(DMRs, chain = chain)) %>% unlist() %>% disjoin() %>% as.data.frame()
        DMRs <- DMRs[,c("seqnames", "start", "end")]
        colnames(DMRs) <- c("chr", "start", "end")
        DMRs$chr <- as.character(DMRs$chr)
        return(DMRs)
}

DMRs_hg18tohg38 <- lapply(DMRs_hg18, liftOverDMRs, chainFile = "Tables/hg18ToHg38.over.chain")
DMRs_hg19tohg38 <- lapply(DMRs_hg19, liftOverDMRs, chainFile = "Tables/hg19ToHg38.over.chain")
DMRs <- append(DMRs_hg18tohg38, DMRs_hg19tohg38) %>% append(DMRs_hg38)

# Get Genes ####
genes <- lapply(DMRs, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_name")
genes <- genes[sort(names(genes), decreasing = TRUE)]

# Get EntrezIDs ####
entrezIDs <- lapply(DMRs, getDMRgeneList, regDomains = regDomains, direction = "all", type = "gene_entrezID")
entrezIDs <- entrezIDs[sort(names(entrezIDs), decreasing = TRUE)]
names(entrezIDs) <- c("VisualCortex_Ellis_mCpH", "VisualCortex_Ellis_mCpG", "TemporalCortex_Wong_Dup15_mCpG", 
                      "TemporalCortex_Wong_ASD_mCpG", "TemporalCortex_Sun_H3K27Ac", "TemporalCortex_LaddAcosta_mCpG", 
                      "Sperm_Feinberg_mCpG", "Placenta_Zhu_mCpG", "PFCneurons_Shulha_H3K4me3", 
                      "PFCneurons_Nardone_mCpG", "PFC_Wong_Dup15_mCpG", "PFC_Wong_ASD_mCpG", "PFC_VogelCiernia_RTT_mCpG",
                      "PFC_VogelCiernia_Dup15_mCpG", "PFC_VogelCiernia_ASD_mCpG", "PFC_Sun_H3K27Ac", "PFC_Nardone_mCpG", 
                      "LCL_Nguyen_mCpG", "CingulateCortex_Nardone_mCpG", "Cerebellum_Wong_Dup15_mCpG", 
                      "Cerebellum_Wong_ASD_mCpG", "Cerebellum_Sun_H3K27Ac", "Cerebellum_LaddAcosta_mCpG", 
                      "Buccal_Berko_mCpG", "Bloodspot_Hannon_mCpG", "Blood_Wong_mCpG", "Blood_Andrews_mCpG")
mapply(write.table, x = entrezIDs, file = paste("Entrez ID Lists/EntrezIDs_", names(entrezIDs), ".txt", sep = ""), 
       MoreArgs = list(sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE))

# Overlap EWAS Genes -----------------------------------------------------------
# Compare Gene Lists ####
allGenes <- unique(regDomains$gene_name) %>% sort()
gom <- newGOM(genes, genome.size = length(allGenes))
gomResults <- cbind(getMatrix(gom, "intersection") %>% melt(),
                    getMatrix(gom, "odds.ratio") %>% melt() %>% .[,"value"],
                    getMatrix(gom, "pval") %>% melt() %>% .[,"value"])
colnames(gomResults) <- c("ListA", "ListB", "Intersection", "OddsRatio", "pValue")
gomResults$ListA_Length <- sapply(as.character(gomResults$ListA), function(x) length(genes[[x]]))
gomResults$ListB_Length <- sapply(as.character(gomResults$ListB), function(x) length(genes[[x]]))
gomResults$PerListA <- gomResults$Intersection * 100 / gomResults$ListA_Length
gomResults$PerListB <- gomResults$Intersection * 100 / gomResults$ListB_Length
gomResults <- subset(gomResults, !as.character(gomResults$ListA) == as.character(gomResults$ListB))
gomResults$qValue <- p.adjust(gomResults$pValue, method = "fdr")
gomResults$log_qValue <- -log10(gomResults$qValue)
gomResults$Significant <- gomResults$qValue < 0.05
gomResults$Significant <- factor(gomResults$Significant, levels = c("TRUE", "FALSE"))
gomResults <- gomResults[,c("ListA", "ListB", "ListA_Length", "ListB_Length", "Intersection", "PerListA", "PerListB", 
                            "OddsRatio", "pValue", "qValue", "log_qValue", "Significant")]
write.table(gomResults, file = "Tables/Autism EWAS Gene Overlap Stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)
table(gomResults$Significant) # 61 / 421 comparisons significantly overlap

plotGOMheatmap(gomResults, type = "log_OddsRatio", file = "Figures/Autism EWAS Gene Overlap log Odds Ratio.png",
               intersect.size = 3, sig.size = 7, axis.text.size = 12, plot.margin = c(1, 5, 1, 1),
               legend.position = c(1.08, 0.855), expand = c(0.027, 0))

# Get Counts for Each Gene ####
geneCounts <- sapply(allGenes, function(x) sapply(genes, function(y) x %in% y)) %>% t() %>% as.data.frame()
geneCounts$Gene <- rownames(geneCounts)
geneCounts$Sum <- rowSums(geneCounts[,colnames(geneCounts)[!colnames(geneCounts) == "Gene"]])
geneCounts <- subset(geneCounts, Sum > 1) # 10108 genes in at least 2 lists
geneCounts$Sum_CpGmeth <- rowSums(geneCounts[,c("LCL_Nguyen", "Blood_Wong", "CingulateCortex_Nardone", "PFC_Nardone", 
                                                "TemporalCortex_LaddAcosta", "Cerebellum_LaddAcosta", "Buccal_Berko", 
                                                "Sperm_Feinberg", "VisualCortex_Ellis", "PFCneurons_Nardone", 
                                                "Bloodspot_Hannon", "Blood_Andrews", "Cerebellum_Wong", "TemporalCortex_Wong", 
                                                "PFC_Wong", "Placenta_Zhu", "PFC_VogelCiernia")])
geneCounts$Sum_Brain <- rowSums(geneCounts[,c("CingulateCortex_Nardone", "PFC_Nardone", "TemporalCortex_LaddAcosta",
                                              "Cerebellum_LaddAcosta", "VisualCortex_Ellis", "PFCneurons_Nardone",
                                              "Cerebellum_Wong", "TemporalCortex_Wong", "PFC_Wong", "PFC_VogelCiernia")])
geneCounts$Sum_Peripheral <- rowSums(geneCounts[,c("LCL_Nguyen", "Blood_Wong", "Buccal_Berko", "Sperm_Feinberg",
                                                   "Bloodspot_Hannon", "Blood_Andrews", "Placenta_Zhu")])
geneCounts$Sum_Blood <- rowSums(geneCounts[,c("Blood_Andrews", "Blood_Wong", "Bloodspot_Hannon", "LCL_Nguyen")])
geneCounts$Sum_Cerebellum <- rowSums(geneCounts[,c("Cerebellum_LaddAcosta", "Cerebellum_Wong")])
geneCounts$Sum_PFC <- rowSums(geneCounts[,c("PFC_Nardone", "PFC_VogelCiernia", "PFC_Wong", "PFCneurons_Nardone")])
geneCounts$Sum_TemporalCortex <- rowSums(geneCounts[,c("TemporalCortex_LaddAcosta", "TemporalCortex_Wong")])
write.table(geneCounts, file = "Tables/Autism EWAS Gene Overlap Counts.txt", sep = "\t", quote = FALSE, row.names = FALSE)

geneCounts_brain <- subset(geneCounts, Sum_Brain == 4)
geneCounts_brain$GeneLists <- apply(geneCounts_brain[,1:22], 1, function(x) paste(colnames(geneCounts_brain[,1:22])[x], collapse = ", "))
geneCounts_brain <- geneCounts_brain[,c("Gene", "Sum", "Sum_CpGmeth", "Sum_Brain", "Sum_Peripheral", "Sum_Blood",
                                        "Sum_Cerebellum", "Sum_PFC", "Sum_TemporalCortex", "GeneLists")]
write.table(geneCounts_brain, file = "Tables/Autism EWAS Top Brain Gene Overlap Counts.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Top Genes ####
table(geneCounts$Sum)
#    2    3    4    5    6    7    8 
# 5389 3040 1268  336   62   10    3 
geneCounts$Gene[geneCounts$Sum %in% c(7, 8)]
# "AGAP1-IT1"    "ALPL"         "DAB2IP"       "GBX2"         "KIF26A"       "KIF26B"       "LINC00676"    "LOC101927815"
# "MIR7850"      "NEURL1-AS1"   "RAP1GAP"      "RAPGEF1"      "SH3PXD2A-AS1"
# 13 genes in 7 or 8 of 22 (34%) differential epigenetic modification lists

table(geneCounts$Sum_CpGmeth)
#    0    1    2    3    4    5    6 
# 1606 4039 3835  568   52    6    2 
geneCounts$Gene[geneCounts$Sum_CpGmeth %in% c(5, 6)]
# "AGAP1-IT1"    "FOXD1"        "GBX2"         "KIF26B"       "LINC00676"    "LOC340090"    "NEURL1-AS1"   "SH3PXD2A-AS1"
# 8 genes in 5 or 6 of 17 (32%) differential CpG methylation lists

table(geneCounts$Sum_Brain)
#    0    1    2    3    4 
# 1760 4352 3641  335   20 
geneCounts$Gene[geneCounts$Sum_Brain == 4]
# "AGAP1-IT1"    "ALPL"         "C14orf180"    "CASZ1"        "DFFA"         "F7"           "GBX2"         "HTRA3"       
# "JAKMIP3"      "KCNG2"        "KIF26A"       "MCF2L-AS1"    "MTSS1L"       "NEURL1-AS1"   "PPT2"         "PPT2-EGFL8"  
# "RAP1GAP"      "RAPGEF1"      "SH3PXD2A-AS1" "SOHLH1"   
# 20 genes in 4 of 10 (40%) brain differential CpG methylation lists

table(geneCounts$Sum_Peripheral)
#    0    1    2    3 
# 9195  884   27    2 
geneCounts$Gene[geneCounts$Sum_Peripheral %in% c(2, 3)]
# "ADAM21P1"  "ADRA2C"    "BEST2"     "CCDC88C"   "CDH20"     "CYP2E1"    "FBXO34"    "FOXD1"     "GPR68"     "HOOK2"    
# "KIF26B"    "LINC00676" "LINC01060" "LINC01262" "LOC340090" "PAX8"      "PAX8-AS1"  "PDHB"      "PPIF"      "PRDM1"    
# "PREP"      "PYGB"      "SDCBP"     "SMYD3"     "SST"       "STPG2"     "SYCE1"     "TMEM14B"   "TRIM2"  
# 29 genes in 2 or 3 of 7  (36%) peripheral tissue differential CpG methylation lists

table(geneCounts$Sum_Blood) 
#    0    1    2 
# 9878  227    3 
geneCounts$Gene[geneCounts$Sum_Blood == 2]
# "PPIF"  "PYGB"  "TRIM2"
# 3 genes in 2 of 3 (67%) blood differential CpG methylation lists

table(geneCounts$Sum_Cerebellum)
#     0     1 
# 10105     3 
# No genes replicated in 2 cerebellum differential CpG methylation lists

table(geneCounts$Sum_PFC)
#    0    1    2    3 
# 4963 4770  366    9 
geneCounts$Gene[geneCounts$Sum_PFC == 3]
# "AGAP1-IT1" "C14orf180" "F7"        "GBX2"      "JAKMIP3"   "KCNG2"     "KIF26A"    "MCF2L-AS1" "MIR6829"  
# 9 genes in 3 of 4 (75%) PFC differential CpG methylation lists

table(geneCounts$Sum_TemporalCortex)
#     0     1 
# 10036    72 
# No genes replicated in 2 temporal cortex differential CpG methylation lists

# Overlap EWAS Regions ---------------------------------------------------------
# Get Overlaps ####
GR_DMRs <- lapply(DMRs, makeGRange, direction = "all")

# Count Number of Intersecting Regions
intersectCounts <- sapply(GR_DMRs, function(x) sapply(GR_DMRs, function(y) intersect(x, y) %>% length()))

# Get Actual Intersecting Regions
intersectRegions <- lapply(GR_DMRs, function(x) lapply(GR_DMRs, function(y) intersect(x, y))) %>%
        lapply(function(x) x[sapply(x, length) > 0]) %>% .[sapply(., length) > 1] %>%
        lapply(function(x) lapply(x, as.data.frame))

# Remove Self Intersects
for(i in 1:length(intersectRegions)){
        intersectRegions[[i]] <- intersectRegions[[i]][!names(intersectRegions[[i]]) == names(intersectRegions)[i]]
}

# Add List Identifiers and Combine
for(i in 1:length(intersectRegions)){
        for(j in 1:length(intersectRegions[[i]])){
                intersectRegions[[i]][[j]]$Region1 <- names(intersectRegions)[i]
                intersectRegions[[i]][[j]]$Region2 <- names(intersectRegions[[i]])[j]
        }
}
intersectRegions <- lapply(intersectRegions, list.rbind) %>% list.rbind()
rownames(intersectRegions) <- 1:nrow(intersectRegions)
intersectRegions <- intersectRegions[,c("seqnames", "start", "end", "Region1", "Region2")]
colnames(intersectRegions)[colnames(intersectRegions) == "seqnames"] <- "chr"
table(intersectRegions$Region1, intersectRegions$Region2)
#' Blood_Andrews: CingulateCortex_Nardone, PFC_Nardone
#' Cerebellum_Sun_H3K27Ac: PFC_Sun_H3K27Ac, TemporalCortex_Sun_H3K27Ac
#' CingulateCortex_Nardone: Blood_Andrews, PFC_Nardone, PFC_Wong, TemporalCortex_Wong
#' PFC_Nardone: Blood_Andrews, CingulateCortex_Nardone, PFC_Wong, TemporalCortex_Wong
#' PFC_Sun_H3K27Ac: Cerebellum_Sun_H3K27Ac, TemporalCortex_Sun_H3K27Ac
#' PFC_Wong: CingulateCortex_Nardone, PFC_Nardone
#' PFCneurons_Shulha_H3K4me3: TemporalCortex_Sun_H3K27Ac
#' TemporalCortex_Sun_H3K27Ac: Cerebellum_Sun_H3K27Ac, PFC_Sun_H3K27Ac, PFCneurons_Shulha_H3K4me3
#' TemporalCortex_Wong: CingulateCortex_Nardone, PFC_Nardone

intersectRegions_top <- subset(intersectRegions, (Region1 == "Blood_Andrews" & Region2 == "CingulateCortex_Nardone") |
                                       (Region1 == "Blood_Andrews" & Region2 == "PFC_Nardone") |
                                       (Region1 == "CingulateCortex_Nardone" & Region2 == "PFC_Wong") |
                                       (Region1 == "CingulateCortex_Nardone" & Region2 == "TemporalCortex_Wong") |
                                       (Region1 == "PFC_Nardone" & Region2 == "PFC_Wong") |
                                       (Region1 == "PFC_Nardone" & Region2 == "TemporalCortex_Wong"))
intersectRegions_top$chr <- as.character(intersectRegions_top$chr) %>% 
        factor(levels = c(paste("chr", 1:22, sep = ""), "chrX", "chrY"))
intersectRegions_top <- intersectRegions_top[order(intersectRegions_top$chr, intersectRegions_top$start),]
intersectRegions_top$DMRid <- paste("DMR", 1:nrow(intersectRegions_top), sep = "_")
intersectRegions_top <- getDMRanno(intersectRegions_top, regDomains = regDomains, 
                                   file = "Tables/Autism EWAS Region Overlap Top Regions.txt")
intersectRegions_topGenes <- getDMRgeneList(intersectRegions_top, regDomains = regDomains, direction = "all", type = "gene_name")
# [1] "ABHD12"       "COL8A2"       "ITGB7"        "ITIH1"        "ITIH3"        "LOC100130238" "LOC101928416" "LOC102724467"
# [9] "MIR6775"      "PRSS37"       "PYGB"         "RBMS2"        "TAS2R5"      

# DMR Gene Autism Gene Overlap Analysis -----------------------------------
# Get Gene List for Negative Control Epigenetic Study ####
# Load Regions
DMRs_hg19 <- read.csv("Tables/Lunnon Cortex DMPs AD hg19.csv", stringsAsFactors = FALSE)

# LiftOver to hg38
DMRs_hg38 <- liftOverDMRs(DMRs_hg19, chainFile = "Tables/hg19ToHg38.over.chain")

# Get EntrezIDs
entrezIDs <- getDMRgeneList(DMRs_hg38, regDomains = regDomains, direction = "all", type = "gene_entrezID")
rm(DMRs_hg19, DMRs_hg38)

# Load Gene Lists ####
# Autism Gene Lists
files <- list.files(c("Entrez ID Lists/Genetic Studies", "Entrez ID Lists/Gene Expression Studies",
                      "Entrez ID Lists/Epigenetic Studies"), full.names = TRUE)
ASDgenes <- sapply(files, read.delim, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
names(ASDgenes) <- gsub(pattern = ".*EntrezIDs_", replacement = "", x = names(ASDgenes)) %>%
        gsub(pattern = ".txt.V1", replacement = "", fixed = TRUE)
ASDgenes$Cortex_Lunnon_AD_mCpG <- entrezIDs
ASDgenes <- ASDgenes[sort(names(ASDgenes))]

# DMR Gene Lists
files <- list.files("Entrez ID Lists/Cord Blood DMRs", full.names = TRUE)
DMRfiles <- files[!grepl("Background", files, fixed = TRUE)]
DMRgenes <- sapply(DMRfiles, read.delim, header = FALSE, sep = "\t", stringsAsFactors = FALSE) %>% lapply(list)
names(DMRgenes) <- gsub(pattern = ".*EntrezIDs_", replacement = "", x = names(DMRgenes)) %>%
        gsub(pattern = "_DMRs.txt.V1", replacement = "", fixed = TRUE)
DMRgenes <- DMRgenes[c("Males_Discovery", "Males_Replication", "Females_Discovery", "Females_Replication")]

# Background Gene Lists
BackgroundFiles <- files[grepl("Background", files, fixed = TRUE)]
BackgroundGenes <- sapply(BackgroundFiles, read.delim, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
names(BackgroundGenes) <- gsub(pattern = ".*EntrezIDs_", replacement = "", x = names(BackgroundGenes)) %>%
        gsub(pattern = "_Background.txt.V1", replacement = "", fixed = TRUE)
BackgroundGenes <- BackgroundGenes[c("Males_Discovery", "Males_Replication", "Females_Discovery", "Females_Replication")]
BackgroundGeneCount <- sapply(BackgroundGenes, length)
rm(files, DMRfiles, BackgroundFiles, BackgroundGenes)

# Run GeneOverlap Stats ####
gom <- mapply(newGOM, gsetA = DMRgenes, genome.size = BackgroundGeneCount, MoreArgs = list(gsetB = ASDgenes))
gomResults <- sapply(gom, getMatrix, name = "intersection") %>% melt()
colnames(gomResults) <- c("ASD_GeneList", "DMR_GeneList", "Intersection")
gomResults$OddsRatio <- sapply(gom, getMatrix, name = "odds.ratio") %>% melt() %>% .[,"value"]
gomResults$pValue <- sapply(gom, getMatrix, name = "pval") %>% melt() %>% .[,"value"]
gomResults$ASD_GeneCount <- sapply(ASDgenes, length)
gomResults$DMR_GeneCount <- sapply(DMRgenes, function(x) length(x[[1]])) %>% rep(each = length(ASDgenes))
gomResults$Per_ASDgenes <- gomResults$Intersection * 100 / gomResults$ASD_GeneCount
gomResults$Per_DMRgenes <- gomResults$Intersection * 100 / gomResults$DMR_GeneCount
gomResults$qValue <- split(gomResults$pValue, f = gomResults$DMR_GeneList) %>% sapply(p.adjust, method = "fdr") %>% 
        melt() %>% .[,"value"] # Correct for number of ASD gene lists
gomResults$Significant <- factor(gomResults$qValue < 0.05, levels = c("TRUE", "FALSE"))
gomResults <- gomResults[,c("DMR_GeneList", "ASD_GeneList", "DMR_GeneCount", "ASD_GeneCount", "Intersection",
                            "Per_DMRgenes", "Per_ASDgenes", "OddsRatio", "pValue", "qValue", "Significant")]
write.table(gomResults, file = "Tables/DMR Gene Autism Gene Overlap Analysis Results.txt", sep = "\t", quote = FALSE, 
            row.names = FALSE)

# Add Intersecting Genes for Supplemental Table
gomResults$Intersection_GeneEntrezIDs <- sapply(gom, getNestedList, name = "intersection") %>% sapply(unlist) %>% sapply(paste, collapse = ", ")
gomResults$Intersection_GeneNames <- sapply(gomResults$Intersection_GeneEntrezIDs, entrezIDs_to_genes, regDomains = regDomains)
gomResults <- gomResults[,c("DMR_GeneList", "ASD_GeneList", "DMR_GeneCount", "ASD_GeneCount", "Intersection", 
                            "Intersection_GeneNames", "Intersection_GeneEntrezIDs", "Per_DMRgenes", "Per_ASDgenes", 
                            "OddsRatio", "pValue", "qValue", "Significant")]
Genetic <- list.files("Entrez ID Lists/Genetic Studies/") %>% str_replace_all(pattern = c("EntrezIDs_" = "", ".txt" = ""))
Epigenetic <- list.files("Entrez ID Lists/Epigenetic Studies/") %>% str_replace_all(pattern = c("EntrezIDs_" = "", ".txt" = ""))
Expression <- list.files("Entrez ID Lists/Gene Expression Studies//") %>% str_replace_all(pattern = c("EntrezIDs_" = "", ".txt" = ""))
gomResults$Category <- NA
gomResults$Category[gomResults$ASD_GeneList %in% Genetic] <- "Genetic"
gomResults$Category[gomResults$ASD_GeneList %in% Epigenetic] <- "Epigenetic"
gomResults$Category[gomResults$ASD_GeneList == "Cortex_Lunnon_AD_mCpG"] <- "Epigenetic"
gomResults$Category[gomResults$ASD_GeneList %in% Expression] <- "Expression"
gomResults$Category <- factor(gomResults$Category, levels = c("Genetic", "Epigenetic", "Expression"))
pattern <- c("Blood_Gilissen_de_novo_CNV_ID" = "ID de novo CNV Gilissen", 
             "Blood_Gilissen_Known_Genes" = "ASD Gilissen",
             "Blood_Iossifov_Genes_With_LGD_Mutations_In_ASD" = "ASD LGD Mutations Iossifov", 
             "Blood_Iossifov_Genes_With_LGD_Mutations_In_ID" = "ID LGD Mutations Iossifov ", 
             "Blood_Kochinke_ID_genes" = "ID Kochinke", "Blood_Sanders_High_Risk_ASD" = "ASD High Risk Sanders", 
             "Grove_GWAS" = "ASD Common Variants Grove", "Lambert" = "Alzheimer's GWAS Lambert", 
             "SFARI_HighConfidence" = "ASD SFARI High Confidence", "SFARI_Hypothesized" = "ASD SFARI Hypothesized", 
             "SFARI_Minimal" = "ASD SFARI Minimal Confidence", "SFARI_StrongCandidate" = "ASD SFARI Strong Confidence", 
             "SFARI_Suggestive" = "ASD SFARI Suggestive", "SFARI_Syndromic" = "ASD SFARI Syndromic", 
             "TheASD_Consortium_GWAS" = "ASD GWAS PGC", "Blood_Andrews_mCpG" = "ASD mCpG Blood Andrews", 
             "Blood_Wong_mCpG" = "ASD mCpG Blood Wong", "Bloodspot_Hannon_mCpG" = "ASD mCpG Bloodspot Hannon",         
             "Buccal_Berko_mCpG" = "ASD mCpG Buccal Berko", "Cerebellum_LaddAcosta_mCpG" = "ASD mCpG Cerebellum Ladd-Acosta", 
             "Cerebellum_Sun_H3K27Ac" = "ASD H3K27Ac Cerebellum Sun", "Cerebellum_Wong_ASD_mCpG" = "ASD mCpG Cerebellum Wong",
             "Cerebellum_Wong_Dup15_mCpG" = "Dup15 mCpG Cerebellum Wong", 
             "CingulateCortex_Nardone_mCpG" = "ASD mCpG Cingulate Cortex Nardone", "LCL_Nguyen_mCpG" = "ASD mCpG LCL Nguyen",
             "PFC_Nardone_mCpG" = "ASD mCpG PFC Nardone", "PFC_Sun_H3K27Ac" = "ASD H3K27Ac PFC Sun", 
             "PFC_VogelCiernia_ASD_mCpG" = "ASD mCpG PFC Vogel Ciernia", 
             "PFC_VogelCiernia_Dup15_mCpG" = "Dup15 mCpG PFC Vogel Ciernia",
             "PFC_VogelCiernia_RTT_mCpG" = "RTT mCpG PFC Vogel Ciernia", "PFC_Wong_ASD_mCpG" = "ASD mCpG PFC Wong", 
             "PFC_Wong_Dup15_mCpG" = "Dup15 mCpG PFC Wong", "PFCneurons_Nardone_mCpG" = "ASD mCpG PFC Neurons Nardone",       
             "PFCneurons_Shulha_H3K4me3" = "ASD H3K4me3 PFC Neurons Shulha", "Placenta_Zhu_mCpG" = "ASD mCpG Placenta Zhu",
             "Sperm_Feinberg_mCpG" = "ASD mCpG Sperm Feinberg", 
             "TemporalCortex_LaddAcosta_mCpG" = "ASD mCpG Temporal Cortex Ladd-Acosta", 
             "TemporalCortex_Sun_H3K27Ac" = "ASD mCpG Temporal Cortex Sun", 
             "TemporalCortex_Wong_ASD_mCpG" = "ASD mCpG Temporal Cortex Wong",  
             "TemporalCortex_Wong_Dup15_mCpG" = "Dup15 mCpG Temporal Cortex Wong", 
             "VisualCortex_Ellis_mCpG" = "ASD mCpG Visual Cortex Ellis", 
             "VisualCortex_Ellis_mCpH" = "ASD mCpH Visual Cortex Ellis", "Blood_Tylee_Array_DEG" = "ASD Blood Tylee",                                        
             "CordBlood_Mordaunt_ASDvsTD_DEG" = "ASD Cord Blood Mordaunt", 
             "CordBlood_Mordaunt_NonTDvsTD_DEG" = "Non-TD Cord Blood Mordaunt", "Cortex_Gupta_ASD_DGE" = "ASD Cortex Gupta",                                         
             "FrontalAndTemporalCortex_Gandal_ASD_DGE_RNAseq" = "ASD Cortex Gandal RNA-seq",               
             "FrontalAndTemporalCortex_Gandal_ASD_DGE" = "ASD Cortex Gandal Array",                      
             "FrontalAndTemporalCortex_Lin_Rett_Vs_Ctrl_DGE" = "RTT Cortex Lin",                
             "FrontalAndTemporalCortex_Parikshak_ASD_DGE" = "ASD Cortex Parikshak",                   
             "FrontalAndTemporalCortex_Parikshak_DifferentiallySplicedGenes" = "ASD Cortex Parikshak Splicing",
             "FrontalAndTemporalCortex_Parikshak_Dup15qVs.Ctrl_Genes" = "Dup15 Cortex Parikshak",       
             "ImprintedGenes_MaternalExpression" = "Imprinted, Maternally Expressed",                            
             "ImprintedGenes_PaternalExpression" = "Imprinted, Paternally Expressed",                            
             "LCL_Tylee_RNAseq_DEG" = "ASD LCL Tylee", "Cortex_Lunnon_AD_mCpG" = "Alzheimer's mCpG Cortex Lunnon", "_" = " ")
gomResults$ASD_GeneList <- as.character(gomResults$ASD_GeneList) %>% str_replace_all(pattern = pattern) %>%
        factor(., levels = sort(unique(. ), decreasing = TRUE))
write.table(gomResults, file = "Tables/DMR Gene Autism Gene Overlap Analysis Results with Genes.txt", sep = "\t", quote = FALSE, 
            row.names = FALSE)

# Plot Heatmap for All Overlaps ####
gomResults$OddsRatio[is.infinite(gomResults$OddsRatio)] <- max(gomResults$OddsRatio[!is.infinite(gomResults$OddsRatio)])
gg <- ggplot(data = gomResults)
gg +
        geom_tile(aes(x = DMR_GeneList, y = ASD_GeneList, fill = OddsRatio)) +
        geom_text(aes(x = DMR_GeneList, y = ASD_GeneList, alpha = Significant), label = "*", color = "white", size = 10,
                  nudge_y = -0.44) +
        facet_grid(rows = vars(Category), scales = "free_y", space = "free_y") +
        scale_fill_gradientn("Odds Ratio", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", 
                             limits = c(0, max(gomResults$OddsRatio[gomResults$ASD_GeneCount >= 5])),
                             breaks = pretty_breaks(n = 3)) +
        scale_alpha_manual(breaks = c("TRUE", "FALSE"), values = c(1, 0), guide = FALSE) +
        scale_x_discrete(expand = c(0.175, 0), labels = function(x) str_replace_all(x, pattern = c("_" = " "))) +
        scale_y_discrete(expand = c(0.025, 0)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              plot.margin = unit(c(0.5, 7, 0.5, 1), "lines"), axis.ticks = element_line(size = 1), 
              panel.background = element_rect(fill = "black"),
              axis.text.x = element_text(size = 16, color = "black", angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_text(size = 13, color = "black"), axis.title = element_blank(), 
              legend.key = element_blank(), legend.position = c(1.25, 0.94), legend.background = element_blank(), 
              legend.key.size = unit(1, "lines"), legend.title = element_text(size = 18), 
              legend.text = element_text(size = 16), strip.background = element_blank())
ggsave("Figures/DMR Gene Autism Gene Overlap Odds Ratio Heatmap All.png", dpi = 600, width = 10, height = 13, units = "in")

# Plot Heatmap for Replicated Overlaps ####
replicatedMales <- gomResults$Significant[gomResults$DMR_GeneList == "Males_Discovery"] == "TRUE" &
        gomResults$Significant[gomResults$DMR_GeneList == "Males_Replication"] == "TRUE"
replicatedFemales <- gomResults$Significant[gomResults$DMR_GeneList == "Females_Discovery"] == "TRUE" &
        gomResults$Significant[gomResults$DMR_GeneList == "Females_Replication"] == "TRUE"
gomResults_rep <- subset(gomResults, replicatedMales | replicatedFemales)
gomResults_rep$ASD_GeneList <- str_replace_all(gomResults_rep$ASD_GeneList, 
                                               pattern = c("ASD Gilissen" = "ASD Genes Gilissen",
                                                           "ASD SFARI Strong Confidence" = "ASD Genes SFARI Strong Confidence",
                                                           "RTT Cortex Lin" = "RTT DEG Cortex Lin")) %>%
        factor(levels = rev(c("ASD H3K27Ac Cerebellum Sun", "ASD H3K27Ac PFC Sun", "ASD mCpG Cingulate Cortex Nardone",  
                              "ASD mCpG PFC Vogel Ciernia", "ASD mCpG Sperm Feinberg", "ASD mCpG Temporal Cortex Sun", 
                              "Dup15 mCpG PFC Vogel Ciernia", "RTT mCpG PFC Vogel Ciernia", 
                              "ASD Genes Gilissen", "ASD Genes SFARI Strong Confidence", "ASD H3K4me3 PFC Neurons Shulha", 
                              "ASD mCpG Bloodspot Hannon", "ASD mCpG LCL Nguyen", "ASD mCpG PFC Nardone", 
                              "ASD mCpG Placenta Zhu", "ASD mCpG Temporal Cortex Wong", "ASD mCpH Visual Cortex Ellis",
                              "Dup15 mCpG Temporal Cortex Wong", "RTT DEG Cortex Lin")))

gg <- ggplot(data = gomResults_rep)
gg +
        geom_tile(aes(x = DMR_GeneList, y = ASD_GeneList, fill = OddsRatio)) +
        geom_text(aes(x = DMR_GeneList, y = ASD_GeneList, alpha = Significant), label = "*", color = "white", size = 12,
                  nudge_y = -0.33) +
        scale_fill_gradientn("Odds Ratio", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", 
                             limits = c(0, max(gomResults_rep$OddsRatio[gomResults_rep$ASD_GeneCount >= 5])),
                             breaks = pretty_breaks(n = 3)) +
        scale_alpha_manual(breaks = c("TRUE", "FALSE"), values = c(1, 0), guide = FALSE) +
        scale_x_discrete(expand = c(0.18, 0), labels = function(x) str_replace_all(x, pattern = c("_" = " "))) +
        scale_y_discrete(expand = c(0.033, 0)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), 
              plot.margin = unit(c(0.5, 8, 0.5, 1), "lines"), axis.ticks = element_line(size = 1), 
              panel.background = element_rect(fill = "black"),
              axis.text.x = element_text(size = 17, color = "black", angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_text(size = 15, color = "black"), axis.title = element_blank(), 
              legend.key = element_blank(), legend.position = c(1.23, 0.89), legend.background = element_blank(), 
              legend.key.size = unit(1, "lines"), legend.title = element_text(size = 18), 
              legend.text = element_text(size = 17))
ggsave("Figures/DMR Gene Autism Gene Overlap Odds Ratio Heatmap Replicated.png", dpi = 600, width = 9, height = 8, 
       units = "in")

# Get Replicated ChrX DMR Genes Overlapping ASD Gene Lists ####
files <- list.files("Entrez ID Lists/Cord Blood ChrX DMRs", full.names = TRUE)
chrX_DMRgenes <- sapply(files, read.delim, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
names(chrX_DMRgenes) <- gsub(pattern = ".*EntrezIDs_", replacement = "", x = names(chrX_DMRgenes)) %>%
        gsub(pattern = "_chrX_DMRs.txt.V1", replacement = "", fixed = TRUE)
chrX_DMRgenes <-chrX_DMRgenes[c("Males_Discovery", "Males_Replication", "Females_Discovery", "Females_Replication")]
chrX_DMRgenes$Males_Replicated <- intersect(chrX_DMRgenes$Males_Discovery, chrX_DMRgenes$Males_Replication)
chrX_DMRgenes$Females_Replicated <- intersect(chrX_DMRgenes$Females_Discovery, chrX_DMRgenes$Females_Replication)
chrX_DMRgenes$MalesFemales_Replicated <- intersect(chrX_DMRgenes$Males_Replicated, chrX_DMRgenes$Females_Replicated)
chrX_DMRgenes$MalesOnly_Replicated <- chrX_DMRgenes$Males_Replicated[!chrX_DMRgenes$Males_Replicated %in% chrX_DMRgenes$Females_Replicated]
chrX_DMRgenes$FemalesOnly_Replicated <- chrX_DMRgenes$Females_Replicated[!chrX_DMRgenes$Females_Replicated %in% chrX_DMRgenes$Males_Replicated]
chrX_DMRgenes <- chrX_DMRgenes[c("MalesFemales_Replicated", "MalesOnly_Replicated", "FemalesOnly_Replicated")]
chrX_Overlaps <- melt(chrX_DMRgenes)
colnames(chrX_Overlaps) <- c("EntrezID", "chrX_DMRs")
chrX_Overlaps$GeneSymbol <- regDomains$gene_name[match(chrX_Overlaps$EntrezID, regDomains$gene_entrezID)]
chrX_Overlaps <- chrX_Overlaps[,c("GeneSymbol", "EntrezID", "chrX_DMRs")]
overlaps <- sapply(ASDgenes, function(x) chrX_Overlaps$EntrezID %in% x)
chrX_Overlaps <- cbind(chrX_Overlaps, overlaps)
Genetic <- Genetic[!Genetic %in% c("Blood_Gilissen_de_novo_CNV_ID", "Blood_Iossifov_Genes_With_LGD_Mutations_In_ID",
                                   "Blood_Kochinke_ID_genes", "Lambert", "SFARI_Hypothesized")]
chrX_Overlaps$GeneticStudies <- rowSums(chrX_Overlaps[,Genetic])
Epigenetic <- Epigenetic[!Epigenetic %in% c("TemporalCortex_Wong_Dup15_mCpG", "Cerebellum_Wong_Dup15_mCpG",
                                            "PFC_VogelCiernia_Dup15_mCpG", "PFC_Wong_Dup15_mCpG", 
                                            "PFC_VogelCiernia_RTT_mCpG")]
chrX_Overlaps$EpigeneticStudies <- rowSums(chrX_Overlaps[,Epigenetic])
Expression <- Expression[!Expression %in% c("CordBlood_Mordaunt_NonTDvsTD_DEG", "FrontalAndTemporalCortex_Lin_Rett_Vs_Ctrl_DGE",
                                            "FrontalAndTemporalCortex_Parikshak_Dup15qVs.Ctrl_Genes", 
                                            "ImprintedGenes_MaternalExpression", "ImprintedGenes_PaternalExpression")]
chrX_Overlaps$ExpressionStudies <- rowSums(chrX_Overlaps[,Expression])
chrX_Overlaps$AnyStudy <- rowSums(chrX_Overlaps[,c("GeneticStudies", "EpigeneticStudies", "ExpressionStudies")])
chrX_Overlaps <- chrX_Overlaps[order(chrX_Overlaps$chrX_DMRs, chrX_Overlaps$GeneSymbol),]
write.csv(chrX_Overlaps, "Tables/ChrX DMR Gene Autism Gene Overlaps.csv", quote = FALSE, row.names = FALSE)

# Get Replicated All DMR Genes Overlapping ASD Gene Lists ####
DMRgenes <- lapply(DMRgenes, unlist)
DMRgenes$Males_Replicated <- intersect(DMRgenes$Males_Discovery, DMRgenes$Males_Replication)
DMRgenes$Females_Replicated <- intersect(DMRgenes$Females_Discovery, DMRgenes$Females_Replication)
DMRgenes$MalesFemales_Replicated <- intersect(DMRgenes$Males_Replicated, DMRgenes$Females_Replicated)
DMRgenes$MalesOnly_Replicated <- DMRgenes$Males_Replicated[!DMRgenes$Males_Replicated %in% DMRgenes$Females_Replicated]
DMRgenes$FemalesOnly_Replicated <- DMRgenes$Females_Replicated[!DMRgenes$Females_Replicated %in% DMRgenes$Males_Replicated]
DMRgenes <- DMRgenes[c("MalesFemales_Replicated", "MalesOnly_Replicated", "FemalesOnly_Replicated")]
Overlaps <- melt(DMRgenes)
colnames(Overlaps) <- c("EntrezID", "DMRs")
Overlaps$GeneSymbol <- regDomains$gene_name[match(Overlaps$EntrezID, regDomains$gene_entrezID)]
Overlaps <- Overlaps[,c("GeneSymbol", "EntrezID", "DMRs")]
overlaps <- sapply(ASDgenes, function(x) Overlaps$EntrezID %in% x)
Overlaps <- cbind(Overlaps, overlaps)
Genetic <- list.files("Entrez ID Lists/Genetic Studies/") %>% str_replace_all(pattern = c("EntrezIDs_" = "", ".txt" = ""))
Epigenetic <- list.files("Entrez ID Lists/Epigenetic Studies/") %>% str_replace_all(pattern = c("EntrezIDs_" = "", ".txt" = ""))
Expression <- list.files("Entrez ID Lists/Gene Expression Studies//") %>% str_replace_all(pattern = c("EntrezIDs_" = "", ".txt" = ""))
Genetic <- Genetic[!Genetic %in% c("Blood_Gilissen_de_novo_CNV_ID", "Blood_Iossifov_Genes_With_LGD_Mutations_In_ID",
                                   "Blood_Kochinke_ID_genes", "Lambert", "SFARI_Hypothesized")]
Overlaps$GeneticStudies <- rowSums(Overlaps[,Genetic])
Epigenetic <- Epigenetic[!Epigenetic %in% c("TemporalCortex_Wong_Dup15_mCpG", "Cerebellum_Wong_Dup15_mCpG",
                                            "PFC_VogelCiernia_Dup15_mCpG", "PFC_Wong_Dup15_mCpG", 
                                            "PFC_VogelCiernia_RTT_mCpG")]
Overlaps$EpigeneticStudies <- rowSums(Overlaps[,Epigenetic])
BrainEpigenetic <- Epigenetic[!Epigenetic %in% c("Blood_Andrews_mCpG", "Buccal_Berko_mCpG", "Placenta_Zhu_mCpG",
                                                 "Blood_Wong_mCpG", "Sperm_Feinberg_mCpG", "Bloodspot_Hannon_mCpG",
                                                 "LCL_Nguyen_mCpG")]
Overlaps$BrainEpigeneticStudies <- rowSums(Overlaps[,BrainEpigenetic])
Expression <- Expression[!Expression %in% c("CordBlood_Mordaunt_NonTDvsTD_DEG", "FrontalAndTemporalCortex_Lin_Rett_Vs_Ctrl_DGE",
                                            "FrontalAndTemporalCortex_Parikshak_Dup15qVs.Ctrl_Genes", 
                                            "ImprintedGenes_MaternalExpression", "ImprintedGenes_PaternalExpression")]
Overlaps$ExpressionStudies <- rowSums(Overlaps[,Expression])
Overlaps$AnyStudy <- rowSums(Overlaps[,c("GeneticStudies", "EpigeneticStudies", "ExpressionStudies")])
Overlaps <- Overlaps[order(Overlaps$DMRs, Overlaps$GeneSymbol),]
write.csv(Overlaps, "Tables/All DMR Gene Autism Gene Overlaps.csv", quote = FALSE, row.names = FALSE)

# Get Expression of Replicated ChrX DMR Genes in Fetal Brain ####
# Data
DMRbrainExp <- chrX_Overlaps[,c("GeneSymbol", "EntrezID", "chrX_DMRs")]
brainExp <- read.csv("Tables/Allen_developmental_transcriptome_expression_matrix.csv", header = FALSE, stringsAsFactors = FALSE)
brainGenes <- read.csv("Tables/Allen_developmental_transcriptome_rows_metadata.csv", header = TRUE, stringsAsFactors = FALSE)
brainSamples <- read.csv("Tables/Allen_developmental_transcriptome_columns_metadata.csv", header = TRUE, stringsAsFactors = FALSE)

# Subset Samples in Brain Expression Matrix
brainExp <- brainExp[,2:ncol(brainExp)]
sampleCols <- brainSamples$column_num[brainSamples$structure_acronym == "DFC" & brainSamples$age == "13 pcw"] # 88  96 107
brainSamples_sub <- brainSamples[sampleCols,]
#     column_num donor_id   donor_name    age gender structure_id structure_acronym                 structure_name
# 88          88    12820 H376.IIIA.50 13 pcw      M        10173               DFC dorsolateral prefrontal cortex
# 96          96    12834 H376.IIIA.51 13 pcw      F        10173               DFC dorsolateral prefrontal cortex
# 107        107    12888 H376.IIIA.52 13 pcw      M        10173               DFC dorsolateral prefrontal cortex
brainExp <- brainExp[,sampleCols]
colnames(brainExp) <- apply(brainSamples_sub[,c("gender", "structure_acronym", "age", "donor_id")], 1, paste, collapse = "_") %>%
        gsub(pattern = " ", replacement = "", fixed = TRUE) %>% paste(., "RPKM", sep = "_")
cor(brainExp$M_DFC_13pcw_12820_RPKM, brainExp$F_DFC_13pcw_12834_RPKM) # 0.9375689 1st male is most similar to female
cor(brainExp$M_DFC_13pcw_12820_RPKM, brainExp$M_DFC_13pcw_12888_RPKM) # 0.7725906 
cor(brainExp$F_DFC_13pcw_12834_RPKM, brainExp$M_DFC_13pcw_12888_RPKM) # 0.8378264 
brainExp <- brainExp[, c("M_DFC_13pcw_12820_RPKM", "F_DFC_13pcw_12834_RPKM")] # Exclude 2nd male
brainExp$EntrezID <- brainGenes$entrez_id

# Add Brain Expression for DMR genes by Entrez ID
DMRbrainExp <- merge(x = DMRbrainExp, y = brainExp, by = "EntrezID", all.x = TRUE, all.y = FALSE, sort = FALSE)
colnames(DMRbrainExp) <- c("EntrezID", "GeneSymbol", "chrX_DMRs", "Male_FetalBrainExp", "Female_FetalBrainExp")

# Fill in Missing by Gene Symbol
missing <- subset(DMRbrainExp, is.na(DMRbrainExp$Male_FetalBrainExp) | is.na(DMRbrainExp$Female_FetalBrainExp))
brainExp$GeneSymbol <- brainGenes$gene_symbol
missing <- merge(x = missing, y = brainExp, by = "GeneSymbol", all.x = TRUE, all.y = FALSE, sort = FALSE)
missing <- subset(missing, !is.na(missing$M_DFC_13pcw_12820_RPKM) & !is.na(missing$F_DFC_13pcw_12834_RPKM))
DMRbrainExp$Male_FetalBrainExp[match(missing$EntrezID.x, DMRbrainExp$EntrezID)] <- missing$M_DFC_13pcw_12820_RPKM
DMRbrainExp$Female_FetalBrainExp[match(missing$EntrezID.x, DMRbrainExp$EntrezID)] <- missing$F_DFC_13pcw_12834_RPKM
table(is.na(DMRbrainExp$Male_FetalBrainExp) | is.na(DMRbrainExp$Female_FetalBrainExp)) # 33 missing

# Write Table
DMRbrainExp <- DMRbrainExp[order(DMRbrainExp$GeneSymbol),]
DMRbrainExp <- DMRbrainExp[,c("GeneSymbol", "EntrezID", "chrX_DMRs", "Male_FetalBrainExp", "Female_FetalBrainExp")]
write.csv(DMRbrainExp, "Tables/ChrX DMR Gene Fetal Brain 13pcw Expression RPKM.csv", quote = FALSE, row.names = FALSE)



