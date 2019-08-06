# Autism Gene Overlap -----------------------------------------------------
# Autism Cord Blood Methylation
# Charles Mordaunt
# 8/5/19

# Packages ####
sapply(c("reshape2", "tidyverse", "GenomicRanges", "annotatr", "GeneOverlap", "rlist", "biomaRt", "rtracklayer"), require, character.only = TRUE)

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

# Write Files
names(DMRs_entrez) <- c("Males_Discovery", "Females_Discovery", "Males_Replication", "Females_Replication")
mapply(write.table, x = DMRs_entrez, file = paste("Tables/EntrezIDs_", names(DMRs_entrez), "_DMRs.txt", sep = ""), 
       MoreArgs = list(sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE))

names(Background_entrez) <- c("Males_Discovery", "Females_Discovery", "Males_Replication", "Females_Replication")
mapply(write.table, x = Background_entrez, file = paste("Tables/EntrezIDs_", names(Background_entrez), "_Background.txt", sep = ""), 
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

# Convert ENSGs to Entrez IDs Overlaps from Vogel Ciernia et al 2019 ---------------
# Read in Files
files <- list.files("ENSG ID Lists")
files <- files[!grepl("CpG", files, fixed = TRUE) & !grepl("CpH", files, fixed = TRUE) & !grepl("H3K", files, fixed = TRUE)]
files <- paste("ENSG ID Lists/", files, sep = "")
ENSG <- lapply(files, read.delim, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
ENSG <- lapply(ENSG, unlist, use.names = TRUE)

# Convert to Entrez IDs
ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
entrezIDs <- lapply(ENSG, getBM, attributes = "entrezgene_id", filters = "ensembl_gene_id", mart = ensembl)
entrezIDs <- lapply(entrezIDs, unlist, use.names = TRUE)

# Write Files
newFiles <- gsub("ENSG ID Lists/ENSGs", replacement = "Entrez ID Lists/EntrezIDs", x = files, fixed = TRUE)
mapply(write.table, x = entrezIDs, file = newFiles, MoreArgs = list(sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE))

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
                      "Tables/Wong Temporal Cortex DMPs hg19.csv", "Tables/Wong BA9 DMPs hg19.csv"),
                    read.csv, header = TRUE, stringsAsFactors = FALSE)
names(DMRs_hg19) <- c("CingulateCortex_Nardone", "PFC_Nardone", "PFCneurons_Shulha_H3K4me3", "TemporalCortex_LaddAcosta",
                      "Cerebellum_LaddAcosta", "Buccal_Berko", "Sperm_Feinberg", "Cerebellum_Sun_H3K27Ac", 
                      "TemporalCortex_Sun_H3K27Ac", "PFC_Sun_H3K27Ac", "VisualCortex_Ellis_CpH", "VisualCortex_Ellis", 
                      "PFCneurons_Nardone", "Bloodspot_Hannon", "Blood_Andrews", "Cerebellum_Wong", "TemporalCortex_Wong", 
                      "PFC_Wong")
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
names(entrezIDs) <- c("VisualCortex_Ellis_mCpH", "VisualCortex_Ellis_mCpG", "TemporalCortex_Wong_mCpG", 
                  "TemporalCortex_Sun_H3K27Ac", "TemporalCortex_LaddAcosta_mCpG", "Sperm_Feinberg_mCpG", "Placenta_Zhu_mCpG",
                  "PFCneurons_Shulha_H3K4me3", "PFCneurons_Nardone_mCpG", "PFC_Wong_mCpG", "PFC_VogelCiernia_RTT_mCpG",
                  "PFC_VogelCiernia_Dup15_mCpG", "PFC_VogelCiernia_ASD_mCpG", "PFC_Sun_H3K27Ac", "PFC_Nardone_mCpG", 
                  "LCL_Nguyen_mCpG", "CingulateCortex_Nardone_mCpG", "Cerebellum_Wong_mCpG", "Cerebellum_Sun_H3K27Ac", 
                  "Cerebellum_LaddAcosta_mCpG", "Buccal_Berko_mCpG", "Bloodspot_Hannon_mCpG", "Blood_Wong_mCpG", 
                  "Blood_Andrews_mCpG")
mapply(write.table, x = entrezIDs, file = paste("Tables/EntrezIDs_", names(entrezIDs), ".txt", sep = ""), 
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

# Gene Lists to Add ####
#' Julia
#' SFARI (Annie's table)
#' Gandal et al 2018 DEGs and DTEs (Annie's table)
#' Tylee et al 2017 Blood
#' Tylee et al 2017 LCL
#' Maternal and Paternal imprinted genes (geneimprint.com)
#' Parikshak 2016 Dup15 > Ctrl
#' Lin et al 2016 Rett DEGs
#' Doan et al 2019
#' Grove et al 2019
#' Other ASD GWAS?

#' Me
#' Rett and Dup15 DMR genes from Annie's paper
#' Parikshak 2016
#' Other Dup15 EWAS 
#' ASD mQTL targets
#' Control gene list?











