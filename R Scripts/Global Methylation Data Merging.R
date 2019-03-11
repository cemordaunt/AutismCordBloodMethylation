# Global Methylation Data Merging ####
# MARBLES WGS
# Charles Mordaunt
# 11/29/18

setwd("~/Desktop/Cord Methylation")
library(dplyr)

# Data ####
samples <- read.delim("MARBLES WGS Sample Database.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
placenta <- read.csv("MARBLES_Glo_92_COI.csv", header = TRUE, stringsAsFactors = FALSE)
marbles1 <- read.delim("MARBLES 1 QC/multiqc_bismark_methextract.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
marbles2 <- read.delim("MARBLES 2 QC/multiqc_bismark_methextract.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
marblesWGBS <- read.delim("MARBLES EARLI WGBS Sample Merged Covariate Database with Demo.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Merge ####
# Placenta
placenta$placenta_CG_meth <- placenta$meth_cpg * 100 / (placenta$meth_cpg + placenta$unmeth_cpg)
placenta <- placenta[,c("COI_ID", "placenta_CG_meth")]
samplesMeth <- merge(x = samples, y = placenta, by = "COI_ID", all.x = TRUE, all.y = FALSE, sort = FALSE)

# Cord
marbles1samples <- strsplit(marbles1$Sample, "_", fixed = TRUE)
marbles1samples <- lapply(marbles1samples, function(x) x[1]) %>% as.character
marbles1$SeqID <- marbles1samples
marbles1$cord_CG_meth <- marbles1$meth_cpg * 100 / (marbles1$meth_cpg + marbles1$unmeth_cpg)
marbles1$cord_platform <- "HiSeq4000"
marbles1$cord_platform[marbles1$SeqID == "JLDS061B"] <- "HiSeq2500"
marbles1 <- marbles1[,c("SeqID", "cord_CG_meth", "cord_platform")]

marbles2samples <- strsplit(marbles2$Sample, "_", fixed = TRUE)
marbles2samples <- lapply(marbles2samples, function(x) x[1]) %>% as.character
marbles2$SeqID <- marbles2samples
marbles2$cord_CG_meth <- marbles2$meth_cpg * 100 / (marbles2$meth_cpg + marbles2$unmeth_cpg)
marbles2$cord_platform <- "HiSeqX10"
marbles2$cord_CHH_meth <- marbles2$meth_chh * 100 / (marbles2$meth_chh + marbles2$unmeth_chh) # JLCM064B has high CH methylation
marbles2 <- marbles2[,c("SeqID", "cord_CG_meth", "cord_platform")]
marbles <- rbind(marbles1, marbles2)

rm(marbles1, marbles2, placenta, samples, marbles1samples, marbles2samples)

marblesWGBS <- marblesWGBS[,c("COI_ID", "Sequencing_ID")]
marbles <- merge(x = marbles, y = marblesWGBS, by.x = "SeqID", by.y = "Sequencing_ID", all.x = TRUE, all.y = FALSE, sort = FALSE)
marbles[is.na(marbles$COI_ID),]

#        SeqID cord_CG_meth cord_platform COI_ID
# 3   JLCM001C     74.84412     HiSeq4000     NA
# 23  JLCM007C     75.92681     HiSeq4000     NA
# 32  JLCM009D     74.79324     HiSeq4000     NA
# 36  JLCM010D     74.73346     HiSeq4000     NA
# 43  JLCM012C     76.25245     HiSeq4000     NA
# 86  JLCM064B     81.17858      HiSeqX10     NA high CH Methylation
# 97  JLCM070A     78.54736      HiSeqX10     NA T cells
# 98  JLCM070B     76.15298      HiSeqX10     NA T cells
# 99  JLCM070C     80.30166      HiSeqX10     NA T cells
# 100 JLCM070D     75.63058      HiSeqX10     NA T cells

marbles$COI_ID[marbles$SeqID == "JLCM001C"] <- 508006
marbles$COI_ID[marbles$SeqID == "JLCM007C"] <- 505706
marbles$COI_ID[marbles$SeqID == "JLCM009D"] <- 508707
marbles$COI_ID[marbles$SeqID == "JLCM010D"] <- 502305
marbles$COI_ID[marbles$SeqID == "JLCM012C"] <- 500906

marbles$COI_ID[duplicated(marbles$COI_ID)] #503505
marbles[marbles$COI_ID == 503505,]
#         SeqID cord_CG_meth cord_platform COI_ID
# 18   JLCM006C     76.65968     HiSeq4000 503505
# 48   JLDS061B     75.34537     HiSeq2500 503505
# 91   JLCM069B     77.04866      HiSeqX10 503505

marbles <- marbles[!marbles$SeqID %in% c("JLCM006C", "JLDS061B"),] # remove duplicates
marbles <- marbles[,c("COI_ID", "cord_CG_meth", "cord_platform")]
samplesMeth <- merge(x = samplesMeth, y = marbles, by = "COI_ID", all.x = TRUE, all.y = FALSE, sort = FALSE)
write.table(samplesMeth, "MARBLES WGS Sample Database with Global Methyl.txt", sep = "\t", quote = FALSE, row.names = FALSE)









