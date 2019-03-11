# EARLI MARBLES Merged Sample Database ####
# Autism Cord Blood Methylation
# Charles Mordaunt
# 10/2/18

# Packages ####
library(measurements)

# Data ####
earli <- read.delim(file = "Samples/EARLI_outcomes.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
marbles <- read.delim(file = "Samples/MARBLES_WGBS_database_Cov.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
merged <- read.delim(file = "Samples/Merged Cord Blood WGBS Database.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
dupCols <- read.csv(file = "Samples/duplicateCols.csv", header = TRUE, stringsAsFactors = FALSE)
cleanDupCols <- read.csv(file = "Samples/cleanDuplicateCols.csv", header = TRUE, stringsAsFactors = FALSE)

# Remove duplicate columns ####
dupColNames <- colnames(dupCols)
table(dupColNames %in% colnames(earli))
dupColNames[!dupColNames %in% colnames(earli)]
# "duplicateCols"  "BMI_cat.x.1"    "BMI_cat.y.1"    "BMI_cat.1"      "BMI_source.x.1" "BMI_source.y.1" "BMI_source.1"  
earli <- earli[,!colnames(earli) %in% dupColNames]

# Add in cleaned columns ####
earli <- merge(x = earli, y = cleanDupCols, by.x = "Subject.ID", by.y = "duplicateCols", all.x = TRUE, all.y = TRUE, sort = FALSE)
table(duplicated(colnames(earli))) # All FALSE
table(grepl(".x", colnames(earli), fixed = TRUE)) # 1 TRUE
table(grepl(".y", colnames(earli), fixed = TRUE)) # All FALSE
colnames(earli)[grepl(".x", colnames(earli), fixed = TRUE)] # "maternal_race.x"
colnames(earli)[grepl("maternal_race", colnames(earli), fixed = TRUE)] # "maternal_race.x"
colnames(earli)[colnames(earli) == "maternal_race.x"] <- "maternal_race"
write.table(earli, "Samples/EARLI_Outcome_DHQ_Supplement_Database.txt", sep = "\t", quote = FALSE, row.names = FALSE)
rm(cleanDupCols, dupCols, dupColNames)

# Add new variables to merged ####
# Variables in marbles and earli, but not merged
marblesCols <- sort(colnames(marbles))
earliCols <- sort(colnames(earli))
marbearli <- intersect(marblesCols, earliCols) # None
write.table(marblesCols, "Samples/MARBLES Covariate Names.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(earliCols, "Samples/EARLI Covariate Names.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
matched <- read.delim("Samples/MARBLES EARLI Matched Covariate Names.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
matchedMarbles <- marbles[,as.character(matched$marbles)]
matchedEarli <- earli[,as.character(matched$earli)]
matchedEarli$height_inches <- conv_unit(matchedEarli$height_inches, from = "inch", to = "cm")
matchedEarli$prepreg_weight <- conv_unit(matchedEarli$prepreg_weight, from = "lbs", to = "kg")
colnames(matchedEarli) <- colnames(matchedMarbles)
matchedDf <- rbind(matchedMarbles, matchedEarli)
filterMerged <- colnames(matchedMarbles)[!colnames(matchedMarbles) == "Cord_Blood_IBC"]
merged <- merged[,!colnames(merged) %in% filterMerged]
merged <- merge(x = merged, y = matchedDf, by = "Cord_Blood_IBC", all.x = TRUE, all.y = FALSE, sort = FALSE)
merged <- merged[,c("Cord_Blood_IBC", "Sequencing_ID", "COI_ID", "Study", "Platform" , "Diagnosis_Alg", "Sex", "Site",
                    "GDM_EQ", "PE_EQ", "SmokeYN_Pregnancy", "Supp_mv_mo_1", "Supp_mv_mo1", "MSLvrTscore36", "MSLfmTscore36",
                    "MSLrlTscore36", "MSLelTscore36", "MSLelcStandard36", "ADOScs", "percent_trimmed", "percent_aligned",
                    "percent_duplicate", "dedup_reads", "C_coverage", "CG_coverage", "percent_chg_meth", "percent_chh_meth",
                    "percent_cpg_meth", "MomAgeYr", "Mat_Height_cm", "Mat_Weight_kg_PrePreg", "Mat_BMI_PrePreg")]              
merged <- merged[!duplicated(merged) == TRUE,]
write.table(merged, file = "Samples/Merged Cord Blood WGBS Database.txt", sep = "\t", quote = FALSE, row.names = FALSE)


