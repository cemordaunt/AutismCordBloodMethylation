# EARLI-MARBLES Covariate Database Merging ####
# ASD Cord Blood Methylation
# Charles Mordaunt
# 10/30/18

# Packages ####
library(sas7bdat)
library(dplyr)
library(measurements)

# Plan ####
# Start with WGBS samples with QC
# Add in demographic data for each study, then merge together

# WGBS Samples ####
samples <- read.delim(file = "Merged Database/Cord Blood WGBS Samples.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
samplesMarbles <- subset(samples, Study == "MARBLES")
samplesEarli <- subset(samples, Study == "EARLI")

# Marbles Covariate Database ####
# Data
marblesCotinine <- read.csv(file = "Merged Database/MARBLES Source Files/CMordaunt_add creat cot 07312018.csv", header = TRUE, stringsAsFactors = FALSE)
marblesCovars1 <- read.csv(file = "Merged Database/MARBLES Source Files/2018-04-12 MARBLES Covariates.csv", header = TRUE, stringsAsFactors = FALSE)
marblesCovars2 <- read.csv(file = "Merged Database/MARBLES Source Files/2018-05-29 Additional MARBLES Covariates.csv", header = TRUE, stringsAsFactors = FALSE)
marblesCovars3 <- read.sas7bdat(file = "Merged Database/MARBLES Source Files/covars_12oct18.sas7bdat")

# Merge Cotinine
marblesCotinine <- marblesCotinine[!is.na(marblesCotinine$COI_ID), c("COI_ID", "Final_Conc_mg_dl", "Conc_in_extract_ng_ml", "Conc_in_urine_ng_ml")]
colnames(marblesCotinine) <- c("COI_ID", "Creatine_mg_dl", "Cotinine_extract_ng_ml", "Cotinine_urine_ng_ml")
samplesMarbles <- merge(x = samplesMarbles, y = marblesCotinine, by = "COI_ID", all.x = TRUE, all.y = FALSE, sort = FALSE)
rm(marblesCotinine)

# Merge Covars Together
marblesCovars1colnames <- as.character(colnames(marblesCovars1)) # 66
marblesCovars2colnames <- as.character(colnames(marblesCovars2)) # 16
marblesCovars3colnames <- as.character(colnames(marblesCovars3)) # 3391

marblesCovars1and2 <- intersect(marblesCovars1colnames, marblesCovars2colnames) # 1 COI_ID
marblesCovars1and3 <- intersect(marblesCovars1colnames, marblesCovars3colnames) # 59 in common
marblesCovars2and3 <- intersect(marblesCovars2colnames, marblesCovars3colnames) # 2 COI_ID, ChildRaceEth

marblesCovars1and3 <- marblesCovars1and3[!marblesCovars1and3 == "COI_ID"]
marblesCovars1 <- marblesCovars1[,!colnames(marblesCovars1) %in% marblesCovars1and3] # Remove variables also in Covars3
marblesCovars2 <- marblesCovars2[,!colnames(marblesCovars2) == "ChildRaceEth"] # Remove variables also in Covars3

marblesCovars <- merge(x = marblesCovars1, y = marblesCovars2, by = "COI_ID", all.x = TRUE, all.y = TRUE, sort = FALSE)
marblesCovars <- merge(x = marblesCovars, y = marblesCovars3, by = "COI_ID", all.x = TRUE, all.y = TRUE, sort = FALSE)
write.table(marblesCovars, file = "Merged Database/MARBLES Source Files/MARBLES Covariates Merged.txt", sep = "\t", quote = FALSE, row.names = FALSE)
rm(marblesCovars1, marblesCovars2, marblesCovars3, marblesCovars1and2, marblesCovars1and3, marblesCovars2and3, marblesCovars1colnames, marblesCovars2colnames, marblesCovars3colnames)

# Subset Covar Table for WGBS samples
marblesCovarsWGBS <- subset(marblesCovars, COI_ID %in% samplesMarbles$COI_ID)
rm(marblesCovars)
nacount <- as.data.frame(t(sapply(marblesCovarsWGBS, function(x){c(table(is.na(x))["TRUE"], table(is.na(x))["FALSE"])})))
colnames(nacount) <- c("TRUE", "FALSE")
nacount$Variable <- rownames(nacount)
rownames(nacount) <- 1:nrow(nacount)
nacount <- nacount[,c("Variable", "TRUE", "FALSE")]
nacount[is.na(nacount) == TRUE] <- 0
lowdata <- nacount[nacount$'TRUE' > 60,"Variable"] # Less then 1/3 of MARBLES samples have data
# [1] "birthheadcir"      "birth_breathing"   "birth_NUCHAL"      "birth_NUCHALtimes" "Birth_Resc"        "Apgar10"          
# [7] "GA_mo_EQ3" 
marblesCovarsWGBS <- marblesCovarsWGBS[,!colnames(marblesCovarsWGBS) %in% lowdata] # remove columns if < 1/3 samples have data
write.table(marblesCovarsWGBS, file = "Merged Database/MARBLES Source Files/MARBLES Covariates Merged WGBS Samples.txt", sep = "\t", quote = FALSE, row.names = FALSE)
rm(nacount, lowdata)

# Subset Covar Table for Selected Variables
marblesSelectedCovarNames <- read.delim(file = "Merged Database/MARBLES Source Files/MARBLES Selected Variables.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
marblesSelectedCovarNames <- as.character(unlist(marblesSelectedCovarNames))
marblesCovarsSelect <- marblesCovarsWGBS[,colnames(marblesCovarsWGBS) %in% marblesSelectedCovarNames]
rm(marblesCovarsWGBS, marblesSelectedCovarNames)

# Remove Duplicate Columns
# "DM1_EQ"      "DM1"                  
# "DM2_EQ"      "DM2" 
# "GDM_EQ"      "GDM"             
# "HTN_EQ"      "HTN"
# "PE_EQ"       "PE" 
# "MMC_fine_EQ" "MMC_fine"
# "MMC_5cat_EQ" "MMC_5cat" 
# "Mat_BMI_PrePreg"     "PrePregBMI_EQ"          
# "ChildRace"   "ChildRaceEth" 
# "MomEdu"      "MomEdu_detail" "MomEdu_6cat"   "MomEdu_7cat"
# "DadEdu_detail"       "DadEdu"        "DadEdu_6cat"   "DadEdu_7cat"
# "N_pregs"     "parity"        "IPI_parity_mos"           

table(marblesCovarsSelect$DM1_EQ == marblesCovarsSelect$DM1) #TRUE 89
table(is.na(marblesCovarsSelect$DM1_EQ)) #2 NA
table(is.na(marblesCovarsSelect$DM1)) #0 NA
table(marblesCovarsSelect$DM2_EQ == marblesCovarsSelect$DM2) #TRUE 89
table(is.na(marblesCovarsSelect$DM2_EQ)) #2 NA
table(is.na(marblesCovarsSelect$DM2)) #0 NA
table(marblesCovarsSelect$GDM_EQ == marblesCovarsSelect$GDM) #TRUE 86 FALSE 1
table(is.na(marblesCovarsSelect$GDM_EQ)) #4 NA
table(is.na(marblesCovarsSelect$GDM)) #0 NA
table(marblesCovarsSelect$HTN_EQ == marblesCovarsSelect$HTN) #TRUE 86 FALSE 3
table(is.na(marblesCovarsSelect$HTN_EQ)) #2 NA
table(is.na(marblesCovarsSelect$HTN)) #0 NA
table(marblesCovarsSelect$PE_EQ == marblesCovarsSelect$PE) #TRUE 85 FALSE 1
table(is.na(marblesCovarsSelect$PE_EQ)) #5 NA
table(is.na(marblesCovarsSelect$PE)) #1 NA
table(marblesCovarsSelect$MMC_fine_EQ == marblesCovarsSelect$MMC_fine) #TRUE 76 FALSE 8
table(is.na(marblesCovarsSelect$MMC_fine_EQ)) #7 NA
table(is.na(marblesCovarsSelect$MMC_fine)) #1 NA
table(marblesCovarsSelect$MMC_5cat_EQ == marblesCovarsSelect$MMC_5cat) #TRUE 79 FALSE 5
table(is.na(marblesCovarsSelect$MMC_5cat_EQ)) #7 NA
table(is.na(marblesCovarsSelect$MMC_5cat)) #1 NA
table(marblesCovarsSelect$PrePregBMI_EQ == marblesCovarsSelect$Mat_BMI_PrePreg) #TRUE 91
table(is.na(marblesCovarsSelect$PrePregBMI_EQ)) #0 NA
table(is.na(marblesCovarsSelect$Mat_BMI_PrePreg)) #0 NA

table(marblesCovarsSelect$ChildRace)
#  1  2  4  9 
# 61  4  8 18 
table(marblesCovarsSelect$ChildRaceEth)
#  1  2  4  9 61 62 71 72 
# 36  3  8  9 14 11  5  5 

table(marblesCovarsSelect$MomEdu)
# 1  2  3  4  5 
# 2  4 41 32 12 
table(marblesCovarsSelect$MomEdu_detail)
# 2  3  4  5  6  7  8 
# 2  4 18 23 32 10  2 
table(marblesCovarsSelect$MomEdu_6cat)
# 1  2  3  4  5  6 
# 2  4 18 23 32 12 
table(marblesCovarsSelect$MomEdu_7cat)
# 2  3  4  5  6  7 
# 2  4 41 32 10  2 

table(marblesCovarsSelect$N_pregs)
# 0  1  2  3  4  5  6  7  8 10 
# 1 13 31 18 12  8  3  1  1  1 
table(marblesCovarsSelect$parity)
# 0  1  2  3  4  5  7 
# 1 34 38  6  6  3  1 
summary(marblesCovarsSelect$IPI_parity_mos)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 3.253  26.407  46.965  50.208  61.906 150.374       3

removeCols <- c("DM1_EQ", "DM2_EQ", "GDM_EQ", "HTN_EQ", "PE_EQ", "MMC_fine_EQ", "MMC_5cat_EQ", "PrePregBMI_EQ")
marblesCovarsSelect <- marblesCovarsSelect[,!colnames(marblesCovarsSelect) %in% removeCols]
rm(removeCols)
marblesDatabase <- merge(x = samplesMarbles, y = marblesCovarsSelect, by = "COI_ID", all.x = TRUE, all.y = TRUE, sort = FALSE)
marblesDatabase <- unique(marblesDatabase)
write.table(marblesDatabase, file = "Merged Database/MARBLES WGBS Sample Covariate Database.txt", sep = "\t", row.names = FALSE,
            quote = FALSE)
rm(marblesCovarsSelect, samplesMarbles)

# EARLI Covariate Database ####
# Data
earliIDs <- read.csv(file = "Merged Database/EARLI Source Files/EARLI_ids.txt", header = TRUE, stringsAsFactors = FALSE)
earliOutcomes <- read.delim(file = "Merged Database/EARLI Source Files/earli_outcomes.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
earliSupp <- read.csv(file = "Merged Database/EARLI Source Files/final_earli_supplements_txt.txt", header = TRUE, stringsAsFactors = FALSE)
earliDHQ1 <- read.csv(file = "Merged Database/EARLI Source Files/dhq1_nutrients_common.txt", header = TRUE, stringsAsFactors = FALSE)
earliDHQ2 <- read.csv(file = "Merged Database/EARLI Source Files/dhq2_nutrients_common.txt", header = TRUE, stringsAsFactors = FALSE)

# Duplicate Columns
earliOutcomesCols <- colnames(earliOutcomes) #25
earliSuppCols <- colnames(earliSupp) #3584
earliDHQ1cols <- colnames(earliDHQ1) #106
earliDHQ2cols <- colnames(earliDHQ2) #106
intersect(earliOutcomesCols, earliSuppCols) %>% length #0
intersect(earliOutcomesCols, earliDHQ1cols) %>% length #0
intersect(earliOutcomesCols, earliDHQ2cols) %>% length #0
intersect(earliSuppCols, earliDHQ1cols) %>% length #16
# [1] "subjectid"          "income"             "maternal_race"      "maternal_ethnicity" "mom_edu"            "maternal_age"      
# [7] "prepreg_weight"     "height_inches"      "BMI"                "BMI_cat"            "BMI_source"         "GA_mateducation"   
# [13] "gest_diab"          "preeclampsia"       "Chron_Diabetes"     "prenatal_smoking"  
intersect(earliSuppCols, earliDHQ2cols) %>% length #16
# [1] "subjectid"          "income"             "maternal_race"      "maternal_ethnicity" "mom_edu"            "maternal_age"      
# [7] "prepreg_weight"     "height_inches"      "BMI"                "BMI_cat"            "BMI_source"         "GA_mateducation"   
# [13] "gest_diab"          "preeclampsia"       "Chron_Diabetes"     "prenatal_smoking"  
intersect(earliDHQ1cols, earliDHQ2cols) %>% length #19
# [1] "subjectid"          "income"             "maternal_race"      "maternal_ethnicity" "mom_edu"            "maternal_age"      
# [7] "prepreg_weight"     "height_inches"      "BMI"                "BMI_cat"            "BMI_source"         "GA_mateducation"   
# [13] "gest_diab"          "preeclampsia"       "Chron_Diabetes"     "prenatal_smoking"   "gestage"            "siteid"            
# [19] "singl"  
intersect(earliSuppCols, earliDHQ1cols) %>% intersect(.,earliDHQ2cols) %>% length #16
# [1] "subjectid"          "income"             "maternal_race"      "maternal_ethnicity" "mom_edu"            "maternal_age"      
# [7] "prepreg_weight"     "height_inches"      "BMI"                "BMI_cat"            "BMI_source"         "GA_mateducation"   
# [13] "gest_diab"          "preeclampsia"       "Chron_Diabetes"     "prenatal_smoking"

# earliOutcome doesn't contain any duplicate columns
# 16 columns are duplicated between earliSupp, earliDHQ1, and earliDHQ2 (including subjectid)
# An additional 3 columns are duplicated between earliDHQ1 and earliDHQ2

# Add earliIDs and earliOutcomes
earliDatabase <- merge(x = samplesEarli, y = earliIDs, by.x = "COI_ID", by.y = "child_id", all.x = TRUE, all.y = FALSE, sort = FALSE)
earliDatabase <- merge(x = earliDatabase, y = earliOutcomes, by.x = "COI_ID", by.y = "Subject_ID", all.x = TRUE, all.y = FALSE, sort = FALSE)
earliDatabase <- earliDatabase[order(earliDatabase$mother_id),]
rm(earliIDs, earliOutcomes, samplesEarli, earliOutcomesCols)

# Subset Supp and DHQs for WGBS Samples
earliSupp <- subset(earliSupp, subjectid %in% earliDatabase$mother_id) #66
earliSupp <- earliSupp[order(earliSupp$subjectid),]
earliDHQ1 <- subset(earliDHQ1, subjectid %in% earliDatabase$mother_id) #53
earliDHQ1 <- earliDHQ1[order(earliDHQ1$subjectid),]
earliDHQ2 <- subset(earliDHQ2, subjectid %in% earliDatabase$mother_id) #33
earliDHQ2 <- earliDHQ2[order(earliDHQ2$subjectid),]

# Remove Duplicate Columns in DHQs
DHQ1andDHQ2cols <- intersect(earliDHQ1cols, earliDHQ2cols) #19
DHQandSuppCols <- intersect(DHQ1andDHQ2cols, earliSuppCols) #16
DHQonlyCols <- DHQ1andDHQ2cols[!DHQ1andDHQ2cols %in% DHQandSuppCols] #3
earliDHQs <- merge(x = earliDHQ1, y = earliDHQ2, by = "subjectid", all.x = TRUE, all.y = TRUE, sort = FALSE) #56
DHQ1andDHQ2cols <- DHQ1andDHQ2cols[!DHQ1andDHQ2cols == "subjectid"]
sapply(DHQ1andDHQ2cols, function(x) c("SameValue" = table(earliDHQs[,paste(x,".x", sep = "")] == earliDHQs[,paste(x,".y", sep = "")])["TRUE"],
                                      "DiffValue" = table(earliDHQs[,paste(x,".x", sep = "")] == earliDHQs[,paste(x,".y", sep = "")])["FALSE"],
                                      "NA1" = table(is.na(earliDHQs[,paste(x,".x", sep = "")]))["TRUE"],
                                      "NA2" = table(is.na(earliDHQs[,paste(x,".y", sep = "")]))["TRUE"],
                                      "NA1and2" = table(is.na(earliDHQs[,paste(x,".x", sep = "")]) & 
                                                                is.na(earliDHQs[,paste(x,".y", sep = "")]))["TRUE"]))
#                income maternal_race maternal_ethnicity mom_edu maternal_age prepreg_weight height_inches BMI BMI_cat BMI_source
# SameValue.TRUE     29            30                 30      29           30             29            29  29      30         30
# DiffValue.NA       NA            NA                 NA      NA           NA             NA            NA  NA      NA         NA
# NA1.TRUE            5             3                  3       4            3              4             4   4       3          3
# NA2.TRUE           24            23                 23      24           23             24            24  24      23         23
# NA1and2.TRUE        2            NA                 NA       1           NA              1             1   1      NA         NA

#                GA_mateducation gest_diab preeclampsia Chron_Diabetes prenatal_smoking gestage siteid singl
# SameValue.TRUE              29        29           29             29               28       1     30    30
# DiffValue.NA                NA        NA           NA             NA               NA      29     NA    NA
# NA1.TRUE                     4         4            4              4               14       3      3     3
# NA2.TRUE                    24        24           24             24               25      23     23    23
# NA1and2.TRUE                 1         1            1              1               11      NA     NA    NA

# Only gestage has different values for DHQ1 and 2, the rest are the same or missing.
# For columns not gestage, merge DHQ1 and 2 values together
# For gestage, keep separate and add DHQ1/2suffix
DHQ1andDHQ2cols <- DHQ1andDHQ2cols[!DHQ1andDHQ2cols == "gestage"]
earliDHQsNoDup <- sapply(DHQ1andDHQ2cols, function(x){
        combined <- earliDHQs[,paste(x, ".x", sep = "")]
        combined[is.na(combined)] <- earliDHQs[is.na(combined),paste(x, ".y", sep = "")]
        return(combined)
})
DHQ1andDHQ2colFilter <- paste(rep(DHQ1andDHQ2cols, each = 2), c(".x", ".y"), sep = "")
earliDHQsDupFiltered <- earliDHQs[,!colnames(earliDHQs) %in% DHQ1andDHQ2colFilter]
earliDHQsFinal <- cbind(earliDHQsDupFiltered, earliDHQsNoDup) #same order
colnames(earliDHQsFinal)[colnames(earliDHQsFinal) == "gestage.x"] <- "DHQ1_gestage"
colnames(earliDHQsFinal)[colnames(earliDHQsFinal) == "gestage.y"] <- "DHQ2_gestage"
rm(earliDHQ1, earliDHQ2, earliDHQs, earliDHQsDupFiltered, earliDHQsNoDup, DHQ1andDHQ2colFilter, DHQ1andDHQ2cols, 
   earliDHQ1cols, earliDHQ2cols, DHQandSuppCols, DHQonlyCols)

# Remove Duplicate Columns between Supp and DHQs
earliDHQsupp <- merge(x = earliDHQsFinal, y = earliSupp, by = "subjectid", all.x = TRUE, all.y = TRUE)
earliDHQsFinalCols <- colnames(earliDHQsFinal) #194
earliDHQandSuppCols <- intersect(earliDHQsFinalCols, earliSuppCols) #16
# [1]  "subjectid"          "income"             "maternal_race"      "maternal_ethnicity" "mom_edu"            "maternal_age"      
# [7]  "prepreg_weight"     "height_inches"      "BMI"                "BMI_cat"            "BMI_source"         "GA_mateducation"   
# [13] "gest_diab"          "preeclampsia"       "Chron_Diabetes"     "prenatal_smoking"
earliDHQandSuppCols <- earliDHQandSuppCols[!earliDHQandSuppCols == "subjectid"]
sapply(earliDHQandSuppCols, function(x) c("SameValue" = table(earliDHQsupp[,paste(x,".x", sep = "")] == earliDHQsupp[,paste(x,".y", sep = "")])["TRUE"],
                                      "DiffValue" = table(earliDHQsupp[,paste(x,".x", sep = "")] == earliDHQsupp[,paste(x,".y", sep = "")])["FALSE"],
                                      "NA1" = table(is.na(earliDHQsupp[,paste(x,".x", sep = "")]))["TRUE"],
                                      "NA2" = table(is.na(earliDHQsupp[,paste(x,".y", sep = "")]))["TRUE"],
                                      "NA1and2" = table(is.na(earliDHQsupp[,paste(x,".x", sep = "")]) & 
                                                                is.na(earliDHQsupp[,paste(x,".y", sep = "")]))["TRUE"]))
#                income maternal_race maternal_ethnicity mom_edu maternal_age prepreg_weight height_inches BMI BMI_cat BMI_source
# SameValue.TRUE     54            56                 56      55           56             55            55   4      56         56
# DiffValue.NA       NA            NA                 NA      NA           NA             NA            NA  51      NA         NA
# NA1.TRUE           12            10                 10      11           10             11            11  11      10         10
# NA2.TRUE            4            NA                 NA       1           NA              1             1   1      NA         NA
# NA1and2.TRUE        4            NA                 NA       1           NA              1             1   1      NA         NA

#                GA_mateducation gest_diab preeclampsia Chron_Diabetes prenatal_smoking
# SameValue.TRUE              55        55           55             55               45
# DiffValue.NA                NA        NA           NA             NA               NA
# NA1.TRUE                    11        11           11             11               21
# NA2.TRUE                     1         1            1              1               15
# NA1and2.TRUE                 1         1            1              1               15

# Only BMI is different, keep separate, rest are NA
earliDHQandSuppCols <- earliDHQandSuppCols[!earliDHQandSuppCols == "BMI"]
earliDHQsSuppNoDup <- sapply(earliDHQandSuppCols, function(x){
        combined <- as.character(earliDHQsupp[,paste(x, ".x", sep = "")])
        combined[is.na(combined)] <- as.character(earliDHQsupp[is.na(combined),paste(x, ".y", sep = "")])
        return(combined)
})
DHQsupColFilter <- paste(rep(earliDHQandSuppCols, each = 2), c(".x", ".y"), sep = "")
earliDHQsuppDupFiltered <- earliDHQsupp[,!colnames(earliDHQsupp) %in% DHQsupColFilter]
earliDHQsuppFinal <- cbind(earliDHQsuppDupFiltered, earliDHQsSuppNoDup) #same order
colnames(earliDHQsuppFinal)[colnames(earliDHQsuppFinal) == "BMI.x"] <- "DHQ_BMI"
colnames(earliDHQsuppFinal)[colnames(earliDHQsuppFinal) == "BMI.y"] <- "Supp_BMI"
rm(earliDHQsFinal, earliDHQsSuppNoDup, earliDHQsupp, earliDHQsuppDupFiltered, earliSupp, DHQsupColFilter, 
   earliDHQandSuppCols, earliDHQsFinalCols, earliSuppCols)

# Remove Variables with Data in < 1/3 of samples
nacount <- as.data.frame(t(sapply(earliDHQsuppFinal, function(x){c(table(is.na(x))["TRUE"], table(is.na(x))["FALSE"])})))
colnames(nacount) <- c("TRUE", "FALSE")
nacount$Variable <- rownames(nacount)
rownames(nacount) <- 1:nrow(nacount)
nacount <- nacount[,c("Variable", "TRUE", "FALSE")]
nacount[is.na(nacount) == TRUE] <- 0
lowdata <- nacount[nacount$'TRUE' > 60,"Variable"] # Less then 1/3 of EARLI samples have data, 128 variables, all supp DailyFreqs
earliDHQsuppFinal <- earliDHQsuppFinal[,!colnames(earliDHQsuppFinal) %in% lowdata] # remove columns if < 1/3 samples have data
write.table(earliDHQsuppFinal, file = "Merged Database/EARLI Source Files/EARLI WGBS Sample All DHQ and Supp Variables.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)
rm(nacount, lowdata)

# Subset for Selected Variables
earliSelectedCovarNames <- read.delim(file = "Merged Database/EARLI Source Files/EARLI Selected Variables.txt", sep = "\t", 
                                      header = FALSE, stringsAsFactors = FALSE) %>% unlist %>% as.character
earliCovarsSelect <- earliDHQsuppFinal[,colnames(earliDHQsuppFinal) %in% earliSelectedCovarNames]
rm(earliDHQsuppFinal, earliSelectedCovarNames)

# Merge Covars with earli Database
earliDatabaseCovars <- merge(x = earliDatabase, y = earliCovarsSelect, by.x = "mother_id", by.y = "subjectid", all.x = TRUE, all.y = FALSE, sort = FALSE)
earliDatabaseCovars$mother_id[duplicated(earliDatabaseCovars$mother_id) == TRUE] # 9004701 9031901
earliDatabaseCovars[earliDatabaseCovars$mother_id %in% c(9004701, 9031901),c("mother_id", "COI_ID", "subject_id_sib")]
#    mother_id  COI_ID subject_id_sib
# 11   9004701 9004707        9004706
# 12   9004701 9004707        9004707
# 58   9031901 9031905        9031902
# 59   9031901 9031905        9031905

table(earliDatabaseCovars[11,] == earliDatabaseCovars[12,]) # 1 diff (subject_id_sib)
table(earliDatabaseCovars[58,] == earliDatabaseCovars[59,]) # 1 diff (subject_id_sib)
table(sapply(earliDatabaseCovars[11,], function(x) is.na(x))) # 2 NA
table(sapply(earliDatabaseCovars[12,], function(x) is.na(x))) # 2 NA
table(sapply(earliDatabaseCovars[58,], function(x) is.na(x))) # 10 NA
table(sapply(earliDatabaseCovars[59,], function(x) is.na(x))) # 10 NA
table(earliDatabaseCovars$COI_ID == earliDatabaseCovars$subject_id_sib) # 2 Diff, just remove and take unique rows

earliDatabaseCovars <- earliDatabaseCovars[,!colnames(earliDatabaseCovars) == "subject_id_sib"]
table(duplicated(earliDatabaseCovars)) # 2 duplicates
earliDatabaseCovars <- unique(earliDatabaseCovars)
rm(earliCovarsSelect, earliDatabase)

# Remove Duplicate Columns
colnames(earliDatabaseCovars)
# "Diagnosis_Alg"       "BSRC_group_Sib_36mos"  "alg_outcome"                
# "Sex"         "Subject_Gender"                    
# "Site.x"      "Site.y"        "siteid"                   
# "family_id.x" "Family_ID" "family_id.y"                  
#                                         
table(earliDatabaseCovars$Diagnosis_Alg == earliDatabaseCovars$BSRC_group_Sib_36mos) #All TRUE
table(earliDatabaseCovars$Diagnosis_Alg == earliDatabaseCovars$alg_outcome) #All FALSE, really same, just 1 0 coded
diagnosis <- data.frame("Diagnosis_Alg" = earliDatabaseCovars$Diagnosis_Alg, "alg_outcome" = earliDatabaseCovars$alg_outcome)
table(earliDatabaseCovars$Sex == earliDatabaseCovars$Subject_Gender) # All TRUE
table(earliDatabaseCovars$Site.x == earliDatabaseCovars$Site.y) # All TRUE
table(earliDatabaseCovars$Site.x == earliDatabaseCovars$siteid) # All TRUE, except when siteid is na
table(earliDatabaseCovars$family_id.x == earliDatabaseCovars$family_id.y) # All TRUE, excep when NA in both
table(earliDatabaseCovars$family_id.x == earliDatabaseCovars$Family_ID) # All TRUE, excep when NA in family_id.x

removeCols <- c("BSRC_group_Sib_36mos", "alg_outcome", "Subject_Gender", "Site.y", "siteid", "family_id.y", "family_id.x")
earliDatabaseCovars <- earliDatabaseCovars[,!colnames(earliDatabaseCovars) %in% removeCols]
earliDatabaseCovars$mother_id[earliDatabaseCovars$COI_ID == 9016702] <- 9016701
# Sample 9016702 is missing mother_id and all DHQ and Supp info, mother (9016701) is not present in DHQ and Supp Databases
write.table(earliDatabaseCovars, file = "Merged Database/EARLI WGBS Sample Covariate Database.txt", sep = "\t", row.names = FALSE,
            quote = FALSE)
rm(diagnosis, removeCols)

# Merge MARBLES and EARLI Databases ####
earliDatabase <- earliDatabaseCovars
rm(earliDatabaseCovars)

# Convert EARLI Variables to MARBLES
earliMomEdu <- earliDatabase$mom_edu %>% as.character
table(earliMomEdu)
# 10 11 12 13 14 15 16 17  3  4  5  9 
#  5  9  2  3 17 16  1  1  1  1  2  5 
earliMomEdu[earliMomEdu == 1] <- 0
earliMomEdu[earliMomEdu %in% c(2:4)] <- 1
earliMomEdu[earliMomEdu %in% c(5:8)] <- 2
earliMomEdu[earliMomEdu == 9] <- 3
earliMomEdu[earliMomEdu %in% c(10,11)] <- 4
earliMomEdu[earliMomEdu %in% c(12,13)] <- 5
earliMomEdu[earliMomEdu == 14] <- 6
earliMomEdu[earliMomEdu == 15] <- 7
earliMomEdu[earliMomEdu %in% c(16,17)] <- 8
table(earliMomEdu)
# 1  2  3  4  5  6  7  8 
# 2  2  5 14  5 17 16  2 
earliDatabase$mom_edu <- earliMomEdu

height_inches <- earliDatabase$height_inches %>% as.character %>% as.numeric
height_cm <- conv_unit(x = height_inches, from = "inch", to = "cm")
earliDatabase$height_inches <- height_cm

weight_lb <- earliDatabase$prepreg_weight %>% as.character %>% as.numeric
weight_kg <- conv_unit(x = weight_lb, from = "lbs", to = "kg")
earliDatabase$prepreg_weight <- weight_kg

earliADOSform <- read.delim(file = "Merged Database/EARLI Source Files/EARLI ADOS Forms.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
earliDatabase <- merge(x = earliDatabase, y = earliADOSform, by.x = "COI_ID", by.y = "Subject_ID", all.x = TRUE, all.y = FALSE, sort = FALSE)
table(earliDatabase$ADOS_dx_form_Sib_36mos)
# ADOS1 Mod1 NW ADOS1 Mod1 SW    ADOS1 Mod2 ADOS2 Mod1 SW    ADOS2 Mod2 
#             1             4            14            19            27 
ADOSconv <- read.delim(file = "Merged Database/EARLI Source Files/ADOScs Conversion Table.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# Used Standardizing table in Gotham et al 2009
earliADOSc <- earliDatabase[,c("COI_ID", "ADOS_dx_form_Sib_36mos", "ADOS_overall_Sib_36mos")]
earliADOSc$ADOSc <- NA
for(i in 1:nrow(earliADOSc)){
        if(earliADOSc$ADOS_dx_form_Sib_36mos[i] %in% c("ADOS1 Mod1 NW", "ADOS2 Mod1 NW")){
                earliADOSc$ADOSc[i] <- ADOSconv$ADOScs[ADOSconv$Mod1NWlow <= earliADOSc$ADOS_overall_Sib_36mos[i] &
                                                               ADOSconv$Mod1NWhigh >= earliADOSc$ADOS_overall_Sib_36mos[i]]
        }
        if(earliADOSc$ADOS_dx_form_Sib_36mos[i] %in% c("ADOS1 Mod1 SW", "ADOS2 Mod1 SW")){
                earliADOSc$ADOSc[i] <- ADOSconv$ADOScs[ADOSconv$Mod1SWlow <= earliADOSc$ADOS_overall_Sib_36mos[i] &
                                                               ADOSconv$Mod1SWhigh >= earliADOSc$ADOS_overall_Sib_36mos[i]]
        }
        if(earliADOSc$ADOS_dx_form_Sib_36mos[i] %in% c("ADOS1 Mod2", "ADOS2 Mod2")){
                earliADOSc$ADOSc[i] <- ADOSconv$ADOScs[ADOSconv$Mod2low <= earliADOSc$ADOS_overall_Sib_36mos[i] &
                                                               ADOSconv$Mod2high >= earliADOSc$ADOS_overall_Sib_36mos[i]]
        }
}

earliDatabase$ADOScs <- earliADOSc$ADOSc

# Convert MARBLES Variables to earli
marblesDiabetes <- as.numeric(marblesDatabase$DM1 | marblesDatabase$DM2)
marblesDatabase$DM1or2 <- marblesDiabetes

marblesMomAgeYr <- marblesDatabase$MomAgeYr %>% floor
marblesDatabase$MomAgeYr <- marblesMomAgeYr

# Match Columns and Merge
varMatch <- read.delim("Merged Database/EARLI MARBLES Variables Matched.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
marblesToMerge <- marblesDatabase[,varMatch$marbles]
earliToMerge <- earliDatabase[,varMatch$earli_matched]
colnames(earliToMerge) <- colnames(marblesToMerge)
mergedDatabase <- rbind(marblesToMerge, earliToMerge)
write.table(mergedDatabase, "Merged Database/MARBLES EARLI WGBS Sample Merged Covariate Database.txt", sep = "\t", quote = FALSE, row.names = FALSE)
rm(earliToMerge, marblesToMerge, samples, varMatch, earliMomEdu, height_cm, height_inches, weight_kg, weight_lb,
   marblesDiabetes, marblesMomAgeYr, i, ADOSconv, earliADOSc, earliADOSform)

# Add in variables from earli_marbles_demographics_all ####
mergedDemo <- read.delim(file = "Merged Database/earli_marbles_demographics_all_BP.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Subset for WGBS Samples
mergedDemoSub <- subset(mergedDemo, earli_ID %in% mergedDatabase$COI_ID | marbles_ID %in% mergedDatabase$COI_ID)
mergedDemoSub$COI_ID <- mergedDemoSub$marbles_ID
mergedDemoSub$COI_ID[is.na(mergedDemoSub$COI_ID)] <- mergedDemoSub$earli_ID[!is.na(mergedDemoSub$earli_ID)]

# Count NA samples by study
nacount <- aggregate(mergedDemoSub, by = list(mergedDemoSub$study), FUN = function(x) table(is.na(x))["TRUE"]) %>% t %>% as.data.frame
colnames(nacount) <- c("EARLI", "MARBLES")
nacount$EARLI <- as.character(nacount$EARLI)
nacount$MARBLES <- as.character(nacount$MARBLES)
nacount[is.na(nacount)] <- 0
nacount$Variable <- rownames(nacount)
rownames(nacount) <- 1:nrow(nacount)
nacount <- nacount[!nacount$Variable == "Group.1", c("Variable", "MARBLES", "EARLI")]
nacount$MARBLES <- as.numeric(nacount$MARBLES)
nacount$EARLI <- as.numeric(nacount$EARLI)
nacount$TotalNA <- rowSums(nacount[,c("MARBLES", "EARLI")])

# Remove Columns Missing in Whole Study
removeCols <- nacount$Variable[nacount$MARBLES == 89 | nacount$EARLI == 65]
removeCols <- c(removeCols, "X", "earli_ID", "marbles_ID")
mergedDemoSub <- mergedDemoSub[,!colnames(mergedDemoSub) %in% removeCols]
colnames(mergedDemoSub)
# [1] "mom_edu"             "mom_birthplace"      "dad_edu"             "dad_age"             "parity"             
# [6] "bw_g"                "subject_gender"      "marital_status"      "home_ownership"      "cotinine_urine_ngml"
# [11] "study"               "COI_ID"  
mergedDemoSub <- mergedDemoSub[,!colnames(mergedDemoSub) %in% c("study", "subject_gender")]

# Merge with Database
mergedDatabaseDemo <- merge(x = mergedDatabase, y = mergedDemoSub, by = "COI_ID", all.x = TRUE, all.y = FALSE)
write.table(mergedDatabaseDemo, file = "Merged Database/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Try to Merge earli_marbles_demographics_all variables missing in one study ####
# Get variables missing in one study
mergedDemoMissing <- subset(mergedDemo, earli_ID %in% mergedDatabase$COI_ID | marbles_ID %in% mergedDatabase$COI_ID)
mergedDemoMissing$COI_ID <- mergedDemoMissing$marbles_ID
mergedDemoMissing$COI_ID[is.na(mergedDemoMissing$COI_ID)] <- mergedDemoMissing$earli_ID[!is.na(mergedDemoMissing$earli_ID)]
mergedDemoMissing <- mergedDemoMissing[,colnames(mergedDemoMissing) %in% c("COI_ID", "study", removeCols)]

# Remove unneeded variables
dontNeed <- c("X", "earli_site", "sibling_gender", "mom_age", "age_at_visit_sib_36mos", "marbles_ID", "curr_pregnant", "Dx_CBE", 
              "MomAgeYr", "Mat_Weight_kg_PrePreg", "Mat_BMI_PrePreg", "HaveEQ6", "GDM_EQ",
              "PrePregBMI_EQ", "IBC", "dupereas", "Family_ID", "STS_Participant_Status.x", "GA_LD", "BMI_prepreg", "Sample_ID",
              "COI_born", "Dx_CBE_collapsed", "SmokeYN_Pregnancy", "who", "earli_ID", "birthweight", "proband_gender",
              "Subject_ID.y", "YOB", "NoFinalDiag", "Dx_alg", "Mat_Height_cm", "PE_EQ", "fw_timepoint", "days_into_trim",
              "birth_NUCHALtimes", "BMI_category", "MomEdu_7cat", "DadEdu_detail", "MomEdu_detail", "DadEdu_6cat", "MomEdu_6cat",
              "DadEdu_7cat", "Normanized.conc..ng.mg.Creatinine.", "Creatinine", "Conc__in_extract__ng_ml_", "paternal_ethnicity",
              "income", "dad_race", "mom_race", "maternal_ethnicity","Mat_Weight_kg_Deliv", "birthheadcir", "Birth_Resc", 
              "MomSmokedRegularly", "SmokeYN_Tri2", "SmokeYN_PrePreg", "OthSmokeYN_Tri2", "OthSmokeYN_PrePreg", "SmokeYN_Tri3", 
              "SmokeYN_Tri1", "OthSmokeYN_Tri3", "MomEverSmoked", "OthSmokeYN_Tri1", "OthSmokeYN_Pregnancy")
mergedDemoMissing <- mergedDemoMissing[,!colnames(mergedDemoMissing) %in% dontNeed]

# Split up variables by missing study
marblesNA <- nacount$Variable[nacount$MARBLES == 89]
earliNA <- nacount$Variable[nacount$EARLI == 65]
mergedDemoMissingMarbles <- mergedDemoMissing[,colnames(mergedDemoMissing) %in% c("COI_ID", "study", marblesNA)]
mergedDemoMissingEarli <- mergedDemoMissing[,colnames(mergedDemoMissing) %in% c("COI_ID", "study", earliNA)]

# Add Values from Marbles database to Marbles Missing Vars and Merge with mergedDatabaseDemo
colnames(mergedDemoMissingMarbles)
# "study"                 "ga_w"                  "final_creatinine_mgdl" "COI_ID"      

final_creatinine_mgdl <- mergedDemoMissingMarbles$final_creatinine_mgdl
final_creatinine_mgdl[mergedDemoMissingMarbles$study == "MAR"] <- mergedDemoMissing$final_creatinine_mg.dl[mergedDemoMissing$study == "MAR"]
mergedDemoMissingMarbles$final_creatinine_mgdl <- final_creatinine_mgdl

ga_w <- data.frame("COI_ID" = mergedDemoMissingMarbles$COI_ID, "study" = mergedDemoMissingMarbles$study,
                   "ga_w" = mergedDemoMissingMarbles$ga_w)
ga_w$marblesGestationalAge <- marblesDatabase$gest_age_deliv[match(ga_w$COI_ID, marblesDatabase$COI_ID)]
ga_w$marbles_ga_w <- round(ga_w$marblesGestationalAge / 7, 1)
ga_w$ga_w[ga_w$study == "MAR"] <- ga_w$marbles_ga_w[ga_w$study == "MAR"]
mergedDemoMissingMarbles$ga_w <- ga_w$ga_w

mergedDatabaseDemo <- merge(x = mergedDatabaseDemo, y = mergedDemoMissingMarbles, by = "COI_ID", all.x = TRUE, all.y = FALSE, sort = FALSE)
mergedDatabaseDemo <- mergedDatabaseDemo[,!colnames(mergedDatabaseDemo) == "study"]

# Reorder Columns
colnames(mergedDatabaseDemo) # 67
# [1] "COI_ID"                    "Sequencing_ID"             "Cord_Blood_IBC"            "Study"                    
# [5] "Platform"                  "Sex"                       "Site"                      "Diagnosis_Alg"            
# [9] "ADOScs"                    "MSLelcStandard36"          "MSLelTscore36"             "MSLfmTscore36"            
# [13] "MSLrlTscore36"             "MSLvrTscore36"             "percent_trimmed"           "percent_aligned"          
# [17] "percent_duplicate"         "dedup_reads"               "C_coverage"                "CG_coverage"              
# [21] "percent_cpg_meth"          "percent_chg_meth"          "percent_chh_meth"          "MomAgeYr"                 
# [25] "MomEdu_detail"             "Mat_Height_cm"             "Mat_Weight_kg_PrePreg"     "Mat_BMI_PrePreg"          
# [29] "DM1or2"                    "GDM"                       "PE"                        "SmokeYN_Pregnancy"        
# [33] "AllEQ_PV_YN_Mo_1"          "AllEQ_PV_YN_Mo_2"          "AllEQ_PV_YN_Mo_3"          "AllEQ_PV_YN_Mo1"          
# [37] "AllEQ_PV_YN_Mo2"           "AllEQ_PV_YN_Mo3"           "AllEQ_PV_YN_Mo4"           "AllEQ_PV_YN_Mo5"          
# [41] "AllEQ_PV_YN_Mo6"           "AllEQ_PV_YN_Mo7"           "AllEQ_PV_YN_Mo8"           "AllEQ_PV_YN_Mo9"          
# [45] "AllEQ_tot_All_FA_mcg_Mo_1" "AllEQ_tot_All_FA_mcg_Mo_2" "AllEQ_tot_All_FA_mcg_Mo_3" "AllEQ_tot_All_FA_mcg_Mo1" 
# [49] "AllEQ_tot_All_FA_mcg_Mo2"  "AllEQ_tot_All_FA_mcg_Mo3"  "AllEQ_tot_All_FA_mcg_Mo4"  "AllEQ_tot_All_FA_mcg_Mo5" 
# [53] "AllEQ_tot_All_FA_mcg_Mo6"  "AllEQ_tot_All_FA_mcg_Mo7"  "AllEQ_tot_All_FA_mcg_Mo8"  "AllEQ_tot_All_FA_mcg_Mo9" 
# [57] "mom_edu"                   "mom_birthplace"            "dad_edu"                   "dad_age"                  
# [61] "parity"                    "bw_g"                      "marital_status"            "home_ownership"           
# [65] "cotinine_urine_ngml"       "ga_w"                      "final_creatinine_mgdl"    

mergedDatabaseDemoCols <- c("COI_ID", "Sequencing_ID", "Cord_Blood_IBC", "Study", "Platform", "Sex", "Site", "Diagnosis_Alg", "ADOScs",
                            "MSLelcStandard36", "MSLelTscore36", "MSLfmTscore36", "MSLrlTscore36", "MSLvrTscore36", "percent_trimmed",
                            "percent_aligned", "percent_duplicate", "dedup_reads", "C_coverage", "CG_coverage","percent_cpg_meth",
                            "percent_chg_meth", "percent_chh_meth",  "ga_w", "bw_g", "MomAgeYr", "MomEdu_detail", "mom_edu", 
                            "mom_birthplace", "Mat_Height_cm", "Mat_Weight_kg_PrePreg", "Mat_BMI_PrePreg", "DM1or2", "GDM", "PE", "parity",
                            "dad_age", "dad_edu", "home_ownership", "marital_status", "SmokeYN_Pregnancy", "cotinine_urine_ngml", 
                            "final_creatinine_mgdl", "AllEQ_PV_YN_Mo_3", "AllEQ_PV_YN_Mo_2", "AllEQ_PV_YN_Mo_1", "AllEQ_PV_YN_Mo1", 
                            "AllEQ_PV_YN_Mo2", "AllEQ_PV_YN_Mo3", "AllEQ_PV_YN_Mo4", "AllEQ_PV_YN_Mo5", "AllEQ_PV_YN_Mo6", 
                            "AllEQ_PV_YN_Mo7", "AllEQ_PV_YN_Mo8", "AllEQ_PV_YN_Mo9", "AllEQ_tot_All_FA_mcg_Mo_3", 
                            "AllEQ_tot_All_FA_mcg_Mo_2", "AllEQ_tot_All_FA_mcg_Mo_1", "AllEQ_tot_All_FA_mcg_Mo1", 
                            "AllEQ_tot_All_FA_mcg_Mo2", "AllEQ_tot_All_FA_mcg_Mo3", "AllEQ_tot_All_FA_mcg_Mo4", "AllEQ_tot_All_FA_mcg_Mo5",
                            "AllEQ_tot_All_FA_mcg_Mo6", "AllEQ_tot_All_FA_mcg_Mo7", "AllEQ_tot_All_FA_mcg_Mo8", "AllEQ_tot_All_FA_mcg_Mo9")
mergedDatabaseDemo <- mergedDatabaseDemo[,mergedDatabaseDemoCols]
mergedDatabaseDemo <- mergedDatabaseDemo[,!colnames(mergedDatabaseDemo) == "MomEdu_detail"]
mergedDatabaseDemo$bw_g <- round(mergedDatabaseDemo$bw_g, 0)

mergedDatabaseDemo2 <- mergedDatabaseDemo
table(sapply(mergedDatabaseDemo2, function(x) table(is.nan(x))["TRUE"]))
#  1  2 
# 10  2 
for(i in 1:ncol(mergedDatabaseDemo2)){
        mergedDatabaseDemo2[,i][is.nan(mergedDatabaseDemo2[,i])] <- NA
}
table(sapply(mergedDatabaseDemo2, function(x) table(is.nan(x))["FALSE"])) # All 66 columns have no NaN

table(sapply(mergedDatabaseDemo2, function(x) table(x == "NaN")["TRUE"]))
# 1 4 
# 1 1 
for(i in 1:ncol(mergedDatabaseDemo2)){
        mergedDatabaseDemo2[,i][mergedDatabaseDemo2[,i] == "NaN"] <- NA
}
table(sapply(mergedDatabaseDemo2, function(x) table(x == "NaN")["FALSE"]))
write.table(mergedDatabaseDemo2, file = "Merged Database/MARBLES EARLI WGBS Sample Merged Covariate Database with Demo.txt", sep = "\t",
            quote = FALSE, row.names = FALSE)

