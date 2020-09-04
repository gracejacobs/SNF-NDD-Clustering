
## Regressing out sex, age, and IQ from cortical thickness values for structural covariance analysis

library(dplyr)
library(haven)

directory <- ("data/Data_with_headers")
CT <-read.csv(file.path(directory, paste("CT_snf.csv", sep="")))
ids <-read.csv(file.path(directory, paste("ids.csv", sep="")))
names(ids)[names(ids)=="x"] <- "ID"

ids <- cbind(ids, CT)

pond <- read_sas("data/Clinical_cog_data/pond_extract_06jun2018.sas7bdat")
pond$SUBJECT <- ifelse(grepl("^88",pond$SUBJECT), paste("sub-0", pond$SUBJECT, sep=""), paste("sub-", pond$SUBJECT, sep=""))
names(pond)[names(pond) == "SUBJECT"] <- "ID"
scq1 <- read.csv("data/Clinical_cog_data/SCQ Data 13 Aug 2018.csv")
scq1$Subject <- ifelse(grepl("^88",scq1$Subject), paste("sub-0", scq1$Subject, sep=""), paste("sub-", scq1$Subject, sep=""))
names(scq1)[names(scq1) == "Subject"] <- "ID"
scq1 <- scq1[,c("ID", "SCQTOT")]
pond <- merge(pond, scq1, by="ID")

pond$DOB_2 <- gsub('\\-', '.', pond$DOB)
pond$POND_DATE_2 <- gsub('\\-', '.', pond$POND_DATE)
pond$DOB_2 <- strptime(pond$DOB_2, format = "%Y.%m.%d")
pond$POND_DATE_2 <- strptime(pond$POND_DATE_2, format = "%Y.%m.%d")

pond$age <- difftime(pond$POND_DATE_2, pond$DOB_2, units = "days")/365
pond$DOB_2 <- NULL
pond$POND_DATE_2 <- NULL
pond$age <- as.numeric(pond$age)

pond$clin_diagnosis <- ifelse(pond$RESEARCH_CONFIRM_DIAG_STD == "1", "ASD", NA)
pond$clin_diagnosis <- ifelse(pond$RESEARCH_CONFIRM_DIAG_STD == "2", "ADHD", pond$clin_diagnosis)
pond$clin_diagnosis <- ifelse(pond$RESEARCH_CONFIRM_DIAG_STD == "3", "OCD", pond$clin_diagnosis)
pond$clin_diagnosis <- ifelse(pond$RESEARCH_CONFIRM_DIAG_STD == "6", "ADHD", pond$clin_diagnosis)
pond$clin_diagnosis <- ifelse(pond$RESEARCH_CONFIRM_DIAG_STD == "10", "CTRL", pond$clin_diagnosis)
pond$clin_diagnosis <- ifelse(pond$RESEARCH_CONFIRM_DIAG_STD == "15", "GAD", pond$clin_diagnosis)

data <- merge(pond, ids, by="ID")

data$Age <- as.numeric(data$age)
data$clin_diagnosis <- as.factor(data$clin_diagnosis)
data$GENDER <- as.factor(data$GENDER)

data$IQ <- data$WASI_FSIQ_4
data$IQ <- ifelse(is.na(data$IQ), data$WASI_II_FSIQ_4, data$IQ)
data$IQ <- ifelse(is.na(data$IQ), data$WISC_IV_FSIQ, data$IQ)
data$IQ <- ifelse(is.na(data$IQ), data$WISC_V_FSIQ, data$IQ)
data$IQ <- ifelse(is.na(data$IQ), data$WASI_II_FSIQ_2, data$IQ)

clusters <- read.csv("data/4clust_groups_k18_0.8_1000perms.csv")
clusters$groups_2 <- ifelse(clusters$groups == "1", "2", NA)
clusters$groups_2 <- ifelse(clusters$groups == "2", "4", clusters$groups_2)
clusters$groups_2 <- ifelse(clusters$groups == "3", "3", clusters$groups_2)
clusters$groups_2 <- ifelse(clusters$groups == "4", "1", clusters$groups_2)
clusters$groups <- clusters$groups_2

data$groups <- clusters$groups
data$groups <- as.factor(data$groups)

clusters <- as.data.frame(clusters$groups)
data <- data[which(!is.na(data$IQ)), ]

listt <- data[which(data$groups=="1"), ]

Cortical <- c("L_bankssts_thickavg",
                    "L_caudalanteriorcingulate_thickavg",
                    "L_caudalmiddlefrontal_thickavg",
                    "L_cuneus_thickavg",
                    "L_entorhinal_thickavg"
                    ,"L_fusiform_thickavg"
                    ,"L_inferiorparietal_thickavg"
                    ,"L_inferiortemporal_thickavg"
                    ,"L_isthmuscingulate_thickavg"
                    ,"L_lateraloccipital_thickavg"
                    ,"L_lateralorbitofrontal_thickavg"
                    ,"L_lingual_thickavg"
                    ,"L_medialorbitofrontal_thickavg"
                    ,"L_middletemporal_thickavg"
                    ,"L_parahippocampal_thickavg"
                    ,"L_paracentral_thickavg"
                    ,"L_parsopercularis_thickavg"
                    ,"L_parsorbitalis_thickavg"
                    ,"L_parstriangularis_thickavg"
                    ,"L_pericalcarine_thickavg"
                    ,"L_postcentral_thickavg"
                    ,"L_posteriorcingulate_thickavg"
                    ,"L_precentral_thickavg"
                    ,"L_precuneus_thickavg"
                    ,"L_rostralanteriorcingulate_thickavg"
                    ,"L_rostralmiddlefrontal_thickavg"
                    ,"L_superiorfrontal_thickavg"
                    ,"L_superiorparietal_thickavg"
                    ,"L_superiortemporal_thickavg"
                    ,"L_supramarginal_thickavg"
                    ,"L_frontalpole_thickavg"
                    ,"L_temporalpole_thickavg"
                    ,"L_transversetemporal_thickavg"
                    ,"L_insula_thickavg"
                    ,"R_bankssts_thickavg"
                    ,"R_caudalanteriorcingulate_thickavg"
                    ,"R_caudalmiddlefrontal_thickavg"
                    ,"R_cuneus_thickavg"
                    ,"R_entorhinal_thickavg"
                    ,"R_fusiform_thickavg"
                    ,"R_inferiorparietal_thickavg"
                    ,"R_inferiortemporal_thickavg"
                    ,"R_isthmuscingulate_thickavg"
                    ,"R_lateraloccipital_thickavg"
                    ,"R_lateralorbitofrontal_thickavg"
                    ,"R_lingual_thickavg"
                    ,"R_medialorbitofrontal_thickavg"
                    ,"R_middletemporal_thickavg"
                    ,"R_parahippocampal_thickavg"
                    ,"R_paracentral_thickavg"
                    ,"R_parsopercularis_thickavg"
                    ,"R_parsorbitalis_thickavg"
                    ,"R_parstriangularis_thickavg"
                    ,"R_pericalcarine_thickavg"
                    ,"R_postcentral_thickavg"
                    ,"R_posteriorcingulate_thickavg"
                    ,"R_precentral_thickavg"
                    ,"R_precuneus_thickavg"
                    ,"R_rostralanteriorcingulate_thickavg"
                    ,"R_rostralmiddlefrontal_thickavg"
                    ,"R_superiorfrontal_thickavg"
                    ,"R_superiorparietal_thickavg"
                    ,"R_superiortemporal_thickavg"
                    ,"R_supramarginal_thickavg"
                    ,"R_frontalpole_thickavg"
                    ,"R_temporalpole_thickavg"
                    ,"R_transversetemporal_thickavg"
                    ,"R_insula_thickavg")

mod_number <- 1
mod_resid <- matrix(nrow=length(listt[,Cortical[1]]), ncol=length(Cortical)*mod_number)


index <- 1
for (a in Cortical) {
	mod_resid[,index] <- resid(lm(data=listt, listt[,a] ~ age + GENDER +IQ, na.action=na.exclude))
	index=index+1
}

colnames(mod_resid) <- c("L_bankssts_thickavg",
                 "L_caudalanteriorcingulate_thickavg",
                 "L_caudalmiddlefrontal_thickavg",
                 "L_cuneus_thickavg",
                 "L_entorhinal_thickavg"
                 ,"L_fusiform_thickavg"
                 ,"L_inferiorparietal_thickavg"
                 ,"L_inferiortemporal_thickavg"
                 ,"L_isthmuscingulate_thickavg"
                 ,"L_lateraloccipital_thickavg"
                 ,"L_lateralorbitofrontal_thickavg"
                 ,"L_lingual_thickavg"
                 ,"L_medialorbitofrontal_thickavg"
                 ,"L_middletemporal_thickavg"
                 ,"L_parahippocampal_thickavg"
                 ,"L_paracentral_thickavg"
                 ,"L_parsopercularis_thickavg"
                 ,"L_parsorbitalis_thickavg"
                 ,"L_parstriangularis_thickavg"
                 ,"L_pericalcarine_thickavg"
                 ,"L_postcentral_thickavg"
                 ,"L_posteriorcingulate_thickavg"
                 ,"L_precentral_thickavg"
                 ,"L_precuneus_thickavg"
                 ,"L_rostralanteriorcingulate_thickavg"
                 ,"L_rostralmiddlefrontal_thickavg"
                 ,"L_superiorfrontal_thickavg"
                 ,"L_superiorparietal_thickavg"
                 ,"L_superiortemporal_thickavg"
                 ,"L_supramarginal_thickavg"
                 ,"L_frontalpole_thickavg"
                 ,"L_temporalpole_thickavg"
                 ,"L_transversetemporal_thickavg"
                 ,"L_insula_thickavg"
                 ,"R_bankssts_thickavg"
                 ,"R_caudalanteriorcingulate_thickavg"
                 ,"R_caudalmiddlefrontal_thickavg"
                 ,"R_cuneus_thickavg"
                 ,"R_entorhinal_thickavg"
                 ,"R_fusiform_thickavg"
                 ,"R_inferiorparietal_thickavg"
                 ,"R_inferiortemporal_thickavg"
                 ,"R_isthmuscingulate_thickavg"
                 ,"R_lateraloccipital_thickavg"
                 ,"R_lateralorbitofrontal_thickavg"
                 ,"R_lingual_thickavg"
                 ,"R_medialorbitofrontal_thickavg"
                 ,"R_middletemporal_thickavg"
                 ,"R_parahippocampal_thickavg"
                 ,"R_paracentral_thickavg"
                 ,"R_parsopercularis_thickavg"
                 ,"R_parsorbitalis_thickavg"
                 ,"R_parstriangularis_thickavg"
                 ,"R_pericalcarine_thickavg"
                 ,"R_postcentral_thickavg"
                 ,"R_posteriorcingulate_thickavg"
                 ,"R_precentral_thickavg"
                 ,"R_precuneus_thickavg"
                 ,"R_rostralanteriorcingulate_thickavg"
                 ,"R_rostralmiddlefrontal_thickavg"
                 ,"R_superiorfrontal_thickavg"
                 ,"R_superiorparietal_thickavg"
                 ,"R_superiortemporal_thickavg"
                 ,"R_supramarginal_thickavg"
                 ,"R_frontalpole_thickavg"
                 ,"R_temporalpole_thickavg"
                 ,"R_transversetemporal_thickavg"
                 ,"R_insula_thickavg")
mod_resid <- cbind(mod_resid, select(listt, ID, age, GENDER, IQ, clin_diagnosis))
mod_resid <- mod_resid[c(69:73, 1:68)]

write.table(mod_resid, file="output/4_Out-of-model_features/Structural_Covariance/group1_resid.csv", row.names = FALSE)

########################
#######################

listt <- data[which(data$groups=="2"), ]

Cortical <- c("L_bankssts_thickavg",
                 "L_caudalanteriorcingulate_thickavg",
                 "L_caudalmiddlefrontal_thickavg",
                 "L_cuneus_thickavg",
                 "L_entorhinal_thickavg"
                 ,"L_fusiform_thickavg"
                 ,"L_inferiorparietal_thickavg"
                 ,"L_inferiortemporal_thickavg"
                 ,"L_isthmuscingulate_thickavg"
                 ,"L_lateraloccipital_thickavg"
                 ,"L_lateralorbitofrontal_thickavg"
                 ,"L_lingual_thickavg"
                 ,"L_medialorbitofrontal_thickavg"
                 ,"L_middletemporal_thickavg"
                 ,"L_parahippocampal_thickavg"
                 ,"L_paracentral_thickavg"
                 ,"L_parsopercularis_thickavg"
                 ,"L_parsorbitalis_thickavg"
                 ,"L_parstriangularis_thickavg"
                 ,"L_pericalcarine_thickavg"
                 ,"L_postcentral_thickavg"
                 ,"L_posteriorcingulate_thickavg"
                 ,"L_precentral_thickavg"
                 ,"L_precuneus_thickavg"
                 ,"L_rostralanteriorcingulate_thickavg"
                 ,"L_rostralmiddlefrontal_thickavg"
                 ,"L_superiorfrontal_thickavg"
                 ,"L_superiorparietal_thickavg"
                 ,"L_superiortemporal_thickavg"
                 ,"L_supramarginal_thickavg"
                 ,"L_frontalpole_thickavg"
                 ,"L_temporalpole_thickavg"
                 ,"L_transversetemporal_thickavg"
                 ,"L_insula_thickavg"
                 ,"R_bankssts_thickavg"
                 ,"R_caudalanteriorcingulate_thickavg"
                 ,"R_caudalmiddlefrontal_thickavg"
                 ,"R_cuneus_thickavg"
                 ,"R_entorhinal_thickavg"
                 ,"R_fusiform_thickavg"
                 ,"R_inferiorparietal_thickavg"
                 ,"R_inferiortemporal_thickavg"
                 ,"R_isthmuscingulate_thickavg"
                 ,"R_lateraloccipital_thickavg"
                 ,"R_lateralorbitofrontal_thickavg"
                 ,"R_lingual_thickavg"
                 ,"R_medialorbitofrontal_thickavg"
                 ,"R_middletemporal_thickavg"
                 ,"R_parahippocampal_thickavg"
                 ,"R_paracentral_thickavg"
                 ,"R_parsopercularis_thickavg"
                 ,"R_parsorbitalis_thickavg"
                 ,"R_parstriangularis_thickavg"
                 ,"R_pericalcarine_thickavg"
                 ,"R_postcentral_thickavg"
                 ,"R_posteriorcingulate_thickavg"
                 ,"R_precentral_thickavg"
                 ,"R_precuneus_thickavg"
                 ,"R_rostralanteriorcingulate_thickavg"
                 ,"R_rostralmiddlefrontal_thickavg"
                 ,"R_superiorfrontal_thickavg"
                 ,"R_superiorparietal_thickavg"
                 ,"R_superiortemporal_thickavg"
                 ,"R_supramarginal_thickavg"
                 ,"R_frontalpole_thickavg"
                 ,"R_temporalpole_thickavg"
                 ,"R_transversetemporal_thickavg"
                 ,"R_insula_thickavg")

mod_number <- 1
mod_resid <- matrix(nrow=length(listt[,Cortical[1]]), ncol=length(Cortical)*mod_number)


index <- 1
for (a in Cortical) {
  mod_resid[,index] <- resid(lm(data=listt, listt[,a] ~ age + GENDER +IQ, na.action=na.exclude))
  index=index+1
}

colnames(mod_resid) <- c("L_bankssts_thickavg",
                         "L_caudalanteriorcingulate_thickavg",
                         "L_caudalmiddlefrontal_thickavg",
                         "L_cuneus_thickavg",
                         "L_entorhinal_thickavg"
                         ,"L_fusiform_thickavg"
                         ,"L_inferiorparietal_thickavg"
                         ,"L_inferiortemporal_thickavg"
                         ,"L_isthmuscingulate_thickavg"
                         ,"L_lateraloccipital_thickavg"
                         ,"L_lateralorbitofrontal_thickavg"
                         ,"L_lingual_thickavg"
                         ,"L_medialorbitofrontal_thickavg"
                         ,"L_middletemporal_thickavg"
                         ,"L_parahippocampal_thickavg"
                         ,"L_paracentral_thickavg"
                         ,"L_parsopercularis_thickavg"
                         ,"L_parsorbitalis_thickavg"
                         ,"L_parstriangularis_thickavg"
                         ,"L_pericalcarine_thickavg"
                         ,"L_postcentral_thickavg"
                         ,"L_posteriorcingulate_thickavg"
                         ,"L_precentral_thickavg"
                         ,"L_precuneus_thickavg"
                         ,"L_rostralanteriorcingulate_thickavg"
                         ,"L_rostralmiddlefrontal_thickavg"
                         ,"L_superiorfrontal_thickavg"
                         ,"L_superiorparietal_thickavg"
                         ,"L_superiortemporal_thickavg"
                         ,"L_supramarginal_thickavg"
                         ,"L_frontalpole_thickavg"
                         ,"L_temporalpole_thickavg"
                         ,"L_transversetemporal_thickavg"
                         ,"L_insula_thickavg"
                         ,"R_bankssts_thickavg"
                         ,"R_caudalanteriorcingulate_thickavg"
                         ,"R_caudalmiddlefrontal_thickavg"
                         ,"R_cuneus_thickavg"
                         ,"R_entorhinal_thickavg"
                         ,"R_fusiform_thickavg"
                         ,"R_inferiorparietal_thickavg"
                         ,"R_inferiortemporal_thickavg"
                         ,"R_isthmuscingulate_thickavg"
                         ,"R_lateraloccipital_thickavg"
                         ,"R_lateralorbitofrontal_thickavg"
                         ,"R_lingual_thickavg"
                         ,"R_medialorbitofrontal_thickavg"
                         ,"R_middletemporal_thickavg"
                         ,"R_parahippocampal_thickavg"
                         ,"R_paracentral_thickavg"
                         ,"R_parsopercularis_thickavg"
                         ,"R_parsorbitalis_thickavg"
                         ,"R_parstriangularis_thickavg"
                         ,"R_pericalcarine_thickavg"
                         ,"R_postcentral_thickavg"
                         ,"R_posteriorcingulate_thickavg"
                         ,"R_precentral_thickavg"
                         ,"R_precuneus_thickavg"
                         ,"R_rostralanteriorcingulate_thickavg"
                         ,"R_rostralmiddlefrontal_thickavg"
                         ,"R_superiorfrontal_thickavg"
                         ,"R_superiorparietal_thickavg"
                         ,"R_superiortemporal_thickavg"
                         ,"R_supramarginal_thickavg"
                         ,"R_frontalpole_thickavg"
                         ,"R_temporalpole_thickavg"
                         ,"R_transversetemporal_thickavg"
                         ,"R_insula_thickavg")
mod_resid <- cbind(mod_resid, select(listt, ID, age, GENDER, IQ, clin_diagnosis))
mod_resid <- mod_resid[c(69:73, 1:68)]

write.table(mod_resid, file="output/4_Out-of-model_features/Structural_Covariance/group2_resid.csv", row.names = FALSE)

##########################
#########################

########################
#######################

listt <- data[which(data$groups=="3"), ]

Cortical <- c("L_bankssts_thickavg",
                 "L_caudalanteriorcingulate_thickavg",
                 "L_caudalmiddlefrontal_thickavg",
                 "L_cuneus_thickavg",
                 "L_entorhinal_thickavg"
                 ,"L_fusiform_thickavg"
                 ,"L_inferiorparietal_thickavg"
                 ,"L_inferiortemporal_thickavg"
                 ,"L_isthmuscingulate_thickavg"
                 ,"L_lateraloccipital_thickavg"
                 ,"L_lateralorbitofrontal_thickavg"
                 ,"L_lingual_thickavg"
                 ,"L_medialorbitofrontal_thickavg"
                 ,"L_middletemporal_thickavg"
                 ,"L_parahippocampal_thickavg"
                 ,"L_paracentral_thickavg"
                 ,"L_parsopercularis_thickavg"
                 ,"L_parsorbitalis_thickavg"
                 ,"L_parstriangularis_thickavg"
                 ,"L_pericalcarine_thickavg"
                 ,"L_postcentral_thickavg"
                 ,"L_posteriorcingulate_thickavg"
                 ,"L_precentral_thickavg"
                 ,"L_precuneus_thickavg"
                 ,"L_rostralanteriorcingulate_thickavg"
                 ,"L_rostralmiddlefrontal_thickavg"
                 ,"L_superiorfrontal_thickavg"
                 ,"L_superiorparietal_thickavg"
                 ,"L_superiortemporal_thickavg"
                 ,"L_supramarginal_thickavg"
                 ,"L_frontalpole_thickavg"
                 ,"L_temporalpole_thickavg"
                 ,"L_transversetemporal_thickavg"
                 ,"L_insula_thickavg"
                 ,"R_bankssts_thickavg"
                 ,"R_caudalanteriorcingulate_thickavg"
                 ,"R_caudalmiddlefrontal_thickavg"
                 ,"R_cuneus_thickavg"
                 ,"R_entorhinal_thickavg"
                 ,"R_fusiform_thickavg"
                 ,"R_inferiorparietal_thickavg"
                 ,"R_inferiortemporal_thickavg"
                 ,"R_isthmuscingulate_thickavg"
                 ,"R_lateraloccipital_thickavg"
                 ,"R_lateralorbitofrontal_thickavg"
                 ,"R_lingual_thickavg"
                 ,"R_medialorbitofrontal_thickavg"
                 ,"R_middletemporal_thickavg"
                 ,"R_parahippocampal_thickavg"
                 ,"R_paracentral_thickavg"
                 ,"R_parsopercularis_thickavg"
                 ,"R_parsorbitalis_thickavg"
                 ,"R_parstriangularis_thickavg"
                 ,"R_pericalcarine_thickavg"
                 ,"R_postcentral_thickavg"
                 ,"R_posteriorcingulate_thickavg"
                 ,"R_precentral_thickavg"
                 ,"R_precuneus_thickavg"
                 ,"R_rostralanteriorcingulate_thickavg"
                 ,"R_rostralmiddlefrontal_thickavg"
                 ,"R_superiorfrontal_thickavg"
                 ,"R_superiorparietal_thickavg"
                 ,"R_superiortemporal_thickavg"
                 ,"R_supramarginal_thickavg"
                 ,"R_frontalpole_thickavg"
                 ,"R_temporalpole_thickavg"
                 ,"R_transversetemporal_thickavg"
                 ,"R_insula_thickavg")

mod_number <- 1
mod_resid <- matrix(nrow=length(listt[,Cortical[1]]), ncol=length(Cortical)*mod_number)


index <- 1
for (a in Cortical) {
  mod_resid[,index] <- resid(lm(data=listt, listt[,a] ~ age + GENDER + IQ, na.action=na.exclude))
  index=index+1
}

colnames(mod_resid) <- c("L_bankssts_thickavg",
                         "L_caudalanteriorcingulate_thickavg",
                         "L_caudalmiddlefrontal_thickavg",
                         "L_cuneus_thickavg",
                         "L_entorhinal_thickavg"
                         ,"L_fusiform_thickavg"
                         ,"L_inferiorparietal_thickavg"
                         ,"L_inferiortemporal_thickavg"
                         ,"L_isthmuscingulate_thickavg"
                         ,"L_lateraloccipital_thickavg"
                         ,"L_lateralorbitofrontal_thickavg"
                         ,"L_lingual_thickavg"
                         ,"L_medialorbitofrontal_thickavg"
                         ,"L_middletemporal_thickavg"
                         ,"L_parahippocampal_thickavg"
                         ,"L_paracentral_thickavg"
                         ,"L_parsopercularis_thickavg"
                         ,"L_parsorbitalis_thickavg"
                         ,"L_parstriangularis_thickavg"
                         ,"L_pericalcarine_thickavg"
                         ,"L_postcentral_thickavg"
                         ,"L_posteriorcingulate_thickavg"
                         ,"L_precentral_thickavg"
                         ,"L_precuneus_thickavg"
                         ,"L_rostralanteriorcingulate_thickavg"
                         ,"L_rostralmiddlefrontal_thickavg"
                         ,"L_superiorfrontal_thickavg"
                         ,"L_superiorparietal_thickavg"
                         ,"L_superiortemporal_thickavg"
                         ,"L_supramarginal_thickavg"
                         ,"L_frontalpole_thickavg"
                         ,"L_temporalpole_thickavg"
                         ,"L_transversetemporal_thickavg"
                         ,"L_insula_thickavg"
                         ,"R_bankssts_thickavg"
                         ,"R_caudalanteriorcingulate_thickavg"
                         ,"R_caudalmiddlefrontal_thickavg"
                         ,"R_cuneus_thickavg"
                         ,"R_entorhinal_thickavg"
                         ,"R_fusiform_thickavg"
                         ,"R_inferiorparietal_thickavg"
                         ,"R_inferiortemporal_thickavg"
                         ,"R_isthmuscingulate_thickavg"
                         ,"R_lateraloccipital_thickavg"
                         ,"R_lateralorbitofrontal_thickavg"
                         ,"R_lingual_thickavg"
                         ,"R_medialorbitofrontal_thickavg"
                         ,"R_middletemporal_thickavg"
                         ,"R_parahippocampal_thickavg"
                         ,"R_paracentral_thickavg"
                         ,"R_parsopercularis_thickavg"
                         ,"R_parsorbitalis_thickavg"
                         ,"R_parstriangularis_thickavg"
                         ,"R_pericalcarine_thickavg"
                         ,"R_postcentral_thickavg"
                         ,"R_posteriorcingulate_thickavg"
                         ,"R_precentral_thickavg"
                         ,"R_precuneus_thickavg"
                         ,"R_rostralanteriorcingulate_thickavg"
                         ,"R_rostralmiddlefrontal_thickavg"
                         ,"R_superiorfrontal_thickavg"
                         ,"R_superiorparietal_thickavg"
                         ,"R_superiortemporal_thickavg"
                         ,"R_supramarginal_thickavg"
                         ,"R_frontalpole_thickavg"
                         ,"R_temporalpole_thickavg"
                         ,"R_transversetemporal_thickavg"
                         ,"R_insula_thickavg")
mod_resid <- cbind(mod_resid, select(listt, ID, age, GENDER, IQ, clin_diagnosis))
mod_resid <- mod_resid[c(69:73, 1:68)]

write.table(mod_resid, file="output/4_Out-of-model_features/Structural_Covariance/group3_resid.csv", row.names = FALSE)

########################
#######################

listt <- data[which(data$groups=="4"), ]

Cortical <- c("L_bankssts_thickavg",
                 "L_caudalanteriorcingulate_thickavg",
                 "L_caudalmiddlefrontal_thickavg",
                 "L_cuneus_thickavg",
                 "L_entorhinal_thickavg"
                 ,"L_fusiform_thickavg"
                 ,"L_inferiorparietal_thickavg"
                 ,"L_inferiortemporal_thickavg"
                 ,"L_isthmuscingulate_thickavg"
                 ,"L_lateraloccipital_thickavg"
                 ,"L_lateralorbitofrontal_thickavg"
                 ,"L_lingual_thickavg"
                 ,"L_medialorbitofrontal_thickavg"
                 ,"L_middletemporal_thickavg"
                 ,"L_parahippocampal_thickavg"
                 ,"L_paracentral_thickavg"
                 ,"L_parsopercularis_thickavg"
                 ,"L_parsorbitalis_thickavg"
                 ,"L_parstriangularis_thickavg"
                 ,"L_pericalcarine_thickavg"
                 ,"L_postcentral_thickavg"
                 ,"L_posteriorcingulate_thickavg"
                 ,"L_precentral_thickavg"
                 ,"L_precuneus_thickavg"
                 ,"L_rostralanteriorcingulate_thickavg"
                 ,"L_rostralmiddlefrontal_thickavg"
                 ,"L_superiorfrontal_thickavg"
                 ,"L_superiorparietal_thickavg"
                 ,"L_superiortemporal_thickavg"
                 ,"L_supramarginal_thickavg"
                 ,"L_frontalpole_thickavg"
                 ,"L_temporalpole_thickavg"
                 ,"L_transversetemporal_thickavg"
                 ,"L_insula_thickavg"
                 ,"R_bankssts_thickavg"
                 ,"R_caudalanteriorcingulate_thickavg"
                 ,"R_caudalmiddlefrontal_thickavg"
                 ,"R_cuneus_thickavg"
                 ,"R_entorhinal_thickavg"
                 ,"R_fusiform_thickavg"
                 ,"R_inferiorparietal_thickavg"
                 ,"R_inferiortemporal_thickavg"
                 ,"R_isthmuscingulate_thickavg"
                 ,"R_lateraloccipital_thickavg"
                 ,"R_lateralorbitofrontal_thickavg"
                 ,"R_lingual_thickavg"
                 ,"R_medialorbitofrontal_thickavg"
                 ,"R_middletemporal_thickavg"
                 ,"R_parahippocampal_thickavg"
                 ,"R_paracentral_thickavg"
                 ,"R_parsopercularis_thickavg"
                 ,"R_parsorbitalis_thickavg"
                 ,"R_parstriangularis_thickavg"
                 ,"R_pericalcarine_thickavg"
                 ,"R_postcentral_thickavg"
                 ,"R_posteriorcingulate_thickavg"
                 ,"R_precentral_thickavg"
                 ,"R_precuneus_thickavg"
                 ,"R_rostralanteriorcingulate_thickavg"
                 ,"R_rostralmiddlefrontal_thickavg"
                 ,"R_superiorfrontal_thickavg"
                 ,"R_superiorparietal_thickavg"
                 ,"R_superiortemporal_thickavg"
                 ,"R_supramarginal_thickavg"
                 ,"R_frontalpole_thickavg"
                 ,"R_temporalpole_thickavg"
                 ,"R_transversetemporal_thickavg"
                 ,"R_insula_thickavg")

mod_number <- 1
mod_resid <- matrix(nrow=length(listt[,Cortical[1]]), ncol=length(Cortical)*mod_number)


index <- 1
for (a in Cortical) {
  mod_resid[,index] <- resid(lm(data=listt, listt[,a] ~ age + GENDER + IQ, na.action=na.exclude))
  index=index+1
}

colnames(mod_resid) <- c("L_bankssts_thickavg",
                         "L_caudalanteriorcingulate_thickavg",
                         "L_caudalmiddlefrontal_thickavg",
                         "L_cuneus_thickavg",
                         "L_entorhinal_thickavg"
                         ,"L_fusiform_thickavg"
                         ,"L_inferiorparietal_thickavg"
                         ,"L_inferiortemporal_thickavg"
                         ,"L_isthmuscingulate_thickavg"
                         ,"L_lateraloccipital_thickavg"
                         ,"L_lateralorbitofrontal_thickavg"
                         ,"L_lingual_thickavg"
                         ,"L_medialorbitofrontal_thickavg"
                         ,"L_middletemporal_thickavg"
                         ,"L_parahippocampal_thickavg"
                         ,"L_paracentral_thickavg"
                         ,"L_parsopercularis_thickavg"
                         ,"L_parsorbitalis_thickavg"
                         ,"L_parstriangularis_thickavg"
                         ,"L_pericalcarine_thickavg"
                         ,"L_postcentral_thickavg"
                         ,"L_posteriorcingulate_thickavg"
                         ,"L_precentral_thickavg"
                         ,"L_precuneus_thickavg"
                         ,"L_rostralanteriorcingulate_thickavg"
                         ,"L_rostralmiddlefrontal_thickavg"
                         ,"L_superiorfrontal_thickavg"
                         ,"L_superiorparietal_thickavg"
                         ,"L_superiortemporal_thickavg"
                         ,"L_supramarginal_thickavg"
                         ,"L_frontalpole_thickavg"
                         ,"L_temporalpole_thickavg"
                         ,"L_transversetemporal_thickavg"
                         ,"L_insula_thickavg"
                         ,"R_bankssts_thickavg"
                         ,"R_caudalanteriorcingulate_thickavg"
                         ,"R_caudalmiddlefrontal_thickavg"
                         ,"R_cuneus_thickavg"
                         ,"R_entorhinal_thickavg"
                         ,"R_fusiform_thickavg"
                         ,"R_inferiorparietal_thickavg"
                         ,"R_inferiortemporal_thickavg"
                         ,"R_isthmuscingulate_thickavg"
                         ,"R_lateraloccipital_thickavg"
                         ,"R_lateralorbitofrontal_thickavg"
                         ,"R_lingual_thickavg"
                         ,"R_medialorbitofrontal_thickavg"
                         ,"R_middletemporal_thickavg"
                         ,"R_parahippocampal_thickavg"
                         ,"R_paracentral_thickavg"
                         ,"R_parsopercularis_thickavg"
                         ,"R_parsorbitalis_thickavg"
                         ,"R_parstriangularis_thickavg"
                         ,"R_pericalcarine_thickavg"
                         ,"R_postcentral_thickavg"
                         ,"R_posteriorcingulate_thickavg"
                         ,"R_precentral_thickavg"
                         ,"R_precuneus_thickavg"
                         ,"R_rostralanteriorcingulate_thickavg"
                         ,"R_rostralmiddlefrontal_thickavg"
                         ,"R_superiorfrontal_thickavg"
                         ,"R_superiorparietal_thickavg"
                         ,"R_superiortemporal_thickavg"
                         ,"R_supramarginal_thickavg"
                         ,"R_frontalpole_thickavg"
                         ,"R_temporalpole_thickavg"
                         ,"R_transversetemporal_thickavg"
                         ,"R_insula_thickavg")
mod_resid <- cbind(mod_resid, select(listt, ID, age, GENDER, IQ, groups))
mod_resid <- mod_resid[c(69:73, 1:68)]

write.table(mod_resid, file="output/4_Out-of-model_features/Structural_Covariance/group4_resid.csv", row.names = FALSE)

####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################

listt <- data[which(data$clin_diagnosis=="ASD"), ]

mod_number <- 1
mod_resid <- matrix(nrow=length(listt[,Cortical[1]]), ncol=length(Cortical)*mod_number)


index <- 1
for (a in Cortical) {
  mod_resid[,index] <- resid(lm(data=listt, listt[,a] ~ age + GENDER + IQ, na.action=na.exclude))
  index=index+1
}

colnames(mod_resid) <- c("L_bankssts_thickavg",
                         "L_caudalanteriorcingulate_thickavg",
                         "L_caudalmiddlefrontal_thickavg",
                         "L_cuneus_thickavg",
                         "L_entorhinal_thickavg"
                         ,"L_fusiform_thickavg"
                         ,"L_inferiorparietal_thickavg"
                         ,"L_inferiortemporal_thickavg"
                         ,"L_isthmuscingulate_thickavg"
                         ,"L_lateraloccipital_thickavg"
                         ,"L_lateralorbitofrontal_thickavg"
                         ,"L_lingual_thickavg"
                         ,"L_medialorbitofrontal_thickavg"
                         ,"L_middletemporal_thickavg"
                         ,"L_parahippocampal_thickavg"
                         ,"L_paracentral_thickavg"
                         ,"L_parsopercularis_thickavg"
                         ,"L_parsorbitalis_thickavg"
                         ,"L_parstriangularis_thickavg"
                         ,"L_pericalcarine_thickavg"
                         ,"L_postcentral_thickavg"
                         ,"L_posteriorcingulate_thickavg"
                         ,"L_precentral_thickavg"
                         ,"L_precuneus_thickavg"
                         ,"L_rostralanteriorcingulate_thickavg"
                         ,"L_rostralmiddlefrontal_thickavg"
                         ,"L_superiorfrontal_thickavg"
                         ,"L_superiorparietal_thickavg"
                         ,"L_superiortemporal_thickavg"
                         ,"L_supramarginal_thickavg"
                         ,"L_frontalpole_thickavg"
                         ,"L_temporalpole_thickavg"
                         ,"L_transversetemporal_thickavg"
                         ,"L_insula_thickavg"
                         ,"R_bankssts_thickavg"
                         ,"R_caudalanteriorcingulate_thickavg"
                         ,"R_caudalmiddlefrontal_thickavg"
                         ,"R_cuneus_thickavg"
                         ,"R_entorhinal_thickavg"
                         ,"R_fusiform_thickavg"
                         ,"R_inferiorparietal_thickavg"
                         ,"R_inferiortemporal_thickavg"
                         ,"R_isthmuscingulate_thickavg"
                         ,"R_lateraloccipital_thickavg"
                         ,"R_lateralorbitofrontal_thickavg"
                         ,"R_lingual_thickavg"
                         ,"R_medialorbitofrontal_thickavg"
                         ,"R_middletemporal_thickavg"
                         ,"R_parahippocampal_thickavg"
                         ,"R_paracentral_thickavg"
                         ,"R_parsopercularis_thickavg"
                         ,"R_parsorbitalis_thickavg"
                         ,"R_parstriangularis_thickavg"
                         ,"R_pericalcarine_thickavg"
                         ,"R_postcentral_thickavg"
                         ,"R_posteriorcingulate_thickavg"
                         ,"R_precentral_thickavg"
                         ,"R_precuneus_thickavg"
                         ,"R_rostralanteriorcingulate_thickavg"
                         ,"R_rostralmiddlefrontal_thickavg"
                         ,"R_superiorfrontal_thickavg"
                         ,"R_superiorparietal_thickavg"
                         ,"R_superiortemporal_thickavg"
                         ,"R_supramarginal_thickavg"
                         ,"R_frontalpole_thickavg"
                         ,"R_temporalpole_thickavg"
                         ,"R_transversetemporal_thickavg"
                         ,"R_insula_thickavg")
mod_resid <- cbind(mod_resid, select(listt, ID, age, GENDER, IQ, groups))
mod_resid <- mod_resid[c(69:73, 1:68)]

write.table(mod_resid, file="output/4_Out-of-model_features/Structural_Covariance/ASD_resid.csv", row.names = FALSE)

listt <- data[which(data$clin_diagnosis=="ADHD"), ]

mod_number <- 1
mod_resid <- matrix(nrow=length(listt[,Cortical[1]]), ncol=length(Cortical)*mod_number)


index <- 1
for (a in Cortical) {
  mod_resid[,index] <- resid(lm(data=listt, listt[,a] ~ age + GENDER + IQ, na.action=na.exclude))
  index=index+1
}

colnames(mod_resid) <- c("L_bankssts_thickavg",
                         "L_caudalanteriorcingulate_thickavg",
                         "L_caudalmiddlefrontal_thickavg",
                         "L_cuneus_thickavg",
                         "L_entorhinal_thickavg"
                         ,"L_fusiform_thickavg"
                         ,"L_inferiorparietal_thickavg"
                         ,"L_inferiortemporal_thickavg"
                         ,"L_isthmuscingulate_thickavg"
                         ,"L_lateraloccipital_thickavg"
                         ,"L_lateralorbitofrontal_thickavg"
                         ,"L_lingual_thickavg"
                         ,"L_medialorbitofrontal_thickavg"
                         ,"L_middletemporal_thickavg"
                         ,"L_parahippocampal_thickavg"
                         ,"L_paracentral_thickavg"
                         ,"L_parsopercularis_thickavg"
                         ,"L_parsorbitalis_thickavg"
                         ,"L_parstriangularis_thickavg"
                         ,"L_pericalcarine_thickavg"
                         ,"L_postcentral_thickavg"
                         ,"L_posteriorcingulate_thickavg"
                         ,"L_precentral_thickavg"
                         ,"L_precuneus_thickavg"
                         ,"L_rostralanteriorcingulate_thickavg"
                         ,"L_rostralmiddlefrontal_thickavg"
                         ,"L_superiorfrontal_thickavg"
                         ,"L_superiorparietal_thickavg"
                         ,"L_superiortemporal_thickavg"
                         ,"L_supramarginal_thickavg"
                         ,"L_frontalpole_thickavg"
                         ,"L_temporalpole_thickavg"
                         ,"L_transversetemporal_thickavg"
                         ,"L_insula_thickavg"
                         ,"R_bankssts_thickavg"
                         ,"R_caudalanteriorcingulate_thickavg"
                         ,"R_caudalmiddlefrontal_thickavg"
                         ,"R_cuneus_thickavg"
                         ,"R_entorhinal_thickavg"
                         ,"R_fusiform_thickavg"
                         ,"R_inferiorparietal_thickavg"
                         ,"R_inferiortemporal_thickavg"
                         ,"R_isthmuscingulate_thickavg"
                         ,"R_lateraloccipital_thickavg"
                         ,"R_lateralorbitofrontal_thickavg"
                         ,"R_lingual_thickavg"
                         ,"R_medialorbitofrontal_thickavg"
                         ,"R_middletemporal_thickavg"
                         ,"R_parahippocampal_thickavg"
                         ,"R_paracentral_thickavg"
                         ,"R_parsopercularis_thickavg"
                         ,"R_parsorbitalis_thickavg"
                         ,"R_parstriangularis_thickavg"
                         ,"R_pericalcarine_thickavg"
                         ,"R_postcentral_thickavg"
                         ,"R_posteriorcingulate_thickavg"
                         ,"R_precentral_thickavg"
                         ,"R_precuneus_thickavg"
                         ,"R_rostralanteriorcingulate_thickavg"
                         ,"R_rostralmiddlefrontal_thickavg"
                         ,"R_superiorfrontal_thickavg"
                         ,"R_superiorparietal_thickavg"
                         ,"R_superiortemporal_thickavg"
                         ,"R_supramarginal_thickavg"
                         ,"R_frontalpole_thickavg"
                         ,"R_temporalpole_thickavg"
                         ,"R_transversetemporal_thickavg"
                         ,"R_insula_thickavg")
mod_resid <- cbind(mod_resid, select(listt, ID, age, GENDER, IQ, groups))
mod_resid <- mod_resid[c(69:73, 1:68)]

write.table(mod_resid, file="output/4_Out-of-model_features/Structural_Covariance/ADHD_resid.csv", row.names = FALSE)

listt <- data[which(data$clin_diagnosis=="OCD"), ]

mod_number <- 1
mod_resid <- matrix(nrow=length(listt[,Cortical[1]]), ncol=length(Cortical)*mod_number)


index <- 1
for (a in Cortical) {
  mod_resid[,index] <- resid(lm(data=listt, listt[,a] ~ age + GENDER + IQ, na.action=na.exclude))
  index=index+1
}

colnames(mod_resid) <- c("L_bankssts_thickavg",
                         "L_caudalanteriorcingulate_thickavg",
                         "L_caudalmiddlefrontal_thickavg",
                         "L_cuneus_thickavg",
                         "L_entorhinal_thickavg"
                         ,"L_fusiform_thickavg"
                         ,"L_inferiorparietal_thickavg"
                         ,"L_inferiortemporal_thickavg"
                         ,"L_isthmuscingulate_thickavg"
                         ,"L_lateraloccipital_thickavg"
                         ,"L_lateralorbitofrontal_thickavg"
                         ,"L_lingual_thickavg"
                         ,"L_medialorbitofrontal_thickavg"
                         ,"L_middletemporal_thickavg"
                         ,"L_parahippocampal_thickavg"
                         ,"L_paracentral_thickavg"
                         ,"L_parsopercularis_thickavg"
                         ,"L_parsorbitalis_thickavg"
                         ,"L_parstriangularis_thickavg"
                         ,"L_pericalcarine_thickavg"
                         ,"L_postcentral_thickavg"
                         ,"L_posteriorcingulate_thickavg"
                         ,"L_precentral_thickavg"
                         ,"L_precuneus_thickavg"
                         ,"L_rostralanteriorcingulate_thickavg"
                         ,"L_rostralmiddlefrontal_thickavg"
                         ,"L_superiorfrontal_thickavg"
                         ,"L_superiorparietal_thickavg"
                         ,"L_superiortemporal_thickavg"
                         ,"L_supramarginal_thickavg"
                         ,"L_frontalpole_thickavg"
                         ,"L_temporalpole_thickavg"
                         ,"L_transversetemporal_thickavg"
                         ,"L_insula_thickavg"
                         ,"R_bankssts_thickavg"
                         ,"R_caudalanteriorcingulate_thickavg"
                         ,"R_caudalmiddlefrontal_thickavg"
                         ,"R_cuneus_thickavg"
                         ,"R_entorhinal_thickavg"
                         ,"R_fusiform_thickavg"
                         ,"R_inferiorparietal_thickavg"
                         ,"R_inferiortemporal_thickavg"
                         ,"R_isthmuscingulate_thickavg"
                         ,"R_lateraloccipital_thickavg"
                         ,"R_lateralorbitofrontal_thickavg"
                         ,"R_lingual_thickavg"
                         ,"R_medialorbitofrontal_thickavg"
                         ,"R_middletemporal_thickavg"
                         ,"R_parahippocampal_thickavg"
                         ,"R_paracentral_thickavg"
                         ,"R_parsopercularis_thickavg"
                         ,"R_parsorbitalis_thickavg"
                         ,"R_parstriangularis_thickavg"
                         ,"R_pericalcarine_thickavg"
                         ,"R_postcentral_thickavg"
                         ,"R_posteriorcingulate_thickavg"
                         ,"R_precentral_thickavg"
                         ,"R_precuneus_thickavg"
                         ,"R_rostralanteriorcingulate_thickavg"
                         ,"R_rostralmiddlefrontal_thickavg"
                         ,"R_superiorfrontal_thickavg"
                         ,"R_superiorparietal_thickavg"
                         ,"R_superiortemporal_thickavg"
                         ,"R_supramarginal_thickavg"
                         ,"R_frontalpole_thickavg"
                         ,"R_temporalpole_thickavg"
                         ,"R_transversetemporal_thickavg"
                         ,"R_insula_thickavg")
mod_resid <- cbind(mod_resid, select(listt, ID, age, GENDER, IQ, groups))
mod_resid <- mod_resid[c(69:73, 1:68)]

write.table(mod_resid, file="output/4_Out-of-model_features/Structural_Covariance/OCD_resid.csv", row.names = FALSE)