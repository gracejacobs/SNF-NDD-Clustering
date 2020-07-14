### Setting up data for the analysis of SNF clusters

clusters <- read.csv("output/1_SNF_analysis/4clust_groups_k18_0.8_1000perms.csv")

clusters$groups_2 <- ifelse(clusters$groups == "1", "2", NA)
clusters$groups_2 <- ifelse(clusters$groups == "2", "4", clusters$groups_2)
clusters$groups_2 <- ifelse(clusters$groups == "3", "3", clusters$groups_2)
clusters$groups_2 <- ifelse(clusters$groups == "4", "1", clusters$groups_2)

clusters$groups <- clusters$groups_2

clusters$X <- NULL
clusters$id <- NULL
names(clusters)[names(clusters)=="V1"] <- "ID"

# data with headers
directory <- ("data/Data_with_headers")
CT <-read.csv(file.path(directory, paste("CT_snf.csv", sep="")))
SA <-read.csv(file.path(directory, paste("SA_snf.csv", sep="")))
Volume <-read.csv(file.path(directory, paste("Volumes_snf.csv", sep="")))
FA <-read.csv(file.path(directory, paste("FA_snf.csv", sep="")))
ids <-read.csv(file.path(directory, paste("ids.csv", sep="")))

names(ids)[names(ids)=="x"] <- "ID"
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

all <- merge(pond, ids, by="ID")

all$clin_diagnosis <- factor(all$clin_diagnosis, labels=c("ADHD", "ASD", "OCD"))
table(all$clin_diagnosis)

all <- merge(clusters, all, by="ID")

all <- cbind(all, CT)
all <- cbind(all, FA)
all <- cbind(all, Volume)
all <- cbind(all, SA)

all$IQ <- all$WASI_FSIQ_4
all$IQ <- ifelse(is.na(all$IQ), all$WASI_II_FSIQ_4, all$IQ)
all$IQ <- ifelse(is.na(all$IQ), all$WISC_IV_FSIQ, all$IQ)
all$IQ <- ifelse(is.na(all$IQ), all$WISC_V_FSIQ, all$IQ)
all$IQ <- ifelse(is.na(all$IQ), all$WASI_II_FSIQ_2, all$IQ)


all$GENDER <- as.factor(all$GENDER)
all$group <- as.factor(all$groups)
all$age <- as.numeric(all$age)
all$IQ <- as.numeric(all$IQ)

clust1 <- all[which(all$group=="1"),]
clust2 <- all[which(all$group=="2"),]
clust3 <- all[which(all$group=="3"),]
clust4 <- all[which(all$group=="4"),]


all$intel <- all$IQ

top_measures <- all %>%
  select(ID, clin_diagnosis, group, group_sex, age, GENDER, R_parstriangularis_thickavg, R_insula_thickavg, ADHD_I_SUB, R_middletemporal_thickavg, L_supramarginal_thickavg, R_inferiortemporal_thickavg, L_middletemporal_thickavg, R_superiorfrontal_thickavg, L_rostralmiddlefrontal_thickavg, L_insula_thickavg, L_fusiform_thickavg, L_parsopercularis_thickavg, 
         R_inferiorparietal_thickavg, R_rostralmiddlefrontal_thickavg, L_inferiorparietal_thickavg, L_inferiortemporal_thickavg, R_precentral_thickavg, L_lateralorbitofrontal_thickavg,
         L_lateraloccipital_thickavg, L_lingual_thickavg, R_supramarginal_thickavg, L_parstriangularis_thickavg, R_fusiform_thickavg, R_lateralorbitofrontal_thickavg, L_postcentral_thickavg, ADHD_HI_SUB, 
         R_lateraloccipital_thickavg, L_superiorfrontal_thickavg, R_medialorbitofrontal_thickavg, L_pericalcarine_thickavg, L_caudalmiddlefrontal_thickavg, R_caudalmiddlefrontal_thickavg, R_postcentral_thickavg, L_parsorbitalis_thickavg, Rpal) %>%
  gather(Measure, value, -ID, -age, -clin_diagnosis, -group, -group_sex, -GENDER)

top_10 <- all %>%
  select(ID, clin_diagnosis, group, group_sex, age, GENDER, R_parstriangularis_surfavg, R_insula_surfavg, R_middletemporal_surfavg, L_supramarginal_surfavg, R_inferiortemporal_surfavg, L_middletemporal_surfavg, R_superiorfrontal_surfavg, L_rostralmiddlefrontal_surfavg, L_insula_surfavg, L_fusiform_surfavg) %>%
  gather(Measure, value, -ID, -age, -clin_diagnosis, -group, -group_sex, -GENDER)

basic_demos <- all[,c("ID", "clin_diagnosis", "intel", "group", "group_sex", "CB68IPTOT", "age", "GENDER", "CB68EPTOT", "TPOCS_TOT", "ADHD_I_SUB", "ADHD_HI_SUB", "RBSALLT", "SCQTOT")]
FA_stacked <- cbind(basic_demos, FA)
CT_stacked <- cbind(basic_demos, CT)
Vol_stacked <- cbind(basic_demos, Volume)

FA_stacked <- FA_stacked %>%
  select(ID, clin_diagnosis, group, intel, group_sex, CB68IPTOT, age, GENDER, CB68EPTOT, TPOCS_TOT, ADHD_I_SUB, ADHD_HI_SUB, RBSALLT, SCQTOT, ends_with(".R"), ends_with(".L")) %>%
  gather(Tract, FA, -ID, -age, -clin_diagnosis,-intel, -group, -group_sex, -GENDER, -CB68EPTOT, -TPOCS_TOT, -ADHD_I_SUB, -ADHD_HI_SUB, -RBSALLT, -SCQTOT, -CB68IPTOT)

CT_stacked <- all %>%
  select(ID, clin_diagnosis, group, group_sex, CB68IPTOT, age, GENDER, CB68EPTOT, TPOCS_TOT, ADHD_I_SUB, ADHD_HI_SUB, RBSALLT, SCQTOT, ends_with("thickavg")) %>%
  gather(Region, thickness, -ID, -age, -clin_diagnosis, -group, -group_sex, -GENDER, -CB68EPTOT, -TPOCS_TOT, -ADHD_I_SUB, -ADHD_HI_SUB, -RBSALLT, -SCQTOT, -CB68IPTOT)

SA_stacked <- all %>%
  select(ID, clin_diagnosis, group, IQ, group_sex, CB68IPTOT, age, GENDER, CB68EPTOT, TPOCS_TOT, ADHD_I_SUB, ADHD_HI_SUB, RBSALLT, SCQTOT, ends_with("surfavg")) %>%
  gather(Region, area, -ID, -age, -IQ, -clin_diagnosis, -group, -group_sex, -GENDER, -CB68EPTOT, -TPOCS_TOT, -ADHD_I_SUB, -ADHD_HI_SUB, -RBSALLT, -SCQTOT, -CB68IPTOT)

Vol_stacked <- Vol_stacked %>%
  select(ID, clin_diagnosis, group, group_sex, CB68IPTOT, age, GENDER, CB68EPTOT, TPOCS_TOT, ADHD_I_SUB, ADHD_HI_SUB, RBSALLT, SCQTOT, starts_with("R"), starts_with("L")) %>%
  gather(Region, volume, -ID, -age, -clin_diagnosis, -group, -group_sex, -GENDER, -CB68EPTOT, -TPOCS_TOT, -ADHD_I_SUB, -ADHD_HI_SUB, -RBSALLT, -SCQTOT, -CB68IPTOT)

all_features <- all %>%
  select(ID, clin_diagnosis, group, group_sex, CB68IPTOT, age, GENDER, CB68EPTOT, TPOCS_TOT, ADHD_I_SUB, ADHD_HI_SUB, RBSALLT, SCQTOT, ends_with("thickavg"), ends_with(".R"), ends_with(".L"), Lthal, Rthal, Lcaud, Rcaud, Lput, Rput, Lpal, Rpal, Lhippo, Rhippo, Lamyg, Ramyg, Laccumb, Raccumb, BCC, SCC, GCC, FX) %>%
  gather(Measure, thickness, -ID, -age, -clin_diagnosis, -group, -group_sex, -GENDER)

symps <- all
symps$CB68IPTOT <- (symps$CB68IPTOT - mean(symps$CB68IPTOT, na.rm = TRUE))/sd(symps$CB68IPTOT, na.rm = TRUE)
symps$CB68EPTOT <- (symps$CB68EPTOT - mean(symps$CB68EPTOT, na.rm = TRUE))/sd(symps$CB68EPTOT, na.rm = TRUE)
symps$TPOCS_TOT <- (symps$TPOCS_TOT - mean(symps$TPOCS_TOT, na.rm = TRUE))/sd(symps$TPOCS_TOT, na.rm = TRUE)
symps$ADHD_I_SUB <- (symps$ADHD_I_SUB - mean(symps$ADHD_I_SUB, na.rm = TRUE))/sd(symps$ADHD_I_SUB, na.rm = TRUE)
symps$ADHD_HI_SUB <- (symps$ADHD_HI_SUB - mean(symps$ADHD_HI_SUB, na.rm = TRUE))/sd(symps$ADHD_HI_SUB, na.rm = TRUE)
symps$RBSALLT <- (symps$RBSALLT - mean(symps$RBSALLT, na.rm = TRUE))/sd(symps$RBSALLT, na.rm = TRUE)
symps$SCQTOT <- (symps$SCQTOT - mean(symps$SCQTOT, na.rm = TRUE))/sd(symps$SCQTOT, na.rm = TRUE)

symptoms_stacked <- symps %>%
  select(ID, clin_diagnosis, intel, group, group_sex, CB68IPTOT, age, GENDER, Rpal, CB68EPTOT, TPOCS_TOT, ADHD_I_SUB, ADHD_HI_SUB, RBSALLT, SCQTOT, L_insula_thickavg, R_superiorfrontal_thickavg, R_rostralmiddlefrontal_thickavg) %>%
  gather(Symptom, severity, -ID,-intel, -age, -clin_diagnosis, -group, -group_sex,-Rpal, -GENDER, -L_insula_thickavg, -R_superiorfrontal_thickavg, -R_rostralmiddlefrontal_thickavg)



