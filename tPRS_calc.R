library(tidyverse)
setwd("/external/rprshnas01/kcni/mding/sc-tprs-mdd")

#Read UKBB phenotype, depression, and mean FA + mean MD datasets. Merge all into one df
abcd_pheno <- read.csv("UKBB_data/UKBB_phenotype_data.csv") #502,413
ukbb_mean_fa_md <- read_csv("UKBB_data/UKBB_Mean_FA_and_Mean_MD_data.csv") #502,413
ukbb_mdd <- read.csv("UKBB_data/UKBB_phq2_MDD_diagnosis_data.csv") #39,462
ukbb_all <- merge(ukbb_pheno, ukbb_mean_fa_md, by = "eid") %>% 
  merge(ukbb_mdd) #n=39,462 

#Create age^2 column
#rename age and assessment center columns
ukbb_all <- ukbb_all %>%
  dplyr::rename(Age = Age.when.attended.assessment.centre_inst2_idx0, 
                site = UK.Biobank.assessment.centre_inst2_idx0)%>%
  mutate(Age2 = Age^2)

#rename the 10 genetic principal components 
colnames(ukbb_all)[636:645] <- c('gPC1','gPC2','gPC3','gPC4','gPC5','gPC6','gPC7','gPC8','gPC9','gPC10')
colnames(ukbb_all[, c(636:645)]) #check new names

#Remove related individuals from the dataset 
#We'll use field 22020 (people used in calculating genetic principal components) as a proxy for kinship coeff (0.20), which retains 1 subject from each related pair
keep_kinship_vars20 <- "Yes"
unrelated_ukbb_all <- ukbb_all %>% filter(Used.in.genetic.principal.components %in% keep_kinship_vars20)
table(unrelated_ukbb_all$Used.in.genetic.principal.components) # check number of participants after filter (31901/39462)

#Function to read tprs files: #Read in tPRSs and isolate eids
read_tprs_file <- function(f_path){
  col_names = c("ID", "tprs") #colnames
  f <- read_delim(paste(f_path), col_names = col_names) %>%
    mutate(`ID` = sub("^sample_", "", ID))
  return(f)
}
#Female tPRSs
Seney_tprsF_gtex <- read_tprs_file("results/Computed_tPRSs/with_GTex/Seney_final_tPRS_F.txt")
Sig_Seney_tprsF_gtex <- read_tprs_file("results/Computed_tPRSs/with_GTex/Sig_Seney_Gtex_final_tPRS_F.txt")
Labonte_tprsF_gtex  <- read_tprs_file("results/Computed_tPRSs/with_GTex/Labonte_final_tPRS_F.txt")
Mansouri_tprsF_gtex  <- read_tprs_file("results/Computed_tPRSs/with_GTex/Mansouri_final_tPRS_F.txt")
MLithwick_tprsF_gtex <- read_tprs_file("results/Computed_tPRSs/with_GTex/MLithwick_final_tPRS_F.txt")
Sig_MLithwick_tprsF_gtex <- read_tprs_file("results/Computed_tPRSs/with_GTex/Sig_MLithwick_Gtex_final_tPRS_F.txt")
Maitra_tprsF_gtex <- read_tprs_file("results/Computed_tPRSs/with_GTex/Maitra_final_tPRS_F.txt")

Seney_tprsF_cmc <- read_tprs_file("results/Computed_tPRSs/with_CMC/Seney_CMC_final_tPRS_F.txt")
Sig_Seney_tprsF_cmc <- read_tprs_file("results/Computed_tPRSs/with_CMC/Sig_Seney_CMC_final_tPRS_F.txt")
Labonte_tprsF_cmc <- read_tprs_file("results/Computed_tPRSs/with_CMC/Labonte_CMC_final_tPRS_F.txt")
Mansouri_tprsF_cmc <- read_tprs_file("results/Computed_tPRSs/with_CMC/Mansouri_CMC_final_tPRS_F.txt")
MLithwick_tprsF_cmc <- read_tprs_file("results/Computed_tPRSs/with_CMC/MLithwick_CMC_final_tPRS_F.txt")
Sig_MLithwick_tprsF_cmc <- read_tprs_file("results/Computed_tPRSs/with_CMC/Sig_MLithwick_CMC_final_tPRS_F.txt")
Maitra_tprsF_cmc <- read_tprs_file("results/Computed_tPRSs/with_CMC/Maitra_CMC_final_tPRS_F.txt")

#Male tPRSs
Seney_tprsM_gtex  <- read_tprs_file("results/Computed_tPRSs/with_GTex/Seney_final_tPRS_M.txt")
Sig_Seney_tprsM_gtex <- read_tprs_file("results/Computed_tPRSs/with_GTex/Sig_Seney_Gtex_final_tPRS_M.txt")
Labonte_tprsM_gtex <- read_tprs_file("results/Computed_tPRSs/with_GTex/Labonte_final_tPRS_M.txt")
Mansouri_tprsM_gtex  <- read_tprs_file("results/Computed_tPRSs/with_GTex/Mansouri_final_tPRS_M.txt")
MLithwick_tprsM_gtex <- read_tprs_file("results/Computed_tPRSs/with_GTex/MLithwick_final_tPRS_M.txt")
Sig_MLithwick_tprsM_gtex <- read_tprs_file("results/Computed_tPRSs/with_GTex/Sig_MLithwick_Gtex_final_tPRS_M.txt")
Maitra_tprsM_gtex <- read_tprs_file("results/Computed_tPRSs/with_GTex/Maitra_final_tPRS_M.txt")

Seney_tprsM_cmc <- read_tprs_file("results/Computed_tPRSs/with_CMC/Seney_CMC_final_tPRS_M.txt")
Sig_Seney_tprsM_cmc <- read_tprs_file("results/Computed_tPRSs/with_CMC/Sig_Seney_CMC_final_tPRS_M.txt")
Labonte_tprsM_cmc <- read_tprs_file("results/Computed_tPRSs/with_CMC/Labonte_CMC_final_tPRS_M.txt")
Mansouri_tprsM_cmc <- read_tprs_file("results/Computed_tPRSs/with_CMC/Mansouri_CMC_final_tPRS_M.txt")
MLithwick_tprsM_cmc <- read_tprs_file("results/Computed_tPRSs/with_CMC/MLithwick_CMC_final_tPRS_M.txt")
Sig_MLithwick_tprsM_cmc <- read_tprs_file("results/Computed_tPRSs/with_CMC/Sig_MLithwick_CMC_final_tPRS_M.txt")
Maitra_tprsM_cmc <- read_tprs_file("results/Computed_tPRSs/with_CMC/Maitra_CMC_final_tPRS_M.txt")


#Merge tPRSs with ukbb df
#Create a named list of the tprs dfs first
all_tprs <- list("Seney_tprsF_gtex" = Seney_tprsF_gtex, "Sig_Seney_tprsF_gtex" = Sig_Seney_tprsF_gtex,
                 "Labonte_tprsF_gtex" = Labonte_tprsF_gtex, "Mansouri_tprsF_gtex" = Mansouri_tprsF_gtex, 
                 "MLithwick_tprsF_gtex" = MLithwick_tprsF_gtex, "Sig_MLithwick_tprsF_gtex" = Sig_MLithwick_tprsF_gtex,
                 "Maitra_tprsF_gtex" = Maitra_tprsF_gtex, "Seney_tprsF_cmc" = Seney_tprsF_cmc, 
                 "Sig_Seney_tprsF_cmc" = Sig_Seney_tprsF_cmc, "Labonte_tprsF_cmc" = Labonte_tprsF_cmc, 
                 "Mansouri_tprsF_cmc" = Mansouri_tprsF_cmc, "MLithwick_tprsF_cmc" = MLithwick_tprsF_cmc,
                 "Sig_MLithwick_tprsF_cmc" = Sig_MLithwick_tprsF_cmc, "Maitra_tprsF_cmc" = Maitra_tprsF_cmc,
                 
                 "Seney_tprsM_gtex" = Seney_tprsM_gtex, "Sig_Seney_tprsM_gtex" = Sig_Seney_tprsM_gtex,
                 "Labonte_tprsM_gtex" = Labonte_tprsM_gtex, "Mansouri_tprsM_gtex" = Mansouri_tprsM_gtex, 
                 "MLithwick_tprsM_gtex" = MLithwick_tprsM_gtex, "Sig_MLithwick_tprsM_gtex" = Sig_MLithwick_tprsM_gtex,
                 "Maitra_tprsM_gtex" = Maitra_tprsM_gtex, "Seney_tprsM_cmc" = Seney_tprsM_cmc,
                 "Sig_Seney_tprsM_cmc" = Sig_Seney_tprsM_cmc, "Labonte_tprsM_cmc" = Labonte_tprsM_cmc,
                 "Mansouri_tprsM_cmc" = Mansouri_tprsM_cmc, "MLithwick_tprsM_cmc" = MLithwick_tprsM_cmc,
                 "Sig_MLithwick_tprsM_cmc" =  Sig_MLithwick_tprsM_cmc, "Maitra_tprsM_cmc" = Maitra_tprsM_cmc)

#Merge all the tPRSs into one df then add their column names
all_tprs_df <- all_tprs %>% purrr::reduce(inner_join, by='ID')
names(all_tprs_df) <- c("eid", names(all_tprs)) #Add "eid" to match UKBB df
#Confirm there are no NAs
sum(is.na(all_tprs_df))
#Check the top of each column to ensure there aren't any duplicated columns
sum(duplicated(unlist(all_tprs_df[1, ])))

#Merge tPRS and gwas-prs with ukbb df
gwas_prs <- read_tsv("UKBB_data/raw_score_GWAS_PRS_depression.tsv")
ukbb_and_scores <- merge(unrelated_ukbb_all, all_tprs_df, by = "eid") %>% #31901 eids merged
  merge(., gwas_prs, by.x = "eid", by.y = "IID")
#We will still filter to retain only those with complete neuroimaging data but export this just in case 
write.csv(ukbb_and_scores, "UKBB_data/Merged_UKBB_phenotype_MDD_w_scores.csv", row.names = F)

#Filter to retain only individuals with complete brain neuroimaging fields at instance 2 (first neuroimaging visit)
#Imaging data: https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=192 (CT,SA,CV,SCV)
ukbb_imaging <- ukbb_and_scores %>% filter(Area.of.lingual..left.hemisphere._inst2_idx0 != "NA") #28,377
colnames(ukbb_imaging) # inspect final dataframe


#Sex stratification
table(ukbb_imaging$Sex) #Check distribution: 14,802 F and 13,575 M
f_ukbb_imaging <- ukbb_imaging %>% filter(Sex == "Female")
m_ukbb_imaging <- ukbb_imaging %>% filter(Sex == "Male")


#Export dfs
write.csv(f_ukbb_imaging, "UKBB_data/Female_Sex_Stratified_UKBB_phenotype_w_scores.csv", row.names = F)
write.csv(m_ukbb_imaging, "UKBB_data/Male_Sex_Stratified_UKBB_phenotype_w_scores.csv", row.names = F)



