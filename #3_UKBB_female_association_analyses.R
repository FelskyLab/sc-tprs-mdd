### This code calculates associations between the tPRSs and 
#surface area, cortical thickness, cortical volume = https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=192
#subcortical volume = https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=1102
#mean FA =
#mean MD =
#depressive symptoms =
#and MDD diagnosis =
#in FEMALE individuals

library(tidyverse)
library(dplyr)
setwd("/external/rprshnas01/netdata_kcni/dflab/team/eoa/tPRS_analysis")
female_df <- read.csv("UKBB_data/Female_Sex_Stratified_UKBB_phenotype_w_scores.csv")

#Create a LM function to test the main effects of tprs on brain phenotype/structure
lm_func <- function(df, tprs_colno, region_colno, phenotype, effect = "main", covariates = NULL){
  df_names <- colnames(df)
  summ <- c() #to return lm summary 
  results <- list() #to return curated results
  
  for (t in tprs_colno){   
    for (r in region_colno){
      tprs <- df[, t]
      region <- df[, r]
      #Create lm formula with the option to add a covariate
      if (is.null(covariates)) {
        formula <- as.formula("region ~ tprs + Age + Age2 + gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + gPC6 + gPC7 + gPC8 + gPC9 + gPC10 + site")
      } else {
        formula <- as.formula(paste("region ~ tprs + Age + Age2 + gPC1 + gPC2 + gPC3 + gPC4 + gPC5 + gPC6 + gPC7 + gPC8 + gPC9 + gPC10 + site +", covariates))
      }
      #Fit the model with the formula
      lm_run <- lm(formula, data=df, na.action=na.omit)
      summ <- summary(lm_run)
      estimate <- summ$coefficients[2,1]
      std_err <- summ$coefficients[2,2] 
      tval <- summ$coefficients[2,3]
      pval <- summ$coefficients[2,4]
      results <- append(results, list(data.frame(score=names(df)[t], 
                                                 `region/variable` =names(df)[r], 
                                                 phenotype=phenotype, 
                                                 effect=effect, 
                                                 tval=round(tval,4),
                                                 pval=round(pval,4), 
                                                 estimate=estimate,
                                                 std_err=std_err,
                                                 n=nrow(df)-sum(is.na(df[,t])), 
                                                 check.names = F)
      )
      )
      summ <- append(summ, summ)
    }
  }

  #Turn the list of dfs into one df
  results <- bind_rows(results)
  #Apply FDR correction to each tprs 
  tprsnames <- df_names[tprs_colno] #tprs col names
  for (name in tprsnames){
    results[results$score==name, 'pFDR']  <- round(p.adjust(
      results[results$score==name, ]$pval, method='fdr'),4
    )
  }
  
  return(list(summary=summ, results=results))
}
                        
#--------**1. BASE MODELS: tPRS-F and tPRS-M on female brain structure**----
#Select female samples
df <- female_df
#Column names
df_colnames <- colnames(df) 
#define the inputs for the lm function. Use either tPRSs computed with either the gtex or CMC reference model
tprsF_cols_gtex <- c(911:917) #female gtex tPRS columns 
tprsF_cols_cmc <- c(918:924) #female cmc tPRS columns

tprsM_cols_gtex <- c(925:931) #male gtex tPRS columns 
tprsM_cols_cmc <- c(932:938) #male cmc tPRS columns

#####SURFACE AREA####
#define surface area regional columns measured at inst2 
area_cols <- df_colnames[grepl("Area", df_colnames)] #Search for "area"
area_cols <- area_cols[grep("inst2", area_cols)] #Search for inst2 among area columns
SA_colno <- which(df_colnames %in% area_cols) #Select column indices matching the SA column names
#Or the indices can be manually entered:
#SA_colno <- c(264,268,272,276,280,284,288,292,296,300,304,308,312,316,320,324,328,332,336,340,344,348,352,356,360,364,368,372,376,380,384,266,270,274,278,282,286,290,294,298,302,306,310,314,318,322,326,330,334,338,342,346,350,354,358,362,366,370,374,378,382,386)
results_SA <- lm_func(df, tprs_colno = tprsF_cols_cmc, region_colno = SA_colno, phenotype = 'surface.area')$results
write.csv(results_SA, "results/Association_results/with_CMC/CMC_UKBB_tprsF_on_female_SA_results_no_GWAS.csv",row.names = F)

results_SA_tprsM_on_F <- lm_func(df, tprs_colno = tprsM_cols_cmc, region_colno = SA_colno, phenotype = 'surface.area')$results
write.csv(results_SA_tprsM_on_F, "results/Association_results/with_CMC/CMC_UKBB_tprsM_on_female_SA_results_no_GWAS.csv",row.names = F)


#####CORTICAL THICKNESS####
#define cortical thickness regional columns measured at inst2 
thickness_cols <- df_colnames[grepl("thickness", df_colnames)] #Search for "thickness"
thickness_cols <- thickness_cols[grep("inst2", thickness_cols)] #Search for inst2 in results
CT_colno <- which(df_colnames %in% thickness_cols) #Select column indices matching CT column names
#Or the indices can be manually:
#CT_colno <- c(388,392,396,400,404,408,412,416,420,424,428,432,436,440,444,448,452,456,460,464,468,472,476,480,484,488,492,496,500,504,508,390,394,398,402,406,410,414,418,422,426,430,434,438,442,446,450,454,458,462,466,470,474,478,482,486,490,494,498,502,506,510)
#apply arguments
results_CT <- lm_func(df, tprs_colno = tprsF_cols_cmc, region_colno = CT_colno, phenotype = 'cortical.thickness')$results
write.csv(results_CT, "results/Association_results/with_CMC/CMC_UKBB_tprsF_on_female_CT_results_no_GWAS.csv",row.names = F)

results_CT_tprsM_on_F <- lm_func(df, tprs_colno = tprsM_cols_cmc, region_colno = CT_colno, phenotype = 'cortical.thickness')$results
write.csv(results_CT_tprsM_on_F, "results/Association_results/with_CMC/CMC_UKBB_tprsM_on_female_CT_results_no_GWAS.csv",row.names = F)


#####CORTICAL VOLUME####
#define cortical volume regional columns measured at inst2 
volume_cols <- df_colnames[grepl("Volume", df_colnames)] #Search for "Volume"
volume_cols <- volume_cols[grep("inst2", volume_cols)] #Search for inst2 in results
CV_colno <- which(df_colnames %in% volume_cols) #Select column indices matching CV column names
#Some of the columns represent subcortical regions so keep only columns below the 700th position (the cortical ones)
CV_colno <- CV_colno[CV_colno < 700]
#Or the indices can be manually entered as done below:
#CV_colno <- c(512,516,520,524,528,532,536,540,544,548,552,556,560,564,568,572,576,582,584,588,592,596,600,604,608,612,616,620,624,628,632,514,518,522,526,530,534,538,542,546,550,554,558,562,566,570,574,578,580, 586,590,594,598,602,606,610,614,618,622,626,630,634)
#apply arguments
results_CV <- lm_func(df, tprs_colno=tprsF_cols_cmc, region_colno = CV_colno, phenotype = 'cortical.volume')$results
write.csv(results_CV, "results/Association_results/with_CMC/CMC_UKBB_tprsF_on_female_CV_results_no_GWAS.csv",row.names = F)

results_CV_tprsM_on_F <- lm_func(df, tprs_colno=tprsM_cols_cmc, region_colno = CV_colno, phenotype = 'cortical.volume')$results
write.csv(results_CV_tprsM_on_F, "results/Association_results/with_CMC/CMC_UKBB_tprsM_on_female_CV_results_no_GWAS.csv",row.names = F)


#####SUBCORTICAL VOLUME####
#eTIV (Volume.ratio.of.BrainSegVol.to.eTIV..whole.brain._inst2_idx0) will be added here
#Define subcortical volume regional columns from previous columns containing 'volume' and 'inst2' 
SCV_colno <- which(df_colnames %in% volume_cols) 
#Keep only columns above the 700th position (the subcortical ones)
SCV_colno <- SCV_colno[SCV_colno > 700]
df_colnames[SCV_colno] #the last three are not scv so remove below
SCV_colno <- SCV_colno[1:14] 
#Or the indices can be manually entered:
# SCV_colno <- c(754,758,762,766,770,774,778,756,760,764,768,772,776,780)

#Run lm with eTIV covariate
results_SCV <- lm_func(df, tprs_colno = tprsF_cols_cmc, region_colno = SCV_colno, phenotype = 'subcortical.volume', covariates = 'Volume.ratio.of.BrainSegVol.to.eTIV..whole.brain._inst2_idx0')
#Inspect covariates for eTIV
results_SCV$summary$coefficients
#res <- results_SCV$results
write.csv(results_SCV$results, "results/Association_results/with_CMC/CMC_UKBB_tprsF_on_female_SCV_results_no_GWAS.csv",row.names = F)

results_SCV_tprsM_on_F <- lm_func(df, tprs_colno = tprsM_cols_cmc, region_colno = SCV_colno, phenotype = 'subcortical.volume', covariates = 'Volume.ratio.of.BrainSegVol.to.eTIV..whole.brain._inst2_idx0')
write.csv(results_SCV_tprsM_on_F$results, "results/Association_results/with_CMC/CMC_UKBB_tprsM_on_female_SCV_results_no_GWAS.csv",row.names = F)


#####MEAN DIFFUSIVITY - MD####
#define MD regional columns 
md_cols <- df_colnames[grepl("Mean.MD", df_colnames)] #Search for mean MD cols
MD_colno <- which(df_colnames %in% md_cols)
#run lm
results_MD <- lm_func(df, tprs_colno = tprsF_cols_cmc, region_colno = MD_colno, phenotype = 'mean.MD')$results
write.csv(results_MD, "results/Association_results/with_CMC/CMC_UKBB_tprsF_on_female_meanMD_results_no_GWAS.csv",row.names = F)

results_MD_tprsM_on_F <- lm_func(df, tprs_colno = tprsM_cols_cmc, region_colno = MD_colno, phenotype = 'mean.MD')$results
write.csv(results_MD_tprsM_on_F, "results/Association_results/with_CMC/CMC_UKBB_tprsM_on_female_meanMD_results_no_GWAS.csv",row.names = F)

#####FRACTIONAL ANISOTROPY####
#define FA regional columns  
fa_cols <- df_colnames[grepl("Mean.FA", df_colnames)] #Search for mean FA cols
FA_colno <- which(df_colnames %in% fa_cols)
results_FA <- lm_func(df, tprs_colno = tprsF_cols_cmc, region_colno = FA_colno, phenotype = 'mean.FA')$results
write.csv(results_FA, "results/Association_results/with_CMC/CMC_UKBB_tprsF_on_female_meanFA_results_no_GWAS.csv",row.names = F)

results_FA_tprsM_on_F <- lm_func(df, tprs_colno = tprsM_cols_cmc, region_colno = FA_colno, phenotype = 'mean.FA')$results
write.csv(results_FA_tprsM_on_F, "results/Association_results/with_CMC/CMC_UKBB_tprsM_on_female_meanFA_results_no_GWAS.csv",row.names = F)

#####DEPRESSIVE SYMPTOMS####
#define DS columns  
symptoms_colno <- c(906,907)
results_DS <- lm_func(df, tprs_colno = tprsF_cols_cmc, region_colno = symptoms_colno, phenotype = 'depressive.symptoms')$results
write.csv(results_DS, "results/Association_results/with_CMC/CMC_UKBB_tprsF_on_female_depressive_symptoms_results_no_GWAS.csv",row.names = F)

results_DS_tprsM_on_F <- lm_func(df, tprs_colno = tprsM_cols_cmc, region_colno = symptoms_colno, phenotype = 'depressive.symptoms')$results
write.csv(results_DS_tprsM_on_F, "results/Association_results/with_CMC/CMC_UKBB_tprsM_on_female_depressive_symptoms_results_no_GWAS.csv",row.names = F)

#####MDD DIAGNOSIS####
#define diagnosis columns: F32 and F33  
diagnosis_colno <- c(908,909)
results_MDD_diag <- lm_func(df, tprs_colno = tprsF_cols_cmc, region_colno = diagnosis_colno, phenotype = 'MDD.diagnosis')$results
write.csv(results_MDD_diag, "results/Association_results/with_CMC/CMC_UKBB_tprsF_on_female_MDD_diagnosis_results_no_GWAS.csv",row.names = F)

results_MDD_diag_tprsM_on_F <- lm_func(df, tprs_colno = tprsM_cols_cmc, region_colno = diagnosis_colno, phenotype = 'MDD.diagnosis')$results
write.csv(results_MDD_diag_tprsM_on_F, "results/Association_results/with_CMC/CMC_UKBB_tprsM_on_female_MDD_diagnosis_results_no_GWAS.csv",row.names = F)



                            
#--------**2. GWAS-PRS MODELS: tPRS-F and tPRS-M on female brain structure**-------------
#GWAS-PRS column will be included as a covariate 
effect <- 'main + GWAS-PRS'
#####SURFACE AREA####
results_SA_PRS <- lm_func(df, tprs_colno = tprsF_cols_cmc, effect, region_colno = SA_colno, phenotype = 'surface.area', covariates = 'PRS')
#Inspect covariates for GWAS-PRS
results_SA_PRS$summary$coefficients
results_SA_PRS$results
write.csv(results_SA_PRS$results, "results/Association_results/with_CMC/CMC_UKBB_tprsF_on_female_SA_results_w_GWAS.csv",row.names = F)

results_SA_tprsM_on_F_w_PRS <- lm_func(df, tprs_colno = tprsM_cols_cmc, effect, region_colno = SA_colno, phenotype = 'surface.area', covariates = 'PRS')
write.csv(results_SA_tprsM_on_F_w_PRS$results, "results/Association_results/with_CMC/CMC_UKBB_tprsM_on_female_SA_results_w_GWAS.csv",row.names = F)

#####CORTICAL THICKNESS####
results_CT_PRS <- lm_func(df, tprs_colno = tprsF_cols_cmc, effect, region_colno = CT_colno, phenotype = 'cortical.thickness', covariates = 'PRS')$results
write.csv(results_CT_PRS, "results/Association_results/with_CMC/CMC_UKBB_tprsF_on_female_CT_results_w_GWAS.csv",row.names = F)

results_CT_tprsM_on_F_w_PRS <- lm_func(df, tprs_colno = tprsM_cols_cmc, effect, region_colno = CT_colno, phenotype = 'cortical.thickness', covariates = 'PRS')$results
write.csv(results_CT_tprsM_on_F_w_PRS, "results/Association_results/with_CMC/CMC_UKBB_tprsM_on_female_CT_results_w_GWAS.csv",row.names = F)

#####CORTICAL VOLUME####
results_CV_PRS <- lm_func(df, tprs_colno = tprsF_cols_cmc, effect, region_colno = CV_colno, phenotype = 'cortical.volume', covariates = 'PRS')$results
write.csv(results_CV_PRS, "results/Association_results/with_CMC/CMC_UKBB_tprsF_on_female_CV_results_w_GWAS.csv",row.names = F)

results_CV_tprsM_on_F_w_PRS <- lm_func(df, tprs_colno = tprsM_cols_cmc, effect, region_colno = CV_colno, phenotype = 'cortical.volume', covariate = 'PRS')$results
write.csv(results_CV_tprsM_on_F_w_PRS, "results/Association_results/with_CMC/CMC_UKBB_tprsM_on_female_CV_results_w_GWAS.csv",row.names = F)

#####SUBCORTICAL VOLUME####
#Add eTIV too: Volume.ratio.of.BrainSegVol.to.eTIV..whole.brain._inst2_idx0
results_SCV_PRS <- lm_func(df, tprs_colno = tprsF_cols_cmc, effect, region_colno = SCV_colno, phenotype = 'subcortical.volume', covariates = 'Volume.ratio.of.BrainSegVol.to.eTIV..whole.brain._inst2_idx0 + PRS')
#Inspect covariates for eTIV and PRS
results_SCV_PRS$summary$coefficients
write.csv(results_SCV_PRS$results, "results/Association_results/with_CMC/CMC_UKBB_tprsF_on_female_SCV_results_w_GWAS.csv",row.names = F)

results_SCV_tprsM_on_F_w_PRS <- lm_func(df, tprs_colno = tprsM_cols_cmc, effect, region_colno = SCV_colno, phenotype = 'subcortical.volume', covariates = 'Volume.ratio.of.BrainSegVol.to.eTIV..whole.brain._inst2_idx0 + PRS')
#Inspect covariates for eTIV and PRS
results_SCV_tprsM_on_F_w_PRS$summary$coefficients
write.csv(results_SCV_tprsM_on_F_w_PRS$results, "results/Association_results/with_CMC/CMC_UKBB_tprsM_on_female_SCV_results_w_GWAS.csv",row.names = F)

#####MEAN DIFFUSIVITY####
results_MD_PRS <- lm_func(df, tprs_colno = tprsF_cols_cmc, effect, region_colno = MD_colno, phenotype = 'mean.MD', covariates = 'PRS')$results
write.csv(results_MD_PRS, "results/Association_results/with_CMC/CMC_UKBB_tprsF_on_female_meanMD_results_w_GWAS.csv",row.names = F)

results_MD_tprsM_on_F_w_PRS <- lm_func(df, tprs_colno = tprsM_cols_cmc, effect, region_colno = MD_colno, phenotype = 'mean.MD', covariates = 'PRS')$results
write.csv(results_MD_tprsM_on_F_w_PRS, "results/Association_results/with_CMC/CMC_UKBB_tprsM_on_female_meanMD_results_w_GWAS.csv",row.names = F)

#####FRACTIONAL ANISOTROPY####
results_FA_PRS <- lm_func(df, tprs_colno = tprsF_cols_cmc, effect, region_colno = FA_colno, phenotype = 'mean.FA', covariates = 'PRS')$results
write.csv(results_FA_PRS, "results/Association_results/with_CMC/CMC_UKBB_tprsF_on_female_meanFA_results_w_GWAS.csv",row.names = F)

results_FA_tprsM_on_F_w_PRS <- lm_func(df, tprs_colno = tprsM_cols_cmc, effect, region_colno = FA_colno, phenotype = 'mean.FA', covariates = 'PRS')$results
write.csv(results_FA_tprsM_on_F_w_PRS, "results/Association_results/with_CMC/CMC_UKBB_tprsM_on_female_meanFA_results_w_GWAS.csv",row.names = F)


#####DEPRESSIVE SYMPTOMS####
#define DS columns  
results_DS_PRS <- lm_func(df, tprs_colno = tprsF_cols_cmc, effect, region_colno = symptoms_colno, phenotype = 'depressive.symptoms', covariates = 'PRS')$results
write.csv(results_DS_PRS, "results/Association_results/with_CMC/CMC_UKBB_tprsF_on_female_depressive_symptoms_results_w_GWAS.csv",row.names = F)

results_DS_tprsM_on_F_w_PRS <- lm_func(df, tprs_colno = tprsM_cols_cmc, effect, region_colno = symptoms_colno, phenotype = 'depressive.symptoms', covariates = 'PRS')$results
write.csv(results_DS_tprsM_on_F_w_PRS, "results/Association_results/with_CMC/CMC_UKBB_tprsM_on_female_depressive_symptoms_results_w_GWAS.csv",row.names = F)

#####MDD DIAGNOSIS####
#define DS columns  
results_MDD_diag_PRS <- lm_func(df, tprs_colno = tprsF_cols_cmc, effect, region_colno = diagnosis_colno, phenotype = 'MDD.diagnosis', covariates = 'PRS')$results
write.csv(results_MDD_diag_PRS, "results/Association_results/with_CMC/CMC_UKBB_tprsF_on_female_MDD_diagnosis_results_W_GWAS.csv",row.names = F)

results_MDD_diag_tprsM_on_F_w_PRS <- lm_func(df, tprs_colno = tprsM_cols_cmc, effect, region_colno = diagnosis_colno, phenotype = 'MDD.diagnosis', covariates = 'PRS')$results
write.csv(results_MDD_diag_tprsM_on_F_w_PRS, "results/Association_results/with_CMC/CMC_UKBB_tprsM_on_female_MDD_diagnosis_results_W_GWAS.csv",row.names = F)

