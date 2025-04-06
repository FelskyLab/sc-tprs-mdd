# required: Wittenberg_to_ensembl.R
# goal: get all covars and calc Wittenberg tPRS
# next: PCA, blood correlations, association analyses


library(readr)
library(readxl)
library(stringr)
setwd("/external/rprshnas01/kcni/mding")
gen_y_pihat <- read_csv("abcd/releases/5/1/core/genetics/gen_y_pihat.csv")
genpcs <- as.data.frame(cbind(gen_y_pihat$src_subject_id,  gen_y_pihat[, 21:30]))
colnames(genpcs)[1] <- "IID"
wnhsp <- read.delim("/external/rprshnas01/kcni/mding/abcd/genomics/genomics_sample/abcd_imputed/non_Hispanic_white_QCd/white_non_Hispanic.txt", header = TRUE, sep = "\t", dec = ".")
wnhsp$IID <- mapply(sub, ".{9}_?", "", wnhsp$IID)
# dep diagnosis not prevalent enough, use anx/dep, with/dep
# mh_p_ksads_dep <- read_csv("abcd/releases/5/0/tabular/core/mental-health/mh_p_ksads_dep.csv")
# mh_y_ksads_dep <- read_csv("abcd/releases/5/0/tabular/core/mental-health/mh_y_ksads_dep.csv")

# want src_subject_id, event_name, site_id_l, interview_age
abcd_y_lt <- read_csv("abcd/releases/5/1/core/abcd-general/abcd_y_lt.csv") 
cov1 <- as.data.frame(cbind(abcd_y_lt$src_subject_id,  abcd_y_lt$eventname, abcd_y_lt$site_id_l, 
                            abcd_y_lt$interview_age))
colnames(cov1) <- c("IID", "event_name", "site", "age")
cov1$age <- as.numeric(cov1$age)
# want src_subject_id, event_name, demo_sex_v2 (1 = male, 2 = female)
gish_p_gi <- read_csv("abcd/releases/5/1/core/gender-identity-sexual-health/gish_p_gi.csv")
cov2 <- as.data.frame(cbind(gish_p_gi$src_subject_id, gish_p_gi$demo_sex_v2))
cov2 <- cov2[which(gish_p_gi$eventname == "baseline_year_1_arm_1"),]
colnames(cov2) <- c("IID", "sex")
cov2$sex <- as.factor(cov2$sex)

#want cbcl_scr_syn_anxdep_t and cbcl_scr_syn_withdep_t
mh_p_cbcl <- read_csv("abcd/releases/5/1/core/mental-health/mh_p_cbcl.csv")
dep <- as.data.frame(cbind(mh_p_cbcl$src_subject_id,  mh_p_cbcl$eventname, mh_p_cbcl$cbcl_scr_syn_anxdep_t, 
                           mh_p_cbcl$cbcl_scr_syn_withdep_t))
colnames(dep) <- c("IID", "event_name", "anxdep_t", "withdep_t")
dep$anxdep_t <- as.numeric(dep$anxdep_t)
dep$withdep_t <- as.numeric(dep$withdep_t)

ph_y_bld <- read_csv("abcd/releases/5/1/core/physical-health/ph_y_bld.csv")
bld <- as.data.frame(cbind(ph_y_bld$src_subject_id,  ph_y_bld$eventname, 
                           ph_y_bld$biospec_blood_wbc_count,
                           ph_y_bld$biospec_blood_neut_abs/ph_y_bld$biospec_blood_lymph_abs,
                           ph_y_bld$biospec_blood_baso_abs, 
                           ph_y_bld$biospec_blood_eos_abs, ph_y_bld$biospec_blood_mono_abs,
                           ph_y_bld$biospec_blood_neut_abs, 
                           ph_y_bld$biospec_blood_imm_gran_abs, 
                           ph_y_bld$biospec_blood_plt_count,
                           ph_y_bld$biospec_blood_plt_count/ph_y_bld$biospec_blood_lymph_abs,
                           ph_y_bld$biospec_blood_cholesterol,
                           ph_y_bld$biospec_blood_hdl_cholesterol,
                           ph_y_bld$biospec_blood_rbc_count,
                           ph_y_bld$biospec_blood_ferritin, 
                           ph_y_bld$biospec_blood_mchc,
                           ph_y_bld$biospec_blood_mpv, ph_y_bld$biospec_blood_nrbc_abs,
                           ph_y_bld$biospec_blood_rdw, 
                           ph_y_bld$biospec_blood_hematocrit))
colnames(bld) <- c("IID", "event_name", "wbc", "nlr", "baso", "eos", "mono", "neut", "imm_gran", "plt",
                   "plr", "chol", "hdl_chol", "rbc", "ferr", "mchc", "mpv", "nrbc", "rdw", "hemcrit")
for (i in 3:ncol(bld)){
  bld[[i]] <- as.numeric(bld[[i]])
}

all <- merge(cov1, cov2, by = c("IID"), sort = FALSE)
all <- merge(all, dep, by = c("IID", "event_name"), sort = FALSE)
all <- merge(all, bld, by = c("IID", "event_name"), sort = FALSE, all = TRUE)
all <- merge(all, genpcs, by = c("IID"), sort = FALSE)
# Due to the observed skewness of these CBCL outcomes, we square-root transformed t-scores 
# prior to modeling. To identify subsets of participants with clinically relevant levels of 
# problematic symptoms and behaviours, we used cut off t-scores of <65 (defined as normal), 
# 65–69 (borderline), and ≥70 (clinically significant) (Fig. 1), as in previous work (Wainberg et al., 2022)

hist(dep$anxdep_t)
hist(dep$withdep_t)

white_all <- merge(wnhsp, all, by = c("IID"), sort = FALSE)
witt_deg<- read_excel("/external/rprshnas01/kcni/mding/sc-tprs-mdd/Wittenberg_deg.xlsx")

library(data.table)
cell_files <- list.dirs(path="~/abcd/workspace/mding/abcd_imputed/with_sc_blood38", full.names=TRUE, recursive=FALSE)
cell_types <- basename(cell_files)
tPRS_uw <- data.frame(IID = unique(white_all$IID))
tPRS_uw[, cell_types] <- NA
tPRS_uw$IID <- unique(white_all$IID)
tPRS_w <- tPRS_uw

#filename <- commandArgs(trailingOnly = TRUE)
tPRS_calc <- function(filename, goi = witt_deg$ensembl_gene_id){
  #filename <- "/nethome/kcni/mding/abcd/workspace/mding/abcd_imputed/with_sc_blood38/transitional_stage_B_cell"
  df <- fread(paste0(filename,"/abcd_pred_", 1, ".txt"), header = TRUE)
  df$IID <- mapply(sub, ".{9}_?", "", df$IID)
  overlap_genes <- intersect(goi, colnames(df))
  df_sub <- df[(df$IID %in%  unique(white_all$IID)), ]
  df_sub <- df_sub[, c("FID", "IID", ..overlap_genes)]
  uw_tprs <- data.frame(FID = df_sub$FID, IID = df_sub$IID)
  w_tprs <- uw_tprs
  g_weights <- vector(length = length(overlap_genes))
  for (i in 1:(ncol(df_sub) - 2)){
    g_weights[i] <- witt_deg$b[which(goi == colnames(df_sub)[i+2])]
  }
  
  for (i in 1:nrow(df_sub)){
    
    uw_tprs$score[i] <- sum(df_sub[i, ][ ,-c(1, 2)])
    w_tprs$score[i] <- sum(df_sub[i, ][ ,-c(1, 2)] * g_weights)
  }
  uw_df <- as.data.frame(uw_tprs$score)
  colnames(uw_df) <- basename(filename)
  w_df <- as.data.frame(w_tprs$score)
  colnames(w_df) <- basename(filename)
  
  return(list(uw_df, w_df))
}


gene_ct_overlap <- function(files){
  genes_by_ct <- vector(mode = "list", length = length(files))
  for (i in 1:length(files)){
    con <- file(paste0(files[i], "/abcd_pred_1.txt"),"r")
    first_line <- readLines(con,n=1)
    close(con)
    first_line <- as.vector(str_split_1(first_line, "\t"))
    overlap_genes <- intersect(witt_deg$ensembl_gene_id, first_line)
    genes_by_ct[[i]] <- overlap_genes
  }
  names(genes_by_ct) <- basename(files)
  return (genes_by_ct)
}

for (f in cell_files){
  a <- tPRS_calc(f)
  tPRS_uw <- cbind(tPRS_uw, a[[1]])
  tPRS_w <- cbind(tPRS_w, a[[2]])
}

colnames(tPRS_uw)[(2:(length(cell_types) + 1))] <- 
  paste0(colnames(tPRS_uw)[2:(length(cell_types) + 1)], "_uw")
colnames(tPRS_w)[(2:(length(cell_types) + 1))] <- 
  paste0(colnames(tPRS_w)[2:(length(cell_types) + 1)], "_w")
og_white_all_tprs <- merge(white_all, tPRS_uw, by = c("IID"), sort = FALSE)

og_white_all_tprs <- merge(og_white_all_tprs, tPRS_w, by = c("IID"), sort = FALSE)

overlaps <- gene_ct_overlap(cell_files)


save(og_white_all_tprs,file="sc-tprs-mdd/og_white_all_tprs_immune.Rda")
save(overlaps,file="sc-tprs-mdd/sc_tPRS_gene_overlaps_immune.Rda")




# #install.packages("UpSetR")
# library(UpSetR)
# upset(fromList(overlaps),mb.ratio = c(0.40, 0.60), nsets = 22, nintersects = 23,
#       text.scale = c(1.3, 1.3, 1, 1, 1.6, 1.4), order.by = "freq")
