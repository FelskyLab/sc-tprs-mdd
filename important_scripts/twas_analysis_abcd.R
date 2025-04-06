# required: load_abcd_data.R
# goal: train/test split TWASs by sex
# - prelim correlations with anx/dep, with/dep
# next: twas_tprs_cbcl_assoc.R

if (! requireNamespace("splitTools", quietly=TRUE)) {
  install.packages("splitTools")
}
library(splitTools)


load("sc-tprs-mdd/og_white_all_tprs_immune.Rda")

# write train/test split data to files=======================

events <- c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1", "2_year_follow_up_y_arm_1",
            "3_year_follow_up_y_arm_1", "4_year_follow_up_y_arm_1")
cell_files <- list.dirs(path="~/abcd/workspace/mding/abcd_imputed/with_sc_blood38", full.names=TRUE, recursive=FALSE)

n_ct <- length(cell_types)
tprs_start_idx <- 19
white_all_tprs <- og_white_all_tprs[which(og_white_all_tprs$event_name == events[1]), -(9:26)]
#white_all_tprs <- white_all_tprs[-(which(white_all_tprs$sex == 3)),]
white_all_tprs$sex <- droplevels(white_all_tprs$sex)
white_all_tprs <- white_all_tprs[-(which(is.na(white_all_tprs$anxdep_t))), ] 
colnames(white_all_tprs)[2] <- "FID"
white_all_tprs$anx_bin <- as.factor(as.integer(white_all_tprs$anxdep_t > 65))
white_all_tprs$with_bin <- as.factor(as.integer(white_all_tprs$withdep_t > 65))
levels(white_all_tprs$sex) <- c("M", "F")
w_uw_corr <- cor(unname(white_all_tprs[,tprs_start_idx:(tprs_start_idx+n_ct-1)]),
                 unname(white_all_tprs[,(tprs_start_idx+n_ct):(tprs_start_idx+2*n_ct-1)]))
corrplot::corrplot(w_uw_corr, method = "square")

onehot <- data.frame(male = as.numeric((white_all_tprs$sex == "M")))
for (site in sort(unique(white_all_tprs$site))[-1]){
  onehot[[site]] <- as.numeric((white_all_tprs$site == site))
  
}
strat_ad <- multi_strata(white_all_tprs[,c("sex", "anxdep_t", "anx_bin")])
strat_wd <- multi_strata(white_all_tprs[,c("sex", "withdep_t", "with_bin")])
inds_ad <- partition(strat_ad, p = c(train = 0.7, test = 0.3), seed = 330, split_into_list = FALSE)
inds_wd <- partition(strat_wd, p = c(train = 0.7, test = 0.3), seed = 331, split_into_list = FALSE)

write.table(cbind(white_all_tprs, onehot)[which(inds_ad == "train"),], 
            file = "/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_data/abcd_ad_train.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")
write.table(cbind(white_all_tprs, onehot)[which(inds_wd == "train"),], 
            file = "/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_data/abcd_wd_train.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")


for (f in cell_files){
  df <- fread(paste0(f,"/abcd_pred_", 1, ".txt"), header = TRUE)
  df$IID <- mapply(sub, ".{9}_?", "", df$IID)
  dfa <- merge(data.frame(IID = white_all_tprs[, c("IID")][inds_ad == "train"]), df, by = c("IID"))
  write.table(dfa, file = paste0(f,"/abcd_pred_ad", "train", ".txt"), 
              quote = FALSE, row.names = FALSE, sep = "\t")
  dfw <- merge(data.frame(IID = white_all_tprs[, c("IID")][inds_wd == "train"]), df, by = c("IID"))
  write.table(dfw, file = paste0(f,"/abcd_pred_wd", "train", ".txt"), 
              quote = FALSE, row.names = FALSE, sep = "\t")
  
}

# by sex
strat_ad <- multi_strata(white_all_tprs[which(white_all_tprs$sex == "M"), c("anxdep_t", "anx_bin")])
strat_wd <- multi_strata(white_all_tprs[which(white_all_tprs$sex == "M"), c("withdep_t", "with_bin")])
inds_ad <- partition(strat_ad, p = c(train = 0.7, test = 0.3), seed = 330, split_into_list = FALSE)
inds_wd <- partition(strat_wd, p = c(train = 0.7, test = 0.3), seed = 331, split_into_list = FALSE)

write.table(cbind(white_all_tprs, onehot)[which(white_all_tprs$sex == "M"),][which(inds_ad == "train"),], 
            file = "/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_data/abcd_ad_train_male.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")
write.table(cbind(white_all_tprs, onehot)[which(white_all_tprs$sex == "M"),][which(inds_wd == "train"),], 
            file = "/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_data/abcd_wd_train_male.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")


for (f in cell_files){
  df <- fread(paste0(f,"/abcd_pred_", 1, ".txt"), header = TRUE)
  df$IID <- mapply(sub, ".{9}_?", "", df$IID)
  dfa <- merge(data.frame(IID = white_all_tprs[which(white_all_tprs$sex == "M"), c("IID")][inds_ad == "train"]), df, by = c("IID"))
  write.table(dfa, file = paste0(f,"/abcd_pred_ad", "train_male", ".txt"), 
              quote = FALSE, row.names = FALSE, sep = "\t")
  dfw <- merge(data.frame(IID = white_all_tprs[which(white_all_tprs$sex == "M"), c("IID")][inds_wd == "train"]), df, by = c("IID"))
  write.table(dfw, file = paste0(f,"/abcd_pred_wd", "train_male", ".txt"), 
              quote = FALSE, row.names = FALSE, sep = "\t")
  
}

strat_ad <- multi_strata(white_all_tprs[which(white_all_tprs$sex == "F"), c("anxdep_t", "anx_bin")])
strat_wd <- multi_strata(white_all_tprs[which(white_all_tprs$sex == "F"), c("withdep_t", "with_bin")])
inds_ad <- partition(strat_ad, p = c(train = 0.7, test = 0.3), seed = 330, split_into_list = FALSE)
inds_wd <- partition(strat_wd, p = c(train = 0.7, test = 0.3), seed = 331, split_into_list = FALSE)

write.table(cbind(white_all_tprs, onehot)[which(white_all_tprs$sex == "F"), ][which(inds_ad == "train"),], 
            file = "/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_data/abcd_ad_train_female.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")
write.table(cbind(white_all_tprs, onehot)[which(white_all_tprs$sex == "F"), ][which(inds_wd == "train"),], 
            file = "/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_data/abcd_wd_train_female.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")


for (f in cell_files){
  df <- fread(paste0(f,"/abcd_pred_", 1, ".txt"), header = TRUE)
  df$IID <- mapply(sub, ".{9}_?", "", df$IID)
  dfa <- merge(data.frame(IID = white_all_tprs[which(white_all_tprs$sex == "F"), c("IID")][inds_ad == "train"]), df, by = c("IID"))
  write.table(dfa, file = paste0(f,"/abcd_pred_ad", "train_female", ".txt"), 
              quote = FALSE, row.names = FALSE, sep = "\t")
  dfw <- merge(data.frame(IID = white_all_tprs[which(white_all_tprs$sex == "F"), c("IID")][inds_wd == "train"]), df, by = c("IID"))
  write.table(dfw, file = paste0(f,"/abcd_pred_wd", "train_female", ".txt"), 
              quote = FALSE, row.names = FALSE, sep = "\t")
  
}

# run abcd_twas_dep.sh, abcd_twas_dep_by_sex.sh======================


# new tPRS w/ TWAS, gene inclusion p < 0.02 ===========================

cell_files <- list.dirs(path="~/abcd/workspace/mding/abcd_imputed/with_sc_blood38", full.names=TRUE, recursive=FALSE)
cell_types <- basename(cell_files)

twas_tPRS_calc <- function(filename, twas_file, p_thresh = 0.02){
  sig_genes <- fread(paste0(filename,"/", twas_file), header = TRUE)
  sig_genes$qval <- qvalue(sig_genes$pvalue)
  sig_genes <- as.data.frame(sig_genes[which(sig_genes$pval < p_thresh), c("gene", "effect")])

  df <- fread(paste0(filename,"/abcd_pred_",1,".txt"), header = TRUE)
  df$IID <- mapply(sub, ".{9}_?", "", df$IID)
  df_sub <- df[(df$IID %in%  unique(white_all_tprs$IID)), ]
  cols = c("FID", "IID", sig_genes$gene[(sig_genes$gene %in% colnames(df_sub))])
  df_sub <- df_sub[,..cols]
  uw_tprs <- data.frame(FID = df_sub$FID, IID = df_sub$IID, score = NA)
  
  w_tprs <- uw_tprs
  uw_tprs$score <- rowSums(df_sub[ ,-c(1, 2)])
  w_tprs$score <- rowSums(data.frame(mapply(`*`,df_sub[,-c(1,2)],sig_genes$effect,SIMPLIFY=FALSE)))
  
  uw_df <- uw_tprs[, c("IID", "score")]
  colnames(uw_df)[2] <- paste0(basename(filename), "_uw")
  w_df <- w_tprs[, c("IID", "score")]
  colnames(w_df)[2] <- paste0(basename(filename), "_w")
  return(list(uw_df, w_df))
}
twas_ad_tPRS_uw <- data.frame(IID = unique(white_all_tprs$IID))
twas_ad_tPRS_uw$IID <- unique(white_all_tprs$IID)
twas_ad_tPRS_w <- twas_ad_tPRS_uw
twas_wd_tPRS_uw <- twas_ad_tPRS_uw
twas_wd_tPRS_w <- twas_ad_tPRS_uw
# sex agnostic ad, wd trps
for (f in cell_files){
  a <- twas_tPRS_calc(f, "abcd_anxdep_assoc_covars.txt")
  twas_ad_tPRS_uw <- merge(twas_ad_tPRS_uw, a[[1]], by = c("IID"), sort = FALSE, all = TRUE)
  twas_ad_tPRS_w <- merge(twas_ad_tPRS_w, a[[2]], by = c("IID"), sort = FALSE, all = TRUE)
}
for (f in cell_files){
  a <- twas_tPRS_calc(f, "abcd_withdep_assoc_covars.txt")
  twas_wd_tPRS_uw <- merge(twas_wd_tPRS_uw, a[[1]], by = c("IID"), sort = FALSE, all = TRUE)
  twas_wd_tPRS_w <- merge(twas_wd_tPRS_w, a[[2]], by = c("IID"), sort = FALSE, all = TRUE)
}

# male tprs
for (f in cell_files){
  a <- twas_tPRS_calc(f, "abcd_anxdep_assoc_covars_male.txt")
  colnames(a[[1]])[2] <- paste0(colnames(a[[1]])[2], "_m")
  colnames(a[[2]])[2] <- paste0(colnames(a[[2]])[2], "_m")
  
  twas_ad_tPRS_uw <- merge(twas_ad_tPRS_uw, a[[1]], by = c("IID"), sort = FALSE, all = TRUE)
  twas_ad_tPRS_w <- merge(twas_ad_tPRS_w, a[[2]], by = c("IID"), sort = FALSE, all = TRUE)
}
for (f in cell_files){
  a <- twas_tPRS_calc(f, "abcd_withdep_assoc_covars_male.txt")
  colnames(a[[1]])[2] <- paste0(colnames(a[[1]])[2], "_m")
  colnames(a[[2]])[2] <- paste0(colnames(a[[2]])[2], "_m")
  
  twas_wd_tPRS_uw <- merge(twas_wd_tPRS_uw, a[[1]], by = c("IID"), sort = FALSE, all = TRUE)
  twas_wd_tPRS_w <- merge(twas_wd_tPRS_w, a[[2]], by = c("IID"), sort = FALSE, all = TRUE)
}

# female tprs
for (f in cell_files){
  a <- twas_tPRS_calc(f, "abcd_anxdep_assoc_covars_female.txt")
  colnames(a[[1]])[2] <- paste0(colnames(a[[1]])[2], "_f")
  colnames(a[[2]])[2] <- paste0(colnames(a[[2]])[2], "_f")
  
  twas_ad_tPRS_uw <- merge(twas_ad_tPRS_uw, a[[1]], by = c("IID"), sort = FALSE, all = TRUE)
  twas_ad_tPRS_w <- merge(twas_ad_tPRS_w, a[[2]], by = c("IID"), sort = FALSE, all = TRUE)
  
}
for (f in cell_files){
  a <- twas_tPRS_calc(f, "abcd_withdep_assoc_covars_female.txt")
  colnames(a[[1]])[2] <- paste0(colnames(a[[1]])[2], "_f")
  colnames(a[[2]])[2] <- paste0(colnames(a[[2]])[2], "_f")
  twas_wd_tPRS_uw <- merge(twas_wd_tPRS_uw, a[[1]], by = c("IID"), sort = FALSE, all = TRUE)
  twas_wd_tPRS_w <- merge(twas_wd_tPRS_w, a[[2]], by = c("IID"), sort = FALSE, all = TRUE)
}


twas_ad_white_all_tprs <- merge(og_white_all_tprs[,-(0:(2*n_ct-1) + 37)], twas_ad_tPRS_uw, by = c("IID"), sort = FALSE)

twas_ad_white_all_tprs <- merge(twas_ad_white_all_tprs, twas_ad_tPRS_w, by = c("IID"), sort = FALSE)

twas_wd_white_all_tprs <- merge(og_white_all_tprs[,-(0:(2*n_ct-1) + 37)], twas_wd_tPRS_uw, by = c("IID"), sort = FALSE)

twas_wd_white_all_tprs <- merge(twas_wd_white_all_tprs, twas_wd_tPRS_w, by = c("IID"), sort = FALSE)

save(twas_ad_white_all_tprs,file="sc-tprs-mdd/twas_ad_white_all_tprs.Rda")
save(twas_wd_white_all_tprs,file="sc-tprs-mdd/twas_wd_white_all_tprs.Rda")
save(twas_ad_tPRS_uw,file="sc-tprs-mdd/twas_ad_tPRS_uw.Rda")
save(twas_ad_tPRS_w,file="sc-tprs-mdd/twas_ad_tPRS_w.Rda")
save(twas_wd_tPRS_uw,file="sc-tprs-mdd/twas_wd_tPRS_uw.Rda")
save(twas_wd_tPRS_w,file="sc-tprs-mdd/twas_wd_tPRS_w.Rda")

# tPRS distr =====================
load("sc-tprs-mdd/twas_ad_tPRS_uw.Rda")
load("sc-tprs-mdd/twas_ad_tPRS_w.Rda")
load("sc-tprs-mdd/twas_wd_tPRS_uw.Rda")
load("sc-tprs-mdd/twas_wd_tPRS_w.Rda")


# general
ggplot(stack( (twas_ad_tPRS_uw[, c(2:8)]), varying = colnames(twas_ad_tPRS_uw)[2:8]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (twas_ad_tPRS_uw[, c(9:15)]), varying = colnames(twas_ad_tPRS_uw)[9:15]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (twas_ad_tPRS_uw[, c(16:23)]), varying = colnames(twas_ad_tPRS_uw)[16:23]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (twas_ad_tPRS_w[, c(2:8)]), varying = colnames(twas_ad_tPRS_w)[2:8]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (twas_ad_tPRS_w[, c(9:15)]), varying = colnames(twas_ad_tPRS_w)[9:15]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (twas_ad_tPRS_w[, c(16:23)]), varying = colnames(twas_ad_tPRS_w)[16:23]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))


ggplot(stack( (twas_wd_tPRS_uw[, c(2:8)]), varying = colnames(twas_wd_tPRS_uw)[2:8]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (twas_wd_tPRS_uw[, c(9:15)]), varying = colnames(twas_wd_tPRS_uw)[9:15]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (twas_wd_tPRS_uw[, c(16:23)]), varying = colnames(twas_wd_tPRS_uw)[16:23]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (twas_wd_tPRS_w[, c(2:8)]), varying = colnames(twas_wd_tPRS_w)[2:8]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (twas_wd_tPRS_w[, c(9:15)]), varying = colnames(twas_wd_tPRS_w)[9:15]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (twas_wd_tPRS_w[, c(16:23)]), varying = colnames(twas_wd_tPRS_w)[16:23]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
# male in male 
ggplot(stack( (twas_wd_tPRS_uw[which(white_all_tprs$sex=="M"), c(24:30)]), varying = colnames(twas_ad_tPRS_uw)[24:30]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (twas_wd_tPRS_uw[which(white_all_tprs$sex=="M"), c(31:37)]), varying = colnames(twas_ad_tPRS_uw)[31:37]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (twas_wd_tPRS_uw[which(white_all_tprs$sex=="M"), c(38:45)]), varying = colnames(twas_ad_tPRS_uw)[38:45]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))

ggplot(stack( (twas_wd_tPRS_w[which(white_all_tprs$sex=="M"), c(24:30)]), varying = colnames(twas_ad_tPRS_w)[24:30]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (twas_wd_tPRS_w[which(white_all_tprs$sex=="M"), c(31:37)]), varying = colnames(twas_ad_tPRS_w)[31:37]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (twas_wd_tPRS_w[which(white_all_tprs$sex=="M"), c(38:45)]), varying = colnames(twas_ad_tPRS_w)[38:45]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
# female in female
ggplot(stack( (twas_wd_tPRS_uw[which(white_all_tprs$sex=="F"), c(46:52)]), varying = colnames(twas_ad_tPRS_uw)[46:52]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (twas_wd_tPRS_uw[which(white_all_tprs$sex=="F"), c(53:59)]), varying = colnames(twas_ad_tPRS_uw)[53:59]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (twas_wd_tPRS_uw[which(white_all_tprs$sex=="F"), c(60:67)]), varying = colnames(twas_ad_tPRS_uw)[60:67]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))

ggplot(stack( (twas_wd_tPRS_w[which(white_all_tprs$sex=="F"), c(46:52)]), varying = colnames(twas_ad_tPRS_w)[46:52]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (twas_wd_tPRS_w[which(white_all_tprs$sex=="F"), c(53:59)]), varying = colnames(twas_ad_tPRS_w)[53:59]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (twas_wd_tPRS_w[which(white_all_tprs$sex=="F"), c(60:67)]), varying = colnames(twas_ad_tPRS_w)[60:67]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))


# male in female, female in male tPRS corr=============

# wittenberg vs twas tPRS
#anxdep
colnames(white_all_tprs) <- gsub("-", "_", colnames(white_all_tprs))
colnames(twas_ad_tPRS_w) <- gsub("-", "_", colnames(twas_ad_tPRS_w))
colnames(twas_wd_tPRS_w) <- gsub("-", "_", colnames(twas_wd_tPRS_w))


tprs_corr <- cor(unname(white_all_tprs[,n_ct+19- 1 + c(3,4,2,5,9,10,12,21,
                                                       14, 22, 15, 17,
                                                       16, 13,8, 1,
                                                       20, 11, 
                                                       18,19,7,6)]),
                 unname(twas_wd_tPRS_w[, 1 + c(3,4,2,5,9,10,12,21,
                                               14, 22, 15, 17,
                                               16, 13,8, 1,
                                               20, 11, 
                                               18,19,7,6)]), 
                 use = "pairwise.complete.obs")
testRes <- corr.test(white_all_tprs[ ,n_ct+19- 1 + c(3,4,2,5,9,10,12,21,
                                                     14, 22, 15, 17,
                                                     16, 13,8, 1,
                                                     20, 11, 
                                                     18,19,7,6)],
                     (twas_wd_tPRS_w[, 1 + c(3,4,2,5,9,10,12,21,
                                                   14, 22, 15, 17,
                                                   16, 13,8, 1,
                                                   20, 11, 
                                                   18,19,7,6)]))

corrplot(tprs_corr, p.mat = testRes$p, sig.level = c(0.01, 0.05, 0.1), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'green', is.cor = FALSE, method = "square")
tprs_corr <- cor(unname(white_all_tprs[which(white_all_tprs$sex == "F"),n_ct+19- 1 + c(3,4,2,5,9,10,12,21,
                                                                                       14, 22, 15, 17,
                                                                                       16, 13,8, 1,
                                                                                       20, 11, 
                                                                                       18,19,7,6)]),
                 unname(twas_wd_tPRS_w[which(white_all_tprs$sex == "F"), 45 + c(3,4,2,5,9,10,12,21,
                                                                          14, 22, 15, 17,
                                                                          16, 13,8, 1,
                                                                          20, 11, 
                                                                          18,19,7,6)]), 
                 use = "pairwise.complete.obs")
testRes <- corr.test(white_all_tprs[which(white_all_tprs$sex == "F"),n_ct+19- 1 + c(3,4,2,5,9,10,12,21,
                                                                                    14, 22, 15, 17,
                                                                                    16, 13,8, 1,
                                                                                    20, 11, 
                                                                                    18,19,7,6)],
                     (twas_wd_tPRS_w[which(white_all_tprs$sex == "F"),45 + c(3,4,2,5,9,10,12,21,
                                                                             14, 22, 15, 17,
                                                                             16, 13,8, 1,
                                                                             20, 11, 
                                                                             18,19,7,6)]))

corrplot(tprs_corr, p.mat = testRes$p, sig.level = c(0.01, 0.05, 0.1), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'green', is.cor = FALSE, method = "square")
tprs_corr <- cor(unname(white_all_tprs[which(white_all_tprs$sex == "M"),n_ct+19- 1 + c(3,4,2,5,9,10,12,21,
                                                                                       14, 22, 15, 17,
                                                                                       16, 13,8, 1,
                                                                                       20, 11, 
                                                                                       18,19,7,6)]),
                 unname(twas_wd_tPRS_w[which(white_all_tprs$sex == "M"),23 + c(3,4,2,5,9,10,12,21,
                                                                         14, 22, 15, 17,
                                                                         16, 13,8, 1,
                                                                         20, 11, 
                                                                         18,19,7,6)]), 
                 use = "pairwise.complete.obs")
testRes <- corr.test(white_all_tprs[which(white_all_tprs$sex == "M"),n_ct+19- 1 + c(3,4,2,5,9,10,12,21,
                                                                                    14, 22, 15, 17,
                                                                                    16, 13,8, 1,
                                                                                    20, 11, 
                                                                                    18,19,7,6)],
                     (twas_wd_tPRS_w[which(white_all_tprs$sex == "M"),23 + c(3,4,2,5,9,10,12,21,
                                                                             14, 22, 15, 17,
                                                                             16, 13,8, 1,
                                                                             20, 11, 
                                                                             18,19,7,6)]))

corrplot(tprs_corr, p.mat = testRes$p, sig.level = c(0.01, 0.05, 0.1), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'green', is.cor = FALSE, method = "square")
# MvF
tprs_corr <- cor(unname(twas_ad_tPRS_w[which(white_all_tprs$sex == "F"),23 + c(3,4,2,5,9,10,12,21,
                                                                                       14, 22, 15, 17,
                                                                                       16, 13,8, 1,
                                                                                       20, 11, 
                                                                                       18,19,7,6)]),
                 unname(twas_ad_tPRS_w[which(white_all_tprs$sex == "F"),45 + c(3,4,2,5,9,10,12,21,
                                                                               14, 22, 15, 17,
                                                                               16, 13,8, 1,
                                                                               20, 11, 
                                                                               18,19,7,6)]), 
                 use = "pairwise.complete.obs")
testRes <- corr.test(twas_ad_tPRS_w[which(white_all_tprs$sex == "F"),23 + c(3,4,2,5,9,10,12,21,
                                                                                    14, 22, 15, 17,
                                                                                    16, 13,8, 1,
                                                                                    20, 11, 
                                                                                    18,19,7,6)],
                     (twas_ad_tPRS_w[which(white_all_tprs$sex == "F"),45 + c(3,4,2,5,9,10,12,21,
                                                                             14, 22, 15, 17,
                                                                             16, 13,8, 1,
                                                                             20, 11, 
                                                                             18,19,7,6)]))

corrplot(tprs_corr, p.mat = testRes$p, sig.level = c(0.01, 0.05, 0.1), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'green', is.cor = FALSE, method = "color", 
         cl.lim=c(-1,1), col=colorRampPalette(c("red", "white","blue"))(200))

tprs_corr <- cor(unname(twas_wd_tPRS_w[which(white_all_tprs$sex == "M"),23 + c(3,4,2,5,9,10,12,21,
                                                                               14, 22, 15, 17,
                                                                               16, 13,8, 1,
                                                                               20, 11, 
                                                                               18,19,7,6)]),
                 unname(twas_wd_tPRS_w[which(white_all_tprs$sex == "M"),45 + c(3,4,2,5,9,10,12,21,
                                                                               14, 22, 15, 17,
                                                                               16, 13,8, 1,
                                                                               20, 11, 
                                                                               18,19,7,6)]), 
                 use = "pairwise.complete.obs")
testRes <- corr.test(twas_wd_tPRS_w[which(white_all_tprs$sex == "M"),23 + c(3,4,2,5,9,10,12,21,
                                                                            14, 22, 15, 17,
                                                                            16, 13,8, 1,
                                                                            20, 11, 
                                                                            18,19,7,6)],
                     (twas_wd_tPRS_w[which(white_all_tprs$sex == "M"),45 + c(3,4,2,5,9,10,12,21,
                                                                             14, 22, 15, 17,
                                                                             16, 13,8, 1,
                                                                             20, 11, 
                                                                             18,19,7,6)]))

corrplot(tprs_corr, p.mat = testRes$p, sig.level = c(0.01, 0.05, 0.1), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'green', is.cor = FALSE, method = "color", 
         cl.lim=c(-1, 0), col=colorRampPalette(c("red","white"))(200))
corrplot(tprs_corr, p.mat = testRes$p, sig.level = c(0.01, 0.05, 0.1), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'green', is.cor = FALSE, method = "color", 
         cl.lim=c(0,1), col=colorRampPalette(c("white","blue"))(200))

#corr mat between dep scores and tPRSs================== 

tprs_corr <- cbind(cbind(cbind(cbind(cbind(cor((twas_ad_tPRS_w[,1 + c(3,4,2,5,9,10,12,21,
                                                     14, 22, 15, 17,
                                                     16, 13,8, 1,
                                                     20, 11, 
                                                     18,19,7,6)]),
                             (white_all_tprs[,c("anxdep_t")]), use = "pairwise.complete.obs"), 
                             cor((twas_wd_tPRS_w[,1 + c(3,4,2,5,9,10,12,21,
                                                        14, 22, 15, 17,
                                                        16, 13,8, 1,
                                                        20, 11, 
                                                        18,19,7,6)]),
                                 (white_all_tprs[,c("withdep_t")]), use = "pairwise.complete.obs")),
                         cor((twas_ad_tPRS_w[which(white_all_tprs$sex == "M"),23 + c(3,4,2,5,9,10,12,21,
                                                                                           14, 22, 15, 17,
                                                                                           16, 13,8, 1,
                                                                                           20, 11, 
                                                                                           18,19,7,6)]), 
                             (white_all_tprs[which(white_all_tprs$sex == "M"),c("anxdep_t")]), use = "pairwise.complete.obs")),
                         cor((twas_wd_tPRS_w[which(white_all_tprs$sex == "M"),23 + c(3,4,2,5,9,10,12,21,
                                                                                     14, 22, 15, 17,
                                                                                     16, 13,8, 1,
                                                                                     20, 11, 
                                                                                     18,19,7,6)]), 
                             (white_all_tprs[which(white_all_tprs$sex == "M"),c("withdep_t")]), use = "pairwise.complete.obs")),
                   cor((twas_ad_tPRS_w[which(white_all_tprs$sex == "F"),45 + c(3,4,2,5,9,10,12,21,
                                                                                     14, 22, 15, 17,
                                                                                     16, 13,8, 1,
                                                                                     20, 11, 
                                                                                     18,19,7,6)]), 
                       (white_all_tprs[which(white_all_tprs$sex == "F"),c("anxdep_t")]), use = "pairwise.complete.obs")), 
                   cor((twas_wd_tPRS_w[which(white_all_tprs$sex == "F"),45 + c(3,4,2,5,9,10,12,21,
                                                                               14, 22, 15, 17,
                                                                               16, 13,8, 1,
                                                                               20, 11, 
                                                                               18,19,7,6)]), 
                       (white_all_tprs[which(white_all_tprs$sex == "F"),c("withdep_t")]), use = "pairwise.complete.obs"))
colnames(tprs_corr) <- c("anx/dep", "with/dep", "anx/dep M", "with/dep M",
                        "anx/dep F", "with/dep F")
testRes <- cbind(cbind(cbind(cbind(cbind(corr.test((twas_ad_tPRS_w[,1 + c(3,4,2,5,9,10,12,21,
                                                                    14, 22, 15, 17,
                                                                    16, 13,8, 1,
                                                                    20, 11, 
                                                                    18,19,7,6)]),
                                             (white_all_tprs[,c("anxdep_t")]))$p, 
                                         corr.test((twas_wd_tPRS_w[,1 + c(3,4,2,5,9,10,12,21,
                                                                    14, 22, 15, 17,
                                                                    16, 13,8, 1,
                                                                    20, 11, 
                                                                    18,19,7,6)]),
                                             (white_all_tprs[,c("withdep_t")]))$p),
                                   corr.test((twas_ad_tPRS_w[which(white_all_tprs$sex == "M"),23 + c(3,4,2,5,9,10,12,21,
                                                                                               14, 22, 15, 17,
                                                                                               16, 13,8, 1,
                                                                                               20, 11, 
                                                                                               18,19,7,6)]), 
                                       (white_all_tprs[which(white_all_tprs$sex == "M"),c("anxdep_t")]))$p),
                             corr.test((twas_wd_tPRS_w[which(white_all_tprs$sex == "M"),23 + c(3,4,2,5,9,10,12,21,
                                                                                         14, 22, 15, 17,
                                                                                         16, 13,8, 1,
                                                                                         20, 11, 
                                                                                         18,19,7,6)]), 
                                 (white_all_tprs[which(white_all_tprs$sex == "M"),c("withdep_t")]))$p),
                       corr.test((twas_ad_tPRS_w[which(white_all_tprs$sex == "F"),45 + c(3,4,2,5,9,10,12,21,
                                                                                   14, 22, 15, 17,
                                                                                   16, 13,8, 1,
                                                                                   20, 11, 
                                                                                   18,19,7,6)]), 
                           (white_all_tprs[which(white_all_tprs$sex == "F"),c("anxdep_t")]))$p), 
                 corr.test((twas_wd_tPRS_w[which(white_all_tprs$sex == "F"),45 + c(3,4,2,5,9,10,12,21,
                                                                             14, 22, 15, 17,
                                                                             16, 13,8, 1,
                                                                             20, 11, 
                                                                             18,19,7,6)]), 
                     (white_all_tprs[which(white_all_tprs$sex == "F"),c("withdep_t")]))$p)
colnames(testRes) <- c("anx/dep", "with/dep", "anx/dep M", "with/dep M",
                        "anx/dep F", "with/dep F")  
  
corrplot(tprs_corr, p.mat = testRes, method = "square", is.cor = F, 
         sig.level = c(0.01, 0.05, 0.1), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'green', cl.ratio = 0.5)
# test data
strat_ad <- multi_strata(white_all_tprs[,c("sex", "anxdep_t", "anx_bin")])
strat_wd <- multi_strata(white_all_tprs[,c("sex", "withdep_t", "with_bin")])
inds_ad_all <- partition(strat_ad, p = c(train = 0.7, test = 0.3), seed = 330, split_into_list = FALSE)
inds_wd_all <- partition(strat_wd, p = c(train = 0.7, test = 0.3), seed = 331, split_into_list = FALSE)

strat_ad <- multi_strata(white_all_tprs[which(white_all_tprs$sex == "M"), c("anxdep_t", "anx_bin")])
strat_wd <- multi_strata(white_all_tprs[which(white_all_tprs$sex == "M"), c("withdep_t", "with_bin")])
inds_ad_m <- partition(strat_ad, p = c(train = 0.7, test = 0.3), seed = 330, split_into_list = FALSE)
inds_wd_m <- partition(strat_wd, p = c(train = 0.7, test = 0.3), seed = 331, split_into_list = FALSE)


strat_ad <- multi_strata(white_all_tprs[which(white_all_tprs$sex == "F"), c("anxdep_t", "anx_bin")])
strat_wd <- multi_strata(white_all_tprs[which(white_all_tprs$sex == "F"), c("withdep_t", "with_bin")])
inds_ad_f <- partition(strat_ad, p = c(train = 0.7, test = 0.3), seed = 330, split_into_list = FALSE)
inds_wd_f <- partition(strat_wd, p = c(train = 0.7, test = 0.3), seed = 331, split_into_list = FALSE)
tprs_corr <- cbind(cbind(cbind(cbind(cbind(cor((twas_ad_tPRS_w[which(inds_ad_all=="train"),1 + c(3,4,2,5,9,10,12,21,
                                                                      14, 22, 15, 17,
                                                                      16, 13,8, 1,
                                                                      20, 11, 
                                                                      18,19,7,6)]),
                                               (white_all_tprs[which(inds_ad_all=="train"),c("anxdep_t")]), use = "pairwise.complete.obs"), 
                                           cor((twas_wd_tPRS_w[which(inds_wd_all=="train"),1 + c(3,4,2,5,9,10,12,21,
                                                                      14, 22, 15, 17,
                                                                      16, 13,8, 1,
                                                                      20, 11, 
                                                                      18,19,7,6)]),
                                               (white_all_tprs[which(inds_wd_all=="train"),c("withdep_t")]), use = "pairwise.complete.obs")),
                                     cor((twas_ad_tPRS_w[which(white_all_tprs$sex == "M"),23 + c(3,4,2,5,9,10,12,21,
                                                                                                 14, 22, 15, 17,
                                                                                                 16, 13,8, 1,
                                                                                                 20, 11, 
                                                                                                 18,19,7,6)][which(inds_ad_m=="train"),]), 
                                         (white_all_tprs[which(white_all_tprs$sex == "M"),c("anxdep_t")][which(inds_ad_m=="train")]), use = "pairwise.complete.obs")),
                               cor((twas_wd_tPRS_w[which(white_all_tprs$sex == "M"),23 + c(3,4,2,5,9,10,12,21,
                                                                                           14, 22, 15, 17,
                                                                                           16, 13,8, 1,
                                                                                           20, 11, 
                                                                                           18,19,7,6)][which(inds_wd_m=="train"),]), 
                                   (white_all_tprs[which(white_all_tprs$sex == "M"),c("withdep_t")][which(inds_wd_m=="train")]), use = "pairwise.complete.obs")),
                         cor((twas_ad_tPRS_w[which(white_all_tprs$sex == "F"),45 + c(3,4,2,5,9,10,12,21,
                                                                                     14, 22, 15, 17,
                                                                                     16, 13,8, 1,
                                                                                     20, 11, 
                                                                                     18,19,7,6)][which(inds_ad_f=="train"),]), 
                             (white_all_tprs[which(white_all_tprs$sex == "F"),c("anxdep_t")][which(inds_ad_f=="train")]), use = "pairwise.complete.obs")), 
                   cor((twas_wd_tPRS_w[which(white_all_tprs$sex == "F"),45 + c(3,4,2,5,9,10,12,21,
                                                                               14, 22, 15, 17,
                                                                               16, 13,8, 1,
                                                                               20, 11, 
                                                                               18,19,7,6)][which(inds_wd_f=="train"),]), 
                       (white_all_tprs[which(white_all_tprs$sex == "F"),c("withdep_t")][which(inds_wd_f=="train")]), use = "pairwise.complete.obs"))
colnames(tprs_corr) <- c("anx/dep", "with/dep", "anx/dep M", "with/dep M",
                         "anx/dep F", "with/dep F")
testRes <- cbind(cbind(cbind(cbind(cbind(corr.test((twas_ad_tPRS_w[which(inds_ad_all=="train"),1 + c(3,4,2,5,9,10,12,21,
                                                                          14, 22, 15, 17,
                                                                          16, 13,8, 1,
                                                                          20, 11, 
                                                                          18,19,7,6)]),
                                                   (white_all_tprs[which(inds_ad_all=="train"),c("anxdep_t")]))$p, 
                                         corr.test((twas_wd_tPRS_w[which(inds_wd_all=="train"),1 + c(3,4,2,5,9,10,12,21,
                                                                          14, 22, 15, 17,
                                                                          16, 13,8, 1,
                                                                          20, 11, 
                                                                          18,19,7,6)]),
                                                   (white_all_tprs[which(inds_wd_all=="train"),c("withdep_t")]))$p),
                                   corr.test((twas_ad_tPRS_w[which(white_all_tprs$sex == "M"),23 + c(3,4,2,5,9,10,12,21,
                                                                                                     14, 22, 15, 17,
                                                                                                     16, 13,8, 1,
                                                                                                     20, 11, 
                                                                                                     18,19,7,6)][which(inds_ad_m=="train"),]), 
                                             (white_all_tprs[which(white_all_tprs$sex == "M"),c("anxdep_t")][which(inds_ad_m=="train")]))$p),
                             corr.test((twas_wd_tPRS_w[which(white_all_tprs$sex == "M"),23 + c(3,4,2,5,9,10,12,21,
                                                                                               14, 22, 15, 17,
                                                                                               16, 13,8, 1,
                                                                                               20, 11, 
                                                                                               18,19,7,6)][which(inds_wd_m=="train"),]), 
                                       (white_all_tprs[which(white_all_tprs$sex == "M"),c("withdep_t")][which(inds_wd_m=="train")]))$p),
                       corr.test((twas_ad_tPRS_w[which(white_all_tprs$sex == "F"),45 + c(3,4,2,5,9,10,12,21,
                                                                                         14, 22, 15, 17,
                                                                                         16, 13,8, 1,
                                                                                         20, 11, 
                                                                                         18,19,7,6)][which(inds_ad_f=="train"),]), 
                                 (white_all_tprs[which(white_all_tprs$sex == "F"),c("anxdep_t")][which(inds_ad_f=="train")]))$p), 
                 corr.test((twas_wd_tPRS_w[which(white_all_tprs$sex == "F"),45 + c(3,4,2,5,9,10,12,21,
                                                                                   14, 22, 15, 17,
                                                                                   16, 13,8, 1,
                                                                                   20, 11, 
                                                                                   18,19,7,6)][which(inds_wd_f=="train"),]), 
                           (white_all_tprs[which(white_all_tprs$sex == "F"),c("withdep_t")][which(inds_wd_f=="train")]))$p)
colnames(testRes) <- c("anx/dep", "with/dep", "anx/dep M", "with/dep M",
                       "anx/dep F", "with/dep F")  

corrplot(tprs_corr, p.mat = testRes, method = "square", is.cor = F, 
         sig.level = c(0.01, 0.05, 0.1), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'green', cl.ratio = 0.5)





