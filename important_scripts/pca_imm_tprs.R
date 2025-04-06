# required: twas_analysis_abcd.R
# goal: PCA on TWAS based tPRSs
# next: TODO investigate meanings of PCs

semantic_order <- c(3,4,2,5,9,10,12,21,
                    14, 22, 15, 17,
                    16, 13,8, 1,
                    20, 11, 
                    18,19,7,6)
semantic_order0 <- semantic_order - 1
nct <- 22
load("/external/rprshnas01/kcni/mding/sc-tprs-mdd/twas_ad_white_all_tprs.Rda")
load("/external/rprshnas01/kcni/mding/sc-tprs-mdd/twas_wd_white_all_tprs.Rda")

scale_twas_ad_tPRS <- scale(twas_ad_tPRS_w[,-1])
scale_twas_wd_tPRS <- scale(twas_wd_tPRS_w[,-1])

pca_twas_tprs_ad <- prcomp((scale_twas_ad_tPRS[,1:22]))
pca_twas_tprs_wd <- prcomp((scale_twas_wd_tPRS[,1:22]))
pca_twas_tprs_ad_m <- prcomp((scale_twas_ad_tPRS[which(white_all_tprs$sex=="M"),23:44]))
pca_twas_tprs_wd_m <- prcomp((scale_twas_wd_tPRS[which(white_all_tprs$sex=="M"),23:44]))
pca_twas_tprs_ad_f <- prcomp((scale_twas_ad_tPRS[which(white_all_tprs$sex=="F"),45:66]))
pca_twas_tprs_wd_f <- prcomp((scale_twas_wd_tPRS[which(white_all_tprs$sex=="F"),45:66]))

pca_twas_tprs_ad_t <- prcomp(t(scale_twas_ad_tPRS[,1:22]))
pca_twas_tprs_wd_t <- prcomp(t(scale_twas_wd_tPRS[,1:22]))
pca_twas_tprs_ad_m_t <- prcomp(t(scale_twas_ad_tPRS[which(white_all_tprs$sex=="M"),23:44]))
pca_twas_tprs_wd_m_t <- prcomp(t(scale_twas_wd_tPRS[which(white_all_tprs$sex=="M"),23:44]))
pca_twas_tprs_ad_f_t <- prcomp(t(scale_twas_ad_tPRS[which(white_all_tprs$sex=="F"),45:66]))
pca_twas_tprs_wd_f_t <- prcomp(t(scale_twas_wd_tPRS[which(white_all_tprs$sex=="F"),45:66]))

plot(pca_twas_tprs_ad)
summary(pca_twas_tprs_ad) #8 pcs to get 80% of var
summary(pca_twas_tprs_wd) #8 pcs to get 80% of var
summary(pca_twas_tprs_ad_m) #2 pcs to get 81% of var
summary(pca_twas_tprs_wd_m) #4 pcs to get 81% of var
summary(pca_twas_tprs_ad_f) #4 pcs to get 80% of var
summary(pca_twas_tprs_wd_f) #4 pcs to get 81% of var

ggplot(stack(as.data.frame(pca_twas_tprs_ad$rotation)[, c(1:6)], varying = colnames(pca_twas_tprs_ad$rotation)[1:6]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack(as.data.frame(pca_twas_tprs_ad$rotation)[, c(7:12)], varying = colnames(pca_twas_tprs_ad$rotation)[7:12]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack(as.data.frame(pca_twas_tprs_ad$rotation)[, c(13:18)], varying = colnames(pca_twas_tprs_ad$rotation)[13:18]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack(as.data.frame(pca_twas_tprs_ad$rotation)[, c(19:22)], varying = colnames(pca_twas_tprs_ad$rotation)[19:22]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
# Plot projections along the components into a scatterplot.
# Axes for points are scaled as values, for vectors as variance
# Default for biplot() is the first and second component.

biplot(pca_twas_tprs_ad)
biplot(pca_twas_tprs_wd)
biplot(pca_twas_tprs_ad_m)
biplot(pca_twas_tprs_ad_f)
biplot(pca_twas_tprs_wd_m)
biplot(pca_twas_tprs_wd_f)
# Examine the actual principal components in a parallel-coordinates

N <- 4
matplot(pca_twas_tprs_ad$rotation[, 1:N],
        type="b", lwd=1,
        xlab = "ABCD participants", ylab="PCs")

matplot(pca_twas_tprs_wd$rotation[, 1:N],
        type="b", lwd=1,
        xlab = "ABCD participants", ylab="PCs")

# A function to find the best correlation:

bestCor <- function(V, dat = scale_tprs[
  which(scale_tprs$event_name=="baseline_year_1_arm_1"),], tprs_start = tprs_start_idx,
  tprs_end = tprs_start_idx + n_ct - 1) {
  # Identify that column in dataframe dat that has the highest correlation
  # or anticorrelation with vector V. Print some information about it.
  # Value: The index of the column (invisibly)
  #
  # ToDo: identify and return n-best correlations
  #
  cors <- numeric(tprs_end - tprs_start + 1)
  for (i in tprs_start:tprs_end) {
    cors[i-tprs_start+1] <- cor(V, dat[ , i]) # correlation between input vector and each
  }                              # column of the feature set in turn
  
  myBest <- max(abs(cors), na.rm = TRUE)
  sel <- abs(cors) == myBest   # Is TRUE for the highest absolute value in the
  # vector of correlation coefficients.
  idx <- which(sel)[1]       # which() finds the index of the TRUE value(s).
  # Pick only the first in case there are ties.
  print(sprintf("\n Highest correlation for index %d (%f): %s ",
                idx,
                cors[idx],
                colnames(dat)[idx+tprs_start-1]))
  return(idx)     # return the index.
}



bestCor(pca_twas_tprs_ad_t$rotation[ , 1], dat = scale_twas_ad_tPRS[,1:22], tprs_start = 1, tprs_end = 22)  # This is PC1 ...
bestCor(pca_twas_tprs_ad_t$rotation[ , 2], dat = scale_twas_ad_tPRS[,1:22], tprs_start = 1, tprs_end = 22)  # PC2
bestCor(pca_twas_tprs_ad_t$rotation[ , 3], dat = scale_twas_ad_tPRS[,1:22], tprs_start = 1, tprs_end = 22)  # etc.
bestCor(pca_twas_tprs_ad_t$rotation[ , 4], dat = scale_twas_ad_tPRS[,1:22], tprs_start = 1, tprs_end = 22)
bestCor(pca_twas_tprs_ad_t$rotation[ , 5], dat = scale_twas_ad_tPRS[,1:22], tprs_start = 1, tprs_end = 22) 
bestCor(pca_twas_tprs_ad_t$rotation[ , 6], dat = scale_twas_ad_tPRS[,1:22], tprs_start = 1, tprs_end = 22)  
bestCor(pca_twas_tprs_ad_t$rotation[ , 7], dat = scale_twas_ad_tPRS[,1:22], tprs_start = 1, tprs_end = 22) 
bestCor(pca_twas_tprs_ad_t$rotation[ , 8], dat = scale_twas_ad_tPRS[,1:22], tprs_start = 1, tprs_end = 22)  

# gamma_delta_T_cell, 0.357156
# transitional_stage_B_cell, 0.380850
# dendritic_cell, -0.345042
# CD4_positive_alpha_beta_T_cell, -0.394845
# platelet, -0.377758
# conventional_dendritic_cell, -0.331595
# memory_B_cell, 0.262334
# plasmacytoid_dendritic_cell, -0.265770

bestCor(pca_twas_tprs_wd_t$rotation[ , 1], dat = scale_twas_wd_tPRS[,1:22], tprs_start = 1, tprs_end = 22)  # This is PC1 ...
bestCor(pca_twas_tprs_wd_t$rotation[ , 2], dat = scale_twas_wd_tPRS[,1:22], tprs_start = 1, tprs_end = 22)  # PC2
bestCor(pca_twas_tprs_wd_t$rotation[ , 3], dat = scale_twas_wd_tPRS[,1:22], tprs_start = 1, tprs_end = 22)  # etc.
bestCor(pca_twas_tprs_wd_t$rotation[ , 4], dat = scale_twas_wd_tPRS[,1:22], tprs_start = 1, tprs_end = 22)
bestCor(pca_twas_tprs_wd_t$rotation[ , 5], dat = scale_twas_wd_tPRS[,1:22], tprs_start = 1, tprs_end = 22) 
bestCor(pca_twas_tprs_wd_t$rotation[ , 6], dat = scale_twas_wd_tPRS[,1:22], tprs_start = 1, tprs_end = 22) 
bestCor(pca_twas_tprs_wd_t$rotation[ , 7], dat = scale_twas_wd_tPRS[,1:22], tprs_start = 1, tprs_end = 22) 
bestCor(pca_twas_tprs_wd_t$rotation[ , 8], dat = scale_twas_wd_tPRS[,1:22], tprs_start = 1, tprs_end = 22)  
# gamma_delta_T_cell, 0.339292
# transitional_stage_B_cell, 0.375717
# platelet, -0.389439
# CD8_positive_alpha_beta_T_cell, 0.351754
# dendritic_cell, -0.349379
# dendritic_cell, -0.273617
# innate_lymphoid_cell, -0.283885
# platelet, 0.270649



bestCor(pca_twas_tprs_ad_m_t$rotation[ , 1], dat = scale_twas_ad_tPRS[which(white_all_tprs$sex=="M"),22+1:22], tprs_start = 1, tprs_end = 22)  # This is PC1 ...
bestCor(pca_twas_tprs_ad_m_t$rotation[ , 2], dat = scale_twas_ad_tPRS[which(white_all_tprs$sex=="M"),22+1:22], tprs_start = 1, tprs_end = 22)  # PC2
bestCor(pca_twas_tprs_ad_m_t$rotation[ , 3], dat = scale_twas_ad_tPRS[which(white_all_tprs$sex=="M"),22+1:22], tprs_start = 1, tprs_end = 22)  # etc.
bestCor(pca_twas_tprs_ad_m_t$rotation[ , 4], dat = scale_twas_ad_tPRS[which(white_all_tprs$sex=="M"),22+1:22], tprs_start = 1, tprs_end = 22)
bestCor(pca_twas_tprs_ad_m_t$rotation[ , 5], dat = scale_twas_ad_tPRS[which(white_all_tprs$sex=="M"),22+1:22], tprs_start = 1, tprs_end = 22) 
bestCor(pca_twas_tprs_ad_m_t$rotation[ , 6], dat = scale_twas_ad_tPRS[which(white_all_tprs$sex=="M"),22+1:22], tprs_start = 1, tprs_end = 22)  
# CD4_positive_alpha_beta_cytotoxic_T_cell, 0.268299
# transitional_stage_B_cell, 0.294500
# CD14_low_CD16_positive_monocyte, -0.289136
# CD8_positive_alpha_beta_T_cell, 0.215968
# conventional_dendritic_cell, -0.307273
# gamma_delta_T_cell, -0.251454


bestCor(pca_twas_tprs_ad_f_t$rotation[ , 1], dat = scale_twas_ad_tPRS[which(white_all_tprs$sex=="F"),44+1:22], tprs_start = 1, tprs_end = 22)  # This is PC1 ...
bestCor(pca_twas_tprs_ad_f_t$rotation[ , 2], dat = scale_twas_ad_tPRS[which(white_all_tprs$sex=="F"),44+1:22], tprs_start = 1, tprs_end = 22)  # PC2
bestCor(pca_twas_tprs_ad_f_t$rotation[ , 3], dat = scale_twas_ad_tPRS[which(white_all_tprs$sex=="F"),44+1:22], tprs_start = 1, tprs_end = 22)  # etc.
bestCor(pca_twas_tprs_ad_f_t$rotation[ , 4], dat = scale_twas_ad_tPRS[which(white_all_tprs$sex=="F"),44+1:22], tprs_start = 1, tprs_end = 22)
bestCor(pca_twas_tprs_ad_f_t$rotation[ , 5], dat = scale_twas_ad_tPRS[which(white_all_tprs$sex=="F"),44+1:22], tprs_start = 1, tprs_end = 22) 
bestCor(pca_twas_tprs_ad_f_t$rotation[ , 6], dat = scale_twas_ad_tPRS[which(white_all_tprs$sex=="F"),44+1:22], tprs_start = 1, tprs_end = 22)  
# CD4_positive_alpha_beta_cytotoxic_T_cell, 0.250126
# transitional_stage_B_cell, -0.264319
# CD4_positive_alpha_beta_T_cell, -0.276261
# CD8_positive_alpha_beta_T_cell, -0.253551
# plasmablast, -0.333986
# double_negative_thymocyte, 0.291798


bestCor(pca_twas_tprs_wd_m_t$rotation[ , 1], dat = scale_twas_wd_tPRS[which(white_all_tprs$sex=="M"),22+1:22], tprs_start = 1, tprs_end = 22)  # This is PC1 ...
bestCor(pca_twas_tprs_wd_m_t$rotation[ , 2], dat = scale_twas_wd_tPRS[which(white_all_tprs$sex=="M"),22+1:22], tprs_start = 1, tprs_end = 22)  # PC2
bestCor(pca_twas_tprs_wd_m_t$rotation[ , 3], dat = scale_twas_wd_tPRS[which(white_all_tprs$sex=="M"),22+1:22], tprs_start = 1, tprs_end = 22)  # etc.
bestCor(pca_twas_tprs_wd_m_t$rotation[ , 4], dat = scale_twas_wd_tPRS[which(white_all_tprs$sex=="M"),22+1:22], tprs_start = 1, tprs_end = 22)
bestCor(pca_twas_tprs_wd_m_t$rotation[ , 5], dat = scale_twas_wd_tPRS[which(white_all_tprs$sex=="M"),22+1:22], tprs_start = 1, tprs_end = 22) 
bestCor(pca_twas_tprs_wd_m_t$rotation[ , 6], dat = scale_twas_wd_tPRS[which(white_all_tprs$sex=="M"),22+1:22], tprs_start = 1, tprs_end = 22)  
# dendritic_cell, -0.286939
# transitional_stage_B_cell, 0.319418
# CD4_positive_alpha_beta_T_cell, -0.325063
# erythrocyte, 0.329100
# innate_lymphoid_cell, -0.261205
# peripheral_blood_mononuclear_cell, 0.363856


bestCor(pca_twas_tprs_wd_f_t$rotation[ , 1], dat = scale_twas_wd_tPRS[which(white_all_tprs$sex=="F"),44+1:22], tprs_start = 1, tprs_end = 22)  # This is PC1 ...
bestCor(pca_twas_tprs_wd_f_t$rotation[ , 2], dat = scale_twas_wd_tPRS[which(white_all_tprs$sex=="F"),44+1:22], tprs_start = 1, tprs_end = 22)  # PC2
bestCor(pca_twas_tprs_wd_f_t$rotation[ , 3], dat = scale_twas_wd_tPRS[which(white_all_tprs$sex=="F"),44+1:22], tprs_start = 1, tprs_end = 22)  # etc.
bestCor(pca_twas_tprs_wd_f_t$rotation[ , 4], dat = scale_twas_wd_tPRS[which(white_all_tprs$sex=="F"),44+1:22], tprs_start = 1, tprs_end = 22)
bestCor(pca_twas_tprs_wd_f_t$rotation[ , 5], dat = scale_twas_wd_tPRS[which(white_all_tprs$sex=="F"),44+1:22], tprs_start = 1, tprs_end = 22) 
bestCor(pca_twas_tprs_wd_f_t$rotation[ , 6], dat = scale_twas_wd_tPRS[which(white_all_tprs$sex=="F"),44+1:22], tprs_start = 1, tprs_end = 22)  
# central_memory_CD8_positive_alpha_beta_T_cell, 0.282942
# naive_B_cell, -0.326724
# memory_B_cell, -0.269131
# platelet, 0.294073
# dendritic_cell, -0.257741
# plasmablast, -0.280891

ggplot(as.data.frame(pca_twas_tprs_ad$rotation), aes(x = PC1, y = PC2)) +
  geom_point(aes(color = white_all_tprs$anxdep_t), alpha = 0.3)+
  scale_color_viridis_c(name = "anx/dep")
ggplot(as.data.frame(pca_twas_tprs_wd$rotation), aes(x = PC1, y = PC2)) +
  geom_point(aes(color = white_all_tprs$withdep_t), alpha = 0.3)+
  scale_color_viridis_c(name = "with/dep")   
ggplot(as.data.frame(pca_twas_tprs_ad_m$rotation), aes(x = PC1, y = PC2)) +
  geom_point(aes(color = white_all_tprs$anxdep_t[which(white_all_tprs$sex == "M")]), alpha = 0.3)+
  scale_color_viridis_c(name = "anx/dep")  
ggplot(as.data.frame(pca_twas_tprs_ad_f$rotation), aes(x = PC1, y = PC2)) +
  geom_point(aes(color = white_all_tprs$withdep_t[which(white_all_tprs$sex == "F")]), alpha = 0.3)+
  scale_color_viridis_c(name = "anx/dep")  

ggplot(as.data.frame(pca_twas_tprs_ad$rotation), aes(x = PC1, y = PC2)) +
  geom_point(aes(color = pc_df$PC1_w), alpha = 0.3)+
  scale_color_viridis_c(name = "PC1 Wittenberg tPRS")

ggplot(as.data.frame(pca_twas_tprs_wd$rotation), aes(x = PC1, y = PC2)) +
  geom_point(aes(color = pc_df$PC1_w), alpha = 0.3)+
  scale_color_viridis_c(name = "PC1 Wittenberg tPRS")   
ggplot(as.data.frame(pca_twas_tprs_ad_m$rotation), aes(x = PC1, y = PC2)) +
  geom_point(aes(color = pc_df$PC1_w[which(white_all_tprs$sex == "M")]), alpha = 0.3)+
  scale_color_viridis_c(name = "PC1 Wittenberg tPRS")  
ggplot(as.data.frame(pca_twas_tprs_ad_f$rotation), aes(x = PC1, y = PC2)) +
  geom_point(aes(color = pc_df$PC1_w[which(white_all_tprs$sex == "F")]), alpha = 0.3)+
  scale_color_viridis_c(name = "PC1 Wittenberg tPRS")



# t-test difference in means between males and females in PCs

diffPC <- function(pcs, dat = white_all_tprs$sex) {
  # identify which PC, if any shows difference in means between 2 groups (usually sex)
  for (i in 1:ncol(pcs)) {
    tt <- t.test(pcs[which(dat == unique(dat)[1]), i],
                 pcs[which(dat == unique(dat)[2]), i])
    if (tt$p.value < 0.05){
      print(sprintf("PC where difference in means observed %d (%f): %f \n",
                    i,
                    tt$p.value, tt$estimate[1] - tt$estimate[2]))
    }
  }                             
}
diffPC(pca_twas_tprs_ad_t$rotation) #PC5 (0.044019): -0.000785
# PC14 (0.028844): -0.000852
# PC17 (0.047160): -0.000774

diffPC(pca_twas_tprs_wd_t$rotation) #PC2 (0.027250): -0.000860

ggplot(as.data.frame(pca_twas_tprs_ad_t$rotation), aes(x = PC14, y = PC15)) +
  geom_point(aes(color = white_all_tprs$sex), alpha = 0.3)+
  scale_color_manual(name = "sex", values = c("blue", "red"))

ggplot(as.data.frame(pca_twas_tprs_wd_t$rotation), aes(x = PC2, y = PC3)) +
  geom_point(aes(color = white_all_tprs$sex), alpha = 0.3)+
  scale_color_manual(name = "sex", values = c("blue", "red"))

# differences by sex in mixed tPRS not very noticeable in PCs

diffPC(pca_twas_tprs_ad_t$rotation, dat = base_bin$anx_bin)
diffPC(pca_twas_tprs_wd_t$rotation, dat = base_bin$with_bin)
ggplot(as.data.frame(pca_twas_tprs_ad_t$rotation), aes(x = PC8, y = PC7)) +
  geom_point(aes(color = base_bin$anx_bin), alpha = 0.3)+
  scale_color_manual(name = "anx/dep low or high", values = c("blue", "red"))

ggplot(as.data.frame(pca_twas_tprs_ad_t$rotation), aes(x = PC15, y = PC14)) +
  geom_point(aes(color = base_bin$anx_bin), alpha = 0.3)+
  scale_color_manual(name = "anx/dep low or high", values = c("blue", "red"))

diffPC(pca_twas_tprs_ad_m_t$rotation, dat = base_bin$anx_bin[base_bin$sex == "1"])
diffPC(pca_twas_tprs_ad_f_t$rotation, dat = base_bin$anx_bin[base_bin$sex == "2"])

ggplot(as.data.frame(pca_twas_tprs_ad_f_t$rotation), aes(x = PC1, y = PC2)) +
  geom_point(aes(color = base_bin$anx_bin[which(base_bin$sex == "2")]), alpha = 0.3)+
  scale_color_manual(name = "anx/dep low or high", values = c("blue", "red"))


diffPC(pca_twas_tprs_wd_m_t$rotation, dat = base_bin$with_bin[base_bin$sex == "1"])
diffPC(pca_twas_tprs_wd_f_t$rotation, dat = base_bin$with_bin[base_bin$sex == "2"])
ggplot(as.data.frame(pca_twas_tprs_wd_f_t$rotation), aes(x = PC12, y = PC10)) +
  geom_point(aes(color = base_bin$with_bin[which(base_bin$sex == "2")]), alpha = 0.3)+
  scale_color_manual(name = "anx/dep low or high", values = c("blue", "red"))

# biplots ===========================
library(ggfortify)
ld_col <- rep("brown2", n_ct)
ld_col[grep("T_cell|thymocyte",cell_types)] <- "darkolivegreen2"
ld_col[grep("B_cell|plasmablast",cell_types)] <- "cyan2"
ld_col[grep("dendritic_cell",cell_types)] <- "darkorchid2"

autoplot(pca_twas_tprs_ad, loadings = TRUE, loadings.colour = ld_col,
         loadings.label.colour = ld_col,
         loadings.label = TRUE, loadings.label.size = 3, x = 2, y = 3) + theme_stata() 
autoplot(pca_twas_tprs_wd, loadings = TRUE,  loadings.colour = ld_col,
         loadings.label.colour = ld_col, loadings.label = TRUE, 
         loadings.label.size = 3, x = 2, y = 3)+ theme_stata()
autoplot(pca_twas_tprs_ad_m, loadings = TRUE,  loadings.colour = ld_col,
         loadings.label.colour = ld_col, loadings.label = TRUE, 
         loadings.label.size = 3)+ theme_stata()
autoplot(pca_twas_tprs_ad_f, loadings = TRUE,  loadings.colour = ld_col,
         loadings.label.colour = ld_col, loadings.label = TRUE, 
         loadings.label.size = 3)+ theme_stata()
autoplot(pca_twas_tprs_wd_m, loadings = TRUE,  loadings.colour = ld_col,
         loadings.label.colour = ld_col, loadings.label = TRUE, 
         loadings.label.size = 3)+ theme_stata()
autoplot(pca_twas_tprs_wd_f, loadings = TRUE,  loadings.colour = ld_col,
         loadings.label.colour = ld_col, loadings.label = TRUE, 
         loadings.label.size = 3)+ theme_stata()

