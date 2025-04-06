# required: load_abcd_data.R
# goal: correlation between wittenberg and TWAS tPRSs and CBC measures
# next: TODO mediation analysis?

ids <- unique(og_white_all_tprs$IID)
imm_tprs  <- as.data.frame(og_white_all_tprs[1, -3])
for (i in 2:length(ids)){
  id <- ids[i]
  sub <- as.data.frame(og_white_all_tprs[which(white_all_tprs$IID == id),])
  for (j in 1:ncol(sub)){
    cn <- colnames(sub)[j]
    if (cn == 'event_name'){
      next
    } else if (cn =='age') {
      imm_tprs[i,'age'] <- min(sub[,'age'])
    } else if (j <= 26 & 9 <= j) {
      val <- mean(sub[, j][!is.na(sub[, j])])
      if (!is.na(val)) {
        imm_tprs[i,cn] <- mean(sub[, j][!is.na(sub[, j])])
      }
    } else {
      imm_tprs[i,cn]<-sub[1, j]
      
    }
    
  }
  
}
colnames(imm_tprs) <- colnames(og_white_all_tprs)[-3]
ct_start_idx <- 36
imm_tprs$sex <- droplevels(imm_tprs$sex)

bld_covars <- imm_tprs[,8:25]
ggplot(stack(bld_covars[, c(1:9)], varying = colnames(cont)[c(1:9)]), aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x') + theme_stata()
ggplot(stack(bld_covars[, c(10:18)], varying = colnames(cont)[c(10:18)]), aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x') + theme_stata()

tprs_corr <- cor(unname(imm_tprs[,n_ct+ct_start_idx- 1 + c(3,4,2,5,9,10,12,21,
                                                           14, 22, 15, 17,
                                                           16, 13,8, 1,
                                                           20, 11, 
                                                           18,19,7,6)]),
                 (imm_tprs[,c(8:25)]), use = "pairwise.complete.obs")
corrplot::corrplot(tprs_corr, method = "color")
corrplot::corrplot(tprs_corr, method = "square", is.corr = FALSE)

# by sex
#male
tprs_corr <- cor(unname(imm_tprs[which(imm_tprs$sex == "1"),n_ct+ct_start_idx- 1 + c(3,4,2,5,9,10,12,21,
                                                                                     14, 22, 15, 17,
                                                                                     16, 13,8, 1,
                                                                                     20, 11, 
                                                                                     18,19,7,6)]),
                 (imm_tprs[which(imm_tprs$sex == "1"),c(8:25)]), use = "pairwise.complete.obs")
corrplot::corrplot(tprs_corr, method = "color")
corrplot::corrplot(tprs_corr, method = "square", is.corr = FALSE)

#female
tprs_corr <- cor(unname(imm_tprs[which(imm_tprs$sex == "2"),n_ct+ct_start_idx- 1 + c(3,4,2,5,9,10,12,21,
                                                                                     14, 22, 15, 17,
                                                                                     16, 13,8, 1,
                                                                                     20, 11, 
                                                                                     18,19,7,6)]),
                 (imm_tprs[which(imm_tprs$sex == "2"),c(8:25)]), use = "pairwise.complete.obs")
corrplot::corrplot(tprs_corr, method = "color")
corrplot::corrplot(tprs_corr, method = "square", is.corr = FALSE)

#================================ corr plots with p-vals
library(psych)
tprs_corr <- cor(unname(imm_tprs[,n_ct+ct_start_idx- 1 + c(3,4,2,5,9,10,12,21,
                                                           14, 22, 15, 17,
                                                           16, 13,8, 1,
                                                           20, 11, 
                                                           18,19,7,6)]),
                 (imm_tprs[,c(8:25)]), use = "pairwise.complete.obs")
testRes <- corr.test((imm_tprs[,n_ct+ct_start_idx- 1 + c(3,4,2,5,9,10,12,21,
                                                         14, 22, 15, 17,
                                                         16, 13,8, 1,
                                                         20, 11, 
                                                         18,19,7,6)]),
                     (imm_tprs[,c(8:25)]))

corrplot(tprs_corr, p.mat = testRes$p, method = "square", is.cor = FALSE, 
         sig.level = c(0.01, 0.05, 0.1), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'green')

#male
tprs_corr <- cor(unname(imm_tprs[which(imm_tprs$sex == "1"),n_ct+ct_start_idx- 1 + c(3,4,2,5,9,10,12,21,
                                                                                     14, 22, 15, 17,
                                                                                     16, 13,8, 1,
                                                                                     20, 11, 
                                                                                     18,19,7,6)]),
                 (imm_tprs[which(imm_tprs$sex == "1"),c(8:25)]), use = "pairwise.complete.obs")
testRes <- corr.test((imm_tprs[which(imm_tprs$sex == "1"),n_ct+ct_start_idx- 1 + c(3,4,2,5,9,10,12,21,
                                                                                   14, 22, 15, 17,
                                                                                   16, 13,8, 1,
                                                                                   20, 11, 
                                                                                   18,19,7,6)]),
                     (imm_tprs[which(imm_tprs$sex == "1"),c(8:25)]))

corrplot(tprs_corr, p.mat = testRes$p, sig.level = c(0.01, 0.05, 0.1), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'green', is.cor = FALSE,  method = "square")

#female
tprs_corr <- cor(unname(imm_tprs[which(imm_tprs$sex == "2"),n_ct+ct_start_idx- 1 + c(3,4,2,5,9,10,12,21,
                                                                                     14, 22, 15, 17,
                                                                                     16, 13,8, 1,
                                                                                     20, 11, 
                                                                                     18,19,7,6)]),
                 (imm_tprs[which(imm_tprs$sex == "2"),c(8:25)]), use = "pairwise.complete.obs")
testRes <- corr.test((imm_tprs[which(imm_tprs$sex == "2"),n_ct+ct_start_idx- 1 + c(3,4,2,5,9,10,12,21,
                                                                                   14, 22, 15, 17,
                                                                                   16, 13,8, 1,
                                                                                   20, 11, 
                                                                                   18,19,7,6)]),
                     (imm_tprs[which(imm_tprs$sex == "2"),c(8:25)]))

corrplot(tprs_corr, p.mat = testRes$p, sig.level = c(0.01, 0.05, 0.1), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'green', is.cor = FALSE, method = "square")

ggplot(as.data.frame(1:n_ct)) +
  geom_text(
    aes(x = 0, y = 1:n_ct, label = rev(paste0(1:n_ct," = ",cell_types[c(3,4,2,5,9,10,12,21,
                                                                        14, 22, 15, 17,
                                                                        16, 13,8, 1,
                                                                        20, 11, 
                                                                        18,19,7,6)]))),
    hjust = 0,
  ) +
  theme_void() + ggtitle("index to cell type tPRS") + theme(plot.title = element_text(hjust = 0.8))

##================================ tPRS vs. blood counts

ggplot(imm_tprs, aes(x=regulatory_T_cell_w, y=baso)) + 
  geom_point() + theme_stata()
ggplot(imm_tprs, aes(x=plasmablast_w, y=baso)) + 
  geom_point() + theme_stata()
ggplot(imm_tprs, aes(x=`CD14-low_CD16-positive_monocyte_w`, y=baso)) + 
  geom_point() + theme_stata()
ggplot(imm_tprs, aes(x=`CD4-positive_alpha-beta_T_cell_w`, y=baso)) + 
  geom_point() + theme_stata()

ggplot(imm_tprs, aes(x=`CD4-positive_alpha-beta_T_cell_w`, y=imm_gran)) + 
  geom_point() + theme_stata()
ggplot(imm_tprs, aes(x=`transitional_stage_B_cell_w`, y=imm_gran)) + 
  geom_point() + theme_stata()
ggplot(imm_tprs, aes(x=`double_negative_thymocyte_w`, y=imm_gran)) + 
  geom_point() + theme_stata()


ggplot(imm_tprs, aes(x=`double_negative_thymocyte_w`, y=rdw)) + 
  geom_point() + theme_stata()
ggplot(imm_tprs, aes(x=`CD14-low_CD16-positive_monocyte_w`, y=mchc)) + 
  geom_point() + theme_stata()
ggplot(imm_tprs, aes(x=`naive_B_cell_w`, y=mchc)) + 
  geom_point() + theme_stata()
#plr  chol rbc  hemcrit

sex_line_grp <- as.numeric(imm_tprs$sex) + 2
ggplot(imm_tprs[!is.na(imm_tprs$plr) & !is.na(imm_tprs$sex),], aes(x = `gamma-delta_T_cell_w`, y = plr, group = sex,
                            color = sex)) +
  scale_color_manual(labels = c("male", "female", 
                                "male linear fit", 
                                "female linear fit"), values = c("darkslategray3", "coral1", 
                                                                 "blue", "brown2"))+
  geom_point(alpha = 0.3) +
  geom_smooth(method=lm, aes(color = factor(sex_line_grp)[!is.na(imm_tprs$plr)& !is.na(imm_tprs$sex)]))+
  theme_stata() + labs(title ="Platelet lymphocyte ratio vs. gamma delta T cell tPRS by sex") + 
  ylab("plr") 


ggplot(imm_tprs[which(imm_tprs$sex=="1"),], aes(x=`gamma-delta_T_cell_w`, y=plr)) + 
  geom_point() + theme_stata()
ggplot(imm_tprs[which(imm_tprs$sex=="2"),], aes(x=`gamma-delta_T_cell_w`, y=plr)) + 
  geom_point() + theme_stata()

ggplot(imm_tprs[!is.na(imm_tprs$rbc) & !is.na(imm_tprs$sex),], 
       aes(x = `gamma-delta_T_cell_w`, y = rbc, group = sex, color = sex)) +
  scale_color_manual(labels = c("male", "female", 
                                "male linear fit", 
                                "female linear fit"), values = c("darkslategray3", "coral1", 
                                                                 "blue", "brown2"))+
  geom_point(alpha = 0.3) +
  geom_smooth(method=lm, aes(color = factor(sex_line_grp)[!is.na(imm_tprs$rbc)& !is.na(imm_tprs$sex)]))+
  theme_stata() + labs(title ="Red blood cell count vs. gamma delta T cell tPRS by sex") + 
  ylab("rbc") 

ggplot(imm_tprs[which(imm_tprs$sex=="1"),], aes(x=`gamma-delta_T_cell_w`, y=rbc)) + 
  geom_point() + theme_stata()
ggplot(imm_tprs[which(imm_tprs$sex=="2"),], aes(x=`gamma-delta_T_cell_w`, y=rbc)) + 
  geom_point() + theme_stata()

ggplot(imm_tprs[!is.na(imm_tprs$hemcrit) & !is.na(imm_tprs$sex),], 
       aes(x = `gamma-delta_T_cell_w`, y = hemcrit, group = sex, color = sex)) +
  scale_color_manual(labels = c("male", "female", 
                                "male linear fit", 
                                "female linear fit"), values = c("darkslategray3", "coral1", 
                                                                 "blue", "brown2"))+
  geom_point(alpha = 0.3) +
  geom_smooth(method=lm, aes(color = factor(sex_line_grp)[!is.na(imm_tprs$hemcrit)& !is.na(imm_tprs$sex)]))+
  theme_stata() + labs(title ="Hematocrit vs. gamma delta T cell tPRS by sex") + 
  ylab("hemcrit") 
ggplot(imm_tprs[which(imm_tprs$sex=="1"),], aes(x=`gamma-delta_T_cell_w`, y=hemcrit)) + 
  geom_point() + theme_stata()
ggplot(imm_tprs[which(imm_tprs$sex=="2"),], aes(x=`gamma-delta_T_cell_w`, y=hemcrit)) + 
  geom_point() + theme_stata()

#=================== blood vs. anx, with

tprs_corr <- cbind(cbind(cor((imm_tprs[,8:25]),
                 (imm_tprs[,c(6:7)]), use = "pairwise.complete.obs"), 
                 cor((imm_tprs[which(imm_tprs$sex == "1"),8:25]), (imm_tprs[which(imm_tprs$sex == "1"),c(6:7)]), use = "pairwise.complete.obs")),
                 cor((imm_tprs[which(imm_tprs$sex == "2"),8:25]), (imm_tprs[which(imm_tprs$sex == "2"),c(6:7)]), use = "pairwise.complete.obs"))
colnames(tprs_corr) = c("anx/dep", "with/dep", "anx/dep M", "with/dep M",
                        "anx/dep F", "with/dep F")
testRes <- cbind(cbind(corr.test((imm_tprs[,8:25]),
                     (imm_tprs[,c(6:7)]))$p, 
                     corr.test((imm_tprs[which(imm_tprs$sex == "1"),8:25]),
                     (imm_tprs[which(imm_tprs$sex == "1"),c(6:7)]))$p), 
                 corr.test((imm_tprs[which(imm_tprs$sex == "2"),8:25]),
                           (imm_tprs[which(imm_tprs$sex == "2"),c(6:7)]))$p)

corrplot(tprs_corr, p.mat = testRes, method = "square", is.cor = F, 
         sig.level = c(0.01, 0.05, 0.1), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'green', cl.ratio = 0.5)


tprs_corr <- cbind(cbind(cor((imm_tprs[,n_ct+ct_start_idx- 1 + c(3,4,2,5,9,10,12,21,
                                                                 14, 22, 15, 17,
                                                                 16, 13,8, 1,
                                                                 20, 11, 
                                                                 18,19,7,6)]),
                             (imm_tprs[,c(6:7)]), use = "pairwise.complete.obs"), 
                         cor((imm_tprs[which(imm_tprs$sex == "1"),n_ct+ct_start_idx- 1 + c(3,4,2,5,9,10,12,21,
                                                                                           14, 22, 15, 17,
                                                                                           16, 13,8, 1,
                                                                                           20, 11, 
                                                                                           18,19,7,6)]), (imm_tprs[which(imm_tprs$sex == "1"),c(6:7)]), use = "pairwise.complete.obs")),
                   cor((imm_tprs[which(imm_tprs$sex == "2"),n_ct+ct_start_idx- 1 + c(3,4,2,5,9,10,12,21,
                                                                                     14, 22, 15, 17,
                                                                                     16, 13,8, 1,
                                                                                     20, 11, 
                                                                                     18,19,7,6)]), (imm_tprs[which(imm_tprs$sex == "2"),c(6:7)]), use = "pairwise.complete.obs"))
colnames(tprs_corr) = c("anx/dep", "with/dep", "anx/dep M", "with/dep M",
                        "anx/dep F", "with/dep F")
testRes <- cbind(cbind(corr.test((imm_tprs[,n_ct+ct_start_idx- 1 + c(3,4,2,5,9,10,12,21,
                                                                     14, 22, 15, 17,
                                                                     16, 13,8, 1,
                                                                     20, 11, 
                                                                     18,19,7,6)]),
                                 (imm_tprs[,c(6:7)]))$p, 
                       corr.test((imm_tprs[which(imm_tprs$sex == "1"),n_ct+ct_start_idx- 1 + c(3,4,2,5,9,10,12,21,
                                                                                               14, 22, 15, 17,
                                                                                               16, 13,8, 1,
                                                                                               20, 11, 
                                                                                               18,19,7,6)]),
                                 (imm_tprs[which(imm_tprs$sex == "1"),c(6:7)]))$p), 
                 corr.test((imm_tprs[which(imm_tprs$sex == "2"),n_ct+ct_start_idx- 1 + c(3,4,2,5,9,10,12,21,
                                                                                         14, 22, 15, 17,
                                                                                         16, 13,8, 1,
                                                                                         20, 11, 
                                                                                         18,19,7,6)]),
                           (imm_tprs[which(imm_tprs$sex == "2"),c(6:7)]))$p)

corrplot(tprs_corr, p.mat = testRes, method = "square", is.cor = F, 
         sig.level = c(0.01, 0.05, 0.1), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'green', cl.ratio = 0.5)


# TWAS tPRSs v blood =================================


tprs_corr <- cor(unname(twas_wd_tPRS_w[match(imm_tprs$IID, twas_wd_tPRS_w$IID),23 + c(3,4,2,5,9,10,12,21,
                                                       14, 22, 15, 17,
                                                       16, 13,8, 1,
                                                       20, 11, 
                                                       18,19,7,6)][which(imm_tprs$sex == "1"),]),
                 (imm_tprs[which(imm_tprs$sex == "1"),8:25]), 
                 use = "pairwise.complete.obs")
testRes <- corr.test(twas_wd_tPRS_w[match(imm_tprs$IID, twas_wd_tPRS_w$IID),23 + c(3,4,2,5,9,10,12,21,
                                                     14, 22, 15, 17,
                                                     16, 13,8, 1,
                                                     20, 11, 
                                                     18,19,7,6)][which(imm_tprs$sex == "1"),],
                     (imm_tprs[which(imm_tprs$sex == "1"),8:25]))

corrplot(tprs_corr, p.mat = testRes$p, sig.level = c(0.01, 0.05, 0.1), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'green', is.cor = FALSE, method = "square")
