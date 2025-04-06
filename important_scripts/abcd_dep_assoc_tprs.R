# required: load_abcd_data.R
# goal: linear regressions with Wittenberg tPRS
# next: glm regressions

load("sc-tprs-mdd/og_white_all_tprs_immune.Rda")
if (! requireNamespace("tidyverse", quietly=TRUE)) {
  install.packages("tidyverse")
}
library(tidyverse)
if (! requireNamespace("vtable", quietly=TRUE)) {
  install.packages("vtable")
}
library(vtable)
packageurl <- "https://github.com/tidyverse/tidyr/archive/refs/tags/v1.2.1.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
library(tidyr)
events <- c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1", "2_year_follow_up_y_arm_1",
            "3_year_follow_up_y_arm_1", "4_year_follow_up_y_arm_1")





#define function to extract overall p-value of model
overall_p <- function(my_model) {
  f <- summary(my_model)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}



#data frames to store linear regressions of depression scores vs. tprs by cell type
reg_tprs <- list(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0", 
                                   "slr_ad1","slr_wd1", "mlr_ad1", "mlr_wd1",
                                   "slr_ad2","slr_wd2", "mlr_ad2", "mlr_wd2",
                                   "slr_ad3","slr_wd3", "mlr_ad3", "mlr_wd3",
                                   "slr_ad4","slr_wd4", "mlr_ad4", "mlr_wd4"))

n_ct <- length(cell_types)
tprs_start_idx <- 19

#try power transform on dep scores, but the transformations are a bit nonsensical
car::powerTransform((og_white_all_tprs[og_white_all_tprs$event_name == events[1],][, 5]))
a <- og_white_all_tprs[c(og_white_all_tprs$event_name == events[1]),][, 7][-c(1026)]
car::powerTransform(a)
b <- og_white_all_tprs[c(og_white_all_tprs$event_name == events[1]),][, 8][-c(1026)]
car::powerTransform(b)

#og_white_all_tprs <- white_all_tprs
white_all_tprs <- og_white_all_tprs[,-(9:26)]
white_all_tprs[, c(7,8)] <- sqrt(og_white_all_tprs[, c(7,8)])
for (event in events){
  white_all_tprs[(white_all_tprs$event_name == event),][, -c(1, 2, 3, 4, 6, 7, 8)] <- 
    scale(white_all_tprs[(white_all_tprs$event_name == event),][, -c(1, 2, 3, 4, 6, 7, 8)])
  
}
white_all_tprs <- white_all_tprs[-(which(white_all_tprs$sex == 3)),] #remove intersex-male
base_white_tprs <- white_all_tprs[(white_all_tprs$event_name == events[1]),]
ggplot(stack( (base_white_tprs[, c(7,8)]), varying = colnames(cont)[7:8]), aes(values)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~ind, scales = 'free_x')

cont <- base_white_tprs[c(white_all_tprs$event_name == events[1]),][, -c(1, 2, 3, 4, 6)]
catg <- base_white_tprs[c(white_all_tprs$event_name == events[1]),][, c(4, 6)]
catg$sex <- as.character(catg$sex)
catg <- catg[, -c(1)]
ggplot(stack(cont, varying = colnames(cont)), aes(values)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~ind, scales = 'free_x')
ggplot(stack( (cont[, c(2,3)]), varying = colnames(cont)[2:3]), aes(values)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~ind, scales = 'free_x')
ggplot(stack(catg, varying = colnames(catg)), aes(values)) + 
  geom_bar() + 
  facet_wrap(~ind, scales = 'free_x')
sumtable(cbind(cont, catg))

tprs_corr <- cor(unname(base_white_all_tprs[tprs_start_idx:(tprs_start_idx - 1 + n_ct)]))
corrplot::corrplot(tprs_corr, method = "square")
tprs_corr <- cor(unname(base_white_all_tprs[(tprs_start_idx + n_ct):length(white_all_tprs)]))
corrplot::corrplot(tprs_corr, method = "square")

adjR2 <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
pval_mod <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
pval_tprs_slope <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
tprs_beta <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
tprs_CI <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
R2 <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))

for (i in tprs_start_idx:ncol(white_all_tprs)){
  for (j in 1:1){
    reg_tprs[[colnames(white_all_tprs)[i]]][[1]] <- lm(
       (base_white_tprs$anxdep_t) ~ base_white_tprs[[i]])
    adjR2[[colnames(white_all_tprs)[i]]][1] <- 
      summary(reg_tprs[[colnames(white_all_tprs)[i]]][[1]])$adj.r.squared
    R2[[colnames(white_all_tprs)[i]]][1] <- 
      summary(reg_tprs[[colnames(white_all_tprs)[i]]][[1]])$r.squared
    pval_mod[[colnames(white_all_tprs)[i]]][1] <- 
      overall_p(reg_tprs[[colnames(white_all_tprs)[i]]][[1]])
    pval_tprs_slope[[colnames(white_all_tprs)[i]]][1] <- 
      summary(reg_tprs[[colnames(white_all_tprs)[i]]][[1]])$coefficients[,4][2]
    tprs_beta[[colnames(white_all_tprs)[i]]][1] <- 
      coef(reg_tprs[[colnames(white_all_tprs)[i]]][[1]])[2]
    tprs_CI[[colnames(white_all_tprs)[i]]][1] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(white_all_tprs)[i]]][[1]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(white_all_tprs)[i]]][[1]]))[2])
    
    
    reg_tprs[[colnames(white_all_tprs)[i]]][[2]] <- lm(
       (base_white_tprs$withdep_t) ~ base_white_tprs[[i]])
    adjR2[[colnames(white_all_tprs)[i]]][2] <- 
      summary(reg_tprs[[colnames(white_all_tprs)[i]]][[2]])$adj.r.squared
    R2[[colnames(white_all_tprs)[i]]][2] <- 
      summary(reg_tprs[[colnames(white_all_tprs)[i]]][[2]])$r.squared
    pval_mod[[colnames(white_all_tprs)[i]]][2] <- 
      overall_p(reg_tprs[[colnames(white_all_tprs)[i]]][[2]])
    pval_tprs_slope[[colnames(white_all_tprs)[i]]][2] <- 
      summary(reg_tprs[[colnames(white_all_tprs)[i]]][[2]])$coefficients[,4][2]
    tprs_beta[[colnames(white_all_tprs)[i]]][2] <- 
      coef(reg_tprs[[colnames(white_all_tprs)[i]]][[2]])[2]
    tprs_CI[[colnames(white_all_tprs)[i]]][2] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(white_all_tprs)[i]]][[2]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(white_all_tprs)[i]]][[2]]))[2])
    
    reg_tprs[[colnames(white_all_tprs)[i]]][[3]] <- lm(
       (base_white_tprs$anxdep_t) ~ 
         base_white_tprs[[i]]  + 
         base_white_tprs[["age"]]  +
         base_white_tprs[["sex"]]  +
         base_white_tprs[["site"]] +
         base_white_tprs[["genetic_pc_1"]]  +
         base_white_tprs[["genetic_pc_2"]]  +
         base_white_tprs[["genetic_pc_3"]]  +
         base_white_tprs[["genetic_pc_4"]]  +
         base_white_tprs[["genetic_pc_5"]]  +
         base_white_tprs[["genetic_pc_6"]]  +
         base_white_tprs[["genetic_pc_7"]]  +
         base_white_tprs[["genetic_pc_8"]]  +
         base_white_tprs[["genetic_pc_9"]]  +
         base_white_tprs[["genetic_pc_10"]] )
    adjR2[[colnames(white_all_tprs)[i]]][3] <- 
      summary(reg_tprs[[colnames(white_all_tprs)[i]]][[3]])$adj.r.squared
    R2[[colnames(white_all_tprs)[i]]][3] <- 
      summary(reg_tprs[[colnames(white_all_tprs)[i]]][[3]])$r.squared
    pval_mod[[colnames(white_all_tprs)[i]]][3] <- 
      overall_p(reg_tprs[[colnames(white_all_tprs)[i]]][[3]])
    pval_tprs_slope[[colnames(white_all_tprs)[i]]][3] <- 
      summary(reg_tprs[[colnames(white_all_tprs)[i]]][[3]])$coefficients[,4][2]
    tprs_beta[[colnames(white_all_tprs)[i]]][3] <- 
      coef(reg_tprs[[colnames(white_all_tprs)[i]]][[3]])[2]
    tprs_CI[[colnames(white_all_tprs)[i]]][3] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(white_all_tprs)[i]]][[3]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(white_all_tprs)[i]]][[3]]))[2])
    
    reg_tprs[[colnames(white_all_tprs)[i]]][[4]] <- lm(
       ( base_white_tprs$withdep_t) ~ 
         base_white_tprs[[i]]  + 
         base_white_tprs[["age"]]  +
         base_white_tprs[["sex"]]  +
         base_white_tprs[["site"]] +
         base_white_tprs[["genetic_pc_1"]]  +
         base_white_tprs[["genetic_pc_2"]]  +
         base_white_tprs[["genetic_pc_3"]]  +
         base_white_tprs[["genetic_pc_4"]]  +
         base_white_tprs[["genetic_pc_5"]]  +
         base_white_tprs[["genetic_pc_6"]]  +
         base_white_tprs[["genetic_pc_7"]]  +
         base_white_tprs[["genetic_pc_8"]]  +
         base_white_tprs[["genetic_pc_9"]]  +
         base_white_tprs[["genetic_pc_10"]] )
    adjR2[[colnames(white_all_tprs)[i]]][4] <- 
      summary(reg_tprs[[colnames( base_white_tprs)[i]]][[4]])$adj.r.squared
    R2[[colnames(white_all_tprs)[i]]][4] <- 
      summary(reg_tprs[[colnames(white_all_tprs)[i]]][[4]])$r.squared
    pval_mod[[colnames(white_all_tprs)[i]]][4] <- 
      overall_p(reg_tprs[[colnames(white_all_tprs)[i]]][[4]])
    pval_tprs_slope[[colnames(white_all_tprs)[i]]][4] <- 
      summary(reg_tprs[[colnames(white_all_tprs)[i]]][[4]])$coefficients[,4][2]
    tprs_beta[[colnames(white_all_tprs)[i]]][4] <- 
      coef(reg_tprs[[colnames(white_all_tprs)[i]]][[4]])[2]
    tprs_CI[[colnames(white_all_tprs)[i]]][4] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(white_all_tprs)[i]]][[4]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(white_all_tprs)[i]]][[4]]))[2])
  }
  
}

pval_ftest <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
for (i in tprs_start_idx:ncol(white_all_tprs)){
  for (j in 1:1){
    red_ad <- lm(
       ( base_white_tprs$anxdep_t ) ~ 1)
    red_wd <- lm(
       ( base_white_tprs$withdep_t ) ~ 1)
    reduced_ad <- lm( ( base_white_tprs$anxdep_t) ~
                        base_white_tprs[["age"]]  +
                        base_white_tprs[["sex"]]  +
                        base_white_tprs[["genetic_pc_1"]]  +
                        base_white_tprs[["genetic_pc_2"]]  +
                        base_white_tprs[["genetic_pc_3"]]  +
                        base_white_tprs[["genetic_pc_4"]]  +
                        base_white_tprs[["genetic_pc_5"]]  +
                        base_white_tprs[["genetic_pc_6"]]  +
                        base_white_tprs[["genetic_pc_7"]]  +
                        base_white_tprs[["genetic_pc_8"]]  +
                        base_white_tprs[["genetic_pc_9"]]  +
                        base_white_tprs[["genetic_pc_10"]] )
    reduced_wd <- lm( ( base_white_tprs$withdep_t) ~ 
                     base_white_tprs[["age"]]  +
                     base_white_tprs[["sex"]]  +
                     base_white_tprs[["genetic_pc_1"]]  +
                     base_white_tprs[["genetic_pc_2"]]  +
                     base_white_tprs[["genetic_pc_3"]]  +
                     base_white_tprs[["genetic_pc_4"]]  +
                     base_white_tprs[["genetic_pc_5"]]  +
                     base_white_tprs[["genetic_pc_6"]]  +
                     base_white_tprs[["genetic_pc_7"]]  +
                     base_white_tprs[["genetic_pc_8"]]  +
                     base_white_tprs[["genetic_pc_9"]]  +
                     base_white_tprs[["genetic_pc_10"]] )
    
    pval_ftest[[colnames(white_all_tprs)[i]]][1] <- anova(reg_tprs[[colnames(white_all_tprs)[i]]][[1]], red_ad)$`Pr(>F)`[2]
    pval_ftest[[colnames(white_all_tprs)[i]]][2] <- anova(reg_tprs[[colnames(white_all_tprs)[i]]][[2]], red_wd)$`Pr(>F)`[2]
    pval_ftest[[colnames(white_all_tprs)[i]]][3] <- anova(reg_tprs[[colnames(white_all_tprs)[i]]][[3]], reduced_ad)$`Pr(>F)`[2]
    pval_ftest[[colnames(white_all_tprs)[i]]][4] <- anova(reg_tprs[[colnames(white_all_tprs)[i]]][[4]], reduced_wd)$`Pr(>F)`[2]
    
    
  }   
}



# try all cell types as covars for fun
ct_ad <- lm(
   ( base_white_tprs$anxdep_t ) ~ 
     `CD14-low_CD16-positive_monocyte_w` + 
     `CD4-positive_alpha-beta_cytotoxic_T_cell_w` +
     `CD4-positive_alpha-beta_T_cell_w` +
     `CD8-positive_alpha-beta_T_cell_w` +
     `central_memory_CD8-positive_alpha-beta_T_cell_w` +
     conventional_dendritic_cell_w +
     dendritic_cell_w +
     double_negative_thymocyte_w +
     `effector_memory_CD4-positive_alpha-beta_T_cell_w` +
     `effector_memory_CD8-positive_alpha-beta_T_cell_w` +
     erythrocyte_w  +
     `gamma-delta_T_cell_w`  +
     innate_lymphoid_cell_w  +
     memory_B_cell_w  +
     naive_B_cell_w  +
     innate_lymphoid_cell_w  +
     peripheral_blood_mononuclear_cell_w +
     plasmablast_w  +
     plasmacytoid_dendritic_cell_w +
     platelet_w  +
     regulatory_T_cell_w +
     transitional_stage_B_cell_w, data = base_white_tprs )
ct_wd <- lm(
   ( base_white_tprs$withdep_t ) ~ 
    `CD14-low_CD16-positive_monocyte_w` + 
    `CD4-positive_alpha-beta_cytotoxic_T_cell_w` +
    `CD4-positive_alpha-beta_T_cell_w` +
    `CD8-positive_alpha-beta_T_cell_w` +
    `central_memory_CD8-positive_alpha-beta_T_cell_w` +
    `conventional_dendritic_cell_w` +
    dendritic_cell_w +
    double_negative_thymocyte_w +
    `effector_memory_CD4-positive_alpha-beta_T_cell_w` +
    `effector_memory_CD8-positive_alpha-beta_T_cell_w` +
    erythrocyte_w  +
    `gamma-delta_T_cell_w`  +
    innate_lymphoid_cell_w  +
    memory_B_cell_w  +
    naive_B_cell_w  +
    innate_lymphoid_cell_w  +
    peripheral_blood_mononuclear_cell_w +
    plasmablast_w  +
    plasmacytoid_dendritic_cell_w +
    platelet_w  +
    regulatory_T_cell_w +
    transitional_stage_B_cell_w, data = base_white_tprs )

#library(car)

# investigate correlation betwee tPRS via VIF
vif(ct_wd)
vif(ct_ad)

# try mixed effects model...need to think about the design a bit more
library(lme4)

mixed_ad <- lmer( (white_all_tprs$anxdep_t) ~  white_all_tprs[[40]] + age + sex + 
                   genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 +
                   genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
                   genetic_pc_9 + genetic_pc_10 + relevel(event_name, ref = "baseline_year_1_arm_1")*white_all_tprs[[40]] + (1|IID) + (1|site), 
                 data = white_all_tprs)
summary(mixed_ad)

pval_mix_time <- data.frame(regs = c("ad0","ad2", "ad3", "ad4", "wd0","wd2", "wd3", "wd4"))

# forest plots for linear regressions================================
library(ggstance)
library(patchwork)
library(scales)
#only unweighted
ad_beta <- c(tprs_beta[1,][, -c(1, (n_ct + 2):(2*n_ct + 1))], tprs_beta[3,][, -c(1,  (n_ct + 2):(2*n_ct + 1))])
wd_beta <- c(tprs_beta[2,][, -c(1,  (n_ct + 2):(2*n_ct + 1))], tprs_beta[4,][, -c(1,  (n_ct + 2):(2*n_ct + 1))])
labs <- c(cell_types, cell_types)
fit_labs <- c(rep("w/o covar", n_ct), rep("w/ covar", n_ct))
ad_lower <- c((ad_beta[1:n_ct] - tprs_CI[1,][, -c(1,  (n_ct + 2):(2*n_ct + 1))]), 
              ad_beta[(n_ct + 1):(2*n_ct)] - tprs_CI[3,][, -c(1,  (n_ct + 2):(2*n_ct + 1))])
wd_lower <- c((wd_beta[1:n_ct] - tprs_CI[2,][, -c(1,  (n_ct + 2):(2*n_ct + 1))]),
              wd_beta[(n_ct + 1):(2*n_ct)] - tprs_CI[4,][, -c(1,  (n_ct + 2):(2*n_ct + 1))])
ad_upper <- c((ad_beta[1:n_ct] + tprs_CI[1,][, -c(1,  (n_ct + 2):(2*n_ct + 1))]), 
              ad_beta[(n_ct + 1):(2*n_ct)] + tprs_CI[3,][, -c(1,  (n_ct + 2):(2*n_ct + 1))])
wd_upper <- c((wd_beta[1:n_ct] + tprs_CI[2,][, -c(1,  (n_ct + 2):(2*n_ct + 1))]),
              wd_beta[(n_ct + 1):(2*n_ct)] + tprs_CI[4,][, -c(1,  (n_ct + 2):(2*n_ct + 1))])

ad_df <- data.frame(betas = as.numeric(unname(ad_beta)), cell_types = unname(labs), 
                        upper = as.numeric(unname(ad_upper)), 
                      lower = as.numeric(unname(ad_lower)), tPRS_type = unname(fit_labs),
                      p_values = paste0("(", round(as.numeric(unname(pval_tprs_slope[1,][, -c(1)]))[1:n_ct], 3),
                                        ", ", round(as.numeric(unname(pval_tprs_slope[3,][, -c(1)]))[1:n_ct], 3),
                                        ")"),
                      R2 = paste0("(", percent(round(as.numeric(unname(R2[1,][, -c(1)]))[1:n_ct], 4)),
                                        ", ", percent(round(as.numeric(unname(R2[3,][, -c(1)]))[1:n_ct], 4)),
                                        ")"))

wd_df <- data.frame(betas = as.numeric(unname(wd_beta)), cell_types = unname(labs), 
                        upper = as.numeric(unname(wd_upper)), 
                        lower = as.numeric(unname(wd_lower)), tPRS_type = unname(fit_labs),
                        p_values = paste0("(", round(as.numeric(unname(pval_tprs_slope[2,][, -c(1)]))[1:n_ct], 3),
                                          ", ", round(as.numeric(unname(pval_tprs_slope[4,][, -c(1)]))[1:n_ct], 3),
                                          ")"),
                        R2 = paste0("(", percent(round(as.numeric(unname(R2[2,][, -c(1)]))[1:n_ct], 4)),
                                    ", ", percent(round(as.numeric(unname(R2[4,][, -c(1)]))[1:n_ct], 4)),
                                                  ")"))

mid_plot <- ggplot(ad_df, aes(x = betas, y = cell_types, xmin = upper, 
                           xmax = lower, group = tPRS_type, color = tPRS_type)) +
  geom_point(size = 4, aes(shape = tPRS_type), position = position_dodgev(height = .5)) +
  geom_pointrange(shape = 1, position = position_dodgev(height = .5)) +
  theme(legend.position = "bottom", axis.text.y = element_text(angle = 0, size = 11)) + 
  ggtitle("Main effect size of tPRSs by cell type with\nanxious/depressed score as outcome") +
  xlab("tPRS main effect sizes") + ylab("cell types")
mid_plot <- mid_plot +
  geom_vline(xintercept = 0, linetype="dashed") + theme_hc()
mid_plot


p_right <- 
  ggplot(ad_df) +
  geom_text(
    aes(x = 0, y = rev(cell_types), label = R2),
    hjust = 0,
  ) + 
  theme_void() + ggtitle( paste0("% variance explained",":\n(w/o covar, w/ covar)")) + theme(plot.title = element_text(hjust = 0.8))
p_right

layout <- c(
  area(t = 0, l = 0, b = 30, r = 9), # left plot, starts at the top of the page (0) and goes 30 units down and 3 units to the right
  area(t = 0, l = 7, b = 30, r = 15) # right most plot starts at top of page, begins where middle plot ends (l=9, and middle plot is r=9), goes to bottom of page (b=30), and extends two units wide (r=11)
)
# final plot arrangement
mid_plot + p_right + plot_layout(design = layout)

mid_plot <- ggplot(wd_df, aes(x = betas, y = cell_types, xmin = upper, 
                                  xmax = lower, group = tPRS_type, color = tPRS_type)) +
  geom_point(size = 4, aes(shape = tPRS_type), position = position_dodgev(height = .5)) +
  geom_pointrange(shape = 1, position = position_dodgev(height = .5)) +
  theme(legend.position = "bottom", axis.text.y = element_text(angle = 0, size = 11)) + 
  ggtitle("Main effect size of tPRSs by cell type with\nwithdrawn/depressed score as outcome") +
  xlab("tPRS main effect sizes") + ylab("cell types")
mid_plot <- mid_plot +
  geom_vline(xintercept = 0, linetype="dashed")+ theme_hc()
mid_plot


p_right <- 
  ggplot(wd_df) +
  geom_text(
    aes(x = 0, y = rev(cell_types), label = R2),
    hjust = 0,
  ) + 
  theme_void() + ggtitle( paste0("% variance explained",":\n(w/o covar, w/ covar)")) + theme(plot.title = element_text(hjust = 0.8))
p_right

# final plot arrangement
mid_plot + p_right + plot_layout(design = layout)

#only weighted
ad_beta <- c(tprs_beta[1,][, -c(1:(n_ct + 1))], tprs_beta[3,][, -c(1:(n_ct + 1))])
wd_beta <- c(tprs_beta[2,][, -c(1:(n_ct + 1))], tprs_beta[4,][, -c(1:(n_ct + 1))])
labs <- c(cell_types, cell_types)
fit_labs <- c(rep("w/o covar", n_ct), rep("w/ covar", n_ct))
ad_lower <- c((ad_beta[1:n_ct] - tprs_CI[1,][, -c(1:(n_ct + 1))]), 
              ad_beta[(n_ct + 1):(2*n_ct)] - tprs_CI[3,][, -c(1:(n_ct + 1))])
wd_lower <- c((wd_beta[1:n_ct] - tprs_CI[2,][, -c(1:(n_ct + 1))]),
              wd_beta[(n_ct + 1):(2*n_ct)] - tprs_CI[4,][, -c(1:(n_ct + 1))])
ad_upper <- c((ad_beta[1:n_ct] + tprs_CI[1,][, -c(1:(n_ct + 1))]), 
              ad_beta[(n_ct + 1):(2*n_ct)] + tprs_CI[3,][, -c(1:(n_ct + 1))])
wd_upper <- c((wd_beta[1:n_ct] + tprs_CI[2,][, -c(1:(n_ct + 1))]),
              wd_beta[(n_ct + 1):(2*n_ct)] + tprs_CI[4,][, -c(1:(n_ct + 1))])

ad_df <- data.frame(betas = as.numeric(unname(ad_beta)), cell_types = unname(labs), 
                    upper = as.numeric(unname(ad_upper)), 
                    lower = as.numeric(unname(ad_lower)), tPRS_type = unname(fit_labs),
                    p_values = paste0("(", round(as.numeric(unname(pval_tprs_slope[1,][, -c(1)]))[(n_ct + 1):(2*n_ct)], 3),
                                      ", ", round(as.numeric(unname(pval_tprs_slope[3,][, -c(1)]))[(n_ct + 1):(2*n_ct)], 3),
                                      ")"),
                    R2 = paste0("(", percent(round(as.numeric(unname(R2[1,][, -c(1)]))[(n_ct + 1):(2*n_ct)], 4)),
                                ", ", percent(round(as.numeric(unname(R2[3,][, -c(1)]))[(n_ct + 1):(2*n_ct)], 4)),
                                ")"))

wd_df <- data.frame(betas = as.numeric(unname(wd_beta)), cell_types = unname(labs), 
                    upper = as.numeric(unname(wd_upper)), 
                    lower = as.numeric(unname(wd_lower)), tPRS_type = unname(fit_labs),
                    p_values = paste0("(", round(as.numeric(unname(pval_tprs_slope[2,][, -c(1)]))[(n_ct + 1):(2*n_ct)], 3),
                                      ", ", round(as.numeric(unname(pval_tprs_slope[4,][, -c(1)]))[(n_ct + 1):(2*n_ct)], 3),
                                      ")"),
                    R2 = paste0("(", percent(round(as.numeric(unname(R2[2,][, -c(1)]))[(n_ct + 1):(2*n_ct)], 4)),
                                ", ", percent(round(as.numeric(unname(R2[4,][, -c(1)]))[(n_ct + 1):(2*n_ct)], 4)),
                                ")"))

mid_plot <- ggplot(ad_df, aes(x = betas, y = cell_types, xmin = upper, 
                              xmax = lower, group = tPRS_type, color = tPRS_type)) +
  geom_point(size = 4, aes(shape = tPRS_type), position = position_dodgev(height = .5)) +
  geom_pointrange(shape = 1, position = position_dodgev(height = .5)) +
  theme(legend.position = "bottom", axis.text.y = element_text(angle = 0, size = 11)) + 
  ggtitle("Main effect size of tPRSs by cell type with\nanxious/depressed score as outcome") +
  xlab("tPRS main effect sizes") + ylab("cell types")
mid_plot <- mid_plot +
  geom_vline(xintercept = 0, linetype="dashed") + theme_hc()
mid_plot


p_right <- 
  ggplot(ad_df) +
  geom_text(
    aes(x = 0, y = rev(cell_types), label = R2),
    hjust = 0,
  ) + 
  theme_void() + ggtitle( paste0("% variance explained",":\n(w/o covar, w/ covar)")) + theme(plot.title = element_text(hjust = 0.8))
p_right

layout <- c(
  area(t = 0, l = 0, b = 30, r = 9), # left plot, starts at the top of the page (0) and goes 30 units down and 3 units to the right
  area(t = 0, l = 7, b = 30, r = 15) # right most plot starts at top of page, begins where middle plot ends (l=9, and middle plot is r=9), goes to bottom of page (b=30), and extends two units wide (r=11)
)
# final plot arrangement
mid_plot + p_right + plot_layout(design = layout)

mid_plot <- ggplot(wd_df, aes(x = betas, y = cell_types, xmin = upper, 
                              xmax = lower, group = tPRS_type, color = tPRS_type)) +
  geom_point(size = 4, aes(shape = tPRS_type), position = position_dodgev(height = .5)) +
  geom_pointrange(shape = 1, position = position_dodgev(height = .5)) +
  theme(legend.position = "bottom", axis.text.y = element_text(angle = 0, size = 11)) + 
  ggtitle("Main effect size of tPRSs by cell type with\nwithdrawn/depressed score as outcome") +
  xlab("tPRS main effect sizes") + ylab("cell types")
mid_plot <- mid_plot +
  geom_vline(xintercept = 0, linetype="dashed")+ theme_hc()
mid_plot


p_right <- 
  ggplot(wd_df) +
  geom_text(
    aes(x = 0, y = rev(cell_types), label = R2),
    hjust = 0,
  ) + 
  theme_void() + ggtitle( paste0("% variance explained",":\n(w/o covar, w/ covar)")) + 
  theme(plot.title = element_text(hjust = 0.8))
p_right

# final plot arrangement
mid_plot + p_right + plot_layout(design = layout)

ct_ind_df <- data.frame(info = paste0(cell_types, " = ", (1:(n_ct))))
p_right <- 
  ggplot(ct_ind_df) +
  geom_text(
    aes(x = 0, y = rev(info), label = info),
    hjust = 1,
  ) + 
  theme_void() + ggtitle("Cell type to index") + 
  theme(plot.title = element_text(hjust = 0.2)) 
  
p_right

#lin reg with sex ixn on tPRS=========================================================
pval_sex_ixn<- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
beta_sex_ixn<- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
CI_sex_ixn <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))

pval_tprs_slope_sx <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
tprs_beta_sx <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
tprs_CI_sx <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
R2_sx <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))

for (i in tprs_start_idx:ncol(white_all_tprs)){
  for (j in 1:1){
    ssex_ixn_ad <- lm( ( base_white_tprs$anxdep_t) ~
                         base_white_tprs[[i]]  + 
                         base_white_tprs[[i]] * base_white_tprs[["sex"]] +
                         base_white_tprs[["sex"]] )
    ssex_ixn_wd <- lm( ( base_white_tprs$withdep_t) ~ 
                         base_white_tprs[[i]]  +
                         base_white_tprs[[i]] * base_white_tprs[["sex"]] +
                         base_white_tprs[["sex"]])
    
    sex_ixn_ad <- lm( ( base_white_tprs$anxdep_t) ~
                        base_white_tprs[[i]]  +
                        base_white_tprs[[i]] * base_white_tprs[["sex"]] +
                        base_white_tprs[["age"]]  +
                        base_white_tprs[["sex"]]  +
                        base_white_tprs[["site"]]  +
                        base_white_tprs[["genetic_pc_1"]]  +
                        base_white_tprs[["genetic_pc_2"]]  +
                        base_white_tprs[["genetic_pc_3"]]  +
                        base_white_tprs[["genetic_pc_4"]]  +
                        base_white_tprs[["genetic_pc_5"]]  +
                        base_white_tprs[["genetic_pc_6"]]  +
                        base_white_tprs[["genetic_pc_7"]]  +
                        base_white_tprs[["genetic_pc_8"]]  +
                        base_white_tprs[["genetic_pc_9"]]  +
                        base_white_tprs[["genetic_pc_10"]] )
    sex_ixn_wd <- lm( ( base_white_tprs$withdep_t) ~ 
                        base_white_tprs[[i]]  +
                        base_white_tprs[[i]] * base_white_tprs[["sex"]] +
                        base_white_tprs[["age"]]  +
                        base_white_tprs[["sex"]]  +
                        base_white_tprs[["site"]]  +
                        base_white_tprs[["genetic_pc_1"]]  +
                        base_white_tprs[["genetic_pc_2"]]  +
                        base_white_tprs[["genetic_pc_3"]]  +
                        base_white_tprs[["genetic_pc_4"]]  +
                        base_white_tprs[["genetic_pc_5"]]  +
                        base_white_tprs[["genetic_pc_6"]]  +
                        base_white_tprs[["genetic_pc_7"]]  +
                        base_white_tprs[["genetic_pc_8"]]  +
                        base_white_tprs[["genetic_pc_9"]]  +
                        base_white_tprs[["genetic_pc_10"]] )
    
    beta_sex_ixn[[colnames(white_all_tprs)[i]]][1] <- coef(ssex_ixn_ad)[3]
    beta_sex_ixn[[colnames(white_all_tprs)[i]]][2] <- coef(ssex_ixn_wd)[3]
    beta_sex_ixn[[colnames(white_all_tprs)[i]]][3] <- coef(sex_ixn_ad)[3]
    beta_sex_ixn[[colnames(white_all_tprs)[i]]][4] <- coef(sex_ixn_wd)[3]
    
    pval_sex_ixn[[colnames(white_all_tprs)[i]]][1] <- summary(ssex_ixn_ad)$coefficients[,4][3]
    pval_sex_ixn[[colnames(white_all_tprs)[i]]][2] <- summary(ssex_ixn_wd)$coefficients[,4][3]
    pval_sex_ixn[[colnames(white_all_tprs)[i]]][3] <- summary(sex_ixn_ad)$coefficients[,4][3]
    pval_sex_ixn[[colnames(white_all_tprs)[i]]][4] <- summary(sex_ixn_wd)$coefficients[,4][3]
    
    CI_sex_ixn[[colnames(white_all_tprs)[i]]][1] <- qt(0.975, df.residual(ssex_ixn_ad))*
      sqrt(diag(vcov(ssex_ixn_ad))[3])
    CI_sex_ixn[[colnames(white_all_tprs)[i]]][2] <- qt(0.975, df.residual(ssex_ixn_wd))*
      sqrt(diag(vcov(ssex_ixn_wd))[3])
    CI_sex_ixn[[colnames(white_all_tprs)[i]]][3] <- qt(0.975, df.residual(sex_ixn_ad))*
      sqrt(diag(vcov(sex_ixn_ad))[3])
    CI_sex_ixn[[colnames(white_all_tprs)[i]]][4] <- qt(0.975, df.residual(sex_ixn_wd))*
      sqrt(diag(vcov(sex_ixn_wd))[3])
    
    pval_tprs_slope_sx[[colnames(white_all_tprs)[i]]][1] <- summary(ssex_ixn_ad)$coefficients[,4][2]
    pval_tprs_slope_sx[[colnames(white_all_tprs)[i]]][2] <- summary(ssex_ixn_wd)$coefficients[,4][2]
    pval_tprs_slope_sx[[colnames(white_all_tprs)[i]]][3] <- summary(sex_ixn_ad)$coefficients[,4][2]
    pval_tprs_slope_sx[[colnames(white_all_tprs)[i]]][4] <- summary(sex_ixn_wd)$coefficients[,4][2]
    
    tprs_beta_sx[[colnames(white_all_tprs)[i]]][1] <- coef(ssex_ixn_ad)[2]
    tprs_beta_sx[[colnames(white_all_tprs)[i]]][2] <- coef(ssex_ixn_wd)[2]
    tprs_beta_sx[[colnames(white_all_tprs)[i]]][3] <- coef(sex_ixn_ad)[2]
    tprs_beta_sx[[colnames(white_all_tprs)[i]]][4] <- coef(sex_ixn_wd)[2]
    
    tprs_CI_sx[[colnames(white_all_tprs)[i]]][1] <- qt(0.975, df.residual(ssex_ixn_ad))*
      sqrt(diag(vcov(ssex_ixn_ad))[2])
    tprs_CI_sx[[colnames(white_all_tprs)[i]]][2] <- qt(0.975, df.residual(ssex_ixn_wd))*
      sqrt(diag(vcov(ssex_ixn_wd))[2])
    tprs_CI_sx[[colnames(white_all_tprs)[i]]][3] <- qt(0.975, df.residual(sex_ixn_ad))*
      sqrt(diag(vcov(sex_ixn_ad))[2])
    tprs_CI_sx[[colnames(white_all_tprs)[i]]][4] <- qt(0.975, df.residual(sex_ixn_wd))*
      sqrt(diag(vcov(sex_ixn_wd))[2])
    
    R2_sx[[colnames(white_all_tprs)[i]]][1] <- summary(ssex_ixn_ad)$r.squared
    R2_sx[[colnames(white_all_tprs)[i]]][2] <- summary(ssex_ixn_wd)$r.squared
    R2_sx[[colnames(white_all_tprs)[i]]][3] <- summary(sex_ixn_ad)$r.squared
    R2_sx[[colnames(white_all_tprs)[i]]][4] <- summary(sex_ixn_wd)$r.squared
    
  }   
}

# sex ixn is significant, look at sex strat lin assoc

# peripheral_blood_mononuclear_cell, effector_memory_CD8-positive_alpha-beta_T_cell,
# conventional_dendritic_cell, double_negative_thymocyte
sex_line_grp <- as.numeric(base_white_tprs$sex) + 2
ggplot(base_white_tprs, aes(x = peripheral_blood_mononuclear_cell_uw, y = anxdep_t, group = sex,
                           color = sex)) +
  scale_color_manual(labels = c("male", "female", 
                                "male linear fit", 
                              "female linear fit"), values = c("darkslategray3", "coral1", 
                                                               "blue", "brown2"))+
  geom_point(alpha = 0.3) +
  geom_smooth(method=lm, aes(color = factor(sex_line_grp)))+
  theme_stata() + labs(title ="Anxious/depressed score vs. PBMC tPRS by sex") + 
  ylab("sqrt anx/dep") +
  geom_hline(yintercept=sqrt(65),  linetype="dashed")

male_biv <- lm (anxdep_t ~ peripheral_blood_mononuclear_cell_uw, 
                data = base_white_tprs[which(base_white_tprs$sex == 1),])
fem_biv <- lm (anxdep_t ~ peripheral_blood_mononuclear_cell_uw, 
               data = base_white_tprs[which(base_white_tprs$sex == 2),])
kable(tidy(male_biv)[, -c(4)])
kable(tidy(fem_biv)[, -c(4)])
ggplot(base_white_tprs, aes(x = peripheral_blood_mononuclear_cell_uw, y = withdep_t, group = sex,
                           color = sex)) +
  scale_color_manual(labels = c("male", "female", 
                                "male linear fit", 
                                "female linear fit"), values = c("darkslategray3", "coral1", 
                                                                 "blue", "brown2"))+
  geom_point(alpha = 0.3) +
  geom_smooth(method=lm, aes(color = factor(sex_line_grp)))+
  theme_stata() + labs(title ="Withdrawn/depressed score vs. PBMC tPRS by sex") +
  ylab("sqrt with/dep") +
  geom_hline(yintercept=sqrt(65),  linetype="dashed")

male_biv <- lm (withdep_t ~ peripheral_blood_mononuclear_cell_uw, 
                data = base_white_tprs[which(base_white_tprs$sex == 1),])
fem_biv <- lm (withdep_t ~ peripheral_blood_mononuclear_cell_uw, 
               data = base_white_tprs[which(base_white_tprs$sex == 2),])
kable(tidy(male_biv)[, -c(4)])
kable(tidy(fem_biv)[, -c(4)])

ggplot(base_white_tprs, aes(x = double_negative_thymocyte_uw, y = anxdep_t, group = sex,
                            color = sex)) +
  scale_color_manual(labels = c("male", "female", 
                                "male linear fit", 
                                "female linear fit"), values = c("darkslategray3", "coral1", 
                                                                 "blue", "brown2"))+
  geom_point(alpha = 0.3) +
  geom_smooth(method=lm, aes(color = factor(sex_line_grp)))+
  theme_stata() + labs(title ="Anxious/depressed score vs. thymocyte tPRS by sex") +
  ylab("sqrt anx/dep") +
  geom_hline(yintercept=sqrt(65),  linetype="dashed")
male_biv <- lm (anxdep_t ~ naive_B_cell_uw, 
                data = base_white_tprs[which(base_white_tprs$sex == 1),])
fem_biv <- lm (anxdep_t ~ naive_B_cell_uw, 
               data = base_white_tprs[which(base_white_tprs$sex == 2),])
kable(tidy(male_biv)[, -c(4)])
kable(tidy(fem_biv)[, -c(4)])

ggplot(base_white_tprs, aes(x = double_negative_thymocyte_uw, y = withdep_t, group = sex,
                            color = sex)) +
  scale_color_manual(labels = c("male", "female", 
                                "male linear fit", 
                                "female linear fit"), values = c("darkslategray3", "coral1", 
                                                                 "blue", "brown2"))+
  geom_point(alpha = 0.3) +
  geom_smooth(method=lm, aes(color = factor(sex_line_grp)))+
  theme_stata() + labs(title ="Withdrawn/depressed score vs. thymocyte tPRS by sex") +
  ylab("sqrt with/dep") +
  geom_hline(yintercept=sqrt(65),  linetype="dashed")

ggplot(base_white_tprs, aes(x = naive_B_cell_uw, y = anxdep_t, group = sex,
                            color = sex)) +
  scale_color_manual(labels = c("male", "female", 
                                "male linear fit", 
                                "female linear fit"), values = c("darkslategray3", "coral1", 
                                                                 "blue", "brown2"))+
  geom_point(alpha = 0.3) +
  geom_smooth(method=lm, aes(color = factor(sex_line_grp)))+
  theme_stata() + labs(title ="Anxious/depressed score vs. naive B cell tPRS by sex") +
  ylab("sqrt anx/dep") +
  geom_hline(yintercept=sqrt(65),  linetype="dashed")

ggplot(base_white_tprs, aes(x = naive_B_cell_uw, y = withdep_t, group = sex,
                            color = sex)) +
  scale_color_manual(labels = c("male", "female", 
                                "male linear fit", 
                                "female linear fit"), values = c("darkslategray3", "coral1", 
                                                                 "blue", "brown2"))+
  geom_point(alpha = 0.3) +
  geom_smooth(method=lm, aes(color = factor(sex_line_grp)))+
  theme_stata() + labs(title ="Withdrawn/depressed score vs. naive B cell tPRS by sex") +
  ylab("sqrt with/dep") +
  geom_hline(yintercept=sqrt(65),  linetype="dashed")

male_biv <- lm (withdep_t ~ naive_B_cell_uw, 
                data = base_white_tprs[which(base_white_tprs$sex == 1),])
fem_biv <- lm (withdep_t ~ naive_B_cell_uw, 
               data = base_white_tprs[which(base_white_tprs$sex == 2),])
kable(tidy(male_biv)[, -c(4)])
kable(tidy(fem_biv)[, -c(4)])


# forest plots for sex ixn linear regressions================================

ad_slr_beta <- tprs_beta_sx[1,][, -c(1)]
wd_slr_beta <- tprs_beta_sx[2,][, -c(1)]
ad_mlr_beta <- tprs_beta_sx[3,][, -c(1)]
wd_mlr_beta <- tprs_beta_sx[4,][, -c(1)]
labs <- c(cell_types, cell_types)
weight_labs <- c(rep("unweighted", n_ct), rep("weighted", n_ct))
ad_slr_lower <- ad_slr_beta - tprs_CI_sx[1,][, -c(1)]
wd_slr_lower <- wd_slr_beta - tprs_CI_sx[2,][, -c(1)]
ad_mlr_lower <- ad_mlr_beta - tprs_CI_sx[3,][, -c(1)]
wd_mlr_lower <- wd_mlr_beta - tprs_CI_sx[4,][, -c(1)]
ad_slr_upper <- ad_slr_beta + tprs_CI_sx[1,][, -c(1)]
wd_slr_upper <- wd_slr_beta + tprs_CI_sx[2,][, -c(1)]
ad_mlr_upper <- ad_mlr_beta + tprs_CI_sx[3,][, -c(1)]
wd_mlr_upper <- wd_mlr_beta + tprs_CI_sx[4,][, -c(1)]

ad_slr_df <- data.frame(betas = as.numeric(unname(ad_slr_beta)), cell_types = unname(labs), 
                        upper = as.numeric(unname(ad_slr_upper)), 
                        lower = as.numeric(unname(ad_slr_lower)), tPRS_type = unname(weight_labs),
                        p_values = paste0("(", round(as.numeric(unname(pval_tprs_slope_sx[1,][, -c(1)]))[1:n_ct], 3),
                                          ", ", round(as.numeric(unname(pval_tprs_slope_sx[1,][, -c(1)]))[(n_ct + 1):(2*n_ct)], 3),
                                          ")"),
                        R2 = paste0("(", formatC(as.numeric(unname(R2_sx[1,][, -c(1)]))[1:n_ct], format = "e", digits = 3),
                                    ", ", formatC(as.numeric(unname(R2_sx[1,][, -c(1)]))[(n_ct + 1):(2*n_ct)], format = "e", digits = 3),
                                    ")"))

wd_slr_df <- data.frame(betas = as.numeric(unname(wd_slr_beta)), cell_types = unname(labs), 
                        upper = as.numeric(unname(wd_slr_upper)), 
                        lower = as.numeric(unname(wd_slr_lower)), tPRS_type = unname(weight_labs),
                        p_values = paste0("(", round(as.numeric(unname(pval_tprs_slope_sx[2,][, -c(1)]))[1:n_ct], 3),
                                          ", ", round(as.numeric(unname(pval_tprs_slope_sx[2,][, -c(1)]))[(n_ct + 1):(2*n_ct)], 3),
                                          ")"),
                        R2 = paste0("(", formatC(as.numeric(unname(R2_sx[2,][, -c(1)]))[1:n_ct], format = "e", digits = 3),
                                    ", ", formatC(as.numeric(unname(R2_sx[2,][, -c(1)]))[(n_ct + 1):(2*n_ct)], format = "e", digits = 3),
                                    ")"))

ad_mlr_df <- data.frame(betas = as.numeric(unname(ad_mlr_beta)), cell_types = unname(labs), 
                        upper = as.numeric(unname(ad_mlr_upper)), 
                        lower = as.numeric(unname(ad_mlr_lower)), tPRS_type = unname(weight_labs),
                        p_values = paste0("(", round(as.numeric(unname(pval_tprs_slope_sx[3,][, -c(1)]))[1:n_ct], 3),
                                          ", ", round(as.numeric(unname(pval_tprs_slope_sx[3,][, -c(1)]))[(n_ct + 1):(2*n_ct)], 3),
                                          ")"),
                        R2 = paste0("(", formatC(as.numeric(unname(R2_sx[3,][, -c(1)]))[1:n_ct], format = "e", digits = 3),
                                    ", ", formatC(as.numeric(unname(R2_sx[3,][, -c(1)]))[(n_ct + 1):(2*n_ct)], format = "e", digits = 3),
                                    ")"))

wd_mlr_df <- data.frame(betas = as.numeric(unname(wd_mlr_beta)), cell_types = unname(labs), 
                        upper = as.numeric(unname(wd_mlr_upper)), 
                        lower = as.numeric(unname(wd_mlr_lower)), tPRS_type = unname(weight_labs),
                        p_values = paste0("(", round(as.numeric(unname(pval_tprs_slope_sx[4,][, -c(1)]))[1:n_ct], 3),
                                          ", ", round(as.numeric(unname(pval_tprs_slope_sx[4,][, -c(1)]))[(n_ct + 1):(2*n_ct)], 3),
                                          ")"),
                        R2 = paste0("(", formatC(as.numeric(unname(R2_sx[4,][, -c(1)]))[1:n_ct], format = "e", digits = 3),
                                    ", ", formatC(as.numeric(unname(R2_sx[4,][, -c(1)]))[(n_ct + 1):(2*n_ct)], format = "e", digits = 3),
                                    ")"))

mid_plot <- ggplot(ad_slr_df, aes(x = betas, y = cell_types, xmin = upper, 
                                  xmax = lower, group = tPRS_type)) +
  geom_point(size = 4, aes(shape = tPRS_type), position = position_dodgev(height = .5)) +
  geom_pointrange(shape = 1, position = position_dodgev(height = .5)) +
  theme(legend.position = "bottom") + ggtitle("Effect sizes by cell type tPRS in \nanxdep ~ tPRS*sex") +
  xlab("effect sizes") + ylab("cell types")
mid_plot <- mid_plot +
  geom_vline(xintercept = 0, linetype="dashed")
mid_plot


p_right <- 
  ggplot(ad_slr_df) +
  geom_text(
    aes(x = 0, y = cell_types, label = R2),
    hjust = 0,
  ) + 
  theme_void() + ggtitle( paste0("R^2",": (unwtd, wtd)")) + theme(plot.title = element_text(hjust = 0.7))
p_right

layout <- c(
  area(t = 0, l = 0, b = 30, r = 9), # left plot, starts at the top of the page (0) and goes 30 units down and 3 units to the right
  area(t = 0, l = 5, b = 30, r = 15) # right most plot starts at top of page, begins where middle plot ends (l=9, and middle plot is r=9), goes to bottom of page (b=30), and extends two units wide (r=11)
)
# final plot arrangement
mid_plot + p_right + plot_layout(design = layout)

mid_plot <- ggplot(wd_slr_df, aes(x = betas, y = cell_types, xmin = upper, 
                                  xmax = lower, group = tPRS_type)) +
  geom_point(size = 4, aes(shape = tPRS_type), position = position_dodgev(height = .5)) +
  geom_pointrange(shape = 1, position = position_dodgev(height = .5)) +
  theme(legend.position = "bottom") + ggtitle("Effect sizes by cell type tPRS in \nwithdep ~ tPRS*sex") +
  xlab("effect sizes") + ylab("cell types")
mid_plot <- mid_plot +
  geom_vline(xintercept = 0, linetype="dashed")
mid_plot


p_right <- 
  ggplot(wd_slr_df) +
  geom_text(
    aes(x = 0, y = (cell_types), label = R2),
    hjust = 0,
  ) + 
  theme_void() + ggtitle( paste0("R^2",": (unwtd, wtd)")) + theme(plot.title = element_text(hjust = 0.7))
p_right

# final plot arrangement
mid_plot + p_right + plot_layout(design = layout)

mid_plot <- ggplot(ad_mlr_df, aes(x = betas, y = cell_types, xmin = upper, 
                                  xmax = lower, group = tPRS_type)) +
  geom_point(size = 4, aes(shape = tPRS_type), position = position_dodgev(height = .5)) +
  geom_pointrange(shape = 1, position = position_dodgev(height = .5)) +
  theme(legend.position = "bottom") + ggtitle("Effect sizes by cell type tPRS in \nanxdep ~ tPRS*sex + covariates") +
  xlab("effect sizes") + ylab("cell types")
mid_plot <- mid_plot +
  geom_vline(xintercept = 0, linetype="dashed")
mid_plot


p_right <- 
  ggplot(ad_mlr_df) +
  geom_text(
    aes(x = 0, y = (cell_types), label = R2),
    hjust = 0,
  ) + 
  theme_void() + ggtitle( paste0("R^2",": (unwtd, wtd)")) + theme(plot.title = element_text(hjust = 0.7))
p_right

# final plot arrangement
mid_plot + p_right + plot_layout(design = layout)

mid_plot <- ggplot(wd_mlr_df, aes(x = betas, y = cell_types, xmin = upper, 
                                  xmax = lower, group = tPRS_type)) +
  geom_point(size = 4, aes(shape = tPRS_type), position = position_dodgev(height = .5)) +
  geom_pointrange(shape = 1, position = position_dodgev(height = .5)) +
  theme(legend.position = "bottom") + ggtitle("Effect sizes by cell type tPRS in \nwithdep ~ tPRS*sex + covariates") +
  xlab("effect sizes") + ylab("cell types")
mid_plot <- mid_plot +
  geom_vline(xintercept = 0, linetype="dashed")
mid_plot


p_right <- 
  ggplot(wd_mlr_df) +
  geom_text(
    aes(x = 0, y = (cell_types), label = R2),
    hjust = 0,
  ) + 
  theme_void() + ggtitle( paste0("R^2",": (unwtd, wtd)")) + theme(plot.title = element_text(hjust = 0.7))
p_right

# final plot arrangement
mid_plot + p_right + plot_layout(design = layout)

save(tprs_beta, file="sc-tprs-mdd/sc_tPRS_betas_immune.Rda")
