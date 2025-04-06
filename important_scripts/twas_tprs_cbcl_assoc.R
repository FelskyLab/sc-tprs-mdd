# required: twas_analysis_abcd.R
# goal: lin, log, gamma regs for TWAS tPRSs and cbcl scores
# - prelim correlations with anx/dep, with/dep
# next: TODO interpret forrest plots, id important cell types 

# lin reg=======================
adjR2 <- data.frame(regs = c("mlr_ad0", "mlr_ad_m","mlr_ad_f", "mlr_wd0", "mlr_wd_m", "mlr_wd_f"))
pval_mod <- data.frame(regs = c("mlr_ad0", "mlr_ad_m","mlr_ad_f", "mlr_wd0", "mlr_wd_m", "mlr_wd_f"))
pval_tprs_slope <- data.frame(regs = c("mlr_ad0", "mlr_ad_m","mlr_ad_f", "mlr_wd0", "mlr_wd_m", "mlr_wd_f"))
tprs_beta <- data.frame(regs = c("mlr_ad0", "mlr_ad_m","mlr_ad_f", "mlr_wd0", "mlr_wd_m", "mlr_wd_f"))
tprs_CI <- data.frame(regs = c("mlr_ad0", "mlr_ad_m","mlr_ad_f", "mlr_wd0", "mlr_wd_m", "mlr_wd_f"))
R2 <- data.frame(regs = c("mlr_ad0", "mlr_ad_m","mlr_ad_f", "mlr_wd0", "mlr_wd_m", "mlr_wd_f"))

twas_ad_tPRS_w <- as.data.frame(cbind(twas_ad_tPRS_w[,1], scale(twas_ad_tPRS_w[,-1])))
twas_wd_tPRS_w <- as.data.frame(cbind(twas_wd_tPRS_w[,1], scale(twas_wd_tPRS_w[,-1])))

for (i in 2:23){
  for (j in 1:1){
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]] <- 
      lm(anx_sqrt ~ twas_ad_tPRS_w[[i]] + age + sex + site + 
           genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
           genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
           genetic_pc_9 + genetic_pc_10, data = base_bin)
    adjR2[[colnames(twas_ad_tPRS_w)[i]]][1] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]])$adj.r.squared
    R2[[colnames(twas_ad_tPRS_w)[i]]][1] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]])$r.squared
    pval_mod[[colnames(twas_ad_tPRS_w)[i]]][1] <- 
      overall_p(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]])
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][1] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][1] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][1] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]]))[2])
    
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]] <- 
      lm(anx_sqrt ~ twas_ad_tPRS_w[which(base_bin$sex == "1"), i+22] + age + site + 
           genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
           genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
           genetic_pc_9 + genetic_pc_10, data = base_bin[which(base_bin$sex == "1"),])
    adjR2[[colnames(twas_ad_tPRS_w)[i]]][2] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]])$adj.r.squared
    R2[[colnames(twas_ad_tPRS_w)[i]]][2] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]])$r.squared
    pval_mod[[colnames(twas_ad_tPRS_w)[i]]][2] <- 
      overall_p(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]])
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][2] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][2] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][2] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]]))[2])
    
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]] <-
      lm(anx_sqrt ~ twas_ad_tPRS_w[which(base_bin$sex == "2"), i+44] + age + site + 
           genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
           genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
           genetic_pc_9 + genetic_pc_10, data = base_bin[which(base_bin$sex == "2"),])
    adjR2[[colnames(twas_ad_tPRS_w)[i]]][3] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]])$adj.r.squared
    R2[[colnames(twas_ad_tPRS_w)[i]]][3] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]])$r.squared
    pval_mod[[colnames(twas_ad_tPRS_w)[i]]][3] <- 
      overall_p(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]])
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][3] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][3] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][3] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]]))[2])
    
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]] <- 
      lm(with_sqrt ~ twas_wd_tPRS_w[[i]] + age + sex + site + 
           genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
           genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
           genetic_pc_9 + genetic_pc_10, data = base_bin)
    adjR2[[colnames(twas_ad_tPRS_w)[i]]][4] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]])$adj.r.squared
    R2[[colnames(twas_ad_tPRS_w)[i]]][4] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]])$r.squared
    pval_mod[[colnames(twas_ad_tPRS_w)[i]]][4] <- 
      overall_p(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]])
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][4] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][4] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][4] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]]))[2])
    
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]] <- 
      lm(with_sqrt ~ twas_wd_tPRS_w[which(base_bin$sex == "1"), i+22] + age + site + 
           genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
           genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
           genetic_pc_9 + genetic_pc_10, data = base_bin[which(base_bin$sex == "1"), ])
    adjR2[[colnames(twas_ad_tPRS_w)[i]]][5] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]])$adj.r.squared
    R2[[colnames(twas_ad_tPRS_w)[i]]][5] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]])$r.squared
    pval_mod[[colnames(twas_ad_tPRS_w)[i]]][5] <- 
      overall_p(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]])
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][5] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][5] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][5] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]]))[2])
    
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]] <- 
      lm(with_sqrt ~ twas_wd_tPRS_w[which(base_bin$sex == "2"), i+44] + age + site + 
           genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
           genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
           genetic_pc_9 + genetic_pc_10, data = base_bin[which(base_bin$sex == "2"), ])
    adjR2[[colnames(twas_ad_tPRS_w)[i]]][6] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]])$adj.r.squared
    R2[[colnames(twas_ad_tPRS_w)[i]]][6] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]])$r.squared
    pval_mod[[colnames(twas_ad_tPRS_w)[i]]][6] <- 
      overall_p(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]])
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][6] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][6] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][6] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]]))[2])
  }
  
}

# forrest plots for lin regs=========================
ad_beta <- unlist(tprs_beta[1:3,-1])
wd_beta <- unlist(tprs_beta[4:6,-1])
labs <- rep(cell_types, each=3)
fit_labs <- rep(c("mixed","male","female"), n_ct)
ad_lower <- ad_beta - unlist(tprs_CI[1:3,-1])
wd_lower <- wd_beta - unlist(tprs_CI[4:6,-1])
ad_upper <- ad_beta + unlist(tprs_CI[1:3,-1])
wd_upper <- wd_beta + unlist(tprs_CI[4:6,-1])

ad_df <- data.frame(betas = as.numeric(unname(ad_beta)), cell_types = unname(labs), 
                    upper = as.numeric(unname(ad_upper)), 
                    lower = as.numeric(unname(ad_lower)), tPRS_type = unname(fit_labs),
                    p_values = rep(paste0("(", round(as.numeric(unname(pval_tprs_slope[1,][, -c(1)]))[1:n_ct], 3),
                                      ", ", round(as.numeric(unname(pval_tprs_slope[2,][, -c(1)]))[1:n_ct], 3),
                                      ",", round(as.numeric(unname(pval_tprs_slope[3,][, -c(1)]))[1:n_ct], 3), ")"), each=3),
                    R2 = rep(paste0("(", percent(round(as.numeric(unname(R2[1,][, -c(1)]))[1:n_ct], 4)),
                                ", ", percent(round(as.numeric(unname(R2[2,][, -c(1)]))[1:n_ct], 4)),
                                ",", percent(round(as.numeric(unname(R2[3,][, -c(1)]))[1:n_ct], 4)), ")"), each=3))

wd_df <- data.frame(betas = as.numeric(unname(wd_beta)), cell_types = unname(labs), 
                    upper = as.numeric(unname(wd_upper)), 
                    lower = as.numeric(unname(wd_lower)), tPRS_type = unname(fit_labs),
                    p_values = rep(paste0("(", round(as.numeric(unname(pval_tprs_slope[4,][, -c(1)]))[1:n_ct], 3),
                                      ", ", round(as.numeric(unname(pval_tprs_slope[5,][, -c(1)]))[1:n_ct], 3),
                                      ",", round(as.numeric(unname(pval_tprs_slope[6,][, -c(1)]))[1:n_ct], 3),")"), each=3),
                    R2 = rep(paste0("(", percent(round(as.numeric(unname(R2[4,][, -c(1)]))[1:n_ct], 4)),
                                ", ", percent(round(as.numeric(unname(R2[5,][, -c(1)]))[1:n_ct], 4)),
                                ",", percent(round(as.numeric(unname(R2[6,][, -c(1)]))[1:n_ct], 4)),")"), each=3))

mid_plot <- ggplot(ad_df, aes(x = betas, y = cell_types, xmin = upper, 
                              xmax = lower, group = tPRS_type, color = tPRS_type)) +
  geom_point(size = 4, aes(shape = tPRS_type), position = position_dodgev(height = .5)) +
  geom_pointrange(shape = 1, position = position_dodgev(height = .5)) +
  theme(legend.position = "bottom", axis.text.y = element_text(angle = 0, size = 11)) + 
  ggtitle("Main effect size of anx/dep tPRSs by cell type\nwith anxious/depressed score as outcome") +
  xlab("tPRS main effect sizes") + ylab("cell types")
mid_plot <- mid_plot +
  geom_vline(xintercept = 0, linetype="dashed")+ theme_hc()

p_right <- 
  ggplot(ad_df) +
  geom_text(
    aes(x = 0, y = (cell_types), label = R2),
    hjust = 0,
  ) + 
  theme_void() + ggtitle( paste0("% variance explained",":\n(mixed, male, female)")) + theme(plot.title = element_text(hjust = 0.8))

# final plot arrangement
mid_plot + p_right + plot_layout(design = layout)

mid_plot <- ggplot(wd_df, aes(x = betas, y = cell_types, xmin = upper, 
                              xmax = lower, group = tPRS_type, color = tPRS_type)) +
  geom_point(size = 4, aes(shape = tPRS_type), position = position_dodgev(height = .5)) +
  geom_pointrange(shape = 1, position = position_dodgev(height = .5)) +
  theme(legend.position = "bottom", axis.text.y = element_text(angle = 0, size = 11)) + 
  ggtitle("Main effect size of with/dep tPRSs by cell type\nwith withdrawn/depressed score as outcome") +
  xlab("tPRS main effect sizes") + ylab("cell types")
mid_plot <- mid_plot +
  geom_vline(xintercept = 0, linetype="dashed")+ theme_hc()

p_right <- 
  ggplot(wd_df) +
  geom_text(
    aes(x = 0, y = (cell_types), label = R2),
    hjust = 0,
  ) + 
  theme_void() + ggtitle( paste0("% variance explained",":\n(mixed, male, female)")) + theme(plot.title = element_text(hjust = 0.8))

# final plot arrangement
mid_plot + p_right + plot_layout(design = layout)

# train lin regs==========================
for (i in 2:23){
  for (j in 1:1){
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]] <- 
      lm(anxdep_t ~ twas_ad_tPRS_w[which(inds_ad_all=="train"), i] + age + sex + site + 
           genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
           genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
           genetic_pc_9 + genetic_pc_10, data = base_bin[which(inds_ad_all=="train"),])
    adjR2[[colnames(twas_ad_tPRS_w)[i]]][1] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]])$adj.r.squared
    R2[[colnames(twas_ad_tPRS_w)[i]]][1] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]])$r.squared
    pval_mod[[colnames(twas_ad_tPRS_w)[i]]][1] <- 
      overall_p(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]])
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][1] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][1] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][1] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]]))[2])
    
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]] <- 
      lm(anxdep_t ~ twas_ad_tPRS_w[which(base_bin$sex == "1"), i+22][which(inds_ad_m=="train")] + age + site + 
           genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
           genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
           genetic_pc_9 + genetic_pc_10, data = base_bin[which(base_bin$sex == "1"),][which(inds_ad_m=="train"),])
    adjR2[[colnames(twas_ad_tPRS_w)[i]]][2] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]])$adj.r.squared
    R2[[colnames(twas_ad_tPRS_w)[i]]][2] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]])$r.squared
    pval_mod[[colnames(twas_ad_tPRS_w)[i]]][2] <- 
      overall_p(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]])
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][2] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][2] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][2] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]]))[2])
    
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]] <-
      lm(anxdep_t ~ twas_ad_tPRS_w[which(base_bin$sex == "2"), i+44][which(inds_ad_f=="train")] + age + site + 
           genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
           genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
           genetic_pc_9 + genetic_pc_10, data = base_bin[which(base_bin$sex == "2"),][which(inds_ad_f=="train"),])
    adjR2[[colnames(twas_ad_tPRS_w)[i]]][3] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]])$adj.r.squared
    R2[[colnames(twas_ad_tPRS_w)[i]]][3] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]])$r.squared
    pval_mod[[colnames(twas_ad_tPRS_w)[i]]][3] <- 
      overall_p(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]])
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][3] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][3] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][3] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]]))[2])
    
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]] <- 
      lm(withdep_t ~ twas_wd_tPRS_w[which(inds_wd_all=="train"), i] + age + sex + site + 
           genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
           genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
           genetic_pc_9 + genetic_pc_10, data = white_all_tprs[which(inds_wd_all=="train"),])
    adjR2[[colnames(twas_ad_tPRS_w)[i]]][4] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]])$adj.r.squared
    R2[[colnames(twas_ad_tPRS_w)[i]]][4] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]])$r.squared
    pval_mod[[colnames(twas_ad_tPRS_w)[i]]][4] <- 
      overall_p(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]])
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][4] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][4] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][4] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]]))[2])
    
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]] <- 
      lm(withdep_t ~ twas_wd_tPRS_w[which(base_bin$sex == "1"), i+22][which(inds_wd_m=="train")] + age + site + 
           genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
           genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
           genetic_pc_9 + genetic_pc_10, data = base_bin[which(base_bin$sex == "1"), ][which(inds_wd_m=="train"),])
    adjR2[[colnames(twas_ad_tPRS_w)[i]]][5] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]])$adj.r.squared
    R2[[colnames(twas_ad_tPRS_w)[i]]][5] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]])$r.squared
    pval_mod[[colnames(twas_ad_tPRS_w)[i]]][5] <- 
      overall_p(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]])
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][5] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][5] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][5] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]]))[2])
    
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]] <- 
      lm(withdep_t ~ twas_wd_tPRS_w[which(base_bin$sex == "2"), i+44][which(inds_wd_f=="train")] + age + site + 
           genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
           genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
           genetic_pc_9 + genetic_pc_10, data = base_bin[which(base_bin$sex == "2"), ][which(inds_wd_f=="train"),])
    adjR2[[colnames(twas_ad_tPRS_w)[i]]][6] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]])$adj.r.squared
    R2[[colnames(twas_ad_tPRS_w)[i]]][6] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]])$r.squared
    pval_mod[[colnames(twas_ad_tPRS_w)[i]]][6] <- 
      overall_p(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]])
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][6] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][6] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][6] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]]))[2])
  }
  
}
ad_beta <- unlist(tprs_beta[1:3,-1])
wd_beta <- unlist(tprs_beta[4:6,-1])
labs <- rep(cell_types, each=3)
fit_labs <- rep(c("mixed","male","female"), n_ct)
ad_lower <- ad_beta - unlist(tprs_CI[1:3,-1])
wd_lower <- wd_beta - unlist(tprs_CI[4:6,-1])
ad_upper <- ad_beta + unlist(tprs_CI[1:3,-1])
wd_upper <- wd_beta + unlist(tprs_CI[4:6,-1])

ad_df <- data.frame(betas = as.numeric(unname(ad_beta)), cell_types = unname(labs), 
                    upper = as.numeric(unname(ad_upper)), 
                    lower = as.numeric(unname(ad_lower)), tPRS_type = unname(fit_labs),
                    p_values = rep(paste0("(", round(as.numeric(unname(pval_tprs_slope[1,][, -c(1)]))[1:n_ct], 3),
                                          ", ", round(as.numeric(unname(pval_tprs_slope[2,][, -c(1)]))[1:n_ct], 3),
                                          ",", round(as.numeric(unname(pval_tprs_slope[3,][, -c(1)]))[1:n_ct], 3), ")"), each=3),
                    R2 = rep(paste0("(", percent(round(as.numeric(unname(R2[1,][, -c(1)]))[1:n_ct], 4)),
                                    ", ", percent(round(as.numeric(unname(R2[2,][, -c(1)]))[1:n_ct], 4)),
                                    ",", percent(round(as.numeric(unname(R2[3,][, -c(1)]))[1:n_ct], 4)), ")"), each=3))

wd_df <- data.frame(betas = as.numeric(unname(wd_beta)), cell_types = unname(labs), 
                    upper = as.numeric(unname(wd_upper)), 
                    lower = as.numeric(unname(wd_lower)), tPRS_type = unname(fit_labs),
                    p_values = rep(paste0("(", round(as.numeric(unname(pval_tprs_slope[4,][, -c(1)]))[1:n_ct], 3),
                                          ", ", round(as.numeric(unname(pval_tprs_slope[5,][, -c(1)]))[1:n_ct], 3),
                                          ",", round(as.numeric(unname(pval_tprs_slope[6,][, -c(1)]))[1:n_ct], 3),")"), each=3),
                    R2 = rep(paste0("(", percent(round(as.numeric(unname(R2[4,][, -c(1)]))[1:n_ct], 4)),
                                    ", ", percent(round(as.numeric(unname(R2[5,][, -c(1)]))[1:n_ct], 4)),
                                    ",", percent(round(as.numeric(unname(R2[6,][, -c(1)]))[1:n_ct], 4)),")"), each=3))

mid_plot <- ggplot(ad_df, aes(x = betas, y = cell_types, xmin = upper, 
                              xmax = lower, group = tPRS_type, color = tPRS_type)) +
  geom_point(size = 4, aes(shape = tPRS_type), position = position_dodgev(height = .5)) +
  geom_pointrange(shape = 1, position = position_dodgev(height = .5)) +
  theme(legend.position = "bottom", axis.text.y = element_text(angle = 0, size = 11)) + 
  ggtitle("Main effect size of anx/dep tPRSs by cell type\nwith anxious/depressed score as outcome\nin training data") +
  xlab("tPRS main effect sizes") + ylab("cell types")
mid_plot <- mid_plot +
  geom_vline(xintercept = 0, linetype="dashed")+ theme_hc()

p_right <- 
  ggplot(ad_df) +
  geom_text(
    aes(x = 0, y = (cell_types), label = R2),
    hjust = 0,
  ) + 
  theme_void() + ggtitle( paste0("% variance explained",":\n(mixed, male, female)")) + theme(plot.title = element_text(hjust = 0.8))

# final plot arrangement
mid_plot + p_right + plot_layout(design = layout)

mid_plot <- ggplot(wd_df, aes(x = betas, y = cell_types, xmin = upper, 
                              xmax = lower, group = tPRS_type, color = tPRS_type)) +
  geom_point(size = 4, aes(shape = tPRS_type), position = position_dodgev(height = .5)) +
  geom_pointrange(shape = 1, position = position_dodgev(height = .5)) +
  theme(legend.position = "bottom", axis.text.y = element_text(angle = 0, size = 11)) + 
  ggtitle("Main effect size of with/dep tPRSs by cell type\nwith withdrawn/depressed score as outcome\nin training data") +
  xlab("tPRS main effect sizes") + ylab("cell types")
mid_plot <- mid_plot +
  geom_vline(xintercept = 0, linetype="dashed")+ theme_hc()

p_right <- 
  ggplot(wd_df) +
  geom_text(
    aes(x = 0, y = (cell_types), label = R2),
    hjust = 0,
  ) + 
  theme_void() + ggtitle( paste0("% variance explained",":\n(mixed, male, female)")) + theme(plot.title = element_text(hjust = 0.8))

# final plot arrangement
mid_plot + p_right + plot_layout(design = layout)

# test lin regs==================================
for (i in 2:23){
  for (j in 1:1){
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]] <- 
      lm(anxdep_t ~ twas_ad_tPRS_w[which(inds_ad_all=="test"), i] + age + sex + site + 
           genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
           genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
           genetic_pc_9 + genetic_pc_10, data = base_bin[which(inds_ad_all=="test"),])
    adjR2[[colnames(twas_ad_tPRS_w)[i]]][1] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]])$adj.r.squared
    R2[[colnames(twas_ad_tPRS_w)[i]]][1] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]])$r.squared
    pval_mod[[colnames(twas_ad_tPRS_w)[i]]][1] <- 
      overall_p(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]])
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][1] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][1] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][1] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]]))[2])
    
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]] <- 
      lm(anxdep_t ~ twas_ad_tPRS_w[which(base_bin$sex == "1"), i+22][which(inds_ad_m=="test")] + age + site + 
           genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
           genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
           genetic_pc_9 + genetic_pc_10, data = base_bin[which(base_bin$sex == "1"),][which(inds_ad_m=="test"),])
    adjR2[[colnames(twas_ad_tPRS_w)[i]]][2] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]])$adj.r.squared
    R2[[colnames(twas_ad_tPRS_w)[i]]][2] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]])$r.squared
    pval_mod[[colnames(twas_ad_tPRS_w)[i]]][2] <- 
      overall_p(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]])
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][2] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][2] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][2] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]]))[2])
    
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]] <-
      lm(anxdep_t ~ twas_ad_tPRS_w[which(base_bin$sex == "2"), i+44][which(inds_ad_f=="test")] + age + site + 
           genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
           genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
           genetic_pc_9 + genetic_pc_10, data = base_bin[which(base_bin$sex == "2"),][which(inds_ad_f=="test"),])
    adjR2[[colnames(twas_ad_tPRS_w)[i]]][3] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]])$adj.r.squared
    R2[[colnames(twas_ad_tPRS_w)[i]]][3] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]])$r.squared
    pval_mod[[colnames(twas_ad_tPRS_w)[i]]][3] <- 
      overall_p(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]])
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][3] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][3] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][3] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]]))[2])
    
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]] <- 
      lm(withdep_t ~ twas_wd_tPRS_w[which(inds_wd_all=="test"), i] + age + sex + site + 
           genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
           genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
           genetic_pc_9 + genetic_pc_10, data = base_bin[which(inds_wd_all=="test"),])
    adjR2[[colnames(twas_ad_tPRS_w)[i]]][4] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]])$adj.r.squared
    R2[[colnames(twas_ad_tPRS_w)[i]]][4] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]])$r.squared
    pval_mod[[colnames(twas_ad_tPRS_w)[i]]][4] <- 
      overall_p(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]])
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][4] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][4] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][4] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]]))[2])
    
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]] <- 
      lm(withdep_t ~ twas_wd_tPRS_w[which(base_bin$sex == "1"), i+22][which(inds_wd_m=="test")] + age + site + 
           genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
           genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
           genetic_pc_9 + genetic_pc_10, data = base_bin[which(base_bin$sex == "1"), ][which(inds_wd_m=="test"),])
    adjR2[[colnames(twas_ad_tPRS_w)[i]]][5] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]])$adj.r.squared
    R2[[colnames(twas_ad_tPRS_w)[i]]][5] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]])$r.squared
    pval_mod[[colnames(twas_ad_tPRS_w)[i]]][5] <- 
      overall_p(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]])
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][5] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][5] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][5] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]]))[2])
    
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]] <- 
      lm(withdep_t ~ twas_wd_tPRS_w[which(base_bin$sex == "2"), i+44][which(inds_wd_f=="test")] + age + site + 
           genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
           genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
           genetic_pc_9 + genetic_pc_10, data = base_bin[which(base_bin$sex == "2"), ][which(inds_wd_f=="test"),])
    adjR2[[colnames(twas_ad_tPRS_w)[i]]][6] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]])$adj.r.squared
    R2[[colnames(twas_ad_tPRS_w)[i]]][6] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]])$r.squared
    pval_mod[[colnames(twas_ad_tPRS_w)[i]]][6] <- 
      overall_p(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]])
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][6] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][6] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][6] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]]))[2])
  }
  
}
ad_beta <- unlist(tprs_beta[1:3,-1])
wd_beta <- unlist(tprs_beta[4:6,-1])
labs <- rep(cell_types, each=3)
fit_labs <- rep(c("mixed","male","female"), n_ct)
ad_lower <- ad_beta - unlist(tprs_CI[1:3,-1])
wd_lower <- wd_beta - unlist(tprs_CI[4:6,-1])
ad_upper <- ad_beta + unlist(tprs_CI[1:3,-1])
wd_upper <- wd_beta + unlist(tprs_CI[4:6,-1])

ad_df <- data.frame(betas = as.numeric(unname(ad_beta)), cell_types = unname(labs), 
                    upper = as.numeric(unname(ad_upper)), 
                    lower = as.numeric(unname(ad_lower)), tPRS_type = unname(fit_labs),
                    p_values = rep(paste0("(", round(as.numeric(unname(pval_tprs_slope[1,][, -c(1)]))[1:n_ct], 3),
                                          ", ", round(as.numeric(unname(pval_tprs_slope[2,][, -c(1)]))[1:n_ct], 3),
                                          ",", round(as.numeric(unname(pval_tprs_slope[3,][, -c(1)]))[1:n_ct], 3), ")"), each=3),
                    R2 = rep(paste0("(", percent(round(as.numeric(unname(R2[1,][, -c(1)]))[1:n_ct], 4)),
                                    ", ", percent(round(as.numeric(unname(R2[2,][, -c(1)]))[1:n_ct], 4)),
                                    ",", percent(round(as.numeric(unname(R2[3,][, -c(1)]))[1:n_ct], 4)), ")"), each=3))

wd_df <- data.frame(betas = as.numeric(unname(wd_beta)), cell_types = unname(labs), 
                    upper = as.numeric(unname(wd_upper)), 
                    lower = as.numeric(unname(wd_lower)), tPRS_type = unname(fit_labs),
                    p_values = rep(paste0("(", round(as.numeric(unname(pval_tprs_slope[4,][, -c(1)]))[1:n_ct], 3),
                                          ", ", round(as.numeric(unname(pval_tprs_slope[5,][, -c(1)]))[1:n_ct], 3),
                                          ",", round(as.numeric(unname(pval_tprs_slope[6,][, -c(1)]))[1:n_ct], 3),")"), each=3),
                    R2 = rep(paste0("(", percent(round(as.numeric(unname(R2[4,][, -c(1)]))[1:n_ct], 4)),
                                    ", ", percent(round(as.numeric(unname(R2[5,][, -c(1)]))[1:n_ct], 4)),
                                    ",", percent(round(as.numeric(unname(R2[6,][, -c(1)]))[1:n_ct], 4)),")"), each=3))

mid_plot <- ggplot(ad_df, aes(x = betas, y = cell_types, xmin = upper, 
                              xmax = lower, group = tPRS_type, color = tPRS_type)) +
  geom_point(size = 4, aes(shape = tPRS_type), position = position_dodgev(height = .5)) +
  geom_pointrange(shape = 1, position = position_dodgev(height = .5)) +
  theme(legend.position = "bottom", axis.text.y = element_text(angle = 0, size = 11)) + 
  ggtitle("Main effect size of anx/dep tPRSs by cell type\nwith anxious/depressed score as outcome\nin test data") +
  xlab("tPRS main effect sizes") + ylab("cell types")
mid_plot <- mid_plot +
  geom_vline(xintercept = 0, linetype="dashed")+ theme_hc()

p_right <- 
  ggplot(ad_df) +
  geom_text(
    aes(x = 0, y = (cell_types), label = R2),
    hjust = 0,
  ) + 
  theme_void() + ggtitle( paste0("% variance explained",":\n(mixed, male, female)")) + theme(plot.title = element_text(hjust = 0.8))

# final plot arrangement
mid_plot + p_right + plot_layout(design = layout)

mid_plot <- ggplot(wd_df, aes(x = betas, y = cell_types, xmin = upper, 
                              xmax = lower, group = tPRS_type, color = tPRS_type)) +
  geom_point(size = 4, aes(shape = tPRS_type), position = position_dodgev(height = .5)) +
  geom_pointrange(shape = 1, position = position_dodgev(height = .5)) +
  theme(legend.position = "bottom", axis.text.y = element_text(angle = 0, size = 11)) + 
  ggtitle("Main effect size of with/dep tPRSs by cell type\nwith withdrawn/depressed score as outcome\nin test data") +
  xlab("tPRS main effect sizes") + ylab("cell types")
mid_plot <- mid_plot +
  geom_vline(xintercept = 0, linetype="dashed")+ theme_hc()

p_right <- 
  ggplot(wd_df) +
  geom_text(
    aes(x = 0, y = (cell_types), label = R2),
    hjust = 0,
  ) + 
  theme_void() + ggtitle( paste0("% variance explained",":\n(mixed, male, female)")) + theme(plot.title = element_text(hjust = 0.8))

# final plot arrangement
mid_plot + p_right + plot_layout(design = layout)


# gamma regs ========================================
pval_tprs_slope <- data.frame(regs = c("mlr_ad0", "mlr_ad_m","mlr_ad_f", "mlr_wd0", "mlr_wd_m", "mlr_wd_f"))
tprs_beta <- data.frame(regs = c("mlr_ad0", "mlr_ad_m","mlr_ad_f", "mlr_wd0", "mlr_wd_m", "mlr_wd_f"))
tprs_CI <- data.frame(regs = c("mlr_ad0", "mlr_ad_m","mlr_ad_f", "mlr_wd0", "mlr_wd_m", "mlr_wd_f"))
R2 <- data.frame(regs = c("mlr_ad0", "mlr_ad_m","mlr_ad_f", "mlr_wd0", "mlr_wd_m", "mlr_wd_f"))

for (i in 2:23){
  for (j in 1:1){
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]] <- 
      glm(anx_sqrt ~ twas_ad_tPRS_w[[i]] + age + sex + site + 
           genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
           genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
           genetic_pc_9 + genetic_pc_10, data = base_bin, family=Gamma(link='log'))
    R2[[colnames(twas_ad_tPRS_w)[i]]][1] <- 
      1 - summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]])$deviance/summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]])$null.deviance
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][1] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][1] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][1] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]]))[2])
    
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]] <- 
      glm(anx_sqrt ~ twas_ad_tPRS_w[which(base_bin$sex == "1"), i+22] + age + site + 
           genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
           genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
           genetic_pc_9 + genetic_pc_10, data = base_bin[which(base_bin$sex == "1"),],
           family=Gamma(link='log'))
    R2[[colnames(twas_ad_tPRS_w)[i]]][2] <- 
      1 - summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]])$deviance/summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]])$null.deviance
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][2] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][2] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][2] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]]))[2])
    
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]] <-
      glm(anx_sqrt ~ twas_ad_tPRS_w[which(base_bin$sex == "2"), i+44] + age + site + 
           genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
           genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
           genetic_pc_9 + genetic_pc_10, data = base_bin[which(base_bin$sex == "2"),],
           family=Gamma(link='log'))
    R2[[colnames(twas_ad_tPRS_w)[i]]][3] <- 
      1 - summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]])$deviance/summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]])$null.deviance
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][3] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][3] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][3] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]]))[2])
    
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]] <- 
      glm(with_sqrt ~ twas_wd_tPRS_w[[i]] + age + sex + site + 
           genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
           genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
           genetic_pc_9 + genetic_pc_10, data = base_bin, family=Gamma(link='log'))
    R2[[colnames(twas_ad_tPRS_w)[i]]][4] <- 
      1 - summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]])$deviance/summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]])$null.deviance
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][4] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][4] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][4] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]]))[2])
    
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]] <- 
      glm(with_sqrt ~ twas_wd_tPRS_w[which(base_bin$sex == "1"), i+22] + age + site + 
           genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
           genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
           genetic_pc_9 + genetic_pc_10, data = base_bin[which(base_bin$sex == "1"), ],
           family=Gamma(link='log'))
    R2[[colnames(twas_ad_tPRS_w)[i]]][5] <- 
      1 - summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]])$deviance/summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]])$null.deviance
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][5] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][5] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][5] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]]))[2])
    
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]] <- 
      glm(with_sqrt ~ twas_wd_tPRS_w[which(base_bin$sex == "2"), i+44] + age + site + 
           genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
           genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
           genetic_pc_9 + genetic_pc_10, data = base_bin[which(base_bin$sex == "2"), ],
           family=Gamma(link='log'))
    R2[[colnames(twas_ad_tPRS_w)[i]]][6] <- 
      1 - summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]])$deviance/summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]])$null.deviance
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][6] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][6] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][6] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]]))[2])
  }
  
}
ad_beta <- exp(unlist(tprs_beta[1:3,-1]))
wd_beta <- exp(unlist(tprs_beta[4:6,-1]))
labs <- rep(cell_types, each=3)
fit_labs <- rep(c("mixed","male","female"), n_ct)
ad_lower <- exp(unlist(tprs_beta[1:3,-1]) - unlist(tprs_CI[1:3,-1]))
wd_lower <- exp(unlist(tprs_beta[4:6,-1]) - unlist(tprs_CI[4:6,-1]))
ad_upper <- exp(unlist(tprs_beta[1:3,-1]) + unlist(tprs_CI[1:3,-1]))
wd_upper <- exp(unlist(tprs_beta[4:6,-1]) + unlist(tprs_CI[4:6,-1]))

ad_df <- data.frame(betas = as.numeric(unname(ad_beta)), cell_types = unname(labs), 
                    upper = as.numeric(unname(ad_upper)), 
                    lower = as.numeric(unname(ad_lower)), tPRS_type = unname(fit_labs),
                    p_values = rep(paste0("(", round(as.numeric(unname(pval_tprs_slope[1,][, -c(1)]))[1:n_ct], 3),
                                          ", ", round(as.numeric(unname(pval_tprs_slope[2,][, -c(1)]))[1:n_ct], 3),
                                          ",", round(as.numeric(unname(pval_tprs_slope[3,][, -c(1)]))[1:n_ct], 3), ")"), each=3),
                    R2 = rep(paste0("(", percent(round(as.numeric(unname(R2[1,][, -c(1)]))[1:n_ct], 4)),
                                    ", ", percent(round(as.numeric(unname(R2[2,][, -c(1)]))[1:n_ct], 4)),
                                    ",", percent(round(as.numeric(unname(R2[3,][, -c(1)]))[1:n_ct], 4)), ")"), each=3))

wd_df <- data.frame(betas = as.numeric(unname(wd_beta)), cell_types = unname(labs), 
                    upper = as.numeric(unname(wd_upper)), 
                    lower = as.numeric(unname(wd_lower)), tPRS_type = unname(fit_labs),
                    p_values = rep(paste0("(", round(as.numeric(unname(pval_tprs_slope[4,][, -c(1)]))[1:n_ct], 3),
                                          ", ", round(as.numeric(unname(pval_tprs_slope[5,][, -c(1)]))[1:n_ct], 3),
                                          ",", round(as.numeric(unname(pval_tprs_slope[6,][, -c(1)]))[1:n_ct], 3),")"), each=3),
                    R2 = rep(paste0("(", percent(round(as.numeric(unname(R2[4,][, -c(1)]))[1:n_ct], 4)),
                                    ", ", percent(round(as.numeric(unname(R2[5,][, -c(1)]))[1:n_ct], 4)),
                                    ",", percent(round(as.numeric(unname(R2[6,][, -c(1)]))[1:n_ct], 4)),")"), each=3))

mid_plot <- ggplot(ad_df, aes(x = betas, y = cell_types, xmin = upper, 
                              xmax = lower, group = tPRS_type, color = tPRS_type)) +
  geom_point(size = 4, aes(shape = tPRS_type), position = position_dodgev(height = .5)) +
  geom_pointrange(shape = 1, position = position_dodgev(height = .5)) +
  theme(legend.position = "bottom", axis.text.y = element_text(angle = 0, size = 11)) + 
  ggtitle("Main rate ratio of anx/dep tPRSs by cell type\nwith anxious/depressed score as outcome") +
  xlab("tPRS main rate ratios") + ylab("cell types")
mid_plot <- mid_plot +
  geom_vline(xintercept = 1, linetype="dashed")+ theme_hc()

p_right <- 
  ggplot(ad_df) +
  geom_text(
    aes(x = 0, y = (cell_types), label = R2),
    hjust = 0,
  ) + 
  theme_void() + ggtitle( paste0("% variance explained",":\n(mixed, male, female)")) + theme(plot.title = element_text(hjust = 0.8))

# final plot arrangement
mid_plot + p_right + plot_layout(design = layout)

mid_plot <- ggplot(wd_df, aes(x = betas, y = cell_types, xmin = upper, 
                              xmax = lower, group = tPRS_type, color = tPRS_type)) +
  geom_point(size = 4, aes(shape = tPRS_type), position = position_dodgev(height = .5)) +
  geom_pointrange(shape = 1, position = position_dodgev(height = .5)) +
  theme(legend.position = "bottom", axis.text.y = element_text(angle = 0, size = 11)) + 
  ggtitle("Main rate ratio of with/dep tPRSs by cell type\nwith withdrawn/depressed score as outcome") +
  xlab("tPRS main rate ratios") + ylab("cell types")
mid_plot <- mid_plot +
  geom_vline(xintercept = 1, linetype="dashed")+ theme_hc()

p_right <- 
  ggplot(wd_df) +
  geom_text(
    aes(x = 0, y = (cell_types), label = R2),
    hjust = 0,
  ) + 
  theme_void() + ggtitle( paste0("% variance explained",":\n(mixed, male, female)")) + theme(plot.title = element_text(hjust = 0.8))

# final plot arrangement
mid_plot + p_right + plot_layout(design = layout)

# logistic regs ==============================
pval_tprs_slope <- data.frame(regs = c("mlr_ad0", "mlr_ad_m","mlr_ad_f", "mlr_wd0", "mlr_wd_m", "mlr_wd_f"))
tprs_beta <- data.frame(regs = c("mlr_ad0", "mlr_ad_m","mlr_ad_f", "mlr_wd0", "mlr_wd_m", "mlr_wd_f"))
tprs_CI <- data.frame(regs = c("mlr_ad0", "mlr_ad_m","mlr_ad_f", "mlr_wd0", "mlr_wd_m", "mlr_wd_f"))
R2 <- data.frame(regs = c("mlr_ad0", "mlr_ad_m","mlr_ad_f", "mlr_wd0", "mlr_wd_m", "mlr_wd_f"))

for (i in 2:23){
  for (j in 1:1){
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]] <- 
      glm(anx_bin ~ twas_ad_tPRS_w[[i]] + age + sex + site + 
            genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
            genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
            genetic_pc_9 + genetic_pc_10, data = base_bin, family=binomial(link='logit'))
    R2[[colnames(twas_ad_tPRS_w)[i]]][1] <- 
      1 - summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]])$deviance/summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]])$null.deviance
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][1] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][1] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][1] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[1]]))[2])
    
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]] <- 
      glm(anx_bin ~ twas_ad_tPRS_w[which(base_bin$sex == "1"), i+22] + age + site + 
            genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
            genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
            genetic_pc_9 + genetic_pc_10, data = base_bin[which(base_bin$sex == "1"),],
          family=binomial(link='logit'))
    R2[[colnames(twas_ad_tPRS_w)[i]]][2] <- 
      1 - summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]])$deviance/summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]])$null.deviance
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][2] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][2] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][2] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[2]]))[2])
    
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]] <-
      glm(anx_bin ~ twas_ad_tPRS_w[which(base_bin$sex == "2"), i+44] + age + site + 
            genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
            genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
            genetic_pc_9 + genetic_pc_10, data = base_bin[which(base_bin$sex == "2"),],
          family=binomial(link='logit'))
    R2[[colnames(twas_ad_tPRS_w)[i]]][3] <- 
      1 - summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]])$deviance/summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]])$null.deviance
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][3] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][3] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][3] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[3]]))[2])
    
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]] <- 
      glm(with_bin ~ twas_wd_tPRS_w[[i]] + age + sex + site + 
            genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
            genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
            genetic_pc_9 + genetic_pc_10, data = base_bin, family=binomial(link='logit'))
    R2[[colnames(twas_ad_tPRS_w)[i]]][4] <- 
      1 - summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]])$deviance/summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]])$null.deviance
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][4] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][4] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][4] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[4]]))[2])
    
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]] <- 
      glm(with_bin ~ twas_wd_tPRS_w[which(base_bin$sex == "1"), i+22] + age + site + 
            genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
            genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
            genetic_pc_9 + genetic_pc_10, data = base_bin[which(base_bin$sex == "1"), ],
          family=binomial(link='logit'))
    R2[[colnames(twas_ad_tPRS_w)[i]]][5] <- 
      1 - summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]])$deviance/summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]])$null.deviance
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][5] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][5] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][5] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[5]]))[2])
    
    reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]] <- 
      glm(with_bin ~ twas_wd_tPRS_w[which(base_bin$sex == "2"), i+44] + age + site + 
            genetic_pc_1 + genetic_pc_2 + genetic_pc_3 + genetic_pc_4 + 
            genetic_pc_5 + genetic_pc_6 + genetic_pc_7 + genetic_pc_8 +
            genetic_pc_9 + genetic_pc_10, data = base_bin[which(base_bin$sex == "2"), ],
          family=binomial(link='logit'))
    R2[[colnames(twas_ad_tPRS_w)[i]]][6] <- 
      1 - summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]])$deviance/summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]])$null.deviance
    pval_tprs_slope[[colnames(twas_ad_tPRS_w)[i]]][6] <- 
      summary(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]])$coefficients[,4][2]
    tprs_beta[[colnames(twas_ad_tPRS_w)[i]]][6] <- 
      coef(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]])[2]
    tprs_CI[[colnames(twas_ad_tPRS_w)[i]]][6] <- 
      qt(0.975, df.residual(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]]))*
      sqrt(diag(vcov(reg_tprs[[colnames(twas_ad_tPRS_w)[i]]][[6]]))[2])
  }
  
}
ad_beta <- exp(unlist(tprs_beta[1:3,-1]))
wd_beta <- exp(unlist(tprs_beta[4:6,-1]))
labs <- rep(cell_types, each=3)
fit_labs <- rep(c("mixed","male","female"), n_ct)
ad_lower <- exp(unlist(tprs_beta[1:3,-1]) - unlist(tprs_CI[1:3,-1]))
wd_lower <- exp(unlist(tprs_beta[4:6,-1]) - unlist(tprs_CI[4:6,-1]))
ad_upper <- exp(unlist(tprs_beta[1:3,-1]) + unlist(tprs_CI[1:3,-1]))
wd_upper <- exp(unlist(tprs_beta[4:6,-1]) + unlist(tprs_CI[4:6,-1]))

ad_df <- data.frame(betas = as.numeric(unname(ad_beta)), cell_types = unname(labs), 
                    upper = as.numeric(unname(ad_upper)), 
                    lower = as.numeric(unname(ad_lower)), tPRS_type = unname(fit_labs),
                    p_values = rep(paste0("(", round(as.numeric(unname(pval_tprs_slope[1,][, -c(1)]))[1:n_ct], 3),
                                          ", ", round(as.numeric(unname(pval_tprs_slope[2,][, -c(1)]))[1:n_ct], 3),
                                          ",", round(as.numeric(unname(pval_tprs_slope[3,][, -c(1)]))[1:n_ct], 3), ")"), each=3),
                    R2 = rep(paste0("(", percent(round(as.numeric(unname(R2[1,][, -c(1)]))[1:n_ct], 4)),
                                    ", ", percent(round(as.numeric(unname(R2[2,][, -c(1)]))[1:n_ct], 4)),
                                    ",", percent(round(as.numeric(unname(R2[3,][, -c(1)]))[1:n_ct], 4)), ")"), each=3))

wd_df <- data.frame(betas = as.numeric(unname(wd_beta)), cell_types = unname(labs), 
                    upper = as.numeric(unname(wd_upper)), 
                    lower = as.numeric(unname(wd_lower)), tPRS_type = unname(fit_labs),
                    p_values = rep(paste0("(", round(as.numeric(unname(pval_tprs_slope[4,][, -c(1)]))[1:n_ct], 3),
                                          ", ", round(as.numeric(unname(pval_tprs_slope[5,][, -c(1)]))[1:n_ct], 3),
                                          ",", round(as.numeric(unname(pval_tprs_slope[6,][, -c(1)]))[1:n_ct], 3),")"), each=3),
                    R2 = rep(paste0("(", percent(round(as.numeric(unname(R2[4,][, -c(1)]))[1:n_ct], 4)),
                                    ", ", percent(round(as.numeric(unname(R2[5,][, -c(1)]))[1:n_ct], 4)),
                                    ",", percent(round(as.numeric(unname(R2[6,][, -c(1)]))[1:n_ct], 4)),")"), each=3))

mid_plot <- ggplot(ad_df, aes(x = betas, y = cell_types, xmin = upper, 
                              xmax = lower, group = tPRS_type, color = tPRS_type)) +
  geom_point(size = 4, aes(shape = tPRS_type), position = position_dodgev(height = .5)) +
  geom_pointrange(shape = 1, position = position_dodgev(height = .5)) +
  theme(legend.position = "bottom", axis.text.y = element_text(angle = 0, size = 11)) + 
  ggtitle("Odds ratio of anx/dep tPRSs by cell type\nwith anxious/depressed score as outcome") +
  xlab("tPRS odds ratios") + ylab("cell types")
mid_plot <- mid_plot +
  geom_vline(xintercept = 1, linetype="dashed")+ theme_hc()

p_right <- 
  ggplot(ad_df) +
  geom_text(
    aes(x = 0, y = (cell_types), label = R2),
    hjust = 0,
  ) + 
  theme_void() + ggtitle( paste0("% variance explained",":\n(mixed, male, female)")) + theme(plot.title = element_text(hjust = 0.8))

# final plot arrangement
mid_plot + p_right + plot_layout(design = layout)

mid_plot <- ggplot(wd_df, aes(x = betas, y = cell_types, xmin = upper, 
                              xmax = lower, group = tPRS_type, color = tPRS_type)) +
  geom_point(size = 4, aes(shape = tPRS_type), position = position_dodgev(height = .5)) +
  geom_pointrange(shape = 1, position = position_dodgev(height = .5)) +
  theme(legend.position = "bottom", axis.text.y = element_text(angle = 0, size = 11)) + 
  ggtitle("Odds ratio of with/dep tPRSs by cell type\nwith withdrawn/depressed score as outcome") +
  xlab("tPRS odds ratios") + ylab("cell types")
mid_plot <- mid_plot +
  geom_vline(xintercept = 1, linetype="dashed")+ theme_hc()

p_right <- 
  ggplot(wd_df) +
  geom_text(
    aes(x = 0, y = (cell_types), label = R2),
    hjust = 0,
  ) + 
  theme_void() + ggtitle( paste0("% variance explained",":\n(mixed, male, female)")) + theme(plot.title = element_text(hjust = 0.8))

# final plot arrangement
mid_plot + p_right + plot_layout(design = layout)



