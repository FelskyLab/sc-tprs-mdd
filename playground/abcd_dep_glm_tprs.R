load("sc-tprs-mdd/og_white_all_tprs_immune.Rda")
if (! requireNamespace("tidyverse", quietly=TRUE)) {
  install.packages("tidyverse")
}
library(tidyverse)
if (! requireNamespace("vtable", quietly=TRUE)) {
  install.packages("vtable")
}
library(vtable)
if (! requireNamespace("tidyr", quietly=TRUE)) {
  packageurl <- "https://github.com/tidyverse/tidyr/archive/refs/tags/v1.2.1.tar.gz"
  install.packages(packageurl, repos=NULL, type="source")
}
if (! requireNamespace("brms", quietly=TRUE)) {
  install.packages("brms")
}
library(brms)
if (! requireNamespace("cmdstanr", quietly=TRUE)) {
  install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
}
library(cmdstanr)
#install_cmdstan(overwrite = TRUE, release_url = "https://github.com/stan-dev/cmdstan/releases/download/v2.9.0/cmdstan-2.9.0.tar.gz")
#
if (! requireNamespace("ggeffects", quietly=TRUE)) {
  install.packages("ggeffects")
}
library(ggeffects)
library(tidyr)
library(ggthemes)

events <- c("baseline_year_1_arm_1", "1_year_follow_up_y_arm_1", "2_year_follow_up_y_arm_1",
            "3_year_follow_up_y_arm_1", "4_year_follow_up_y_arm_1")

reg_tprs <- list(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0", 
                          "slr_ad1","slr_wd1", "mlr_ad1", "mlr_wd1",
                          "slr_ad2","slr_wd2", "mlr_ad2", "mlr_wd2",
                          "slr_ad3","slr_wd3", "mlr_ad3", "mlr_wd3",
                          "slr_ad4","slr_wd4", "mlr_ad4", "mlr_wd4"))
n_ct <- length(cell_types)
tprs_start_idx <- 19
#og_white_all_tprs <- white_all_tprs
bin_all_tprs <- og_white_all_tprs[,-(9:26)]
colnames(bin_all_tprs) <- gsub("-", "_", colnames(bin_all_tprs))
bin_all_tprs$anx_bin <- as.factor(as.integer(bin_all_tprs$anxdep_t > 65))
bin_all_tprs$with_bin <- as.factor(as.integer(bin_all_tprs$withdep_t > 65))
bin_all_tprs$anx_sqrt <- sqrt(bin_all_tprs$anxdep_t)
bin_all_tprs$with_sqrt <- sqrt(bin_all_tprs$withdep_t)
#bin_all_tprs <- bin_all_tprs[-(which(bin_all_tprs$sex==3)), ] # remove intersex-male
bin_all_tprs <- bin_all_tprs[-(which(is.na(bin_all_tprs$anx_bin))), ] # remove rows with NA in depression scores
bin_all_tprs <- bin_all_tprs[-(which(bin_all_tprs$IID == "NDAR_INVJHJDGEFN")), ] # remove participant where their baseline
# depression score is NA
bin_all_tprs$sex <- droplevels(bin_all_tprs$sex)
scale_tprs <- bin_all_tprs
for (event in events){
  scale_tprs[(scale_tprs$event_name == event),][, -c(1, 2, 3, 4, 6, 7, 8, 
                                                     (tprs_start_idx+2*n_ct), (tprs_start_idx+2*n_ct+1), 
                                                     (tprs_start_idx+2*n_ct+2), (tprs_start_idx+2*n_ct+3))] <- 
    scale(scale_tprs[(scale_tprs$event_name == event),][, -c(1, 2, 3, 4, 6, 7, 8, 
                                                             (tprs_start_idx+2*n_ct), (tprs_start_idx+2*n_ct+1), 
                                                             (tprs_start_idx+2*n_ct+2), (tprs_start_idx+2*n_ct+3))])
  
}

cont <- bin_all_tprs[c(bin_all_tprs$event_name == events[1]),][, -c(1, 2, 3, 4, 6, 
                                                                    (tprs_start_idx+2*n_ct), (tprs_start_idx+2*n_ct+1))]
catg <- bin_all_tprs[c(bin_all_tprs$event_name == events[1]),][, c(4, 6, 
                                                                   (tprs_start_idx+2*n_ct), (tprs_start_idx+2*n_ct+1))]

levels(catg$sex) <- c("male", "female")
catg$site <- gsub("site", "", catg$site)
catg$anx_bin <- as.character(catg$anx_bin)
catg$with_bin <- as.character(catg$with_bin)

#have 22 diff ct
ggplot(stack(cont[, c(1,4:13)], varying = colnames(cont)[c(1,5:14)]), aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x') + theme_stata()

ggplot(stack( (cont[, c(14:20)]), varying = colnames(cont)[14:20]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (cont[, c(21:27)]), varying = colnames(cont)[21:27]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x') + theme_stata() +
  theme(text=element_text(size=9))
ggplot(stack( (cont[, c(28:35)]), varying = colnames(cont)[28:35]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x') + theme_stata()+
  theme(text=element_text(size=9))
ggplot(stack( (cont[, c(n_ct+14:20 )]), varying = colnames(cont)[n_ct+14:20]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')
ggplot(stack( (cont[, c(n_ct+21:27)]), varying = colnames(cont)[n_ct+21:27]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')
ggplot(stack( (cont[, c(n_ct+28:34)]), varying = colnames(cont)[n_ct+28:34]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')

ggplot(stack(catg[,c(1,2)], varying = colnames(catg)[1:2]), aes(values)) + 
  geom_bar() + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata()

ggplot(stack( (cont[, c(2, 3, (ncol(cont) - 1),ncol(cont))]), varying = colnames(cont)[c(2, 3, (ncol(cont) - 1),ncol(cont))]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata()


library(ggbeeswarm)
library(viridis)
ggplot(cbind(stack( (cont[, c(26:34)]), varying = colnames(cont)[26:34]), 
             anxdep = rep(catg$anx_bin, 9),
             sex = rep(catg$sex, 9)), 
       aes(x = anxdep, y = values, fill = anxdep, color = sex)) + 
  geom_violin(position="dodge", alpha=0.5, color = 'grey') +
  scale_fill_viridis(discrete=T, name="") +
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))

ggplot(as.data.frame(cbind(tPRS = cont[, c(24)], 
             anxdep = rep(catg$anx_bin, 1),
             sex = rep(catg$sex, 1))), 
       aes(x = anxdep, y = tPRS, fill = anxdep, color = sex)) + 
  scale_fill_viridis(discrete=T, name="") +
  geom_beeswarm(alpha = 0.3) +
  geom_violin(position="dodge", alpha=0.5, color = NA) +
  theme_stata() + 
  theme(text=element_text(size=9))

ggplot(stack( (cont[, c(21:27)]), varying = colnames(cont)[21:27]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x') + theme_stata() +
  theme(text=element_text(size=9))
ggplot(stack( (cont[, c(28:35)]), varying = colnames(cont)[28:35]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x') + theme_stata()+
  theme(text=element_text(size=9))


df1_long <- melt(cbind(catg$sex,(cont[, c(2, 3)])), id.vars=c("catg$sex"))
ggplot(df1_long, aes(x=`catg$sex`,y=value,fill=`catg$sex`)) + 
  geom_boxplot() + facet_wrap(~variable)+ theme_stata() + labs(fill='sex') 
df1_long <- melt(cbind(catg$sex,(cont[, c((ncol(cont) - 1),ncol(cont))])), id.vars=c("catg$sex"))
ggplot(df1_long, aes(x=`catg$sex`,y=value,fill=`catg$sex`)) + 
  geom_boxplot() + facet_wrap(~variable)+ theme_stata() + labs(fill='sex') 
a <- do.call("cbind", list(catg$sex,cont$memory_B_cell_uw, 
                           cont$double_negative_thymocyte_uw, cont$innate_lymphoid_cell_uw))
a <- as.data.frame(a)
a[, -c(1)] <- sapply(a[, -c(1)], as.numeric)
colnames(a) <- c("sex", "memory B cell uw", "double negative thymocyte uw", "innate lymphoid cell uw")
df1_long <- melt(a, 
                 id.vars=c("sex"))
ggplot(df1_long, aes(x=`sex`,y=value,fill=sex)) + 
  geom_boxplot() + facet_wrap(~variable)+ theme_stata() + labs(fill='sex') 

base_bin <- scale_tprs[c(bin_all_tprs$event_name == events[1]),]

gam_anx_sq <- glm(anx_sqrt ~ plasmablast_uw*sex + age  + genetic_pc_1 + genetic_pc_2 +
                    genetic_pc_3 + genetic_pc_4 + genetic_pc_5 + genetic_pc_6 + genetic_pc_7 +
                    genetic_pc_8 + genetic_pc_9 + genetic_pc_10, data=base_bin, family=Gamma(link='log'))
summary(gam_anx_sq)

gam_with_sq <- glm(with_sqrt ~ plasmablast_uw*sex + age  + genetic_pc_1 + genetic_pc_2 +
                     genetic_pc_3 + genetic_pc_4 + genetic_pc_5 + genetic_pc_6 + genetic_pc_7 +
                     genetic_pc_8 + genetic_pc_9 + genetic_pc_10, data=base_bin, family=Gamma(link='log'))
summary(gam_with_sq)$coefficients

log_anx <- glm(anx_bin ~ plasmablast_uw*sex + age  + genetic_pc_1 + genetic_pc_2 +
                 genetic_pc_3 + genetic_pc_4 + genetic_pc_5 + genetic_pc_6 + genetic_pc_7 +
                 genetic_pc_8 + genetic_pc_9 + genetic_pc_10, data=base_bin, family = binomial(link = "logit"))
summary(log_anx)
log_with <- glm(with_bin ~ plasmablast_uw*sex + age  + genetic_pc_1 + genetic_pc_2 +
                  genetic_pc_3 + genetic_pc_4 + genetic_pc_5 + genetic_pc_6 + genetic_pc_7 +
                  genetic_pc_8 + genetic_pc_9 + genetic_pc_10, data=base_bin, family = binomial(link = "logit"))
summary(log_with)


library('survival')
wei_anx_sq <- survreg(Surv(anx_sqrt) ~ plasmablast_uw*sex + age  + genetic_pc_1 + genetic_pc_2 +
                        genetic_pc_3 + genetic_pc_4 + genetic_pc_5 + genetic_pc_6 + genetic_pc_7 +
                        genetic_pc_8 + genetic_pc_9 + genetic_pc_10, data=base_bin, dist='weibull')
summary(wei_anx_sq)




#adjR2 <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
#pval_mod <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
pval_tprs_slope_gam_sx <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
tprs_beta_gam_sx <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
tprs_CI_gam_sx <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
R2_gam_sx <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))

pval_sex_ixn_gam<- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
beta_sex_ixn_gam<- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
CI_sex_ixn_gam <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
for (i in tprs_start_idx:(tprs_start_idx - 1 + 2*n_ct)){
  for (j in 1:1){
    ad1 <- glm(anx_sqrt ~ base_bin[,i]*sex, data=base_bin, family=Gamma(link='log'))
    
    pval_tprs_slope_gam_sx[[colnames(base_bin)[i]]][1] <- 
      summary(ad1)$coefficients[,4][2]
    tprs_beta_gam_sx[[colnames(base_bin)[i]]][1] <- 
      coef(ad1)[2]
    tprs_CI_gam_sx[[colnames(base_bin)[i]]][1] <- 
      qt(0.975, df.residual(ad1))*
      sqrt(diag(vcov(ad1))[2])
    pval_sex_ixn_gam[[colnames(base_bin)[i]]][1] <-summary(ad1)$coefficients[,4][4]
    beta_sex_ixn_gam[[colnames(base_bin)[i]]][1] <-coef(ad1)[4]
    CI_sex_ixn_gam[[colnames(base_bin)[i]]][1] <-
      qt(0.975, df.residual(ad1))*sqrt(diag(vcov(ad1))[4])
    R2_gam_sx[[colnames(base_bin)[i]]][1] <- 1 - summary(ad1)$deviance/summary(ad1)$null.deviance
    
    wd1 <- glm(with_sqrt ~ base_bin[,i]*sex, data = base_bin, family=Gamma(link='log'))
    
    pval_tprs_slope_gam_sx[[colnames(base_bin)[i]]][2] <- 
      summary(wd1)$coefficients[,4][2]
    tprs_beta_gam_sx[[colnames(base_bin)[i]]][2] <- 
      coef(wd1)[2]
    tprs_CI_gam_sx[[colnames(base_bin)[i]]][2] <- 
      qt(0.975, df.residual(wd1))*
      sqrt(diag(vcov(wd1))[2])
    pval_sex_ixn_gam[[colnames(base_bin)[i]]][2] <-summary(wd1)$coefficients[,4][4]
    beta_sex_ixn_gam[[colnames(base_bin)[i]]][2] <-coef(wd1)[4]
    CI_sex_ixn_gam[[colnames(base_bin)[i]]][2] <-
      qt(0.975, df.residual(ad1))*sqrt(diag(vcov(ad1))[4])
    R2_gam_sx[[colnames(base_bin)[i]]][2] <- 1 - summary(wd1)$deviance/summary(wd1)$null.deviance
    
    ad1 <- glm(anx_sqrt ~ base_bin[,i]*sex + age + site + genetic_pc_1 + genetic_pc_2 +
                 genetic_pc_3 + genetic_pc_4 + genetic_pc_5 + genetic_pc_6 + genetic_pc_7 +
                 genetic_pc_8 + genetic_pc_9 + genetic_pc_10, data=base_bin, family=Gamma(link='log'))
    
    pval_tprs_slope_gam_sx[[colnames(base_bin)[i]]][3] <- 
      summary(ad1)$coefficients[,4][2]
    tprs_beta_gam_sx[[colnames(base_bin)[i]]][3] <- 
      coef(ad1)[2]
    tprs_CI_gam_sx[[colnames(base_bin)[i]]][3] <- 
      qt(0.975, df.residual(ad1))*
      sqrt(diag(vcov(ad1))[2])
    pval_sex_ixn_gam[[colnames(base_bin)[i]]][3] <- 
      summary(ad1)$coefficients[,4][15]
    beta_sex_ixn_gam[[colnames(base_bin)[i]]][3] <- 
      coef(ad1)[15]
    CI_sex_ixn_gam[[colnames(base_bin)[i]]][3] <- 
      qt(0.975, df.residual(ad1))*
      sqrt(diag(vcov(ad1))[15])
    R2_gam_sx[[colnames(base_bin)[i]]][3] <- 1 - summary(ad1)$deviance/summary(ad1)$null.deviance
    
    wd1 <- glm(with_sqrt ~ base_bin[,i]*sex + age + site  + genetic_pc_1 + genetic_pc_2 +
                 genetic_pc_3 + genetic_pc_4 + genetic_pc_5 + genetic_pc_6 + genetic_pc_7 +
                 genetic_pc_8 + genetic_pc_9 + genetic_pc_10, data=base_bin, family=Gamma(link='log'))
    
    pval_tprs_slope_gam_sx[[colnames(base_bin)[i]]][4] <- 
      summary(wd1)$coefficients[,4][2]
    tprs_beta_gam_sx[[colnames(base_bin)[i]]][4] <- 
      coef(wd1)[2]
    tprs_CI_gam_sx[[colnames(base_bin)[i]]][4] <- 
      qt(0.975, df.residual(wd1))*
      sqrt(diag(vcov(wd1))[2])
    pval_sex_ixn_gam[[colnames(base_bin)[i]]][4] <- 
      summary(wd1)$coefficients[,4][15]
    beta_sex_ixn_gam[[colnames(base_bin)[i]]][4] <- 
      coef(wd1)[15]
    CI_sex_ixn_gam[[colnames(base_bin)[i]]][4] <- 
      qt(0.975, df.residual(wd1))*
      sqrt(diag(vcov(wd1))[15])
    R2_gam_sx[[colnames(base_bin)[i]]][4] <- 1 - summary(wd1)$deviance/summary(wd1)$null.deviance
  }
  
}
pval_tprs_slope_gam <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
tprs_beta_gam <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
tprs_CI_gam <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
R2_gam <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
for (i in tprs_start_idx:(tprs_start_idx - 1 + 2*n_ct)){
  for (j in 1:1){
    ad1 <- glm(anx_sqrt ~ base_bin[,i], data=base_bin, family=Gamma(link='log'))
    
    pval_tprs_slope_gam[[colnames(base_bin)[i]]][1] <- 
      summary(ad1)$coefficients[,4][2]
    tprs_beta_gam[[colnames(base_bin)[i]]][1] <- 
      coef(ad1)[2]
    tprs_CI_gam[[colnames(base_bin)[i]]][1] <- 
      qt(0.975, df.residual(ad1))*
      sqrt(diag(vcov(ad1))[2])
    R2_gam[[colnames(base_bin)[i]]][1] <- 1 - summary(ad1)$deviance/summary(ad1)$null.deviance
    
    wd1 <- glm(with_sqrt ~ base_bin[,i], data = base_bin, family=Gamma(link='log'))
    
    pval_tprs_slope_gam[[colnames(base_bin)[i]]][2] <- 
      summary(wd1)$coefficients[,4][2]
    tprs_beta_gam[[colnames(base_bin)[i]]][2] <- 
      coef(wd1)[2]
    tprs_CI_gam[[colnames(base_bin)[i]]][2] <- 
      qt(0.975, df.residual(wd1))*
      sqrt(diag(vcov(wd1))[2])
    R2_gam[[colnames(base_bin)[i]]][2] <- 1 - summary(wd1)$deviance/summary(wd1)$null.deviance
    
    ad1 <- glm(anx_sqrt ~ base_bin[,i] + sex + age + site + genetic_pc_1 + genetic_pc_2 +
                 genetic_pc_3 + genetic_pc_4 + genetic_pc_5 + genetic_pc_6 + genetic_pc_7 +
                 genetic_pc_8 + genetic_pc_9 + genetic_pc_10, data=base_bin, family=Gamma(link='log'))
    
    pval_tprs_slope_gam[[colnames(base_bin)[i]]][3] <- 
      summary(ad1)$coefficients[,4][2]
    tprs_beta_gam[[colnames(base_bin)[i]]][3] <- 
      coef(ad1)[2]
    tprs_CI_gam[[colnames(base_bin)[i]]][3] <- 
      qt(0.975, df.residual(ad1))*
      sqrt(diag(vcov(ad1))[2])
    R2_gam[[colnames(base_bin)[i]]][3] <- 1 - summary(ad1)$deviance/summary(ad1)$null.deviance
    
    wd1 <- glm(with_sqrt ~ base_bin[,i] + sex + age + site + genetic_pc_1 + genetic_pc_2 +
                 genetic_pc_3 + genetic_pc_4 + genetic_pc_5 + genetic_pc_6 + genetic_pc_7 +
                 genetic_pc_8 + genetic_pc_9 + genetic_pc_10, data=base_bin, family=Gamma(link='log'))
    
    pval_tprs_slope_gam[[colnames(base_bin)[i]]][4] <- 
      summary(wd1)$coefficients[,4][2]
    tprs_beta_gam[[colnames(base_bin)[i]]][4] <- 
      coef(wd1)[2]
    tprs_CI_gam[[colnames(base_bin)[i]]][4] <- 
      qt(0.975, df.residual(wd1))*
      sqrt(diag(vcov(wd1))[2])
    R2_gam[[colnames(base_bin)[i]]][4] <- 1 - summary(wd1)$deviance/summary(wd1)$null.deviance
  }
  
}

pval_tprs_slope_log_sx <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
tprs_beta_log_sx <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
tprs_CI_log_sx <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
R2_log_sx <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))

pval_sex_ixn_log<- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
beta_sex_ixn_log<- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
CI_sex_ixn_log <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
for (i in tprs_start_idx:(tprs_start_idx - 1 + 2*n_ct)){
  for (j in 1:1){
    ad1 <- glm(anx_bin ~ base_bin[,i]*sex, data=base_bin, family=binomial(link='logit'))
    
    pval_tprs_slope_log_sx[[colnames(base_bin)[i]]][1] <- 
      summary(ad1)$coefficients[,4][2]
    tprs_beta_log_sx[[colnames(base_bin)[i]]][1] <- 
      coef(ad1)[2]
    tprs_CI_log_sx[[colnames(base_bin)[i]]][1] <- 
      qt(0.975, df.residual(ad1))*
      sqrt(diag(vcov(ad1))[2])
    pval_sex_ixn_log[[colnames(base_bin)[i]]][1] <-summary(ad1)$coefficients[,4][4]
    beta_sex_ixn_log[[colnames(base_bin)[i]]][1] <-coef(ad1)[4]
    CI_sex_ixn_log[[colnames(base_bin)[i]]][1] <-
      qt(0.975, df.residual(ad1))*
      sqrt(diag(vcov(ad1))[4])
    R2_log_sx[[colnames(base_bin)[i]]][1] <- 1 - summary(ad1)$deviance/summary(ad1)$null.deviance
    
    wd1 <- glm(with_bin ~ base_bin[,i]*sex, data = base_bin, family=binomial(link='logit'))
    
    pval_tprs_slope_log_sx[[colnames(base_bin)[i]]][2] <- 
      summary(wd1)$coefficients[,4][2]
    tprs_beta_log_sx[[colnames(base_bin)[i]]][2] <- 
      coef(wd1)[2]
    tprs_CI_log_sx[[colnames(base_bin)[i]]][2] <- 
      qt(0.975, df.residual(wd1))*
      sqrt(diag(vcov(wd1))[2])
    pval_sex_ixn_log[[colnames(base_bin)[i]]][2] <-summary(wd1)$coefficients[,4][4]
    beta_sex_ixn_log[[colnames(base_bin)[i]]][2] <-coef(wd1)[4]
    CI_sex_ixn_log[[colnames(base_bin)[i]]][2] <-
      qt(0.975, df.residual(wd1))*
      sqrt(diag(vcov(ad1))[4])
    R2_log_sx[[colnames(base_bin)[i]]][2] <- 1 - summary(wd1)$deviance/summary(wd1)$null.deviance
    
    ad1 <- glm(anx_bin ~ base_bin[,i]*sex + age + site + genetic_pc_1 + genetic_pc_2 +
                 genetic_pc_3 + genetic_pc_4 + genetic_pc_5 + genetic_pc_6 + genetic_pc_7 +
                 genetic_pc_8 + genetic_pc_9 + genetic_pc_10, data=base_bin, family=binomial(link='logit'))
    
    pval_tprs_slope_log_sx[[colnames(base_bin)[i]]][3] <- 
      summary(ad1)$coefficients[,4][2]
    tprs_beta_log_sx[[colnames(base_bin)[i]]][3] <- 
      coef(ad1)[2]
    tprs_CI_log_sx[[colnames(base_bin)[i]]][3] <- 
      qt(0.975, df.residual(ad1))*
      sqrt(diag(vcov(ad1))[2])
    pval_sex_ixn_log[[colnames(base_bin)[i]]][3] <- 
      summary(ad1)$coefficients[,4][15]
    beta_sex_ixn_log[[colnames(base_bin)[i]]][3] <- 
      coef(ad1)[15]
    CI_sex_ixn_log[[colnames(base_bin)[i]]][3] <- 
      qt(0.975, df.residual(ad1))*
      sqrt(diag(vcov(ad1))[15])
    R2_log_sx[[colnames(base_bin)[i]]][3] <- 1 - summary(ad1)$deviance/summary(ad1)$null.deviance
    
    wd1 <- glm(with_bin ~ base_bin[,i]*sex + age + site + genetic_pc_1 + genetic_pc_2 +
                 genetic_pc_3 + genetic_pc_4 + genetic_pc_5 + genetic_pc_6 + genetic_pc_7 +
                 genetic_pc_8 + genetic_pc_9 + genetic_pc_10, data=base_bin, family=binomial(link='logit'))
    
    pval_tprs_slope_log_sx[[colnames(base_bin)[i]]][4] <- 
      summary(wd1)$coefficients[,4][2]
    tprs_beta_log_sx[[colnames(base_bin)[i]]][4] <- 
      coef(wd1)[2]
    tprs_CI_log_sx[[colnames(base_bin)[i]]][4] <- 
      qt(0.975, df.residual(wd1))*
      sqrt(diag(vcov(wd1))[2])
    pval_sex_ixn_log[[colnames(base_bin)[i]]][4] <- 
      summary(wd1)$coefficients[,4][15]
    beta_sex_ixn_log[[colnames(base_bin)[i]]][4] <- 
      coef(wd1)[15]
    CI_sex_ixn_log[[colnames(base_bin)[i]]][4] <- 
      qt(0.975, df.residual(wd1))*
      sqrt(diag(vcov(wd1))[15])
    R2_log_sx[[colnames(base_bin)[i]]][4] <- 1 - summary(wd1)$deviance/summary(wd1)$null.deviance
  }
  
}
pval_tprs_slope_log <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
tprs_beta_log <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
tprs_CI_log <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
R2_log <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
for (i in tprs_start_idx:(tprs_start_idx - 1 + 2*n_ct)){
  for (j in 1:1){
    ad1 <- glm(anx_bin ~ base_bin[,i], data=base_bin, family=binomial(link='logit'))
    
    pval_tprs_slope_log[[colnames(base_bin)[i]]][1] <- 
      summary(ad1)$coefficients[,4][2]
    tprs_beta_log[[colnames(base_bin)[i]]][1] <- 
      coef(ad1)[2]
    tprs_CI_log[[colnames(base_bin)[i]]][1] <- 
      qt(0.975, df.residual(ad1))*
      sqrt(diag(vcov(ad1))[2])
    R2_log[[colnames(base_bin)[i]]][1] <- 1 - summary(ad1)$deviance/summary(ad1)$null.deviance
    
    wd1 <- glm(with_bin ~ base_bin[,i], data = base_bin, family=binomial(link='logit'))
    
    pval_tprs_slope_log[[colnames(base_bin)[i]]][2] <- 
      summary(wd1)$coefficients[,4][2]
    tprs_beta_log[[colnames(base_bin)[i]]][2] <- 
      coef(wd1)[2]
    tprs_CI_log[[colnames(base_bin)[i]]][2] <- 
      qt(0.975, df.residual(wd1))*
      sqrt(diag(vcov(wd1))[2])
    R2_log[[colnames(base_bin)[i]]][2] <- 1 - summary(wd1)$deviance/summary(wd1)$null.deviance
    
    ad1 <- glm(anx_bin ~ base_bin[,i] + sex + age + site + genetic_pc_1 + genetic_pc_2 +
                 genetic_pc_3 + genetic_pc_4 + genetic_pc_5 + genetic_pc_6 + genetic_pc_7 +
                 genetic_pc_8 + genetic_pc_9 + genetic_pc_10, data=base_bin, family=binomial(link='logit'))
    
    pval_tprs_slope_log[[colnames(base_bin)[i]]][3] <- 
      summary(ad1)$coefficients[,4][2]
    tprs_beta_log[[colnames(base_bin)[i]]][3] <- 
      coef(ad1)[2]
    tprs_CI_log[[colnames(base_bin)[i]]][3] <- 
      qt(0.975, df.residual(ad1))*
      sqrt(diag(vcov(ad1))[2])
    R2_log[[colnames(base_bin)[i]]][3] <- 1 - summary(ad1)$deviance/summary(ad1)$null.deviance
    
    wd1 <- glm(with_bin ~ base_bin[,i] + sex + age + site + genetic_pc_1 + genetic_pc_2 +
                 genetic_pc_3 + genetic_pc_4 + genetic_pc_5 + genetic_pc_6 + genetic_pc_7 +
                 genetic_pc_8 + genetic_pc_9 + genetic_pc_10, data=base_bin, family=binomial(link='logit'))
    
    pval_tprs_slope_log[[colnames(base_bin)[i]]][4] <- 
      summary(wd1)$coefficients[,4][2]
    tprs_beta_log[[colnames(base_bin)[i]]][4] <- 
      coef(wd1)[2]
    tprs_CI_log[[colnames(base_bin)[i]]][4] <- 
      qt(0.975, df.residual(wd1))*
      sqrt(diag(vcov(wd1))[2])
    R2_log[[colnames(base_bin)[i]]][4] <- 1 - summary(wd1)$deviance/summary(wd1)$null.deviance
  }
}
#for logistic reg

library(ggstance)
library(patchwork)
library(scales)
#only unweighted
ad_beta <- c(tprs_beta_log[1,][, -c(1, (n_ct + 2):(2*n_ct + 1))], tprs_beta_log[3,][, -c(1,  (n_ct + 2):(2*n_ct + 1))])
wd_beta <- c(tprs_beta_log[2,][, -c(1,  (n_ct + 2):(2*n_ct + 1))], tprs_beta_log[4,][, -c(1,  (n_ct + 2):(2*n_ct + 1))])
labs <- c(cell_types, cell_types)
fit_labs <- c(rep("w/o covar", n_ct), rep("w/ covar", n_ct))
ad_lower <- c((ad_beta[1:n_ct] - tprs_CI_log[1,][, -c(1,  (n_ct + 2):(2*n_ct + 1))]), 
              ad_beta[(n_ct + 1):(2*n_ct)] - tprs_CI_log[3,][, -c(1,  (n_ct + 2):(2*n_ct + 1))])
wd_lower <- c((wd_beta[1:n_ct] - tprs_CI_log[2,][, -c(1,  (n_ct + 2):(2*n_ct + 1))]),
              wd_beta[(n_ct + 1):(2*n_ct)] - tprs_CI_log[4,][, -c(1,  (n_ct + 2):(2*n_ct + 1))])
ad_upper <- c((ad_beta[1:n_ct] + tprs_CI_log[1,][, -c(1,  (n_ct + 2):(2*n_ct + 1))]), 
              ad_beta[(n_ct + 1):(2*n_ct)] + tprs_CI_log[3,][, -c(1,  (n_ct + 2):(2*n_ct + 1))])
wd_upper <- c((wd_beta[1:n_ct] + tprs_CI_log[2,][, -c(1,  (n_ct + 2):(2*n_ct + 1))]),
              wd_beta[(n_ct + 1):(2*n_ct)] + tprs_CI_log[4,][, -c(1,  (n_ct + 2):(2*n_ct + 1))])

ad_df <- data.frame(betas = as.numeric(unname(ad_beta)), cell_types = unname(labs), 
                    upper = as.numeric(unname(ad_upper)), 
                    lower = as.numeric(unname(ad_lower)), tPRS_type = unname(fit_labs),
                    p_values = paste0("(", round(as.numeric(unname(pval_tprs_slope_log[1,][, -c(1)]))[1:n_ct], 3),
                                      ", ", round(as.numeric(unname(pval_tprs_slope_log[3,][, -c(1)]))[1:n_ct], 3),
                                      ")"),
                    R2 = paste0("(", percent(round(as.numeric(unname(R2_log[1,][, -c(1)]))[1:n_ct], 4)),
                                ", ", percent(round(as.numeric(unname(R2_log[3,][, -c(1)]))[1:n_ct], 4)),
                                ")"))

wd_df <- data.frame(betas = as.numeric(unname(wd_beta)), cell_types = unname(labs), 
                    upper = as.numeric(unname(wd_upper)), 
                    lower = as.numeric(unname(wd_lower)), tPRS_type = unname(fit_labs),
                    p_values = paste0("(", round(as.numeric(unname(pval_tprs_slope_log[2,][, -c(1)]))[1:n_ct], 3),
                                      ", ", round(as.numeric(unname(pval_tprs_slope_log[4,][, -c(1)]))[1:n_ct], 3),
                                      ")"),
                    R2 = paste0("(", percent(round(as.numeric(unname(R2_log[2,][, -c(1)]))[1:n_ct], 4)),
                                ", ", percent(round(as.numeric(unname(R2_log[4,][, -c(1)]))[1:n_ct], 4)),
                                ")"))

mid_plot <- ggplot(ad_df, aes(x = betas, y = cell_types, xmin = upper, 
                              xmax = lower, group = tPRS_type, color = tPRS_type)) +
  geom_point(size = 4, aes(shape = tPRS_type), position = position_dodgev(height = .5)) +
  geom_pointrange(shape = 1, position = position_dodgev(height = .5)) +
  theme(legend.position = "bottom", axis.text.y = element_text(angle = 0, size = 11)) + 
  ggtitle("Main log-odds ratio by unweighted cell type\ntPRS in anxdep ~ tPRS") +
  xlab("tPRS main log-odds ratio") + ylab("cell types")
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
  ggtitle("Main log-odds ratio by unweighted cell type\n tPRS in withdep ~ tPRS") +
  xlab("tPRS main log-odds ratio") + ylab("cell types")
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
ad_beta <- c(tprs_beta_log[1,][, -c(1:(n_ct + 1))], tprs_beta_log[3,][, -c(1:(n_ct + 1))])
wd_beta <- c(tprs_beta_log[2,][, -c(1:(n_ct + 1))], tprs_beta_log[4,][, -c(1:(n_ct + 1))])
labs <- c(cell_types, cell_types)
fit_labs <- c(rep("w/o covar", n_ct), rep("w/ covar", n_ct))
ad_lower <- c((ad_beta[1:n_ct] - tprs_CI_log[1,][, -c(1:(n_ct + 1))]), 
              ad_beta[(n_ct + 1):(2*n_ct)] - tprs_CI_log[3,][, -c(1:(n_ct + 1))])
wd_lower <- c((wd_beta[1:n_ct] - tprs_CI_log[2,][, -c(1:(n_ct + 1))]),
              wd_beta[(n_ct + 1):(2*n_ct)] - tprs_CI_log[4,][, -c(1:(n_ct + 1))])
ad_upper <- c((ad_beta[1:n_ct] + tprs_CI_log[1,][, -c(1:(n_ct + 1))]), 
              ad_beta[(n_ct + 1):(2*n_ct)] + tprs_CI_log[3,][, -c(1:(n_ct + 1))])
wd_upper <- c((wd_beta[1:n_ct] + tprs_CI_log[2,][, -c(1:(n_ct + 1))]),
              wd_beta[(n_ct + 1):(2*n_ct)] + tprs_CI_log[4,][, -c(1:(n_ct + 1))])

ad_df <- data.frame(betas = as.numeric(unname(ad_beta)), cell_types = unname(labs), 
                    upper = as.numeric(unname(ad_upper)), 
                    lower = as.numeric(unname(ad_lower)), tPRS_type = unname(fit_labs),
                    p_values = paste0("(", round(as.numeric(unname(pval_tprs_slope_log[1,][, -c(1)]))[(n_ct + 1):(2*n_ct)], 3),
                                      ", ", round(as.numeric(unname(pval_tprs_slope_log[3,][, -c(1)]))[(n_ct + 1):(2*n_ct)], 3),
                                      ")"),
                    R2 = paste0("(", percent(round(as.numeric(unname(R2_log[1,][, -c(1)]))[(n_ct + 1):(2*n_ct)], 4)),
                                ", ", percent(round(as.numeric(unname(R2_log[3,][, -c(1)]))[(n_ct + 1):(2*n_ct)], 4)),
                                ")"))

wd_df <- data.frame(betas = as.numeric(unname(wd_beta)), cell_types = unname(labs), 
                    upper = as.numeric(unname(wd_upper)), 
                    lower = as.numeric(unname(wd_lower)), tPRS_type = unname(fit_labs),
                    p_values = paste0("(", round(as.numeric(unname(pval_tprs_slope_log[2,][, -c(1)]))[(n_ct + 1):(2*n_ct)], 3),
                                      ", ", round(as.numeric(unname(pval_tprs_slope_log[4,][, -c(1)]))[(n_ct + 1):(2*n_ct)], 3),
                                      ")"),
                    R2 = paste0("(", percent(round(as.numeric(unname(R2_log[2,][, -c(1)]))[(n_ct + 1):(2*n_ct)], 4)),
                                ", ", percent(round(as.numeric(unname(R2_log[4,][, -c(1)]))[(n_ct + 1):(2*n_ct)], 4)),
                                ")"))

mid_plot <- ggplot(ad_df, aes(x = betas, y = cell_types, xmin = upper, 
                              xmax = lower, group = tPRS_type, color = tPRS_type)) +
  geom_point(size = 4, aes(shape = tPRS_type), position = position_dodgev(height = .5)) +
  geom_pointrange(shape = 1, position = position_dodgev(height = .5)) +
  theme(legend.position = "bottom", axis.text.y = element_text(angle = 0, size = 11)) + 
  ggtitle("Main log-odds ratio by weighted cell type\ntPRS in anxdep ~ tPRS") +
  xlab("tPRS main log-odds ratio") + ylab("cell types")
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
  ggtitle("Main log-odds ratio by weighted cell type\ntPRS in withdep ~ tPRS") +
  xlab("tPRS main log-odds ratio") + ylab("cell types")
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


wd1 <- glm(with_bin ~ dendritic_cell_w, data = base_bin, family=binomial(link='logit'))
plot(ggpredict(wd1,"dendritic_cell_w"))

newdat <- data.frame(dendritic_cell_w=seq(min(base_bin$dendritic_cell_w), 
                                          max(base_bin$dendritic_cell_w),len=100))
newdat$binarized_withdep <- predict(wd1, newdata=newdat, type="response")
binarized
plot(~dendritic_cell_w, data=base_bin, col="red4")
lines(binarized_withdep~dendritic_cell_w, data=newdat, col="green4", lwd=2)

# fit zibeta model
base_bin$prop_anx <- (base_bin$anxdep_t - min(base_bin$anxdep_t))
base_bin$prop_anx <-   base_bin$prop_anx/(max(base_bin$prop_anx) + 0.001)
base_bin$prop_with <- (base_bin$withdep_t - min(base_bin$withdep_t))
base_bin$prop_with <-   base_bin$prop_with/(max(base_bin$prop_with) + 0.001)
#colnames(base_bin) <- gsub("-", "_", colnames(base_bin))
pval_tprs_slope_zib <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
tprs_beta_zib <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
tprs_CI_zib <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))
R2_zib <- data.frame(regs = c("slr_ad0","slr_wd0", "mlr_ad0", "mlr_wd0"))

#set_cmdstan_path(path = "/external/rprshnas01/kcni/mding/R/x86_64-conda-linux-gnu-library/4.0/cmdstanr")

priors <- c(
  set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
  set_prior("normal(0, 1)", class = "b"),
  set_prior("logistic(0, 1)", class = "Intercept", dpar = "zi")
)
ad1 <- brm(
  bf(prop_anx~ dendritic_cell_w, phi ~ 1, zi ~ 1),
  data = base_bin,
  family = zero_inflated_beta(),
  chains = 4, iter = 2000, warmup = 1000,
  cores = 4, seed = 1009,
  #prior = priors,
  #backend = "cmdstanr"
)
ad2 <- brm(
  bf(prop_anx~ dendritic_cell_w + sex + age + site + genetic_pc_1 + genetic_pc_2 +
       genetic_pc_3 + genetic_pc_4 + genetic_pc_5 + genetic_pc_6 + genetic_pc_7 +
       genetic_pc_8 + genetic_pc_9 + genetic_pc_10, phi ~ 1, zi ~ 1),
  data = base_bin,
  family = zero_inflated_beta(),
  chains = 4, iter = 2000, warmup = 1000,
  cores = 4, seed = 1009,
  #prior = priors,
  #backend = "cmdstanr"
)

ad3 <- zibr(
  logistic_cov = base_bin[,c(tprs_start_idx, 5, 9:18)], 
  beta_cov = base_bin[,c(tprs_start_idx, 5, 9:18)], Y = base_bin$prop_anx,
  subject_ind = base_bin$IID, time_ind = as.data.frame(rep.int(1, nrow(base_bin)))
)
get_prior(
  bf(prop_anx~ dendritic_cell_w + sex + age + site + genetic_pc_1 + genetic_pc_2 +
       genetic_pc_3 + genetic_pc_4 + genetic_pc_5 + genetic_pc_6 + genetic_pc_7 +
       genetic_pc_8 + genetic_pc_9 + genetic_pc_10, phi ~ 1, zi ~ 1),
  data = base_bin,
  family = zero_inflated_beta()
)

get_prior(
  bf(prop_anx~ dendritic_cell_w, phi ~ 1, zi ~ 1),
  data = base_bin,
  family = zero_inflated_beta()
)
for (i in tprs_start_idx:(tprs_start_idx - 1 + 2*n_ct)){
  for (j in 1:1){
    ad1 <- brm(
      bf(brmsformula(paste0("prop_anx~", colnames(base_bin)[i])),phi ~ 1, zi ~ 1),
      data = base_bin,
      family = zero_inflated_beta(),
      chains = 4, iter = 2000, warmup = 1000,
      cores = 4, seed = 1009,
      #prior = priors,
      #backend = "cmdstanr",
      #file = "model_beta_zi_int_only"
    )
    
    pval_tprs_slope_log[[colnames(base_bin)[i]]][1] <- 
      summary(ad1)$coefficients[,4][2]
    tprs_beta_log[[colnames(base_bin)[i]]][1] <- 
      coef(ad1)[2]
    tprs_CI_log[[colnames(base_bin)[i]]][1] <- 
      qt(0.975, df.residual(ad1))*
      sqrt(diag(vcov(ad1))[2])
    R2_log[[colnames(base_bin)[i]]][1] <- 1 - summary(ad1)$deviance/summary(ad1)$null.deviance
    
    wd1 <- glm(with_bin ~ base_bin[,i], data = base_bin, family=binomial(link='logit'))
    
    pval_tprs_slope_log[[colnames(base_bin)[i]]][2] <- 
      summary(wd1)$coefficients[,4][2]
    tprs_beta_log[[colnames(base_bin)[i]]][2] <- 
      coef(wd1)[2]
    tprs_CI_log[[colnames(base_bin)[i]]][2] <- 
      qt(0.975, df.residual(wd1))*
      sqrt(diag(vcov(wd1))[2])
    R2_log[[colnames(base_bin)[i]]][2] <- 1 - summary(wd1)$deviance/summary(wd1)$null.deviance
    
    ad1 <- glm(anx_bin ~ base_bin[,i] + sex + age + site + genetic_pc_1 + genetic_pc_2 +
                 genetic_pc_3 + genetic_pc_4 + genetic_pc_5 + genetic_pc_6 + genetic_pc_7 +
                 genetic_pc_8 + genetic_pc_9 + genetic_pc_10, data=base_bin, family=binomial(link='logit'))
    
    pval_tprs_slope_log[[colnames(base_bin)[i]]][3] <- 
      summary(ad1)$coefficients[,4][2]
    tprs_beta_log[[colnames(base_bin)[i]]][3] <- 
      coef(ad1)[2]
    tprs_CI_log[[colnames(base_bin)[i]]][3] <- 
      qt(0.975, df.residual(ad1))*
      sqrt(diag(vcov(ad1))[2])
    R2_log[[colnames(base_bin)[i]]][3] <- 1 - summary(ad1)$deviance/summary(ad1)$null.deviance
    
    wd1 <- glm(with_bin ~ base_bin[,i] + sex + age + site + genetic_pc_1 + genetic_pc_2 +
                 genetic_pc_3 + genetic_pc_4 + genetic_pc_5 + genetic_pc_6 + genetic_pc_7 +
                 genetic_pc_8 + genetic_pc_9 + genetic_pc_10, data=base_bin, family=binomial(link='logit'))
    
    pval_tprs_slope_log[[colnames(base_bin)[i]]][4] <- 
      summary(wd1)$coefficients[,4][2]
    tprs_beta_log[[colnames(base_bin)[i]]][4] <- 
      coef(wd1)[2]
    tprs_CI_log[[colnames(base_bin)[i]]][4] <- 
      qt(0.975, df.residual(wd1))*
      sqrt(diag(vcov(wd1))[2])
    R2_log[[colnames(base_bin)[i]]][4] <- 1 - summary(wd1)$deviance/summary(wd1)$null.deviance
  }
}


