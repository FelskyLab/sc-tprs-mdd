

for (i in 1:22){
  df <- read.delim(paste0("/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_imputed/with_CMC/abcd_pred_", 5, ".txt"), header = TRUE, sep = "\t", dec = ".")
  ABCD_CMC_5 <- as.numeric(as.vector(unlist(df[, -c(1,2)])))
  hist(ABCD_CMC_5, breaks = 100)
  qqnorm(ABCD_CMC_5)
  qqline(ABCD_CMC_5)
  
}

for (i in 1:22){
  df <- read.delim(paste0("with_sc_blood38/CD4-positive_alpha-beta_T_cell/abcd_pred_", 1, ".txt"), header = TRUE, sep = "\t", dec = ".")
  ABCD_platelet_1 <- as.numeric(as.vector(unlist(df[, -c(1,2)])))
  hist(ABCD_platelet_1, breaks = 100)
  qqnorm(ABCD_platelet_1)
  qqline()
}
#install.packages("data.table")
setwd("/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_imputed")
library(data.table)
df <- fread(paste0("with_sc_blood38/CD4-positive_alpha-beta_T_cell/abcd_pred_", 1, ".txt"), header = TRUE)
CD4_positive_alpha_beta_T_cell <- as.numeric(as.vector(unlist(df[, -c(1,2)][(sample(1:(nrow(df)-1), size = 2000)),])))
CD4_positive_alpha_beta_T_cell <- CD4_positive_alpha_beta_T_cell[(sample(1:(length(CD4_positive_alpha_beta_T_cell)), size = 2000))]
png("CD4-positive_alpha-beta_T_cellhist.png")
hist(CD4_positive_alpha_beta_T_cell, breaks = 100)
dev.off()
png("CD4-positive_alpha-beta_T_cellqqnorm.png")
qqnorm(CD4_positive_alpha_beta_T_cell)
qqline(CD4_positive_alpha_beta_T_cell)
dev.off()