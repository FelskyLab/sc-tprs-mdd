
setwd("/external/rprshnas01/kcni/mding/sc-tprs-mdd")
for (i in 1:22){
  df <- read.delim(paste0("abcd_imputed/with_CMC/abcd_pred_", 5, ".txt"), header = TRUE, sep = "\t", dec = ".")
  ABCD_CMC_5 <- as.numeric(as.vector(unlist(df[, -c(1,2)])))
  hist(ABCD_CMC_5, breaks = 100)
  
}

for (i in 1:22){
  df <- read.delim(paste0("abcd_imputed/with_sc_blood/platelet/abcd_pred_", 1, ".txt"), header = TRUE, sep = "\t", dec = ".")
  ABCD_platelet_1 <- as.numeric(as.vector(unlist(df[, -c(1,2)])))
  hist(ABCD_platelet_1, breaks = 100)
  qqnorm(ABCD_platelet_1)
  qqline()
}

df <- read.delim(paste0("abcd_imputed/with_sc_blood/platelet/abcd_pred_", 1, ".txt"), header = TRUE, sep = "\t", dec = ".")
ABCD_platelet_1 <- as.numeric(as.vector(unlist(df[, -c(1,2)])))
hist(ABCD_platelet_1, breaks = 100)
qqnorm(ABCD_platelet_1)
qqline()
df <- read.delim(paste0("abcd_imputed/with_sc_blood/peripheral_blood_mononuclear_cell/abcd_pred_", 21, ".txt"), header = TRUE, sep = "\t", dec = ".")
ABCD_PBMC_21 <- as.numeric(as.vector(unlist(df[, -c(1,2)])))
qqnorm(ABCD_PBMC_21)
qqline()