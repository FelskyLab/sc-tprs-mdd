
load("sc-tprs-mdd/ct_gene_diff_uw.Rda")
load("sc-tprs-mdd/model_wgt_snp_df.Rda")
lin_files <- list.files("/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_data/gwas/only_model", 
                        pattern = "*.linear", full.names = TRUE)
snps_no_covar <- data.frame()
snps_w_covar <- data.frame()
for (f in lin_files){
  x <- read.delim(f, header = TRUE)
  x <- cbind(x, ct = rep(gsub("no_covar.|.glm.linear|w_covar.","",basename(f)), 
                      times = nrow(x)))
  if (grepl("no_covar", f)){
    snps_no_covar <- rbind(snps_no_covar, x)
    next
  }
  snps_w_covar <- rbind(snps_w_covar, x)
  
}
table(snps_no_covar[, "ID"])
table(snps_w_covar[, "ID"])

norm_ct_idx <- c(1,3,6,8,10,12,13,17)
bimod_ct_idx <- c(2,4,5,7,9,11,14,15,16,18,19,20,21,22)
names(ct_gene_diff_uw) <- gsub("-", "_", names(ct_gene_diff_uw))

norm_ct <- names(ct_gene_diff_uw)[norm_ct_idx]
bimod_ct <- names(ct_gene_diff_uw)[bimod_ct_idx]
diff_snps <- setdiff(unique(snps_no_covar[which(snps_no_covar$ct %in% bimod_ct), "ID"]),
                     unique(snps_no_covar[which(snps_no_covar$ct %in% norm_ct), "ID"]))
sim_snps <- intersect(unique(snps_no_covar[which(snps_no_covar$ct %in% bimod_ct), "ID"]),
                    unique(snps_no_covar[which(snps_no_covar$ct %in% norm_ct), "ID"]))
all_diff <- snps_no_covar[which(snps_no_covar$ID %in% diff_snps),]
all_sim <- snps_no_covar[which(snps_no_covar$ID %in% sim_snps),]
diff_snps_norm <- setdiff(unique(snps_no_covar[which(snps_no_covar$ct %in% norm_ct), "ID"]),
                          unique(snps_no_covar[which(snps_no_covar$ct %in% bimod_ct), "ID"]))
all_diff_norm <- snps_no_covar[which(snps_no_covar$ID %in% diff_snps_norm),]

sort(table(all_diff$ID))
sort(table(all_diff_norm$ID))
sort(table(all_sim$ID))
# all snps common between bimod and normal cts have same effect size + direction
# mean(all_diff_norm$BETA)
# [1] 0.01873607
# mean(all_sim$BETA)
# [1] -0.04661741
# mean(all_diff$BETA)
# [1] -0.04160039
# mean(snps_no_covar[which(snps_no_covar$ct %in% norm_ct), "BETA"])
# [1] -0.01549587
# mean(snps_no_covar[which(snps_no_covar$ct %in% bimod_ct), "BETA"])
# [1] -0.043747

snps_no_covar[which(snps_no_covar$ID == "chr19_19610913_G_A_b38"),]
snps_no_covar[which(snps_no_covar$ID == "chr20_35331517_C_A_b38"),]


# chr4_80235944_T_C_b38 common among B cells
# chr10_86964703_G_A_b38 common among dendritic cell and plasmacytoid dendritic cell (not conventional dendritic cell)

diff_snps <- setdiff(unique(snps_w_covar[which(snps_w_covar$ct %in% bimod_ct), "ID"]),
                     unique(snps_w_covar[which(snps_w_covar$ct %in% norm_ct), "ID"]))
sim_snps <- intersect(unique(snps_w_covar[which(snps_w_covar$ct %in% bimod_ct), "ID"]),
                      unique(snps_w_covar[which(snps_w_covar$ct %in% norm_ct), "ID"]))
all_diff <- snps_w_covar[which(snps_w_covar$ID %in% diff_snps),]
all_sim <- snps_w_covar[which(snps_w_covar$ID %in% sim_snps),]
diff_snps_norm <- setdiff(unique(snps_w_covar[which(snps_w_covar$ct %in% norm_ct), "ID"]),
                          unique(snps_w_covar[which(snps_w_covar$ct %in% bimod_ct), "ID"]))
all_diff_norm <- snps_w_covar[which(snps_w_covar$ID %in% diff_snps_norm),]

sort(table(all_diff$ID))
sort(table(all_diff_norm$ID))
sort(table(all_sim$ID))

table(all_diff[which(all_diff$ct %in% c("CD4-positive_alpha-beta_cytotoxic_T_cell_uw", 
                                        "central_memory_CD8-positive_alpha-beta_T_cell_uw",
                                        "erythrocyte_uw",
                                        "dendritic_cell_uw")),"ID"])

