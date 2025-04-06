# required: abcd_twas_dep.sh,  abcd_twas_dep_by_sex.sh
# goal: correlation among wittenberg and TWAS tPRSs, network presentation of overlapping genes among TWAS tPRSs
# next: twas_analysis_abcd.R

# compare TWAS tPRS to Wittenberg and other TWAS tPRSs============================


twas_gwas_betas_covars_ad <- data.frame(gene = witt_deg$ensembl_gene_id, wittenberg = witt_deg$b)

twas_gwas_betas_covars_wd <- twas_gwas_betas_covars_ad
for (f in cell_files[c(3,4,2,5,9,10,12,21,
                       14, 22, 15, 17,
                       16, 13,8, 1,
                       20, 11, 
                       18,19,7,6)]){
  ad <- fread(paste0(f,"/abcd_anxdep_assoc_covars.txt"), header = TRUE)
  wd <- fread(paste0(f,"/abcd_withdep_assoc_covars.txt"), header = TRUE)
  twas_gwas_betas_covars_ad <- merge(twas_gwas_betas_covars_ad, ad[which(ad$pvalue < 0.02), c("gene", "effect")], 
                                     by = c("gene"), sort = FALSE, all = TRUE)
  twas_gwas_betas_covars_wd <- merge(twas_gwas_betas_covars_wd, wd[which(wd$pvalue < 0.02), c("gene", "effect")], 
                                     by = c("gene"), sort = FALSE, all = TRUE)
}

corr_ad <- cor(twas_gwas_betas_covars_ad[, -1][,c(1, 1+c(3,4,2,5,9,10,12,21,
                                                         14, 22, 15, 17,
                                                         16, 13,8, 1,
                                                         20, 11, 
                                                         18,19,7,6))], use = "pairwise.complete.obs")
testRes <- corr.test(twas_gwas_betas_covars_ad[, -1][,c(1, 1+c(3,4,2,5,9,10,12,21,
                                                               14, 22, 15, 17,
                                                               16, 13,8, 1,
                                                               20, 11, 
                                                               18,19,7,6))])
colnames(corr_ad)[-(1)] <- 1:n_ct
row.names(corr_ad)[-(1)] <- 1:n_ct
corrplot(corr_ad, p.mat = testRes$p, sig.level = c(0.01, 0.05, 0.1), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'green', is.cor = FALSE, method = "square")
n_mat <- matrix(data=0,nrow=n_ct+1,ncol=n_ct+1)
for (i in 1:(n_ct+1)){
  for (j in 1:(n_ct+1)){
    if (i <= j){
      n_mat[i,j] <- length(intersect(which(!is.na(twas_gwas_betas_covars_ad[,i+1])), 
                                     which(!is.na(twas_gwas_betas_covars_ad[,j+1]))))
    }
  }
}
for (i in 1:(n_ct+1)){
  for (j in 1:(n_ct+1)){
    if (i <= j){
      corr_ad[i,j] <- n_mat[i,j]/max(n_mat, na.rm = T)
    }
  }
}
lower <- corr_ad[-1,][,-1]
colnames(lower) <- 1:n_ct
rownames(lower) <- 1:n_ct
upper <- n_mat[-1,][,-1]
colnames(upper) <- 1:n_ct
rownames(upper) <- 1:n_ct
corrplot(lower, p.mat = testRes$p[-1,][,-1], sig.level=c(0.01, 0.05, 0.1), pch.cex = 0.9, 
         type='lower', insig = 'label_sig', pch.col = 'green',tl.pos = "tl",
         is.cor = T, method = "square", diag = FALSE, cl.ratio = 0.1)
corrplot(upper, method = "color", addCoef.col = "black",
        type='upper', add = T, number.cex=0.5,
        col = COL1("YlGn"), is.corr = F, 
        number.digits = 0, tl.pos = "t", cl.ratio = 0.1, cl.pos = "n") 


corr_wd <- cor(twas_gwas_betas_covars_wd[,-1][,c(1, 1+c(3,4,2,5,9,10,12,21,
                                                        14, 22, 15, 17,
                                                        16, 13,8, 1,
                                                        20, 11, 
                                                        18,19,7,6))], use = "pairwise.complete.obs")
testRes <- corr.test(twas_gwas_betas_covars_wd[,-1][,c(1, 1+c(3,4,2,5,9,10,12,21,
                                                              14, 22, 15, 17,
                                                              16, 13,8, 1,
                                                              20, 11, 
                                                              18,19,7,6))])
colnames(corr_wd)[-(1)] <- 1:n_ct
row.names(corr_wd)[-(1)] <- 1:n_ct
corrplot(corr_wd, p.mat = testRes$p, sig.level = c(0.01, 0.05, 0.1), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'green', is.cor = FALSE, method = "square")

n_mat <- matrix(data=0,nrow=n_ct+1,ncol=n_ct+1)
for (i in 1:(n_ct+1)){
  for (j in 1:(n_ct+1)){
    if (i <= j){
      n_mat[i,j] <- length(intersect(which(!is.na(twas_gwas_betas_covars_wd[,i+1])), 
                                     which(!is.na(twas_gwas_betas_covars_wd[,j+1]))))
    }
  }
}
for (i in 1:(n_ct+1)){
  for (j in 1:(n_ct+1)){
    if (i <= j){
      corr_ad[i,j] <- n_mat[i,j]/max(n_mat, na.rm = T)
    }
  }
}
lower <- corr_wd[-1,][,-1]
colnames(lower) <- 1:n_ct
rownames(lower) <- 1:n_ct
upper <- n_mat[-1,][,-1]
colnames(upper) <- 1:n_ct
rownames(upper) <- 1:n_ct
corrplot(lower, p.mat = testRes$p[-1,][,-1], sig.level=c(0.01, 0.05, 0.1), pch.cex = 0.9, 
         type='lower', insig = 'label_sig', pch.col = 'green',tl.pos = "tl",
         is.cor = T, method = "square", diag = FALSE, cl.ratio = 0.1)
corrplot(upper, method = "color", addCoef.col = "black",
         type='upper', add = T, number.cex=0.5,
         col = COL1("YlGn"), is.corr = F, 
         number.digits = 0, tl.pos = "t", cl.ratio = 0.1, cl.pos = "n")


twas_genes_covars_ad <- data.frame(gene = witt_deg$ensembl_gene_id, effect = witt_deg$b, 
                                   cell_type = rep("wittenberg", nrow(witt_deg)))
twas_genes_covars_wd <- twas_genes_covars_ad

p_thresh <- 0.02
for (f in cell_files){
  ad <- fread(paste0(f,"/abcd_anxdep_assoc_covars.txt"), header = TRUE)
  wd <- fread(paste0(f,"/abcd_withdep_assoc_covars.txt"), header = TRUE)
  ad <- ad[which(ad$pvalue < p_thresh), c("gene", "effect")]
  ad$cell_type <- rep(basename(f), nrow(ad))
  wd <- wd[which(wd$pvalue < p_thresh), c("gene", "effect")]
  wd$cell_type <- rep(basename(f), nrow(wd))
  twas_genes_covars_ad <- rbind(twas_genes_covars_ad, ad)
  twas_genes_covars_wd <- rbind(twas_genes_covars_wd, wd)
  
}
sort(table(twas_genes_covars_ad$gene)[twas_genes_covars_ad$gene[which(twas_genes_covars_ad$cell_type == "wittenberg")]])

sort(table(twas_genes_covars_wd$gene)[twas_genes_covars_wd$gene[which(twas_genes_covars_wd$cell_type == "wittenberg")]])


sort(intersect(twas_genes_covars_ad$gene[-which(twas_genes_covars_ad$cell_type == "wittenberg")], 
               twas_genes_covars_wd$gene[-which(twas_genes_covars_wd$cell_type == "wittenberg")]))

common_1517_ad <- twas_gwas_betas_covars_ad[which(!is.na(twas_gwas_betas_covars_ad$`5`) & !is.na(twas_gwas_betas_covars_ad$`2`)), c("gene", "1", "2")]
common_1517_wd <- twas_gwas_betas_covars_wd[which(!is.na(twas_gwas_betas_covars_wd$`1`) & !is.na(twas_gwas_betas_covars_wd$`2`)), c("gene", "1", "2")]
nrow(common_1517_ad)
nrow(common_1517_wd)

#with males ====
twas_gwas_betas_male_ad <- data.frame(gene = witt_deg$ensembl_gene_id, wittenberg = witt_deg$b)

twas_gwas_betas_male_wd <- twas_gwas_betas_male_ad
for (f in cell_files[c(3,4,2,5,9,10,12,21,
                       14, 22, 15, 17,
                       16, 13,8, 1,
                       20, 11, 
                       18,19,7,6)]){
  ad <- fread(paste0(f,"/abcd_anxdep_assoc_covars_male.txt"), header = TRUE)
  wd <- fread(paste0(f,"/abcd_withdep_assoc_covars_male.txt"), header = TRUE)
  twas_gwas_betas_male_ad <- merge(twas_gwas_betas_male_ad, ad[which(ad$pvalue < 0.02), c("gene", "effect")], 
                                   by = c("gene"), sort = FALSE, all = TRUE)
  twas_gwas_betas_male_wd <- merge(twas_gwas_betas_male_wd, wd[which(wd$pvalue < 0.02), c("gene", "effect")], 
                                   by = c("gene"), sort = FALSE, all = TRUE)
}



corr_ad <- cor(twas_gwas_betas_male_ad[, -1][,c(1, 1+c(3,4,2,5,9,10,12,21,
                                                       14, 22, 15, 17,
                                                       16, 13,8, 1,
                                                       20, 11, 
                                                       18,19,7,6))], use = "pairwise.complete.obs")
testRes <- corr.test(twas_gwas_betas_male_ad[, -1][,c(1, 1+c(3,4,2,5,9,10,12,21,
                                                             14, 22, 15, 17,
                                                             16, 13,8, 1,
                                                             20, 11, 
                                                             18,19,7,6))])
colnames(corr_ad)[-(1)] <- 1:n_ct
row.names(corr_ad)[-(1)] <- 1:n_ct
corrplot(corr_ad, p.mat = testRes$p, sig.level = c(0.01, 0.05, 0.1), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'green', is.cor = FALSE, method = "square")

n_mat <- matrix(data=0,nrow=n_ct+1,ncol=n_ct+1)
for (i in 1:(n_ct+1)){
  for (j in 1:(n_ct+1)){
    if (i <= j){
      n_mat[i,j] <- length(intersect(which(!is.na(twas_gwas_betas_male_ad[,i+1])), 
                                     which(!is.na(twas_gwas_betas_male_ad[,j+1]))))
    }
  }
}
for (i in 1:(n_ct+1)){
  for (j in 1:(n_ct+1)){
    if (i <= j){
      corr_ad[i,j] <- n_mat[i,j]/max(n_mat, na.rm = T)
    }
  }
}
lower <- corr_ad[-1,][,-1]
colnames(lower) <- 1:n_ct
rownames(lower) <- 1:n_ct
upper <- n_mat[-1,][,-1]
colnames(upper) <- 1:n_ct
rownames(upper) <- 1:n_ct
corrplot(lower, p.mat = testRes$p[-1,][,-1], sig.level=c(0.01, 0.05, 0.1), pch.cex = 0.9, 
         type='lower', insig = 'label_sig', pch.col = 'green',tl.pos = "tl",
         is.cor = T, method = "square", diag = FALSE, cl.ratio = 0.1)
corrplot(upper, method = "color", addCoef.col = "black",
         type='upper', add = T, number.cex=0.5,
         col = COL1("YlGn"), is.corr = F, 
         number.digits = 0, tl.pos = "t", cl.ratio = 0.1, cl.pos = "n")

corr_wd <- cor(twas_gwas_betas_male_wd[,-1][,c(1, 1+c(3,4,2,5,9,10,12,21,
                                                      14, 22, 15, 17,
                                                      16, 13,8, 1,
                                                      20, 11, 
                                                      18,19,7,6))], use = "pairwise.complete.obs")
testRes <- corr.test(twas_gwas_betas_male_wd[,-1][,c(1, 1+c(3,4,2,5,9,10,12,21,
                                                            14, 22, 15, 17,
                                                            16, 13,8, 1,
                                                            20, 11, 
                                                            18,19,7,6))])
colnames(corr_wd)[-(1)] <- 1:n_ct
row.names(corr_wd)[-(1)] <- 1:n_ct
corrplot(corr_wd, p.mat = testRes$p, sig.level = c(0.01, 0.05, 0.1), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'green', is.cor = FALSE, method = "square")
n_mat <- matrix(data=0,nrow=n_ct+1,ncol=n_ct+1)
for (i in 1:(n_ct+1)){
  for (j in 1:(n_ct+1)){
    if (i <= j){
      n_mat[i,j] <- length(intersect(which(!is.na(twas_gwas_betas_male_wd[,i+1])), 
                                     which(!is.na(twas_gwas_betas_male_wd[,j+1]))))
    }
  }
}
for (i in 1:(n_ct+1)){
  for (j in 1:(n_ct+1)){
    if (i <= j){
      corr_ad[i,j] <- n_mat[i,j]/max(n_mat, na.rm = T)
    }
  }
}
lower <- corr_wd[-1,][,-1]
colnames(lower) <- 1:n_ct
rownames(lower) <- 1:n_ct
upper <- n_mat[-1,][,-1]
colnames(upper) <- 1:n_ct
rownames(upper) <- 1:n_ct
corrplot(lower, p.mat = testRes$p[-1,][,-1], sig.level=c(0.01, 0.05, 0.1), pch.cex = 0.9, 
         type='lower', insig = 'label_sig', pch.col = 'green',tl.pos = "tl",
         is.cor = T, method = "square", diag = FALSE, cl.ratio = 0.1)
corrplot(upper, method = "color", addCoef.col = "black",
         type='upper', add = T, number.cex=0.5,
         col = COL1("YlGn"), is.corr = F, 
         number.digits = 0, tl.pos = "t", cl.ratio = 0.1, cl.pos = "n")



twas_genes_male_ad <- data.frame(gene = witt_deg$ensembl_gene_id, effect = witt_deg$b, 
                                 cell_type = rep("wittenberg", nrow(witt_deg)))
twas_genes_male_wd <- twas_genes_male_ad
p_thresh <- 0.02
for (f in cell_files){
  ad <- fread(paste0(f,"/abcd_anxdep_assoc_covars_male.txt"), header = TRUE)
  wd <- fread(paste0(f,"/abcd_withdep_assoc_covars_male.txt"), header = TRUE)
  ad <- ad[which(ad$pvalue < p_thresh), c("gene", "effect")]
  ad$cell_type <- rep(basename(f), nrow(ad))
  wd <- wd[which(wd$pvalue < p_thresh), c("gene", "effect")]
  wd$cell_type <- rep(basename(f), nrow(wd))
  twas_genes_male_ad <- rbind(twas_genes_male_ad, ad)
  twas_genes_male_wd <- rbind(twas_genes_male_wd, wd)
  
}
write.csv(twas_genes_covars_wd$gene,"~/sc-tprs-mdd/wd_combined_twas_genes.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)

sort(table(twas_genes_male_ad$gene)[twas_genes_male_ad$gene[which(twas_genes_male_ad$cell_type == "wittenberg")]])

sort(table(twas_genes_male_wd$gene)[twas_genes_male_wd$gene[which(twas_genes_male_wd$cell_type == "wittenberg")]])


sort(intersect(twas_genes_male_ad$gene[-which(twas_genes_male_ad$cell_type == "wittenberg")], 
               twas_genes_male_wd$gene[-which(twas_genes_male_wd$cell_type == "wittenberg")]))

# with females ==========
twas_gwas_betas_female_ad <- data.frame(gene = witt_deg$ensembl_gene_id, wittenberg = witt_deg$b)

twas_gwas_betas_female_wd <- twas_gwas_betas_female_ad
for (f in cell_files[c(3,4,2,5,9,10,12,21,
                       14, 22, 15, 17,
                       16, 13,8, 1,
                       20, 11, 
                       18,19,7,6)]){
  ad <- fread(paste0(f,"/abcd_anxdep_assoc_covars_female.txt"), header = TRUE)
  wd <- fread(paste0(f,"/abcd_withdep_assoc_covars_female.txt"), header = TRUE)
  twas_gwas_betas_female_ad <- merge(twas_gwas_betas_female_ad, ad[which(ad$pvalue < 0.02), c("gene", "effect")], 
                                     by = c("gene"), sort = FALSE, all = TRUE)
  twas_gwas_betas_female_wd <- merge(twas_gwas_betas_female_wd, wd[which(wd$pvalue < 0.02), c("gene", "effect")], 
                                     by = c("gene"), sort = FALSE, all = TRUE)
}

corr_ad <- cor(twas_gwas_betas_female_ad[, -1][,c(1, 1+c(3,4,2,5,9,10,12,21,
                                                         14, 22, 15, 17,
                                                         16, 13,8, 1,
                                                         20, 11, 
                                                         18,19,7,6))], use = "pairwise.complete.obs")
testRes <- corr.test(twas_gwas_betas_female_ad[, -1][,c(1, 1+c(3,4,2,5,9,10,12,21,
                                                               14, 22, 15, 17,
                                                               16, 13,8, 1,
                                                               20, 11, 
                                                               18,19,7,6))])
colnames(corr_ad)[-(1)] <- 1:n_ct
row.names(corr_ad)[-(1)] <- 1:n_ct
corrplot(corr_ad, p.mat = testRes$p, sig.level = c(0.01, 0.05, 0.1), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'green', is.cor = FALSE, method = "square")
n_mat <- matrix(data=0,nrow=n_ct+1,ncol=n_ct+1)
for (i in 1:(n_ct+1)){
  for (j in 1:(n_ct+1)){
    if (i <= j){
      n_mat[i,j] <- length(intersect(which(!is.na(twas_gwas_betas_female_ad[,i+1])), 
                                     which(!is.na(twas_gwas_betas_female_ad[,j+1]))))
    }
  }
}
for (i in 1:(n_ct+1)){
  for (j in 1:(n_ct+1)){
    if (i <= j){
      corr_ad[i,j] <- n_mat[i,j]/max(n_mat, na.rm = T)
    }
  }
}
lower <- corr_ad[-1,][,-1]
colnames(lower) <- 1:n_ct
rownames(lower) <- 1:n_ct
upper <- n_mat[-1,][,-1]
colnames(upper) <- 1:n_ct
rownames(upper) <- 1:n_ct
corrplot(lower, p.mat = testRes$p[-1,][,-1], sig.level=c(0.01, 0.05, 0.1), pch.cex = 0.9, 
         type='lower', insig = 'label_sig', pch.col = 'green',tl.pos = "tl",
         is.cor = T, method = "square", diag = FALSE, cl.ratio = 0.1)
corrplot(upper, method = "color", addCoef.col = "black",
         type='upper', add = T, number.cex=0.5,
         col = COL1("YlGn"), is.corr = F, 
         number.digits = 0, tl.pos = "t", cl.ratio = 0.1, cl.pos = "n")
corr_wd <- cor(twas_gwas_betas_female_wd[,-1][,c(1, 1+c(3,4,2,5,9,10,12,21,
                                                        14, 22, 15, 17,
                                                        16, 13,8, 1,
                                                        20, 11, 
                                                        18,19,7,6))], use = "pairwise.complete.obs")
testRes <- corr.test(twas_gwas_betas_female_wd[,-1][,c(1, 1+c(3,4,2,5,9,10,12,21,
                                                              14, 22, 15, 17,
                                                              16, 13,8, 1,
                                                              20, 11, 
                                                              18,19,7,6))])
colnames(corr_wd)[-(1)] <- 1:n_ct
row.names(corr_wd)[-(1)] <- 1:n_ct
corrplot(corr_wd, p.mat = testRes$p, sig.level = c(0.01, 0.05, 0.1), pch.cex = 0.9,
         insig = 'label_sig', pch.col = 'green', is.cor = FALSE, method = "square")

n_mat <- matrix(data=0,nrow=n_ct+1,ncol=n_ct+1)
for (i in 1:(n_ct+1)){
  for (j in 1:(n_ct+1)){
    if (i <= j){
      n_mat[i,j] <- length(intersect(which(!is.na(twas_gwas_betas_female_wd[,i+1])), 
                                     which(!is.na(twas_gwas_betas_female_wd[,j+1]))))
    }
  }
}
for (i in 1:(n_ct+1)){
  for (j in 1:(n_ct+1)){
    if (i <= j){
      corr_ad[i,j] <- n_mat[i,j]/max(n_mat, na.rm = T)
    }
  }
}
lower <- corr_wd[-1,][,-1]
colnames(lower) <- 1:n_ct
rownames(lower) <- 1:n_ct
upper <- n_mat[-1,][,-1]
colnames(upper) <- 1:n_ct
rownames(upper) <- 1:n_ct
corrplot(lower, p.mat = testRes$p[-1,][,-1], sig.level=c(0.01, 0.05, 0.1), pch.cex = 0.9, 
         type='lower', insig = 'label_sig', pch.col = 'green',tl.pos = "tl",
         is.cor = T, method = "square", diag = FALSE, cl.ratio = 0.1)
corrplot(upper, method = "color", addCoef.col = "black",
         type='upper', add = T, number.cex=0.5,
         col = COL1("YlGn"), is.corr = F, 
         number.digits = 0, tl.pos = "t", cl.ratio = 0.1, cl.pos = "n")


twas_genes_female_ad <- data.frame(gene = witt_deg$ensembl_gene_id, effect = witt_deg$b, 
                                   cell_type = rep("wittenberg", nrow(witt_deg)))
twas_genes_female_wd <- twas_genes_female_ad
p_thresh <- 0.02
for (f in cell_files){
  ad <- fread(paste0(f,"/abcd_anxdep_assoc_covars_female.txt"), header = TRUE)
  wd <- fread(paste0(f,"/abcd_withdep_assoc_covars_female.txt"), header = TRUE)
  ad <- ad[which(ad$pvalue < p_thresh), c("gene", "effect")]
  ad$cell_type <- rep(basename(f), nrow(ad))
  wd <- wd[which(wd$pvalue < p_thresh), c("gene", "effect")]
  wd$cell_type <- rep(basename(f), nrow(wd))
  twas_genes_female_ad <- rbind(twas_genes_female_ad, ad)
  twas_genes_female_wd <- rbind(twas_genes_female_wd, wd)
  
}

sort(table(twas_genes_female_ad$gene)[twas_genes_female_ad$gene[-which(twas_genes_female_ad$cell_type == "wittenberg")]])

sort(table(twas_genes_female_wd$gene)[twas_genes_female_wd$gene[-which(twas_genes_female_wd$cell_type == "wittenberg")]])


sort(intersect(twas_genes_female_ad$gene[-which(twas_genes_female_ad$cell_type == "wittenberg")], 
               twas_genes_female_wd$gene[-which(twas_genes_female_wd$cell_type == "wittenberg")]))


# TWAS sc-tPRS shared genes among cell types=================
# graph/network form

# graph looks too busy :( need to workshop
library(igraph)
library(visNetwork)

tprs_corr <- cor((twas_ad_tPRS_w[,2:23]),
                 use = "pairwise.complete.obs")
testRes <- corr.test(twas_ad_tPRS_w[,2:23])

mixed_ad_to_mixed_ad <- data.frame(from = numeric(sum(1:(n_ct-1))), to = numeric(sum(1:(n_ct-1))), 
                             weight = numeric(sum(1:(n_ct-1))), color = numeric(sum(1:(n_ct-1))))
color.gradient <- function(x, colors=c("purple","cyan","yellow"), colsteps=1000) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}
edge_corr_colour <- color.gradient(c(tprs_corr))
idx <- 0
for (i in 1:n_ct){
  for (j in 2:n_ct){
    if (i < j){
      idx <- idx + 1
      mixed_ad_to_mixed_ad$from[idx] <- cell_types[i]
      mixed_ad_to_mixed_ad$to[idx] <- cell_types[j]
      mixed_ad_to_mixed_ad$weight[idx] <- length(intersect(twas_genes_covars_ad$gene[which(twas_genes_covars_ad$cell_type ==  cell_types[i])], 
                                                                twas_genes_covars_ad$gene[which(twas_genes_covars_ad$cell_type ==  cell_types[j])]))
      if (testRes$p[i,j] < 0.05){
         mixed_ad_to_mixed_ad$color[idx] <- edge_corr_colour[(i-1)*n_ct+(j)]
      } else{
        mixed_ad_to_mixed_ad$color[idx] <- "#888888"
      }
      mixed_ad_to_mixed_ad$title[idx] <- paste0("# genes in common:",mixed_ad_to_mixed_ad$weight[idx],
                                                ", r = ", c(tprs_corr)[(i-1)*n_ct+(j)],", p = ", testRes$p[i,j])
      
    }
  }
}
mixed_ad_to_mixed_ad$label <- mixed_ad_to_mixed_ad$weight
mixed_ad_to_mixed_ad$value <- mixed_ad_to_mixed_ad$weight

g <- graph_from_data_frame(mixed_ad_to_mixed_ad, directed = FALSE)
visn <- toVisNetworkData(g)
## copy column "weight" to new column "value" in list "edges"
visn$edges$value <- visn$edges$weight

visNetwork(visn$nodes, visn$edges) %>%
  visIgraphLayout(layout = "layout_in_circle") 
# Plot the weighted graph


final_network <- visNetwork(visn$nodes, mixed_ad_to_mixed_ad, width = "100%", height = 2000) %>%
  visIgraphLayout() %>%
  visNodes(
    shape = "dot",
    color = list(
      background = "#0085AF",
      border = "#013848",
      highlight = "#FF8000"
    ),
    shadow = list(enabled = F, size = 10)
  ) %>%
  visEdges(
    shadow = T,
    color = list(highlight = "#C62F4B")
  ) %>%
  visLayout(randomSeed = 330)
visSave(final_network, paste0("mixed_ad_twas_tprs.html"), selfcontained = TRUE, background = "white")


f_ad_to_m_ad <- data.frame(from = numeric(n_ct*n_ct), to = numeric(n_ct*n_ct), 
                                   weight = numeric(n_ct*n_ct), color = numeric(n_ct*n_ct))

idx <- 0
for (i in 1:n_ct){
  for (j in 1:n_ct){
      idx <- idx + 1
      f_ad_to_m_ad$from[idx] <- paste0(which(semantic_order == i),"_f")
      f_ad_to_m_ad$to[idx] <- paste0(which(semantic_order == j),"_m")
      f_ad_to_m_ad$weight[idx] <- length(intersect(twas_genes_female_ad$gene[which(twas_genes_female_ad$cell_type ==  cell_types[i])], 
                                                           twas_genes_male_ad$gene[which(twas_genes_male_ad$cell_type ==  cell_types[j])]))
      
  }
}
f_ad_to_m_ad$label <- f_ad_to_m_ad$weight
f_ad_to_m_ad$value <- f_ad_to_m_ad$weight
f_ad_to_m_ad$title <- f_ad_to_m_ad$weight
f_ad_to_m_ad$color  <- color.gradient(c(f_ad_to_m_ad$weight))

g <- graph_from_data_frame(f_ad_to_m_ad, directed = FALSE)
visn <- toVisNetworkData(g)
## copy column "weight" to new column "value" in list "edges"
visn$edges$value <- visn$edges$weight

visNetwork(visn$nodes, visn$edges) %>%
  visIgraphLayout(layout = "layout_in_circle") 
# Plot the weighted graph


final_network <- visNetwork(visn$nodes, f_ad_to_m_ad, width = "100%", height = 2000) %>%
  visIgraphLayout(layout = "layout_with_sugiyama") %>%
  visNodes(
    shape = "dot",
    color = list(
      background = "#0085AF",
      border = "#013848",
      highlight = "#FF8000"
    ),
    shadow = list(enabled = F, size = 10)
  ) %>%
  visEdges(
    shadow = T,
    color = list(color = "#DDDDDD", highlight = "#C62F4B")
  ) %>%
  visLayout(randomSeed = 330)
visSave(final_network, paste0("fm_ad_twas_tprs.html"), selfcontained = TRUE, background = "white")

f_wd_to_m_wd <- data.frame(from = numeric(n_ct*n_ct), to = numeric(n_ct*n_ct), 
                           weight = numeric(n_ct*n_ct), color = numeric(n_ct*n_ct))

idx <- 0
for (i in 1:n_ct){
  for (j in 1:n_ct){
    idx <- idx + 1
    f_wd_to_m_wd$from[idx] <- paste0(which(semantic_order == i),"_f")
    f_wd_to_m_wd$to[idx] <- paste0(which(semantic_order == j),"_m")
    f_wd_to_m_wd$weight[idx] <- length(intersect(twas_genes_female_wd$gene[which(twas_genes_female_wd$cell_type ==  cell_types[i])], 
                                                 twas_genes_male_wd$gene[which(twas_genes_male_wd$cell_type ==  cell_types[j])]))
    
  }
}
f_wd_to_m_wd$label <- f_wd_to_m_wd$weight
f_wd_to_m_wd$value <- f_wd_to_m_wd$weight
f_wd_to_m_wd$title <- f_wd_to_m_wd$weight
f_wd_to_m_wd$color  <- color.gradient(c(f_wd_to_m_wd$weight))

g <- graph_from_data_frame(f_wd_to_m_wd, directed = FALSE)
visn <- toVisNetworkData(g)
## copy column "weight" to new column "value" in list "edges"
visn$edges$value <- visn$edges$weight

visNetwork(visn$nodes, visn$edges) %>%
  visIgraphLayout(layout = "layout_in_circle") 
# Plot the weighted graph


final_network <- visNetwork(visn$nodes, f_wd_to_m_wd, width = "100%", height = 2000) %>%
  visIgraphLayout(layout = "layout_with_sugiyama") %>%
  visNodes(
    shape = "dot",
    color = list(
      background = "#0085AF",
      border = "#013848",
      highlight = "#FF8000"
    ),
    shadow = list(enabled = F, size = 10)
  ) %>%
  visEdges(
    shadow = T,
    color = list(color = "#DDDDDD", highlight = "#C62F4B")
  ) %>%
  visLayout(randomSeed = 330)
visSave(final_network, paste0("fm_wd_twas_tprs.html"), selfcontained = TRUE, background = "white")
