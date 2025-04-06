# required: load_abcd_data.R
# goal: explore why wittenberg tPRSs are bimodal
# - is it due to a few highly influential SNPs or not?
# next: dead end... no key, distribution defining SNPs found


# histograms of every model's weight distr

load("sc-tprs-mdd/sc_tPRS_betas_immune.Rda")

sort_betas <- cbind(tprs_beta[, 2:(n_ct + 1)][,rev(order(tprs_beta[3,2:(n_ct + 1)]))],
                    tprs_beta[, 2:(n_ct + 1)+n_ct][,rev(order(tprs_beta[3,2:(n_ct + 1)+n_ct]))])

# by beta > 0

ggplot(stack( (scale_tprs[(colnames(sort_betas)[2:8 - 1 ])]), varying = colnames(sort_betas)[2:8 -1]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (scale_tprs[(colnames(sort_betas)[9:15 -1])]), varying = colnames(sort_betas)[9:15-1]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (scale_tprs[(colnames(sort_betas)[16:(n_ct+1)-1])]), varying = colnames(sort_betas)[16:(n_ct+1)-1]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))

normal_range <- c(-1, 1)
low_range <- c(-10,-2)
range_tprs<- scale_tprs[, c(19:(18+2*n_ct))]
range_tprs[,1:ncol(range_tprs)] <- NA
range_tprs[scale_tprs[, c(19:(18+2*n_ct))] > normal_range[1] & scale_tprs[, c(19:(18+2*n_ct))] < normal_range[2]] <- "normal"
range_tprs[scale_tprs[, c(19:(18+2*n_ct))] > low_range[1] & scale_tprs[, c(19:(18+2*n_ct))] < low_range[2]] <- "low"
if (!requireNamespace("RSQLite", quietly = TRUE)){
  install.packages("RSQLite")
}
library(RSQLite)
if (!requireNamespace("DBI", quietly = TRUE)){
  install.packages("DBI")
}
library(DBI)


model_files <- list.files(path="abcd/workspace/mding/immune_cell_models", pattern="*.db", full.names=TRUE, recursive=FALSE)
models <- list()
read_model <- function(filename){
  fn <- sub('\\.db$', '', basename(filename))
  sqlite.driver <- dbDriver("SQLite")
  db <- dbConnect(sqlite.driver, dbname = filename)
  
  model <- dbReadTable(db,"weights")
  dbDisconnect(db)
  #format chr and position
  #model_chr <- as.numeric(str_remove(str_extract(model$varID, "chr[0-9]+"), "chr"))
  #model_pos <- as.numeric(str_remove(str_extract(model$varID, "_[0-9]+"), "_"))
  
  models[[fn]] <<- data.frame(SNP = model$varID, 
                              ZSCORE = scale(model$weight), gene = model$gene)
}

sapply(model_files, read_model)
if (!requireNamespace("qqman", quietly = TRUE)){
  install.packages("qqman")
}
library(qqman)
normal_idx <- c(1,6,8,17)
mid_idx <- c(3, 13,16)
low_idx <- c(15,19,21)

for (i in 1:n_ct){
  print(paste(names(models)[i], length(unique(models[[i]]$SNP)), length(unique(models[[i]]$gene))))
}


manhattan(models[[normal_idx[1]]], ylim = c(-0.015, 0.015), logp = FALSE, p = "P", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(names(models)[normal_idx[1]]," model weights"))
manhattan(models[[normal_idx[2]]], ylim = c(-0.015, 0.015), logp = FALSE, p = "P", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(names(models)[normal_idx[2]]," model weights"))
manhattan(models[[normal_idx[3]]], ylim = c(-0.015, 0.015), logp = FALSE, p = "P", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(names(models)[normal_idx[3]]," model weights"))
manhattan(models[[normal_idx[4]]], ylim = c(-0.015, 0.015), logp = FALSE, p = "P", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(names(models)[normal_idx[4]]," model weights"))

ggplot(stack( (og_white_all_tprs[(colnames(tprs_beta)[1+c(normal_idx[1], normal_idx[1]+n_ct)])]),
              varying = colnames(tprs_beta)[c(1+normal_idx[1], normal_idx[1]+n_ct)]), 
       aes(values)) +
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (og_white_all_tprs[(colnames(tprs_beta)[1+c(normal_idx[2], normal_idx[2]+n_ct)])]),
              varying = colnames(tprs_beta)[1+c(normal_idx[2], normal_idx[2]+n_ct)]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (og_white_all_tprs[(colnames(tprs_beta)[1+c(normal_idx[3], normal_idx[3]+n_ct)])]),
              varying = colnames(tprs_beta)[1+c(normal_idx[3], normal_idx[3]+n_ct)]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (og_white_all_tprs[(colnames(tprs_beta)[1+c(normal_idx[4], normal_idx[4]+n_ct)])]),
              varying = colnames(tprs_beta)[1+c(normal_idx[4], normal_idx[4]+n_ct)]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
manhattan(models[[mid_idx[1]]], ylim = c(-0.015, 0.015), logp = FALSE, p = "P", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(names(models)[mid_idx[1]]," model weights"))
manhattan(models[[mid_idx[2]]], ylim = c(-0.015, 0.015), logp = FALSE, p = "P", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(names(models)[mid_idx[2]]," model weights"))
manhattan(models[[mid_idx[3]]], ylim = c(-0.015, 0.015), logp = FALSE, p = "P", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(names(models)[mid_idx[3]]," model weights"))

ggplot(stack( (og_white_all_tprs[(colnames(tprs_beta)[1+c(mid_idx[1], mid_idx[1]+n_ct)])]),
              varying = colnames(tprs_beta)[c(1+mid_idx[1], mid_idx[1]+n_ct)]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (og_white_all_tprs[(colnames(tprs_beta)[1+c(mid_idx[2], mid_idx[2]+n_ct)])]),
              varying = colnames(tprs_beta)[1+c(mid_idx[2], mid_idx[2]+n_ct)]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (og_white_all_tprs[(colnames(tprs_beta)[1+c(mid_idx[3], mid_idx[3]+n_ct)])]),
              varying = colnames(tprs_beta)[1+c(mid_idx[3], mid_idx[3]+n_ct)]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
manhattan(models[[low_idx[1]]], ylim = c(-0.015, 0.015), logp = FALSE, p = "P", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(names(models)[low_idx[1]]," model weights"))
manhattan(models[[low_idx[2]]], ylim = c(-0.015, 0.015), logp = FALSE, p = "P", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(names(models)[low_idx[2]]," model weights"))
manhattan(models[[low_idx[3]]], ylim = c(-0.015, 0.015), logp = FALSE, p = "P", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(names(models)[low_idx[3]]," model weights"))

ggplot(stack( (og_white_all_tprs[(colnames(tprs_beta)[1+c(low_idx[1], low_idx[1]+n_ct)])]),
              varying = colnames(tprs_beta)[c(1+low_idx[1], low_idx[1]+n_ct)]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (og_white_all_tprs[(colnames(tprs_beta)[1+c(low_idx[2], low_idx[2]+n_ct)])]),
              varying = colnames(tprs_beta)[1+c(low_idx[2], low_idx[2]+n_ct)]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (og_white_all_tprs[(colnames(tprs_beta)[1+c(low_idx[3], low_idx[3]+n_ct)])]),
              varying = colnames(tprs_beta)[1+c(low_idx[3], low_idx[3]+n_ct)]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (scale_tprs[(colnames(sort_betas)[2:8 - 1 ])]), varying = colnames(sort_betas)[2:8 -1]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (scale_tprs[(colnames(sort_betas)[9:15 -1])]), varying = colnames(sort_betas)[9:15-1]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (scale_tprs[(colnames(sort_betas)[16:(n_ct+1)-1])]), varying = colnames(sort_betas)[16:(n_ct+1)-1]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))


tPRS_gene_diff <- function(filename, tprss){
  #filename <- "/nethome/kcni/mding/abcd/workspace/mding/abcd_imputed/with_sc_blood38/transitional_stage_B_cell"
  df <- fread(paste0(filename,"/abcd_pred_", 1, ".txt"), header = TRUE)
  df$IID <- mapply(sub, ".{9}_?", "", df$IID)
  overlap_genes <- intersect(witt_deg$ensembl_gene_id, colnames(df))
  df_sub <- df[(df$IID %in%  unique(tprss$IID)), ]
  df_sub <- df_sub[, c("FID", "IID", ..overlap_genes)]
  df_sub <- merge(df_sub, 
                  tprss[,c("IID", paste0(basename(filename), "_", "uw"), "PC1_uw")], 
                  by = c("IID"), sort = FALSE)
  g_diff <- vector()
  for (k in 3:(ncol(df_sub)-2)){
    a <- ncol(df_sub) - 1
    df_sub[,k] <- scale(df_sub[,..k])
    if (!is.na(df_sub[,..k][1,])){
      tt <- t.test(df_sub[(which(df_sub$PC1_uw > -0.02)),..k],
                   df_sub[(which(df_sub$PC1_uw < -0.02)),..k])
      if (!is.na(tt$p.value) & tt$p.value <= 0.05 ){
        g_diff[colnames(df_sub)[k]] <- tt$estimate[1] - tt$estimate[2]
      }
    }
  }
  return(g_diff)
}

ct_gene_diff_uw <- list()
for(c in cell_files){
  ct_gene_diff_uw[[paste0(basename(c), "_uw")]] <- tPRS_gene_diff(c, pc_df)
  
}

for (c in names(ct_gene_diff_uw)){
  ct_gene_diff_uw[[c]] <- ct_gene_diff_uw[[c]][order(abs(ct_gene_diff_uw[[c]]), decreasing = TRUE)]
}
lapply(ct_gene_diff_uw, head)
# install.packages("vroom")
# library(vroom)
vcf_id <- fread("/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_data/abcd_snps.txt", 
                header = TRUE, skip = 30)

wgt_snp_df <- data.frame(ct = c(), snp = c(), wgt = c(), CHR = c(), BP = c())
for(mod in names(models)){
  genes <- names(ct_gene_diff_uw[[paste0(mod,"_uw")]])
  idx <- which(models[[mod]]$gene %in% genes)
  idx2 <- which(models[[mod]]$SNP[idx] %in% vcf_id$ID)
  #print(idx2)
  df <- data.frame(ct = rep((paste0(mod,"_uw")), length(idx2)), 
                   SNP = models[[mod]]$SNP[idx][idx2], wgt = models[[mod]]$P[idx][idx2],
                   quant = ecdf(abs(models[[mod]]$P))(abs(models[[mod]]$P[idx][idx2])))
  model_chr <- as.numeric(str_remove(str_extract(df$SNP, "chr[0-9]+"), "chr"))
  model_pos <- as.numeric(str_remove(str_extract(df$SNP, "_[0-9]+"), "_"))
  df$CHR <- model_chr
  df$BP <- model_pos
  wgt_snp_df <- rbind(wgt_snp_df, df)
}

manhattan(wgt_snp_df, ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = "Aggregate model weights for SNPs that control expression\nlevel of influential genes")


balance_gene_expr <- unlist(lapply(ct_gene_diff_uw, sum))
balance_gene_expr <- balance_gene_expr[order(abs(balance_gene_expr), decreasing = TRUE)]
sort(table(wgt_snp_df$SNP),decreasing=TRUE)[1:52]
sort(table(names(unlist(unname(ct_gene_diff_uw)))),decreasing=TRUE)[1:59]
fileConn<-file("infl_genes.txt")
writeLines((names(unlist(unname(ct_gene_diff_uw)))), fileConn)
close(fileConn)
ct_labs <- unique(wgt_snp_df$ct)
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[1]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[1], " model weights for SNPs that control expression\nlevel of influential genes"),
          highlight = c("chr12_47902703_A_G_b38"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[2]),], ylim = c(-0.001, 0.001), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, 
          main = paste0(ct_labs[2], " model weights for SNPs that control expression\nlevel of influential genes"),
          highlight = c("chr20_35291403_C_T_b38"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[3]),], ylim = c(-0.0005, 0.0005), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[3], " model weights for SNPs that control expression\nlevel of influential genes"),
          highlight = c("chr3_14978842_A_G_b38"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[4]),], ylim = c(-0.001, 0.001), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[4], " model weights for SNPs that control expression\nlevel of influential genes"),
          highlight = c("chr15_89091211_A_G_b38", "chr20_35291675_A_G_b38", "chr11_75398824_T_C_b38"), xlim = c(98000000, 150000000))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[5]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[5], " model weights for SNPs that control expression\nlevel of influential genes"),
          highlight = c("chr7_90656015_G_A_b38") , xlim = c(140000000, 190000000))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[6]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[6], " model weights for SNPs that control expression\nlevel of influential genes"),
          highlight =c("chr12_47902703_A_G_b38"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[7]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[7], " model weights for SNPs that control expression\nlevel of influential genes"),
          highlight = c("chr12_47902703_A_G_b38", "chr7_90656015_G_A_b38"), xlim = c(194000000, 400000000))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[8]),], ylim = c(-0.001, 0.001), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[8], " model weights for SNPs that control expression\nlevel of influential genes"),
          xlim = c(220000000, 270000000))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[9]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[9], " model weights for SNPs that control expression\nlevel of influential genes"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[10]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[10], " model weights for SNPs that control expression\nlevel of influential genes"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[11]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[11], " model weights for SNPs that control expression\nlevel of influential genes"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[12]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[12], " model weights for SNPs that control expression\nlevel of influential genes"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[13]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[13], " model weights for SNPs that control expression\nlevel of influential genes"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[14]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[14], " model weights for SNPs that control expression\nlevel of influential genes"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[15]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[15], " model weights for SNPs that control expression\nlevel of influential genes"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[16]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[16], " model weights for SNPs that control expression\nlevel of influential genes"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[17]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[17], " model weights for SNPs that control expression\nlevel of influential genes"),
          highlight = c("chr12_47902703_A_G_b38"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[18]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[18], " model weights for SNPs that control expression\nlevel of influential genes"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[19]),], ylim = c(-0.001, 0.001), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[19], " model weights for SNPs that control expression\nlevel of influential genes"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[20]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[20], " model weights for SNPs that control expression\nlevel of influential genes"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[21]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[21], " model weights for SNPs that control expression\nlevel of influential genes"))
manhattan(wgt_snp_df[which(wgt_snp_df$ct == ct_labs[22]),], ylim = c(-0.01, 0.01), logp = FALSE, p = "wgt", ylab = "Weight", genomewideline = FALSE,
          suggestiveline = FALSE, main = paste0(ct_labs[22], " model weights for SNPs that control expression\nlevel of influential genes"),
          highlight = c("chr7_90656015_G_A_b38"))
iid <- data.frame(IID = scale_tprs[,2])
colnames(scale_tprs) <- gsub("-", "_", colnames(scale_tprs))
write.table(cbind(iid, scale_tprs)[(scale_tprs$event_name == events[1]),][, c(1,20:41)], 
            "/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_data/tprs_pheno.txt",
            sep="\t",row.names=FALSE, quote = F)

write.table(cbind(iid, scale_tprs)[(scale_tprs$event_name == events[1]),][, c(1,10:19)], 
            "/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_data/tprs_pheno_genpc.txt",
            sep="\t",row.names=FALSE, quote = F)

save(wgt_snp_df,file="sc-tprs-mdd/model_wgt_snp_df.Rda")
save(ct_gene_diff_uw,file="sc-tprs-mdd/ct_gene_diff_uw.Rda")

for (ct in names(ct_gene_diff_uw)){
  write.table(wgt_snp_df[which(wgt_snp_df$ct == ct), "SNP"], 
              paste0("/external/rprshnas01/kcni/mding/abcd/workspace/mding/immune_cell_models/snps/", ct,"_from_gene.txt"),
              sep="\t",row.names=FALSE, quote = F)
  
}

# run impute_many_expression_mod.sh

# recalc tPRSs with candidate SNPs removed=============================
tPRS_uw_nr <- tPRS_uw[,1:2]
tPRS_w_nr <- tPRS_w[,1:2]
rib_genes <- read.delim("/external/rprshnas01/kcni/mding/sc-tprs-mdd/ribo_genes_ensembl.txt", sep = " ")
witt_no_rib <- witt_deg$ensembl_gene_id[-which(witt_deg$ensembl_gene_id %in% rib_genes$ensembl_gene_id)]
for (f in cell_files[10:length(cell_files)]){
  a <- tPRS_calc(f, goi = witt_no_rib)
  tPRS_uw_nr <- cbind(tPRS_uw_nr, a[[1]])
  tPRS_w_nr <- cbind(tPRS_w_nr, a[[2]])
}


ggplot(stack( (tPRS_uw_nr[, c(2:8)]), varying = colnames(tPRS_uw_nr)[2:8]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (tPRS_uw_nr[, c(9:15)]), varying = colnames(tPRS_uw_nr)[9:15]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x') + theme_stata() +
  theme(text=element_text(size=9))
ggplot(stack( (tPRS_uw[, c(16:23)]), varying = colnames(tPRS_uw)[16:23]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x') + theme_stata()+
  theme(text=element_text(size=9))

ggplot(stack( (tPRS_uw_nr[, c(2:8)]), varying = colnames(cont)[2:8]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x')+ theme_stata() + 
  theme(text=element_text(size=9))
ggplot(stack( (tPRS_uw_nr[, c(9:15)]), varying = colnames(cont)[9:15]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x') + theme_stata() +
  theme(text=element_text(size=9))
ggplot(stack( (tPRS_uw_nr[, c(16:24)]), varying = colnames(cont)[16:23]), 
       aes(values)) + 
  geom_histogram(bins = 20) + 
  facet_wrap(~ind, scales = 'free_x') + theme_stata()+
  theme(text=element_text(size=9))



cell_files_mod <- list.dirs(path="~/abcd/workspace/mding/abcd_imputed/modified_models", 
                            full.names=TRUE, recursive=FALSE)
tPRS_uw <- data.frame(IID = unique(white_all$IID))
tPRS_uw[, cell_types] <- NA
tPRS_uw$IID <- unique(white_all$IID)
tPRS_w <- tPRS_uw




a <- tPRS_calc("/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_imputed/modified_models/chr12_47902703_A_G_b38conventional_dendritic_cell")
tPRS_uw <- a[[1]]
tPRS_w <- a[[2]]
a <- tPRS_calc("/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_imputed/modified_models/chr12_47902703_A_G_b38dendritic_cell")
tPRS_uw <- a[[1]]
tPRS_w <- a[[2]]
ggplot(tPRS_uw, 
       aes(chr12_47902703_A_G_b38dendritic_cell)) + 
  geom_histogram(bins = 20) + theme_stata() + 
  theme(text=element_text(size=9))

a <- tPRS_calc("/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_imputed/modified_models/chr7_90656015_G_A_b38transitional_stage_B_cell")
tPRS_uw <- a[[1]]
tPRS_w <- a[[2]]
ggplot(tPRS_uw, 
       aes(chr7_90656015_G_A_b38transitional_stage_B_cell)) + 
  geom_histogram(bins = 20) + theme_stata() + 
  theme(text=element_text(size=9))
a <- tPRS_calc("/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_imputed/modified_models/chr15_89091211_A_G_b38_chr20_35291675_A_G_b38_chr11_75398824_T_C_b38CD8-positive_alpha-beta_T_cell")
tPRS_uw <- a[[1]]
tPRS_w <- a[[2]]
ggplot(tPRS_uw, 
       aes(`chr15_89091211_A_G_b38_chr20_35291675_A_G_b38_chr11_75398824_T_C_b38CD8-positive_alpha-beta_T_cell`)) + 
  geom_histogram(bins = 20) + theme_stata() + 
  theme(text=element_text(size=9))
# do plink gwas
a <- tPRS_calc("/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_imputed/modified_models/chr10_86964703_G_A_b38dendritic_cell")
tPRS_uw <- a[[1]]
tPRS_w <- a[[2]]
ggplot(tPRS_uw, 
       aes(chr10_86964703_G_A_b38dendritic_cell)) + 
  geom_histogram(bins = 20) + theme_stata() + 
  theme(text=element_text(size=9))

a <- tPRS_calc("/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_imputed/modified_models/chr10_86964703_G_A_b38plasmacytoid_dendritic_cell")
tPRS_uw <- a[[1]]
tPRS_w <- a[[2]]
ggplot(tPRS_uw, 
       aes(chr10_86964703_G_A_b38plasmacytoid_dendritic_cell)) + 
  geom_histogram(bins = 20) + theme_stata() + 
  theme(text=element_text(size=9))

# remove all snps from gwas common with infl genes
a <- tPRS_calc("/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_imputed/modified_models/chr16_4085757_T_C_b38_chr19_19610913_G_A_b38CD4-positive_alpha-beta_cytotoxic_T_cell")
tPRS_uw <- a[[1]]
tPRS_w <- a[[2]]
ggplot(tPRS_uw, 
       aes(`chr16_4085757_T_C_b38_chr19_19610913_G_A_b38CD4-positive_alpha-beta_cytotoxic_T_cell`)) + 
  geom_histogram(bins = 20) + theme_stata() + 
  theme(text=element_text(size=9))

a <- tPRS_calc("/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_imputed/modified_models/chr19_19610913_G_A_b38_chr20_7104633_T_C_b38CD8-positive_alpha-beta_T_cell")
tPRS_uw <- a[[1]]
tPRS_w <- a[[2]]
ggplot(tPRS_uw, 
       aes(`chr19_19610913_G_A_b38_chr20_7104633_T_C_b38CD8-positive_alpha-beta_T_cell`)) + 
  geom_histogram(bins = 20) + theme_stata() + 
  theme(text=element_text(size=9))

a <- tPRS_calc("/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_imputed/modified_models/chr7_91127367_C_T_b38central_memory_CD8-positive_alpha-beta_T_cell")
tPRS_uw <- a[[1]]
tPRS_w <- a[[2]]
ggplot(tPRS_uw, 
       aes(`chr7_91127367_C_T_b38central_memory_CD8-positive_alpha-beta_T_cell`)) + 
  geom_histogram(bins = 20) + theme_stata() + 
  theme(text=element_text(size=9))

tPRS_calc("/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_imputed/modified_models/chr10_86964703_G_A_b38_chr15_54307877_C_T_b38_chr20_35283507_T_C_b38_chr21_38474772_G_A_b38_chr21_38506662_C_T_b38dendritic_cell")
ggplot(tPRS_uw, 
       aes(chr10_86964703_G_A_b38_chr15_54307877_C_T_b38_chr20_35283507_T_C_b38_chr21_38474772_G_A_b38_chr21_38506662_C_T_b38dendritic_cell)) + 
  geom_histogram(bins = 20) + theme_stata() + 
  theme(text=element_text(size=9))

tPRS_calc("/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_imputed/modified_models/chr10_15762487_G_A_b38_chr19_19610913_G_A_b38_chr20_6773599_G_A_b38effector_memory_CD4-positive_alpha-beta_T_cell")
ggplot(tPRS_uw, 
       aes(`chr10_15762487_G_A_b38_chr19_19610913_G_A_b38_chr20_6773599_G_A_b38effector_memory_CD4-positive_alpha-beta_T_cell`)) + 
  geom_histogram(bins = 20) + theme_stata() + 
  theme(text=element_text(size=9))

tPRS_calc("/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_imputed/modified_models/chr1_181698147_A_G_b38_chr6_17979622_G_A_b38_chr6_38198735_G_A_b38_chr20_6773599_G_A_b38erythrocyte")
ggplot(tPRS_uw, 
       aes(chr1_181698147_A_G_b38_chr6_17979622_G_A_b38_chr6_38198735_G_A_b38_chr20_6773599_G_A_b38erythrocyte)) + 
  geom_histogram(bins = 20) + theme_stata() + 
  theme(text=element_text(size=9))

tPRS_calc("/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_imputed/modified_models/chr4_80235944_T_C_b38_chr6_17406677_C_T_b38_chr9_76621037_G_A_b38_chr9_77455995_T_C_b38_chr19_19610913_G_A_b38_chr20_35284920_C_T_b38memory_B_cell")
ggplot(tPRS_uw, 
       aes(chr4_80235944_T_C_b38_chr6_17406677_C_T_b38_chr9_76621037_G_A_b38_chr9_77455995_T_C_b38_chr19_19610913_G_A_b38_chr20_35284920_C_T_b38memory_B_cell)) + 
  geom_histogram(bins = 20) + theme_stata() + 
  theme(text=element_text(size=9))

tPRS_calc("/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_imputed/modified_models/chr4_80235944_T_C_b38_chr19_19610913_G_A_b38naive_B_cell")
ggplot(tPRS_uw, 
       aes(chr4_80235944_T_C_b38_chr19_19610913_G_A_b38naive_B_cell)) + 
  geom_histogram(bins = 20) + theme_stata() + 
  theme(text=element_text(size=9))

tPRS_calc("/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_imputed/modified_models/chr6_37622174_T_C_b38_chr6_38198735_G_A_b38_chr14_75105562_G_A_b38natural_killer_cell")
ggplot(tPRS_uw, 
       aes(chr6_37622174_T_C_b38_chr6_38198735_G_A_b38_chr14_75105562_G_A_b38natural_killer_cell)) + 
  geom_histogram(bins = 20) + theme_stata() + 
  theme(text=element_text(size=9))

tPRS_calc("/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_imputed/modified_models/chr1_181698147_A_G_b38_chr19_19610913_G_A_b38plasmablast")
ggplot(tPRS_uw, 
       aes(chr1_181698147_A_G_b38_chr19_19610913_G_A_b38plasmablast)) + 
  geom_histogram(bins = 20) + theme_stata() + 
  theme(text=element_text(size=9))

tPRS_calc("/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_imputed/modified_models/chr10_15762487_G_A_b38_chr10_86964703_G_A_b38plasmacytoid_dendritic_cell")
ggplot(tPRS_uw, 
       aes(chr10_15762487_G_A_b38_chr10_86964703_G_A_b38plasmacytoid_dendritic_cell)) + 
  geom_histogram(bins = 20) + theme_stata() + 
  theme(text=element_text(size=9))

tPRS_calc("/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_imputed/modified_models/chr9_34388113_T_C_b38_chr20_36180056_T_G_b38platelet")
ggplot(tPRS_uw, 
       aes(chr9_34388113_T_C_b38_chr20_36180056_T_G_b38platelet)) + 
  geom_histogram(bins = 20) + theme_stata() + 
  theme(text=element_text(size=9))

tPRS_calc("/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_imputed/modified_models/chr3_14447862_G_T_b38_chr19_19610913_G_A_b38regulatory_T_cell")
ggplot(tPRS_uw, 
       aes(chr3_14447862_G_T_b38_chr19_19610913_G_A_b38regulatory_T_cell)) + 
  geom_histogram(bins = 20) + theme_stata() + 
  theme(text=element_text(size=9))

tPRS_calc("/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_imputed/modified_models/chr4_80235944_T_C_b38_chr19_19610913_G_A_b38transitional_stage_B_cell")
ggplot(tPRS_uw, 
       aes(chr4_80235944_T_C_b38_chr19_19610913_G_A_b38transitional_stage_B_cell)) + 
  geom_histogram(bins = 20) + theme_stata() + 
  theme(text=element_text(size=9))

# remove all snps from gwas for bimodal high/low categorization ======================

# run tprs_gwas.sh
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



tPRS_calc("/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_imputed/modified_models/chr11_64702052_G_A_b38transitional_stage_B_cell")
ggplot(tPRS_uw, 
       aes(chr11_64702052_G_A_b38transitional_stage_B_cell)) + 
  geom_histogram(bins = 20) + theme_stata() + 
  theme(text=element_text(size=9))