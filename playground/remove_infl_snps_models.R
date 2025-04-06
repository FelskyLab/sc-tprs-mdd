library(RSQLite)
if (!requireNamespace("DBI", quietly = TRUE)){
  install.packages("DBI")
}
library(DBI)
get_snps <- function(filename){
  sqlite.driver <- dbDriver("SQLite")
  db <- dbConnect(sqlite.driver, dbname = filename)
  model <- dbReadTable(db,"weights")
  dbDisconnect(db)
  
  out_filename <- paste0(dirname(filename),"/snps/",
                         sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(filename)), ".txt")
  print(out_filename)
  write.table(model$rsid, out_filename, sep = "\t", row.names = FALSE, quote = FALSE)
  
  #rm(snps, model)
  
}
lapply(model_files[1:22], get_snps)

remove_snps <- function(filename, snps){
  sqlite.driver <- dbDriver("SQLite")
  db <- dbConnect(sqlite.driver, dbname = filename)
  extra <- dbReadTable(db,"extra")
  model <- dbReadTable(db,"weights")
  dbDisconnect(db)
  print(model[which(model$rsid %in% snps),])
  model <- model[-which(model$rsid %in% snps),]
  
  out_filename <- paste0(dirname(filename),"/modified_models/",
                         paste(snps, sep = '', collapse = '_'), basename(filename))
  con <- dbConnect(SQLite(), out_filename)
  dbWriteTable(con, "weights", model)
  dbWriteTable(con, "extra", extra)
  dbDisconnect(con)
  
  #rm(snps, model)
  
}

get_removed_genes <- function(filename, snps){
  sqlite.driver <- dbDriver("SQLite")
  db <- dbConnect(sqlite.driver, dbname = filename)
  extra <- dbReadTable(db,"extra")
  model <- dbReadTable(db,"weights")
  dbDisconnect(db)
  x <- model[which(model$rsid %in% snps),]
  cell <- paste0(gsub(pattern = "-", 
                   replacement = "_", sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(f))), "_uw")
 return(cbind(x, data.frame(ct = rep(cell, nrow(x)))))
}

#load in *.db file
#filename <- "dendritic_cell.db"
#filename <- commandArgs(trailingOnly = TRUE)
#model_files <- list.files(path=commandArgs(trailingOnly = TRUE), pattern="*.db", full.names=TRUE, recursive=FALSE)
model_files <- list.files(path="abcd/workspace/mding/immune_cell_models", pattern="*.db", full.names=TRUE, recursive=FALSE)
# lapply(model_files[1:22], remove_snps)
remove_snps(model_files[22], c("chr7_90656015_G_A_b38"))
remove_snps(model_files[6], c("chr12_47902703_A_G_b38"))
remove_snps(model_files[17], c("chr12_47902703_A_G_b38"))
remove_snps(model_files[4], 
            c("chr15_89091211_A_G_b38", "chr20_35291675_A_G_b38", "chr11_75398824_T_C_b38"))
remove_snps(model_files[4], 
            c("chr15_89091211_A_G_b38"))
remove_snps(model_files[7], c("chr12_47902703_A_G_b38"))
remove_snps(model_files[7], c("chr7_90656015_G_A_b38"))
remove_snps(model_files[5], c("chr7_90656015_G_A_b38"))
for (i in c(9,11,14,15,18,19,20,21,22)){
  remove_snps(model_files[i], c("chr11_133146316_-_C_b38"))
  
}

remove_snps(model_files[7], c("chr10_86964703_G_A_b38"))
remove_snps(model_files[19], c("chr10_86964703_G_A_b38"))

remove_snps(model_files[2], c("chr19_19610913_G_A_b38"))
remove_snps(model_files[9], c("chr19_19610913_G_A_b38"))
remove_snps(model_files[14], c("chr19_19610913_G_A_b38"))
remove_snps(model_files[15], c("chr19_19610913_G_A_b38"))
remove_snps(model_files[22], c("chr19_19610913_G_A_b38"))

for (f in model_files[bimod_ct_idx]){
  x <- all_diff[
    which(all_diff$ct == paste0(gsub(pattern = "-", 
            replacement = "_", sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(f))), "_uw") & 
            all_diff$BETA < 0),"ID"]
  remove_snps(f, x)
  
}

rem_genes <- data.frame()
for (f in model_files[bimod_ct_idx]){
  x <- all_diff[
    which(all_diff$ct == paste0(gsub(pattern = "-", 
                                     replacement = "_", sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(f))), "_uw") & 
            all_diff$BETA < 0),"ID"]
  rem_genes <- rbind(rem_genes, get_removed_genes(f, x))
  
}
write.table(rem_genes$gene, "/external/rprshnas01/kcni/mding/abcd/workspace/mding/abcd_data/rem_genes.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE)