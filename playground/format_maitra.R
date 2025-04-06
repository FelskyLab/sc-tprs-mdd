setwd("/external/rprshnas01/kcni/mding/sc-tprs-mdd")

library(readxl)
cnames <- as.list(read_excel("sc-tprs-mdd/degs_567.xlsx", sheet = "SupplementaryData5", 
                     range = cell_rows(3)))
degm <- read_excel("sc-tprs-mdd/degs_567.xlsx", sheet = "SupplementaryData5", 
                   range = cell_rows(c(3,NA)))
degf <- read_excel("sc-tprs-mdd/degs_567.xlsx", sheet = "SupplementaryData6", 
                   range = cell_rows(c(3,NA)))

library("AnnotationDbi")
library("org.Hs.eg.db")
degm$ensembl_gene_id <- mapIds(org.Hs.eg.db,
                  keys=degm$gene, 
                  column="ENSEMBL",
                  keytype="SYMBOL",
                  multiVals="first")
degf$ensembl_gene_id <- mapIds(org.Hs.eg.db,
                               keys=degf$gene, 
                               column="ENSEMBL",
                               keytype="SYMBOL",
                               multiVals="first")

# remove.packages("tidyselect")
# install.packages("tidyselect")
library(tidyverse)
library(biomaRt)
#install("BiocManager")
library(BiocManager)
#install("BiocFileCache", force = TRUE)
# install.packages("devtools")\
# remove.packages("dbplyr")
# devtools::install_version("dbplyr", version = "2.3.4")

mart <- useMart('ENSEMBL_MART_ENSEMBL', host = 'useast.ensembl.org')
mart <- useDataset('hsapiens_gene_ensembl', mart)

annotLookup <- getBM(mart = mart, attributes = c('hgnc_symbol', 
                                                 'ensembl_gene_id'), 
                     uniqueRows = TRUE)

degm$ensembl_gene_id2 <- annotLookup$ensembl_gene_id[match(unlist(degm$gene), 
                                                           annotLookup$hgnc_symbol)]
degf$ensembl_gene_id2 <- annotLookup$ensembl_gene_id[match(unlist(degf$gene), 
                                                           annotLookup$hgnc_symbol)]

degf$ensembl_gene_id[is.na(degf$ensembl_gene_id)] <- 
  degf$ensembl_gene_id2[is.na(degf$ensembl_gene_id)] 
degm$ensembl_gene_id[is.na(degm$ensembl_gene_id)] <- 
  degm$ensembl_gene_id2[is.na(degm$ensembl_gene_id)] 
degf <- subset(degf, select = -c(ensembl_gene_id2))
degm <- subset(degm, select = -c(ensembl_gene_id2))
degf$gene[is.na(degf$ensembl_gene_id)]
degm$gene[is.na(degm$ensembl_gene_id)]

# TODO manually search up NA ensembl id's
#manualf <- c("ENSG00000265737", "")


library(openxlsx)
cell_types <- c("InN", "ExN", "Oli", "OPC", "Ast", "End", "Mic")

degf_ct <- list()
for (ct in cell_types){
  degf_ct[[ct]] <- degf[grep(degf$cluster_id, pattern = ct),]
  write.xlsx(file = paste0("deg_by_sex_ct/deg_f_", ct, ".xlsx"), x = degf_ct[[ct]])
}

degm_ct <- list()
for (ct in cell_types){
  degm_ct[[ct]] <- degm[grep(degm$cluster_id, pattern = ct),]
  write.xlsx(file = paste0("deg_by_sex_ct/deg_m_", ct, ".xlsx"), x = degm_ct[[ct]])
}

