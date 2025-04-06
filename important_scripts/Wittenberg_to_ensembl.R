#!/external/rprshnas01/kcni/mding/R-4.4.2/bin/Rscript

# goal: turn gene names from Wittenberg et al. 2020 to ensembl format
# next: load_abcd_data.R

setwd("/external/rprshnas01/kcni/mding/sc-tprs-mdd")

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install(version = "3.20", ask = FALSE)
}
library(BiocManager)
if (!requireNamespace("readxl", quietly = TRUE)){
  BiocManager::install("readxl")
}
library(readxl)
cnames <- as.list(read_excel("sc-tprs-mdd/gene_lists/Wittenberg2020.xlsx", sheet = "S9 SMD Meta-Analysis", 
                     range = cell_rows(1)))
deg <- read_excel("sc-tprs-mdd/gene_lists/Wittenberg2020.xlsx", sheet = "S9 SMD Meta-Analysis", 
                   range = cell_rows(c(3,NA)))

if (!requireNamespace("AnnotationDbi", quietly = TRUE)){
  BiocManager::install("AnnotationDbi")
}
library("AnnotationDbi")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)){
  BiocManager::install("org.Hs.eg.db")
}
library("org.Hs.eg.db")
deg$ensembl_gene_id <- mapIds(org.Hs.eg.db,
                  keys=deg$gene, 
                  column="ENSEMBL",
                  keytype="SYMBOL",
                  multiVals="first")

if (!requireNamespace("tidyverse", quietly = TRUE)){
  BiocManager::install("tidyverse")
}
library(tidyverse)
if (!requireNamespace("biomaRt", quietly = TRUE)){
  BiocManager::install("biomaRt")
}
library(biomaRt)

#install("BiocFileCache", force = TRUE)
# install.packages("devtools")\
# remove.packages("dbplyr")
# devtools::install_version("dbplyr", version = "2.3.4")

mart <- useMart('ENSEMBL_MART_ENSEMBL', host = 'useast.ensembl.org')
mart <- useDataset('hsapiens_gene_ensembl', mart)

annotLookup <- getBM(mart = mart, attributes = c('hgnc_symbol', 
                                                 'ensembl_gene_id'), 
                     uniqueRows = TRUE)

deg$ensembl_gene_id2 <- annotLookup$ensembl_gene_id[match(unlist(deg$gene), 
                                                           annotLookup$hgnc_symbol)]

deg$ensembl_gene_id[is.na(deg$ensembl_gene_id)] <- 
  deg$ensembl_gene_id2[is.na(deg$ensembl_gene_id)] 
deg <- subset(deg, select = -c(ensembl_gene_id2))
deg$gene[is.na(deg$ensembl_gene_id)]

# TODO manually search up NA ensembl id's
#manualf <- c("ENSG00000265737", "")

write.xlsx(file = "gene_lists/Wittenberg_deg.xlsx", x = deg)
