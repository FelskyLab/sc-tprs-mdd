#!/external/rprshnas01/kcni/mding/R-4.4.2/bin/Rscript
#setwd("abcd/workspace/mding")
.libPaths( c("/external/rprshnas01/kcni/mding/R-4.4.2/library", .libPaths()) )
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install(version = "3.20", ask = FALSE)
}

if (!requireNamespace("GenomicFeatures", quietly = TRUE)){
  BiocManager::install("GenomicFeatures", ask = FALSE)
}
if (!requireNamespace("AnnotationHub", quietly = TRUE)){
  BiocManager::install("AnnotationHub", ask = FALSE)
}
if (!requireNamespace("BSgenome", quietly = TRUE)){
  BiocManager::install("BSgenome", ask = FALSE)
}
library(BSgenome)
if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)){
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", ask = FALSE)
}
# if (!requireNamespace("SNPlocs", quietly = TRUE)){
#   BiocManager::install("SNPlocs", ask = FALSE)
# }
# library("SNPlocs")
if (!requireNamespace("SNPlocs.Hsapiens.dbSNP155.GRCh38", quietly = TRUE)){
  BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38", ask = FALSE)
}
library(BSgenome.Hsapiens.UCSC.hg38)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)

if (!requireNamespace("stringr", quietly = TRUE)){
  install.packages("stringr")
}
library(stringr)
if (!requireNamespace("stringi", quietly = TRUE)){
  install.packages("stringi")
}
library(stringi)
if (!requireNamespace("AcidGenomes", quietly = TRUE)){
  install.packages(pkgs = "AcidGenomes",
                   repos = c("https://r.acidgenomics.com",
                             BiocManager::repositories()),
                   dependencies = TRUE)
}
library(AcidGenomes)
if (!requireNamespace("RSQLite", quietly = TRUE)){
  install.packages("RSQLite")
}
library(RSQLite)
if (!requireNamespace("DBI", quietly = TRUE)){
  install.packages("DBI")
}
library(DBI)
#some functions

ranges2label <- function(ranges, genes, label) {
  
  gnames <- findOverlaps(ranges, genes)
  gnames <- as(gnames, "List")
  gnames <- extractList(genes@elementMetadata@listData[[label]], gnames)
  gnames <- unique(gnames)
  gnames <- unstrsplit(gnames, ";")
  Rle(gnames)
}

#library(AnnotationHub)
#BiocManager::install("ensembldb")
#library(ensembldb)
#reference Ensembl gene names
# gff <- makeGRangesFromEnsembl(organism = "Homo sapiens", level = "genes",
#                               genomeBuild = "GRCh38")

pos_2_rsid <- function(filename){
  sqlite.driver <- dbDriver("SQLite")
  db <- dbConnect(sqlite.driver, dbname = filename)
  
  model <- dbReadTable(db,"weights")
  dbDisconnect(db)
  
  # #investigate the duplicate values
  # dups <- model[(duplicated(model$varID) | duplicated(model$varID, fromLast=TRUE)),]
  # dups <- dups[(order(dups$varID)),]
  
  #format chr and position
  model_chr <- str_remove(str_extract(model$varID, "chr[0-9]+"), "chr")
  model_pos <- str_remove(str_extract(model$varID, "_[0-9]+"), "_")
  
  coords <- paste0(model_chr, ":",
                   model_pos, "-",
                   model_pos)
  
  
  #create object with model data to query with
  snps <- GRanges(coords,
                  alleles_as_ambig = mergeIUPACLetters(paste0(model$ref_allele, model$eff_allele)),
                  gene = model$gene, varID = model$varID,
                  ref_allele = model$ref_allele, eff_allele = model$eff_allele,
                  weight = model$weight)
  #align snps with Ensembl genes
  #snps$geneSymbols <- ranges2label(snps, gff, "geneId")
  #
  # #Find indices where the gene doesn't match the chr:position
  # N <- length(snps@elementMetadata@listData[["gene"]])
  # matching_mask <- vector(length = N)
  #
  # #longest step, try to optimize
  # matching_mask <- mapply(function(snp, gene) +(stri_detect_fixed(snp, gene)),
  #                         as.vector(snps$geneSymbols), snps@elementMetadata@listData[["gene"]])
  #
  # #for (i in 1:N){
  # #  pBar(i, N)
  # #  matching_mask[i] <- grepl(snps@elementMetadata@listData[["gene"]][i], snps$geneSymbols[i], fixed = TRUE)
  #
  # #}
  # sum(matching_mask)/length(snps)
  # #remove these incongruent rows
  # msk_snps <- snps[(which(matching_mask==1)),]
  # msk_model <- model[matching_mask,]
  
  #genome with rsid annotations
  bsgenome <- BSgenome.Hsapiens.UCSC.hg38
  seqlevelsStyle(bsgenome) <- "NCBI"
  #get the rsids matching snps we have. Takes a loooonnnnngggggg time
  rsid <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38,
                         snps, genome=bsgenome)
  
  # rsid_names <- findOverlaps(msk_snps, rsid)
  # rsid_names <- as(rsid_names, "List")
  # rsid_names <- extractList(rsid@elementMetadata@listData[["RefSNP_id"]], rsid_names)
  # rsid_names <- unique(rsid_names)
  # rsid_names <- unstrsplit(rsid_names, ";")
  # msk_snps$rsid <- Rle(rsid_names)
  snps$rsid <- ranges2label(snps, rsid, "RefSNP_id")
  
  new_model <- as.data.frame(snps@elementMetadata@listData)
  new_model <- new_model[,-c(1, 7)]
  
  out_filename <- paste0(dirname(filename),"/modified_models/mod_",basename(filename))
  con <- dbConnect(SQLite(), out_filename)
  dbWriteTable(con, "weights", new_model)
  dbDisconnect(con)
  
  rm(rsid, snps, model)
  
}

#load in *.db file
#filename <- "dendritic_cell.db"
#filename <- commandArgs(trailingOnly = TRUE)
#model_files <- list.files(path=commandArgs(trailingOnly = TRUE), pattern="*.db", full.names=TRUE, recursive=FALSE)
model_files <- list.files(path="abcd/workspace/mding/immune_cell_models", pattern="*.db", full.names=TRUE, recursive=FALSE)
lapply(model_files[1:29], pos_2_rsid)

add_extra_table <- function(filename_og, filename_mod){
  sqlite.driver <- dbDriver("SQLite")
  db <- dbConnect(sqlite.driver, dbname = filename_og)
  
  extra <- dbReadTable(db,"extra")
  dbDisconnect(db)
  
  con <- dbConnect(SQLite(), filename_mod)
  dbWriteTable(con, "extra", extra)
  dbDisconnect(con)
}
mod_model_files <- list.files(path="abcd/workspace/mding/immune_cell_models/modified_models", pattern="*.db", full.names=TRUE, recursive=FALSE)
mapply(add_extra_table, model_files, mod_model_files)
