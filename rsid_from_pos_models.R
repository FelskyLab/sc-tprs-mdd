#setwd("C:/Users/mding/Desktop/BCB330/sc-tprs-mdd")
if (!require("BiocManager", quietly = TRUE)){
   install.packages("BiocManager")
}
BiocManager::install(version = "3.20", ask = FALSE)
BiocManager::install(c("GenomicFeatures", "AnnotationHub",
                       "BSgenome.Hsapiens.UCSC.hg38",
                       "SNPlocs.Hsapiens.dbSNP155.GRCh38"), force = FALSE,
                     ask = FALSE)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)

if (!require("stringr", quietly = TRUE)){
  install.packages("stringr")
}
library(stringr)
if (!require("stringi", quietly = TRUE)){
  install.packages("stringi")
}
library(stringr)
if (!require("AcidGenomes", quietly = TRUE)){
  install.packages(pkgs = "AcidGenomes",
                   repos = c("https://r.acidgenomics.com",
                             BiocManager::repositories()),
                  dependencies = TRUE)
}
library(AcidGenomes)
if (!require("RSQLite", quietly = TRUE)){
  install.packages("RSQLite")
}
library(RSQLite)

#some functions
cat("  Defining pBar() ...\n")
pBar <- function(i, l, nCh = 50) {
  # Draw a progress bar in the console
  # i: the current iteration
  # l: the total number of iterations
  # nCh: width of the progress bar
  ticks <- round(seq(1, l-1, length.out = nCh))
  if (i < l) {
    if (any(i == ticks)) {
      p <- which(i == ticks)[1]  # use only first, in case there are ties
      p1 <- paste(rep("#", p), collapse = "")
      p2 <- paste(rep("-", nCh - p), collapse = "")
      cat(sprintf("\r|%s%s|", p1, p2))
      flush.console()
    }
  }
  else { # done
    cat("\n")
  }
}

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



#load in *.db file
#filename <- "dendritic_cell.db"
filename <- commandArgs(trailingOnly = TRUE)
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

#reference Ensembl gene names
gff <- makeGRangesFromEnsembl(organism = "Homo sapiens", level = "genes",
                              genomeBuild = "GRCh38")
#create object with model data to query with
snps <- GRanges(coords,
  alleles_as_ambig = mergeIUPACLetters(paste0(model$ref_allele, model$eff_allele)),
  gene = model$gene, varID = model$varID,
  ref_allele = model$ref_allele, eff_allele = model$eff_allele,
  weight = model$weight)
#align snps with Ensembl genes
snps$geneSymbols <- ranges2label(snps, gff, "geneId")

#Find indices where the gene doesn't match the chr:position
N <- length(snps@elementMetadata@listData[["gene"]])
matching_mask <- vector(length = N)

#longest step, try to optimize
matching_mask <- mapply(function(snp, gene) +(stri_detect_fixed(snp, gene)),
                        as.vector(snps$geneSymbols), snps@elementMetadata@listData[["gene"]])

#for (i in 1:N){
#  pBar(i, N)
#  matching_mask[i] <- grepl(snps@elementMetadata@listData[["gene"]][i], snps$geneSymbols[i], fixed = TRUE)

#}
sum(matching_mask)/length(snps)
#remove these incongruent rows
msk_snps <- snps[matching_mask,]
msk_model <- model[matching_mask,]

#genome with rsid annotations
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
seqlevelsStyle(bsgenome) <- "NCBI"
#get the rsids matching snps we have
rsid <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38,
                       msk_snps, genome=bsgenome)

# rsid_names <- findOverlaps(msk_snps, rsid)
# rsid_names <- as(rsid_names, "List")
# rsid_names <- extractList(rsid@elementMetadata@listData[["RefSNP_id"]], rsid_names)
# rsid_names <- unique(rsid_names)
# rsid_names <- unstrsplit(rsid_names, ";")
# msk_snps$rsid <- Rle(rsid_names)
msk_snps$rsid <- ranges2label(msk_snps, rsid, "RefSNP_id")

new_model <- as.data.frame(msk_snps@elementMetadata@listData)
new_model <- new_model[,-c(1, 7)]

out_filename <- paste0(dirname(filename),"/modified_models/mod_",basename(filename))
con <- dbConnect(SQLite(), out_filename)
dbWriteTable(con, "weights", df)
dbDisconnect(con)
# rsid_list <- list(chr = rep(as.numeric(rsid@seqnames@values), times = rsid@seqnames@lengths),
#                       pos = rsid@ranges@pos,
#                       ref_allele = rsid@elementMetadata@listData[["ref_allele"]],
#                       alt_allele = as.vector(rsid@elementMetadata@listData[["alt_alleles"]]),
#                       rsid = rsid@elementMetadata@listData[["RefSNP_id"]])
# rsid_table <- new.env()
# len <- length(rsid_list[[1]])
# for (i in 1:len){
#   for (j in 1:length(rsid_list[[4]][i])){
#     str <- paste0(rsid_list[[1]][i], rsid_list[[2]][i],
#                   rsid_list[[3]][i],
#                   rsid_list[[4]][[i]][j])
#     if (exists(str, rsid_table)){
#       rsid_table[[str]] <- c(rsid_table[[str]], rsid_list[[5]][i])
#       next
#     }
#     rsid_table[[str]] <-rsid_list[[5]][i]
#
#   }
#
# }
# len <- nrow(msk_model)
# to_remove <- c()
# for (i in 1:len){
#   str <- paste0(msk_model$chr[i],
#                 msk_model$pos[i],
#                 msk_model$ref_allele[i],
#                 msk_model$eff_allele[i])
#   value <- rsid_table[[str]]
#   if(is.null(value)){
#     to_remove <- c(to_remove, i)
#     next
#   }
#   if (length(value)>1){
#     msk_model$rsid[i] <- value[1]
#     for(j in 2:length(value)){
#       new_row <- msk_model[i,]
#       new_row$rsid <- value[j]
#       msk_model[(nrow(msk_model)+1),] <- new_row
#
#     }
#     next
#   }
#   msk_model$rsid[i] <- value
# }
# msk_model[(1:100),]
# rem_msk_model<- msk_model[-to_remove,]
#
# for(j in 2:1){
#   print(j)
# }

# library(biomaRt)
# ## Use the default ENSEMBL Variation Mart & Human dataset
# snpMart <- useEnsembl(biomart = "snps",
#                      dataset = "hsapiens_snp")
#
# ## Create an example set of coordinates as a dataframe
# SNP_M <- data.frame(CHR = c(1,1), START = c(10020, 10039), END = c(10020, 10039))
#
#
#
# ## Submit the query
# getBM(attributes = c('refsnp_id', 'chr_name', 'chrom_start', 'chrom_end', 'allele'),
#       filters = c('chromosomal_region'),
#       values = coords[10:100],
#       mart = snpMart)

install.packages("installr")
library(installr)
updateR()