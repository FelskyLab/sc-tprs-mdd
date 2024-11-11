setwd("C:/Users/mding/Desktop/BCB330/sc-tprs-mdd")
#install.packages("RSQLite")
library(RSQLite)

filename <- "transitional_stage_B_cell.db"
sqlite.driver <- dbDriver("SQLite")
db <- dbConnect(sqlite.driver,
                dbname = filename)

## Some operations
dbListTables(db)
trans_b_model <- dbReadTable(db,"weights")
dbDisconnect()

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.20")
#
# BiocManager::install(c("GenomicFeatures", "AnnotationDbi",
#                        "BSgenome.Hsapiens.UCSC.hg38", "SNPlocs.Hsapiens.dbSNP155.GRCh38"))

library(BSgenome.Hsapiens.UCSC.hg38)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
genome <- BSgenome.Hsapiens.UCSC.hg38
seqlevelsStyle(genome) <- "NCBI"
snps <- GRanges(paste0(str_remove(str_extract(trans_b_model$varID, "chr[0-9]+"), "chr"), ":",
                   str_remove(str_extract(trans_b_model$varID, "_[0-9]+"), "_"), "-",
                   as.numeric(str_remove(str_extract(trans_b_model$varID, "_[0-9]+"), "_")) + 1))

snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38, snps, genome=genome)

library(biomaRt)
## Use the default ENSEMBL Variation Mart & Human dataset
snpMart <- useEnsembl(biomart = "snps",
                     dataset = "hsapiens_snp")

## Create an example set of coordinates as a dataframe
SNP_M <- data.frame(CHR = c(1,1), START = c(10020, 10039), END = c(10020, 10039))

## Combine these into the format chr:start:end
## It's important to include the end even if it's a single base,
## otherwise it searches to the end of the chromosome
coords <- apply(SNP_M, 1, paste, collapse = ":")
coords
#> [1] "1:10020:10020" "1:10039:10039"

## Submit the query
getBM(attributes = c('refsnp_id', 'chr_name', 'chrom_start', 'chrom_end', 'allele'),
      filters = c('chromosomal_region'),
      values = coords,
      mart = snpMart)
