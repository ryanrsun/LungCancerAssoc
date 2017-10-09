#' Hg19 Gene table using mostly Ensembl annotations
#'
#' See make_ensembl_gene_list.R for information on how we created this table.
#' There are a small number (17) of genes with N/A chromosome - that's just how
#' they came to us from refgene (and were not in ensembl).
#' A '2' in the Notes column means it is a duplicate entry of another gene with the
#' same HGNC name - sometimes refgene puts the same gene on multiple chromosomes
#' or in multiple bands along the same chromosome.  There will be another entry with
#' the same HGNC name and a '1' in the Notes column, this is the primary entry.
#'
#' @docType data
#'
#' @usage data(ensembl_hg19_oct17)
#'
#' @format A data table with the transcription start/stop point for >48,000 genes.
#'
#' @keywords datasets
#'
#' @examples
#' data(ensembl_hg19_oct17)
#' FGFR2_row <- which(ensembl_hg19_oct17$HGNC_name == 'FGFR2')
#' print(ensembl_hg19_oct17[FGFR2_row, ])
"ensembl_hg19_oct17"
