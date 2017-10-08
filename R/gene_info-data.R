#' Hg19 Gene table using mostly Ensembl annotations
#'
#' See make_ensembl_gene_list.R for information on how we created this table.
#'
#' @docType data
#'
#' @usage data(gene_info)
#'
#' @format A data table with the transcription start/stop point for >48,000 genes.
#'
#' @keywords datasets
#'
#' @examples
#' data(gene_info)
#' FGFR2_row <- which(gene_info$HGNC_name == 'FGFR2')
#' print(gene_info[FGFR2_row, ])
"gene_info"
