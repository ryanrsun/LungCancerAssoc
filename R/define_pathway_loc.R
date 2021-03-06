#' define_pathway_loc.R
#'
#' Pathways are most often defined by sets of genes, while genotype data is at the SNP-level
#' and is identified by chromosome (CHR) and base pair (BP). To link Pathways -- SNPs, we need
#' this function to tell us the genetic positions defined by each gene (and thus each pathway).
#' Incorporates some error checking and checkpointing messages.
#'
#' @param gene_tab A data frame that has columns Gene, CHR, BP for the annotation corresponding
#' to the genotype data.
#' @param pathway_genes A vector with the name of all the genes in the pathway of interest.
#' @param gene_sig_tab A table holding the p-value of each gene. Should only have the columns 'Gene' and 'P_value'.
#' If specified, then the returned pathway_info will be sorted in increasing order of p_value. If null, then no sorting.
#' @param checkpoint A boolean variable - TRUE if you want this function to print some diagnostic messages.
#'
#' @return A data frame containing columns for Gene, CHR, Start, End for all genes in your pathawy.
#'
#' @keywords internal
#' @export
#'
#' @examples
#' new_gene_tab <- data.frame(Gene=paste('Gene', 1:5, sep=''), CHR=1, Start=1:5, End=10:15)
#' define_gene_location(gene_tab=new_gene_tab, pathway_genes=c('Gene1', 'Gene2'), checkpoint=FALSE)

define_pathway_loc <- function(gene_tab, pathway_genes, gene_sig_tab=NULL, checkpoint) {

    # Do we even have any genes?
    na_cells <- which(is.na(pathway_genes) | pathway_genes == '')
    if (length(na_cells) > 0) {
        pathway_genes <- pathway_genes[-na_cells]
    }
    if (length(pathway_genes) == 0) {
        return (1)
    }

    # Build the output table
    pathway_info <- data.frame(Gene=pathway_genes, CHR=NA, Start=NA, End=NA, Diff=NA)
    for (gene_it in 1:length(pathway_genes))
    {
        gene_name <- as.character(pathway_genes[gene_it])

        # Find it
        gene_tab_row <- which(gene_tab$Gene == gene_name)

        # Can't find it? Leave NA
        if (length(gene_tab_row) == 0) {
            next
        }

        # Is it a duplicate?
        if (length(gene_tab_row) > 1) {
            primary_row <- which(gene_tab$Gene == gene_name & gene_tab$Notes == 1)
            if (length(primary_row) == 1) {
                gene_tab_row <- primary_row
            } else {
                gene_tab_row <- gene_tab_row[1]
                if (checkpoint) {cat('Duplicate entries found for: ', gene_name, '\n')}
            }
        }

        # Record
        pathway_info$CHR[gene_it] <- gene_tab$CHR[gene_tab_row]
        pathway_info$Start[gene_it] <-  gene_tab$Start[gene_tab_row]
        pathway_info$End[gene_it] <-  gene_tab$End[gene_tab_row]
        pathway_info$Diff[gene_it] <-  gene_tab$End[gene_tab_row] -  gene_tab$Start[gene_tab_row]
    }

    # Remove bad data or if on sex chromosome.
    missing_rows <- which(is.na(pathway_info$CHR) | is.na(pathway_info$Start) | is.na(pathway_info$End) |
                              pathway_info$CHR >= 23)
    # No found genes?
    if (length(missing_rows) == nrow(pathway_info)) {
        return(-1)
    }

    # Remove missing genes
    if (length(missing_rows) > 0) {
        missing_genes <- as.character(pathway_info$Gene[missing_rows])
        pathway_info <- pathway_info[-missing_rows, ]
    } else {
        missing_genes <- ''
    }

    # If gene significance table supplied, order the genes by increasing p_value
    if (!is.null(gene_sig_tab)) {
        pathway_info <- merge(x=pathway_info, y=gene_sig_tab, by='Gene', all.x=TRUE)
        NA_pvalues <- which(is.na(pathway_info$P_value))
        if (length(NA_pvalues) > 0) {
            pathway_info$P_value[NA_pvalues] <- 1
        }

        # Order
        pathway_info <- pathway_info[order(pathway_info$P_value, decreasing=FALSE), ]
    }

    # Diagnostic information
    if (checkpoint == TRUE) {
        cat('Missing ', length(missing_rows), ' genes: ', missing_genes, '\n')
        cat('Using ', nrow(pathway_info), ' out of ', length(pathway_genes), ' : ', as.character(pathway_info$Gene), '\n')
    }

    # Return
    return(pathway_info)
}
