#' query_1000G.R
#'
#' Get 1000G subject SNPs at specific locations on the genome. This is a little less disk-space
#' intensive than using 'wget' to download entire chromosome files. Only use for one region
#' at a time (region is one contiguous block on one chromosome).
#' NOTE: You must have the vcf_to_ped_covert.pl script (provided on 1000G website) in your working
#' directory!
#'
#' @param CHR The chromosome where your region is located.
#' @param start_bp The starting BP of the region.
#' @param end_bp The ending BP of the region.
#' @param gene_name Name of region you are querying, for error messages.
#' @param buffer Add additional distance to the start and end of the region, possibly to capture nearby
#' effects such as cis-eQTLs.
#' @param pop_vec A vector containing strings describing the populations you want to query. See 1000G
#' websites for the population codes.
#' @param Snum Used in cluster job submission scripts to organize jobs. Must be numeric.
#' @param aID Used in cluster job submission scripts to organize jobs. Must be numeric.
#' @param checkpoint Would you like the function to print diagnostic/error messages?
#'
#' @return Return 0 for success and 1 for error.
#'
#' @export
#'
#' @examples
#' query_1000G(CHR=1, start_bp=100000, end_bp=200000, buffer=10000, pop_vec='GBR')

query_1000G <- function(CHR, start_bp, end_bp, gene_name, buffer, pop_vec, Snum, aID, checkpoint=T)
{
    # Make sure inputs take the correct format
    if ( !is.numeric(c(CHR, start_bp, end_bp, buffer)) ) {
        stop('CHR, start_bp, end_bp, buffer must all be numeric')
    }

    if ( !is.character(pop_vec) ) {
        stop('pop_vec must be character strings.')
    }

    # Some room around the gene to allow for other influence structures
    start_bp <- start_bp - buffer
    end_bp <- end_bp + buffer
    region_string <- paste(CHR, ':', start_bp, '-', end_bp, sep='')

    # Data and panel link
    chr_vcf_link <- paste('ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr', CHR, 	'.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz', sep='')
    panel_link <- 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel'

    # Which populations
    pop_string <- paste('-population', pop_vec[1], sep=' ')
    if (length(pop_vec) > 1)
    {
        for (pop_it in 2:length(pop_vec))
        {
            pop_string <- paste(pop_string, '-population', pop_vec[pop_it], sep=' ')
        }
    }

    # Query
    res <- system2(command='perl', args=c('vcf_to_ped_convert.pl', '-vcf', chr_vcf_link, '-sample_panel_file',
                                          panel_link, '-region', region_string, pop_string, '-base_format', 'letter'), wait=TRUE)

    # Did query work?
    if (res != 0) {
        # This is not necessarily a bad thing, but we generally *should* be able to query all genes in a dense reference dataset.
        if (checkpoint) {
            cat('Error querying ', gene_name, '\n')
        }
        return (1)
    } else {
        # Immediately append the Snum and aID in case another job is working with the same gene
        init_ped_name <- paste(CHR, '_', start_bp, '-', end_bp, '.ped', sep='')
        init_info_name <- paste(CHR, '_', start_bp, '-', end_bp, '.info', sep='')
        new_ped_name <- paste('S', Snum, '_aID', aID, '_', CHR, '_', start_bp, '-', end_bp, '.ped', sep='')
        new_info_name <- paste('S', Snum, '_aID', aID, '_', CHR, '_', start_bp, '-', end_bp, '.info', sep='')
        system2(command='mv', args=c(init_ped_name, new_ped_name))
        system2(command='mv', args=c(init_info_name, new_info_name))

        return (0)
    }
}
