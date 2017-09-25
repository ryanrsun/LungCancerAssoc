#' prune_snps.R
#'
#' We can only do set-based tests with SNPs that are in both the ped file (to estimate correlation)
#' and in the summary statistics file (obviously because we need the summary statistic). This function
#' tells us which SNPs are indeed in both. Use only with one region (contiguous length on one chromosome) at a time.
#' Make sure you have the PLINK binary in your working directory!
#'
#' @param Snum The pruned file will be name S[Snum]_[aID].prune.out. The Snum and aID parameters
#' allow you to ensure that your pruned files do not overwrite each other if you have many jobs running at once.
#' @param aID The pruned file will be name S[Snum]_[aID].prune.out. The Snum and aID parameters
#' allow you to ensure that your pruned files do not overwrite each other if you have many jobs running at once.
#' @param fname_root The root of our .ped and .info files downloaded from 1000G, used to remove them at the end.
#' @param prune_R2 Prune all SNPs that have pairwise correlation greater than this threshold.
#' @param temp_Gmat The genotypes of the SNPs, will remove columns in accordance with PLINK pruning.
#' @param temp_Gmat_record Info on the SNPs, will remove rows in accordance with PLINK pruning.
#' @param checkpoint A boolean, if TRUE, print out diagnostic/error messages.
#'
#' @return A list with the elements temp_Gmat (containing the genotypes at each qualifying SNP) and
#' temp_Gmat_record (containing the info on SNPs), or 1 if nothing to return.
#'
#' @export
#'
#' @examples

prune_snps <- function(Snum, aID, fname_root, prune_R2, temp_Gmat, temp_Gmat_record, checkpoint)
{
    # Pruned outfile file name.
    prune_root <- paste('S', Snum, '_', aID, sep='')
    prune_out_name <- paste(prune_root, '.prune.out', sep='')

    # Get prune list in plink
    system2(command='./plink', args=c('--noweb', '--file', fname_root, '--indep-pairwise',
                                      20, 5, prune_R2, '--silent', '--out', prune_root), wait=T)
    system2(command='sleep', args=c(1))

    # Read the list of SNPs to keep, prune
    snps_to_prune <- tryCatch(read.table(prune_out_name), warning=function(w) w,
                              error=function(e) e)

    # If we have SNPs to prune
    if (class(snps_to_prune)[1] == 'data.frame')
    {
        # Prune temp_Gmat and temp_Gmat_record
        record_rows_to_prune <- which(temp_Gmat_record$RS %in% snps_to_prune[,1])
        temp_Gmat <- as.matrix(temp_Gmat[, -record_rows_to_prune])
        temp_Gmat_record <- temp_Gmat_record[-record_rows_to_prune, ]

        # Check temp_Gmat for SNP column consistency
        temp_Gmat_freqs <- apply(temp_Gmat, 2, mean)/2
        temp_record_freqs <- as.numeric(temp_Gmat_record$ped_A2_freq)
        discrep_freqs <- which(abs(temp_Gmat_freqs-temp_record_freqs)>10^(-10))
        if (length(discrep_freqs) > 0) {
            if (checkpoint) {
                cat( 'After pruning, possible discrepencies in ', length(discrep_freqs), 'freqs \n')
            } else {
                stop('After pruning, possible discrepencies in ', length(discrep_freqs), 'freqs')
            }
        }

        if (checkpoint) {
            cat('After pruning at ', prune_R2, ' we have ', nrow(temp_Gmat_record), ' SNPs left \n')
        }

    } else {
        if (checkpoint) {
            cat('Nothing to prune for', fname_root, '\n')
        }
    }

    #  Remove prune files
    rm_prune_name <- paste('S', Snum, '_', aID, '*', sep='')
    system2(command='rm', args=rm_prune_name)
    system2(command='sleep', args=c(1))

    return (list(temp_Gmat=temp_Gmat, temp_Gmat_record=temp_Gmat_record))
}
