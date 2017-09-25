#' run_tests.R
#'
#' Run GHC, GBJ, and minP on the summmary statistics.
#'
#' @param snp_mat This is the n*p matrix of genotypes from your reference panel. Hopefully the
#' temp_Gmat variable from prune_snps().
#' @param combined_Gmat_record The snp information data frame, hopefully the one output from prune_snps().
#' @param run_GHC Boolean, if true run both GBJ and GHC, otherwise just GBJ.
#' @param P_mat Projection matrix for estimation of covariance matrix.
#' @param checkpoint Boolean, true if you would like to print diagnostic information.
#'
#' @return A list with the elements GBJ, GBJ_p, GBJ_err, GHC, GHC_p, GHC_err, and num_large_Z (number of test
#' statistics greater than 8.2 in magnitude, which R cannot deal with well).
#'
#' @export
#'
#' @examples

run_tests <- function(snp_mat, combined_Gmat_record, P_mat=P_mat, run_GHC=FALSE, checkpoint=TRUE)
{
    # Get the correlations
    cor_mat <- matrix(data=NA, nrow=ncol(snp_mat), ncol=ncol(snp_mat))
    denominators <- rep(NA, ncol(snp_mat))
    for (i in 1:ncol(snp_mat))
    {
        temp_G <- snp_mat[, i]
        denominators[i] <- sqrt(t(temp_G) %*% P_mat %*% temp_G)
    }
    for (i in 2:ncol(snp_mat))
    {
        for (j in 1:(i-1))
        {
            temp_G1 <- snp_mat[, i]
            temp_G2 <- snp_mat[, j]
            numerator <- t(temp_G1) %*% P_mat %*% temp_G2
            sig_hat <- numerator / (denominators[i] * denominators[j])
            cor_mat[i, j] <- sig_hat
            cor_mat[j, i] <- sig_hat
        }
    }

    if (checkpoint) {
            cat(ncol(snp_mat), ' SNPs, now running GBJ \n')
    }

    pvalues <- as.numeric(combined_Gmat_record$Pvalue)
    abs_Z <- abs(qnorm(1-pvalues/2))
    d <- length(abs_Z)
    num_large_Z <- which(abs_Z > 8.19)
    if( length(num_large_Z) > 0) {
        num_large_Z <- length(num_large_Z)
        abs_Z[which(abs_Z > 8.19)] <- 8.19
    } else {
        num_large_Z <- 0
    }

    # Run GBJ
    GBJ_output <- GBJ(test_stats=abs_Z, cor_mat=cor_mat)

    # Run GHC
    if (run_GHC) {
        if (checkpoint) {
            cat('Done GBJ, running GHC. \n')
        }
        GHC_output <- GHC(test_stats=abs_Z, cor_mat=cor_mat)
        return( list(gbj=GBJ_output$GBJ, GBJ_p=GBJ_output$GBJ_pvalue, GBJ_err=GBJ_output$err_code,
                     ghc=GHC_output$GHC, GHC_p=GHC_output$GHC_pvalue, GHC_err=GHC_output$err_code,
                     num_large_Z=num_large_Z) )
    } else {
        return( list(gbj=GBJ_output$GBJ, GBJ_p=GBJ_output$GBJ_pvalue, GBJ_err=GBJ_output$err_code,
                     ghc=NA, GHC_p=NA, GHC_err=NA,
                     num_large_Z=num_large_Z) )
    }
}

