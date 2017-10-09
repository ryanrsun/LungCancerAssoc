#' match_ped_summary.R
#'
#' We can only do set-based tests with SNPs that are in both the ped file (to estimate correlation)
#' and in the summary statistics file (obviously because we need the summary statistic). This function
#' tells us which SNPs are indeed in both. Use only with one region (contiguous length on one chromosome) at a time.
#'
#' @param SS_fname_root Root of the summary statistic filename. The full filename should be [SS_root][CHR].txt.
#' This file should have column headers with the names 'CHR' and 'P-value' and 'BP' and 'RS'.
#' @param fname_root Root of the downloaded 1000G files. Used to delete the files. Leave as NULL if you
#' don't want to delete the downloaded files.
#' @param ped_file A standard PLINK ped file, hopefully cleaned from clean_1000G_raw.
#' @param map_file The standard .map file downloaded from 1000G, hopefully cleaned from clean_1000G_raw.
#' @param CHR The chromosome of the region.
#' @param start_bp The starting BP of the region.
#' @param end_bp The ending BP of the region.
#' @param gene_name A name given to the region (often a gene); used for printing error messages.
#' @param threshold_1000G Only use 1000G SNPs which pass this MAF? Rare alleles may be too unstable
#' for estimating correlations.
#' @param checkpoint A boolean, if TRUE, print out diagnostic/error messages.
#'
#' @return A list with the elements temp_Gmat (containing the genotypes at each qualifying SNP) and
#' temp_Gmat_record (containing the info on SNPs), or 1 if nothing to return.
#'
#' @export
#'
#' @examples
#'

match_ped_summary <- function(SS_fname_root, fname_root, ped_file, map_file, CHR, start_bp, end_bp, gene_name, threshold_1000G, checkpoint)
{

    # Open summary stats file
    SS_fname <- paste(SS_fname_root, CHR, '.txt', sep='')

    # Got rid of the tryCatch here because we were getting some loss of accuracy messages,
    # so make sure these files all exist.
    SS_file <- suppressWarnings(fread(SS_fname, header=T, showProgress=FALSE))

    # Do we have any summary statistics in our region?
    SS_rows <- which(SS_file$Chr==CHR & SS_file$BP>=start_bp & SS_file$BP<=end_bp)
    if (length(SS_rows) == 0) {
        # If not, remove the .ped and .info files.
        rm_dl_name <- paste(fname_root, '*', sep='')
        system2(command='rm', args=rm_dl_name)
        system2(command='sleep', args=c(1))

        if (checkpoint) {
            cat('No summary statistics in the region of ', gene_name, 'moving on\n')
        }
        return (1)
    } else {
        # Cut the entire summary stat file down to just those rows in our region.
        temp_SS <- SS_file[SS_rows,]
    }

    # Filter the SNPs we will ultimately use in testing.
    # Filter by MAF in summary stat sample and also 1000G GBR sample.
    # Keep track of the ones we use so that we can perform LD pruning later.
    temp_Gmat <- matrix(data=NA, nrow=nrow(ped_file), ncol=nrow(map_file))
    temp_Gmat_record <- matrix(data=NA, nrow=nrow(map_file), ncol=8); temp_Gmat_record <- data.frame(temp_Gmat_record)
    colnames(temp_Gmat_record) <- c('RS', 'BP', 'ped_A1', 'ped_A2', 'SS_minorA', 'SS_majorA', 'Pvalue', 'ped_A2_freq')
    have_SS <- rep(NA, nrow(map_file))
    for (snp_it in 1:nrow(map_file))
    {
        # Can we find the SNP in the summary stats file?
        temp_rs <- as.character(map_file$RS[snp_it])
        summary_row <- which(temp_SS$RS == temp_rs)

        if ( length(summary_row) == 1)								# Exists in summary stat file
        {
            # We do not take into account the reference allele here - but it will be necessary
            # if we begin to combine p-values for a single SNP across multiple studies.

            ped_column_1 <- 5 + snp_it*2
            ped_column_2 <- 6 + snp_it*2

            # Everything should have been cleaned and biallelic,
            a1 <- as.character(ped_file[1, ped_column_1])
            freq1 <- ( length(which(ped_file[,ped_column_1] == a1)) +
                           length(which(ped_file[,ped_column_2] == a1)) ) / (2*nrow(ped_file))
            freq2 <- ( length(which(ped_file[,ped_column_1] != a1)) +
                           length(which(ped_file[,ped_column_2] != a1)) ) / (2*nrow(ped_file))

            # Find a2
            if (length(which(ped_file[, ped_column_1] != a1)) > 0) {
                a2_ind <- which(ped_file[, ped_column_1] != a1)[1]
                a2 <- as.character(ped_file[a2_ind, ped_column_1])
            } else {
                a2_ind <- which(ped_file[, ped_column_2] != a1)[1]
                a2 <- as.character(ped_file[a2_ind, ped_column_2])
            }

            # Make sure the minor allele passes the MAF threshold.
            if (freq1 > threshold_1000G & freq2 > threshold_1000G)
            {
                # Capture the genotype information
                temp_snp1 <- rep(0, nrow(ped_file))
                temp_snp2 <- rep(0, nrow(ped_file))
                side1 <- which(ped_file[,ped_column_1] != a1)
                side2 <- which(ped_file[,ped_column_2] != a1)
                temp_snp1[side1] <- 1
                temp_snp2[side2] <- 1
                temp_Gmat[, snp_it] <- temp_snp1 + temp_snp2

                # Record the SNP information
                temp_Gmat_record$RS[snp_it] <- temp_rs
                temp_Gmat_record$BP[snp_it] <- temp_SS$BP[summary_row]
                temp_Gmat_record$ped_A1[snp_it] <- a1
                temp_Gmat_record$ped_A2[snp_it] <- a2
                temp_Gmat_record$ped_A2_freq[snp_it] <- freq2
                temp_Gmat_record$Pvalue[snp_it] <- temp_SS$'P-value'[summary_row]
                #temp_Gmat_record$SS_majorA[snp_it]  <-
                #temp_Gmat_record$SS_minorA[snp_it]  <-
                have_SS[snp_it] <- 1
            }
        }
    }

    # It's possible we have none of these SNPs in the summary stats after filtering.
    if (sum(have_SS, na.rm=TRUE) == 0)
    {
        if (checkpoint) {
            cat('No summary statistics after filtering in the region of ', gene_name, 'moving on\n')
        }
        rm_dl_name <- paste(fname_root, '*', sep='')
        system2(command='rm', args=rm_dl_name)
        system2(command='sleep', args=c(1))
        return (1)
    } else {
        # Print how many qualifying SNPs
        if (checkpoint) {
            cat('We have ', sum(have_SS, na.rm=TRUE), ' out of ', nrow(map_file), ' for ',
                gene_name, '\n')
        }
    }

    # Trim temporary G matrix and record matrix
    temp_Gmat <- as.matrix(temp_Gmat[, which(have_SS==1)])
    temp_Gmat_record <- temp_Gmat_record[which(have_SS==1), ]

    # Make new ped_file and map_has_SS using only the SNPs that were found in summary stats file
    map_has_SS <- which(have_SS==1)
    ped_has_SS <- sort( c(1:6, (map_has_SS*2+5), (map_has_SS*2+6)), decreasing=FALSE )
    map_file <- map_file[map_has_SS, ]
    ped_file <- ped_file[, ped_has_SS]

    # Check our ped_file still has column integrity (if more than one left)
    if (ncol(temp_Gmat) > 1)
    {
        old_ped_freqs <- cbind(as.numeric(temp_Gmat_record$ped_A2_freq), 1-as.numeric(temp_Gmat_record$ped_A2_freq))
        old_ped_freqs <- apply(old_ped_freqs, 1, min)
        for (snp_it in 1:nrow(map_file))
        {
            temp_a <- as.character(ped_file[1, snp_it*2+5])
            total_a <- length(which(ped_file[, snp_it*2+5]==temp_a)) +
                length(which(ped_file[, snp_it*2+6]==temp_a))
            temp_freq <- total_a/(2*nrow(ped_file))
            temp_freq <- min(temp_freq, 1-temp_freq)
            temp_dif <- temp_freq - old_ped_freqs[snp_it]
            if (abs(temp_dif)>10^(-10)) {
                if (checkpoint) {
                    cat('Before pruning, something wrong with snp ', snp_it, '\n')
                } else {
                    stop('Matching of ped and SS file results in MAF discrepancy when recording genotypes.')
                }
            }
        }
    }

    # Write the new ped file (for pruning)
    ped_fname <- paste(fname_root, '.ped', sep='')
    write.table(ped_file, ped_fname, append=F, quote=F, row.names=F, col.names=F)

    # Write the map file (for pruning)
    map_fname <- paste(fname_root, '.map', sep='')
    write.table(map_file, file=map_fname, append=F, quote=F, row.names=F, col.names=F)

    return (list(temp_Gmat=temp_Gmat, temp_Gmat_record=temp_Gmat_record))
}

