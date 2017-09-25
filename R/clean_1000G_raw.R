#' clean_1000G_raw.R
#'
#' After querying data from 1000G, clean the raw data (mostly remove non-biallelic SNPs).
#'
#' @param fname_root The root (everything before .ped) of the downloaded 1000G file.
#' @param gene_name A name (often a gene name, but not always) for the region of SNPs you have downloaded,
#' only to print with error messages.
#' @param start_bp The starting position of the range we want.
#' @param end_bp The ending position of the range we want.
#' @param checkpoint A boolean, if TRUE, print out diagnostic/error messages.
#'
#' @return A list with the elements info_file and ped_file, containing the (clean) genotype info and counts.
#'
#' @export
#'
#' @examples
#' query_1000G(CHR=1, start_bp=100000, end_bp=200000, buffer=10000, pop_vec='GBR')

clean_1000G_raw <- function(fname_root, gene_name, start_bp, end_bp, checkpoint)
{
    # Open the .ped and .info files downloaded from 1000G
    ped_fname <- paste(fname_root, '.ped', sep='')
    ped_file <- tryCatch(read.table(ped_fname, header=F), warning=function(w) w, error=function(e) e)
    map_fname <- paste(fname_root, '.map', sep='')
    map_file <- tryCatch(read.table(map_fname, header=F), warning=function(w) w, error=function(e) e)

    # If can't open, quit
    if (class(map_file)[1] == 'simpleWarning' | class(ped_file)[1] == 'simpleWarning' |
        class(map_file)[1] == 'simpleError' | class(ped_file)[1] == 'simpleError') {
        system2(command='rm', args=map_fname)
        system2(command='rm', args=ped_fname)

        # Print error?
        if (checkpoint) {
            cat('No data for ', fname_root, ' so skipping.\n')
        }
        return (1)
    }

    # Sort the ped file IDs
    ped_file <- ped_file[order(ped_file[,1], decreasing=FALSE), ]

    # Truncate down if there are extra snps outside our range
    colnames(map_file) <- c('Chr', 'RS', 'Other', 'BP')
    keep_rows <- which(map_file$BP >= start_bp & map_file$BP <= end_bp)

    # If nothing left after truncating, move on
    if (length(keep_rows) == 0) {
        if (checkpoint) {
            cat('No data for ', fname_root, ' so skipping.\n')
        }
        return (1)
    }

    # Still SNPs, truncate
    map_file <- map_file[keep_rows, ]
    ped_keep_col1 <- c(1:6, keep_rows*2 + 5)
    ped_keep_col2 <- c(keep_rows*2 + 6)
    keep_col <- sort(unique(c(ped_keep_col1, ped_keep_col2)), decreasing=FALSE)
    ped_file <- ped_file[, keep_col]

    # Sometimes if a column is all 'T', R reads it as logical 'TRUE'
    ped_col_classes <- lapply(ped_file, class)
    bad_ped_col <- which(ped_col_classes == 'logical')
    if (length(bad_ped_col) > 0)
    {
        ped_file[ ,bad_ped_col] <- 'T'
    }

    # Check for SNPs that are not biallelic, remove them
    num_check <- (ncol(ped_file)-6)/2
    not_biallelic_info <- NULL
    not_biallelic_ped <- NULL
    for (biallelic_it in 1:num_check)
    {
        col_1 <- biallelic_it*2+5
        col_2 <- biallelic_it*2+6

        # Should only be two levels (bases) for each SNP
        all_bases <- c(levels(ped_file[, col_1]), levels(ped_file[, col_2]))
        all_bases <- unique(all_bases)

        # If it's not biallelic, add these columns to the 'bad' list.
        if (length(all_bases) != 2) {
            if (is.null(not_biallelic_info))
            {
                not_biallelic_info <- biallelic_it
                not_biallelic_ped <- c(col_1, col_2)
            } else{
                not_biallelic_info  <- c(not_biallelic_info, biallelic_it)
                not_biallelic_ped <- c(not_biallelic_ped, col_1, col_2)
            }
        }
    }

    # Remove if any bad SNPs
    if (length(not_biallelic_info) > 0)
    {
        map_file <- map_file[-not_biallelic_info, ]
        ped_file <- ped_file[, -not_biallelic_ped]
        if (checkpoint) {
            cat ('Removed ', length(not_biallelic_info), ' not biallelic SNPs from ',
                 gene_name, '\n')
        }
    } else {
        if (checkpoint) {
            cat ('Removed no SNPs from ', gene_name, '\n')
        }
    }

    return (list(map_file=map_file, ped_file=ped_file))
}
