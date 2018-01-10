# check_inputs.R

#' Check inputs to run_pathwayanal() for errors.
#'
#' @param pathways_tab A data.frame of pathways defined by genes in the pathway.
#' First column name should be 'Pathway_name', second should be 'Pathway_description', all others should
#' be 'Gene1', 'Gene2', etc. Use NA to fill blanks.
#' @param pathways_tab_fname The name of a file formatted in the manner described by pathways_tab.
#' You only need to specify either pathways_tab or pathways_tab_fname.
#' @param gene_tab A data.frame which defines the location of each gene in the genome.
#' Should have column headings including at least 'Gene', 'CHR', 'Start', 'End'.
#' @param gene_tab_fname The name of a file formatted in the manner described by gene_tab.
#' You only need to specify either gene_tab or gene_tab_fname.
#' @param SS_file A data.frame holding all the summary statistics. Should have column headings
#' including at least 'CHR', 'BP', 'RS' and 'P-value'.
#' @param SS_fname_root The root name of a file formatted in the manner described by SS_file.
#' If you use this option it is assume you have separated your summary statistics by chromosome
#' into files name [SS_fname_root][1].txt, [SS_fname_root][2].txt, etc.
#' You only need to specify either SS_file or SS_fname_root.
#' @param evecs_tab Data.frame of eigenvectors for correlation estimation.  Should have
#' column headings 'Subject', 'EV1', 'EV2', and so on.
#' @param evecs_tab_fname The name of a file formatted in the manner described by evecs_tab.
#' You only need to specify either evecs_tab or evecs_tab_fname.
#' @param gene_sig_tab Data.frame with two columns - 'Gene' and 'P_value'. Holds the significance
#' of each single gene.
#' @param gene_sig_tab_fname The name of a file formatted in the manner described by gene_sig_tab.
#' You only need to specify at most one of gene_sig_tab or gene_sig_tab_fname.  If you don't specify
#' one, then the leave_k_out procedure will not return sensible results, but you may not need that function.
#' @param num_PCs_use Number of PCs to use.
#' @param pathways_per_job How many pathways to test in one call of the run_pathwayanal() function.
#' Will only be used if you also specify aID to determine which part of the pathway_tab to use.
#' @param gene_buffer A buffer region added to the Start and End of each gene region to capture,
#' for example, possible cis-eQTL effects.
#' @param threshold_1000G The minimum MAF needed for a reference panel SNP before we trust it
#' to be used in covariance estimation.
#' @param snp_limit If the pathway has more than this many SNPs, rerun the function to prune
#' more aggressively before testing. Recommended value of 1000, do not set above 2000 or numerical
#' stability will suffer greatly.
#' @param hard_snp_limit If after pruning limit we still haven't gone under snp_limit, can slightly
#' raise the threshold and see if calculation is stable enough.
#' @param prune_factor If the pathway has more than snp_limit SNPs, then multiply the current
#' pruning level by this factor and rerun.
#' @param prune_limit If the pruning factor is less than this amount, stop pruning and move on.
#' @param prune_to_start If true, begin pruning at prune_factor, otherwise don't prune on first run.
#' @param run_GHC Boolean, if true then test with both GBJ and GHC, if false just GBJ.
#' @param out_name_root Root of output filename. If Snum and aID are specified then output name will
#' be [out_name_root]_S[Snum]_[aID].txt.
#' @param refsnp_dir Directory holding reference panel genotypes.
#' @param input_dir Directory holding summary statistics, pathway table, gene table, eigenvectors, PLINK binary.
#' @param output_dir Directory to save output file.
#' @param Snum Used in cluster job submission scripts to organize jobs. Defaults to 1, must be numeric.
#' @param aID Used in cluster job submission scripts to organize jobs. Defaults to 1, must be numeric.
#' @param checkpoint Boolean, if true, print out diagnostic messages.
#'
#' @return A list containing pathways_tab, gene_tab, and the projection matrix for covariance estimation.
#'
#' @keywords internal
#' @export
#'
#' @examples
#'
#'
check_inputs <- function(pathways_tab=NULL, pathways_tab_fname=NULL,
                         gene_tab=NULL, gene_tab_fname=NULL,
                         SS_file=NULL, SS_fname_root=NULL,
                         evecs_tab=NULL, evecs_tab_fname=NULL,
                         gene_sig_tab=NULL, gene_sig_tab_fname=NULL, num_PCs_use=5,
                         pathways_per_job=10, gene_buffer=5000, threshold_1000G=0.03,
                         snp_limit=1000, hard_snp_limit=1500, prune_factor=0.5, prune_limit=0.0625,
                         prune_to_start=TRUE, run_GHC=FALSE, out_name_root='pathway_anal',
                         refsnp_dir=NULL, input_dir=NULL, output_dir=NULL,
                         Snum=1, aID=1, checkpoint=TRUE) {

    # Make sure our numeric arguments are actually numeric
    if (class(c(pathways_per_job, gene_buffer, threshold_1000G, prune_factor,
                prune_limit, snp_limit, aID)) != 'numeric') {
        stop('You have specified a non-numeric argument for a numeric parameter.')
    }
    if (length(which(is.na(c(pathways_per_job, gene_buffer, threshold_1000G, prune_factor,
                prune_limit, snp_limit, aID)))) > 0 ) {
        stop('You have specified a non-numeric argument for a numeric parameter.')
    }


    # Make sure logical arguments are logical
    if (class(c(checkpoint, prune_to_start, run_GHC)) != 'logical') {
        stop('You have specified a non-Boolean argument for a Boolean parameter.')
    }
    if (length(which(is.na(c(checkpoint, prune_to_start, run_GHC)))) > 0) {
        stop('You have specified a non-Boolean argument for a Boolean parameter.')
    }

    # Try to set input and output and refsnp directory
    pwd <- getwd()
    if (!is.null(output_dir)) {
        temp_success <- tryCatch(setwd(output_dir), warning=function(w) w, error=function(e) e)
        if (class(temp_success)[1] %in% c('simpleError', 'simpleWarning')) {
            stop('Could not set output directory.')
        }
    }
    if (!is.null(refsnp_dir)) {
        temp_success <- tryCatch(setwd(refsnp_dir), warning=function(w) w, error=function(e) e)
        if (class(temp_success)[1] %in% c('simpleError', 'simpleWarning')) {
            stop('Could not set refsnp directory.')
        }
    } else {
        stop('Must set refsnp directory.')
    }
    # Try input last so we can revert back to it.
    if (!is.null(input_dir)) {
        temp_success <- tryCatch(setwd(input_dir), warning=function(w) w, error=function(e) e)
        if (class(temp_success)[1] %in% c('simpleError', 'simpleWarning')) {
            stop('Could not set working directory.')
        }
    } else {
        stop('Must set input directory.')
    }

    # Make the projection matrix for correlation estimation
    # Use the fname only if an actual table not provided.
    # If neither fname nor table provided, use the one in LungCancerAssoc/data
    if (is.null(evecs_tab) & is.null(evecs_tab_fname)) {
        data(PCs_1000G_Euro)
        evecs_tab <- PCs_1000G_Euro
    } else if (is.null(evecs_tab)) {
        evecs_tab <- tryCatch(read.table(evecs_tab_fname, header=T),
                              warning=function(w) w, error=function(e) e)
        if (class(evecs_tab)[1] %in% c('simpleWarning', 'simpleError')) {
            stop('Problem opening eigenvectors table')
        }
    }
    # Sort the eigenvectors by ID
    evecs_tab <- evecs_tab[order(evecs_tab$ID, decreasing=FALSE), ]

    # Build the projection matrix
    evec_cols_to_use <- paste('PC', 1:num_PCs_use, sep='')
    X_mat <- as.matrix( subset(evecs_tab, select=evec_cols_to_use) )
    X_mat <- cbind(1, X_mat)
    W_mat <- diag(x=1, nrow=nrow(evecs_tab), ncol=nrow(evecs_tab))
    P_mat <- tryCatch(W_mat - X_mat %*% solve(t(X_mat) %*% X_mat) %*% t(X_mat),
                      warning=function(w) w, error=function(e) e)
    if (class(P_mat)[1] %in% c('simpleWarning', 'simpleError')) {
        stop('Problem creating projection matrix from eigenvectors.')
    }

    # Make sure we can open the pathways table.
    # Use the fname only if an actual table not provided.
    # If neither fname nor table provided, load from data/ folder.
    if ( is.null(pathways_tab) & is.null(pathways_tab_fname) ) {
        data(bader_apr1_pathways)
        pathways_tab <- bader_apr1_pathways
    } else {
        if (is.null(pathways_tab)) {
            pathways_tab <- tryCatch(read.table(pathways_tab_fname, header=T),
                                 warning=function(w) w, error=function(e) e)
            if (class(pathways_tab)[1] != 'data.frame') {
                stop('Problem opening pathways list.')
            }
        }
    }

    # Make sure it has the correct column names
    correct_colnames <- c('Pathway_name', 'Pathway_description', paste('Gene', 1:(ncol(pathways_tab)-2), sep=''))
    if ( length(which(correct_colnames == colnames(pathways_tab))) != ncol(pathways_tab)) {
        stop('Your pathway table does not appear to have the correct column headings.')
    }

    # Cut the pathway list
    if ( !is.null(aID) ) {
        start_row <- (aID-1)*pathways_per_job + 1
        end_row <- aID*pathways_per_job
        if (start_row > nrow(pathways_tab)) {stop()}
        if (end_row > nrow(pathways_tab)) {end_row <- nrow(pathways_tab)}
        pathways_tab <- pathways_tab[start_row:end_row, ]
    }

    # Check a gene table provided.
    # Use the fname only if table not provided.
    # If neither fname not table provided, use the one in LungCancerAssoc/data
    if (is.null(gene_tab) & is.null(gene_tab_fname)) {
        data(ensembl_refgene_hg19_20180109)
        gene_tab <- ensembl_refgene_hg19_20180109
        colnames(gene_tab) <- c('Transcript_ID', 'CHR', 'Strand', 'Start', 'End',
                                'cdsStart', 'cdsEnd', 'num_exons', 'Gene', 'Gene_ID', 'Notes')

        # Just get rid of the duplicates now, why keep them?
        gene_tab <- gene_tab[which(gene_tab$Notes <= 1), ]
        if (nrow(gene_tab) != 49132) {stop('Did you change the default gene table recently?')}
    } else {
        if (is.null(gene_tab)) {
            gene_tab <- tryCatch(data.table::fread(gene_tab_fname, data.table=F, header=T),
                            warning=function(w) w, error=function(e) e)

            if (class(gene_tab)[1] != 'data.frame') {
                stop('Problem opening gene table.')
            }
        }
    }

    # Check gene table for correct column names
    if ( length(which(c('Gene', 'CHR', 'Start', 'End') %in% colnames(gene_tab))) != 4) {
        stop('Your gene table must have columns "Gene", "CHR", "Start", "End".')
    }

    # Make sure the gene table CHR/Start/End are numeric!
    gene_tab$CHR <- as.numeric(as.character(gene_tab$CHR))
    gene_tab$Start <- as.numeric(as.character(gene_tab$Start))
    gene_tab$End <- as.numeric(as.character(gene_tab$End))

    # Print out where FGFR2 is.
    # Should be chr 10 123237844..123357972
    cat('Checking gene table for FGFR2... \n')
    FGFR2_row <- which(gene_tab$Gene == 'FGFR2')
    if (length(FGFR2_row) == 0) {
        cat('Warning: No FGFR2 found.')
    } else if (length(FGFR2_row >= 1)) {
        cat('FGFR2 found at CHR: ', gene_tab$CHR[FGFR2_row], '\n',
            'Start: ', gene_tab$Start[FGFR2_row], '\n',
            'End: ', gene_tab$End[FGFR2_row], '\n')
    }

    # Check a summary statistics file provided
    # Use SS_root only if SS_file not provided.
    if (is.null(SS_fname_root) & is.null(SS_file)) {
        stop('You need to provide a summary statistics file')
    } else if (is.null(SS_file)) {
        cat('Attempting to open Chr 1 summary statistics... \n')
        SS_file <- tryCatch(data.table::fread(paste(SS_fname_root, '1.txt', sep=''), header=T),
                                 warning=function(w) w, error=function(e) e)
        if (class(SS_file)[1] %in% c('simpleWarning', 'simpleError')) {
            stop('Could not read Chr 1 summary statistics.')
        }
    }

    # Check the SS file for must have columns
    must_have_columns <- c('Chr', 'BP', 'P-value', 'RS')
    if (length(which(must_have_columns %in% colnames(SS_file))) != length(must_have_columns)) {
        stop('Your SS_file does not have all the required columns: ', must_have_columns)
    }
    if (class(SS_file)[1] != 'data.table' & class(SS_file)[1] != 'data.frame') {
        stop('SS_file must be a data.frame')
    }

    # Check gene_sig_tab_fname
    if (!is.null(gene_sig_tab_fname)) {
        if (class(gene_sig_tab_fname) != 'character') {
            stop('You have specified an invalid gene_sig_tab_fname name.')
        } else {
            gene_sig_tab <- tryCatch(data.table::fread(gene_sig_tab_fname, header=T),
                                warning=function(w) w, error=function(e) e)
            if (class(gene_sig_tab)[1] %in% c('simpleWarning', 'simpleError')) {
                stop('Could not read gene_sig_tab.')
            }
        }
    }

    # Now check that gene_sig_tab has the correct columns
    if (!is.null(gene_sig_tab)) {
       if (ncol(gene_sig_tab) != 2 | colnames(gene_sig_tab) != c('Gene', 'P_value')) {
           stop('gene_sig_tab should have exactly two columns, Gene and P_value')
       }
    }

    # Check output name
    if (class(out_name_root) != 'character') {
        stop('You have specified an invalid output name.')
    }

    return(list(gene_tab=gene_tab, pathways_tab=pathways_tab, P_mat=P_mat, gene_sig_tab=gene_sig_tab))
}
