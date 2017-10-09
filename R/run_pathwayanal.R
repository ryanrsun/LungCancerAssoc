# run_pathwayanal.R

#' High-level function to control pathway analysis.
#'
#' #' @param pathways_tab A data.frame of pathways defined by genes in the pathway.
#' First column name should be 'Pathway_name', second should be 'Pathway_description', all others should
#' be 'Gene1', 'Gene2', etc. Use NA to fill blanks.
#' @param pathways_tab_fname The name of a file formatted in the manner described by pathways_tab.
#' You only need to specify either pathways_tab or pathways_tab_fname.
#' @param gene_tab A data.frame which defines the location of each gene in the genome.
#' Should have column headings including at least 'Gene', 'CHR', 'Start', 'End'.
#' @param gene_tab_fname The name of a file formatted in the manner described by gene_tab.
#' You only need to specify either gene_tab or gene_tab_fname.
#' @param SS_file A data.frame holding all the summary statistics. Should have column headings
#' including at least 'RS' and 'P-value'.
#' @param SS_fname_root The root name of a file formatted in the manner described by SS_file.
#' If you use this option it is assume you have separated your summary statistics by chromosome
#' into files name [SS_fname_root][1].txt, [SS_fname_root][2].txt, etc.
#' You only need to specify either SS_file or SS_fname_root.
#' @param evecs_tab Data.frame of eigenvectors for correlation estimation.  Should have
#' column headings 'Subject', 'EV1', 'EV2', and so on.
#' @param evecs_tab_fname The name of a file formatted in the manner described by evecs_tab.
#' You only need to specify either evecs_tab or evecs_tab_fname.
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
#' @param Snum Used in cluster job submission scripts to organize jobs.
#' @param aID Used in cluster job submission scripts to organize jobs.
#' @param checkpoint Boolean, if true, print out diagnostic messages.
#'
#' @export
#'
#' @examples
#'

run_pathwayanal <- function(pathways_tab=NULL, pathways_tab_fname=NULL,
                            gene_tab=NULL, gene_tab_fname=NULL,
                            SS_file=NULL, SS_fname_root=NULL,
                            evecs_tab=NULL, evecs_tab_fname=NULL, num_PCs_use=5,
                            pathways_per_job=10, gene_buffer=5000, threshold_1000G=0.03,
                            prune_factor=0.5, prune_limit=0.0625, snp_limit=1000, hard_snp_limit=1500,
                            prune_to_start=TRUE, run_GHC=FALSE, out_name_root='pathway_anal',
                            refsnp_dir=NULL, input_dir=NULL, output_dir=NULL,
                            Snum=1, aID=1, checkpoint=TRUE) {

    # Error checking for inputs
    checked_files <- check_inputs(pathways_tab, pathways_tab_fname,
                                  gene_tab, gene_tab_fname,
                                  SS_file, SS_fname_root,
                                  evecs_tab, evecs_tab_fname,
                                  gene_sig_tab=NULL, gene_sig_tab_fname=NULL, num_PCs_use,
                                  pathways_per_job, gene_buffer, threshold_1000G,
                                  prune_factor, prune_limit, snp_limit, hard_snp_limit,
                                  prune_to_start, run_GHC, out_name_root,
                                  refsnp_dir, input_dir, output_dir,
                                  Snum, aID, checkpoint)
    P_mat <- checked_files$P_mat
    gene_tab <-  checked_files$gene_tab
    pathways_tab <- checked_files$pathways_tab

    # May not use GHC columns
    pathway_results <- data.frame(Pathway_name=pathways_tab$Pathway_name,
                                  Pathway_description=pathways_tab$Pathway_description,
                                  num_snps=NA, num_pathway_genes=NA, num_genes_used=NA,
                                  prune_R2=NA, num_large_Z=NA, gbj=NA, GBJ_p=NA, GBJ_err=NA,
                                  ghc=NA, GHC_p=NA, GHC_err=NA)
    pathway_results$num_pathway_genes <- apply(pathways_tab, 1, function(x) {length(which(!is.na(x) & x != ''))-2})

    # Loop through pathway table
    for (pathway_it in 1:nrow(pathways_tab))
    {
        # Start with pruning?
        ifelse(prune_to_start, prune_R2 <- prune_factor, prune_R2 <- 1)

        # Diagnostics
        if (checkpoint) {
            cat('\n Starting pathway ', pathway_it, 'of ', nrow(pathways_tab), '\n')
            cat('Pathway name is: ', as.character(pathways_tab$Pathway_name[pathway_it]), '\n')
            cat('Pathway description is: ', as.character(pathways_tab$Pathway_description[pathway_it]), '\n')
        }

        # Find the locations of all our genes
        # First have to remove NAs from each row, send in a cleaned gene vector
        # Subset
        subset_colnames <- paste('Gene', 1:(ncol(pathways_tab)-2), sep='')
        subset_row <- rep(FALSE, nrow(pathways_tab)); subset_row[pathway_it] <- TRUE
        path_genes <- subset(pathways_tab, subset=subset_row, select=subset_colnames)

        # Clean
        bad_elements <- which(is.na(path_genes) | path_genes=="")
        if (length(bad_elements) == length(path_genes)) {
            pathway_results$num_snps[pathway_it] <- 0
            pathway_results$num_genes_used[pathway_it] <- 0
            next
        }
        if (length(bad_elements) > 0) {path_genes <- path_genes[-bad_elements]}
        path_genes <- unlist(path_genes)

        # Get locations
        pathway_info <- define_pathway_loc(gene_tab=gene_tab, pathway_genes=path_genes,
                                           checkpoint=checkpoint)

        # Sometimes no genes
        if (class(pathway_info) == 'numeric') {
            pathway_results$num_snps[pathway_it] <- 0
            pathway_results$num_genes_used[pathway_it] <- 0
            next
        } else {
            pathway_results$num_genes_used[pathway_it] <- nrow(pathway_info)
        }

        # Tells us to keep looping if too many SNPs
        pathway_done <- FALSE

        # Return here if pathway has too many SNPs and we need to prune more
        while (pathway_done == FALSE) {

            # Build these matrices gene-by-gene
            combined_Gmat <- NULL
            combined_Gmat_record <- NULL

            # Loop through each gene
            for (gene_it in 1:nrow(pathway_info))
            {
                # Initial stats
                gene_name <- as.character(pathway_info$Gene[gene_it])
                CHR <- pathway_info$CHR[gene_it]
                start_bp <- pathway_info$Start[gene_it] - gene_buffer
                end_bp <- pathway_info$End[gene_it] + gene_buffer

                # Move a copy of the gene to our working directory
                fname_root <- paste('S', Snum, '_aID', aID, sep='')
                copy_from <- paste(refsnp_dir, '/', gene_name, '.ped', sep='')
                copy_to <- paste(input_dir, '/', fname_root, '.ped', sep='')
                cp_success <- system2(command='cp', args=c(copy_from, copy_to), wait=TRUE)
                # Not in one of reference panel or gene table
                if (cp_success) {next}
                copy_from <- paste(refsnp_dir, '/', gene_name, '.map', sep='')
                copy_to <- paste(input_dir, '/', fname_root, '.map', sep='')
                cp_success <- system2(command='cp', args=c(copy_from, copy_to), wait=TRUE)
                if (cp_success) {next}

                # Sometimes copying takes a while
                system2(command='sleep', args=c(2))

                # Clean data
                init_data_list <- clean_1000G_raw(fname_root=fname_root, gene_name=gene_name,
                                                  start_bp=start_bp, end_bp=end_bp,
                                                  checkpoint=checkpoint)

                # This means no data on that gene.
                if (class(init_data_list) == 'numeric') {next}
                # Have data, cleaned it.
                ped_file <- init_data_list$ped_file
                map_file <- init_data_list$map_file

                # Keep only the SNPs that are in summary stats file
                matched_data_list <- match_ped_summary(SS_fname_root=SS_fname_root, fname_root=fname_root, ped_file=ped_file,
                                                       map_file=map_file, CHR=CHR, start_bp=start_bp, end_bp=end_bp,
                                                       gene_name=gene_name, threshold_1000G=threshold_1000G, checkpoint=checkpoint)

                # No summary statistics in the region
                if (class(matched_data_list) == 'numeric') {next}

                # Found some SNP in both SS and ref data, recorded the genotype data and SS info.
                temp_Gmat <- matched_data_list$temp_Gmat
                temp_Gmat_record <- matched_data_list$temp_Gmat_record

                # LD pruning in PLINK
                pruned_list <- prune_snps(Snum=Snum, aID=aID, fname_root=fname_root, prune_R2=prune_R2,
                                          temp_Gmat=temp_Gmat, temp_Gmat_record=temp_Gmat_record,
                                          checkpoint=checkpoint)
                # Append every gene's results
                if (is.null(combined_Gmat))
                {
                    combined_Gmat <- pruned_list$temp_Gmat
                    combined_Gmat_record <- pruned_list$temp_Gmat_record
                } else {
                    combined_Gmat <- cbind(combined_Gmat, pruned_list$temp_Gmat)
                    combined_Gmat_record <- rbind(combined_Gmat_record, pruned_list$temp_Gmat_record)
                }

                # Checkpoint
                cat('Running total of SNPs: ', ncol(combined_Gmat), '\n')
                cat('Done with ', gene_name, ' ', gene_it, '/', nrow(pathway_info), '\n')

                # Remove the queried files?
                rm_dl_name <- paste(fname_root, '*', sep='')
                system2(command='rm', args=rm_dl_name)
                system2(command='sleep', args=c(1))
            }

            ####################################################
            # Went through all the genes in the pathway

            # No SNPs from the pathway genes are in the SS file?
            if (is.null(combined_Gmat))
            {
                pathway_results$num_snps[pathway_it] <- 0
                pathway_results$num_genes_used[pathway_it] <- nrow(pathway_info)
                # Go to next pathway
                pathway_done <- TRUE
                next
            }

            # Too many, redo with more pruning
            if (ncol(combined_Gmat) > snp_limit | ncol(combined_Gmat) < 2)
            {
                if (ncol(combined_Gmat) < 2) {
                    pathway_results$num_snps[pathway_it] <- ncol(combined_Gmat)
                    pathway_results$num_genes_used[pathway_it] <- nrow(pathway_info)

                    # Go to next pathway
                    pathway_done <- TRUE
                    cat ('Not enough snps, done with this pathway \n')
                    next
                }

                # Too many SNPs
                # Prune more and try again
                prune_R2 <- prune_R2 * prune_factor
                if (prune_R2 >= prune_limit)
                {
                    cat ('Too many SNPs, now prune at ', prune_R2, ' and redo \n')
                    next
                }

                # Under the prune limit - try the calculation or just give up?
                if (ncol(combined_Gmat) > hard_snp_limit) 	{
                    pathway_results$num_snps[pathway_it] <- ncol(combined_Gmat)
                    pathway_results$num_genes_used[pathway_it] <- nrow(pathway_info)
                    pathway_done <- TRUE
                    cat ('Too many SNPs and prune_R2 already too low, giving up \n')
                    next
                } else {
                    # If not under the hard snp limit, give it a try
                    cat ('prune_R2 too low, but going to give it a try with ', ncol(combined_Gmat), ' snps \n')
                }
            }

            # Run the tests
            all_pvalues <- run_tests(snp_mat=combined_Gmat, combined_Gmat_record=combined_Gmat_record,
                                     P_mat=P_mat, run_GHC=run_GHC, checkpoint=checkpoint)

            # Record
            pathway_results$gbj[pathway_it] <- all_pvalues$gbj
            pathway_results$GBJ_p[pathway_it] <- all_pvalues$GBJ_p
            pathway_results$GBJ_err[pathway_it] <- all_pvalues$GBJ_err
            pathway_results$num_snps[pathway_it] <- ncol(combined_Gmat)
            pathway_results$num_genes_used[pathway_it] <- nrow(pathway_info)
            pathway_results$num_large_Z[pathway_it] <- all_pvalues$num_large_Z

            if (run_GHC) {
                pathway_results$ghc[pathway_it] <- all_pvalues$ghc
                pathway_results$GHC_p[pathway_it] <- all_pvalues$GHC_p
                pathway_results$GHC_err[pathway_it] <- all_pvalues$GHC_err
            }

            pathway_results$prune_R2[pathway_it] <- prune_R2
            pathway_done <- TRUE
            cat('Done with this pathway \n')
        }
    }

    # Write
    setwd(output_dir)
    out_name <- paste(out_name_root, '_S', Snum, '_', aID, '.txt', sep='')
    write.table(pathway_results, file=out_name, append=F, quote=F, row.names=F, col.names=T, sep='\t')
}




