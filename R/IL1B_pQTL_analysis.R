# Xihong wants to run a separate analysis on 26 individual SNPs
# that are pQTLs of IL1B.

IL1B_snps_analysis <- function() {
    # The pQTLs
    target_snps <- c('rs7588285', 'rs13402561', 'rs4889294', 'rs115242021', 'rs79186011', 'rs185085918',
                     'rs185637423', 'rs181572669', 'rs2255932', 'rs1545747', 'rs74694425', 'rs4786740', 'rs62015704',
                     'rs201997403', 'rs184938669', 'rs9898641', 'rs1942793', 'rs143319329', 'rs2149803', 'rs6559536',
                     'rs141211125', 'rs183059183', 'rs139079409', 'rs188550518', 'rs188572687', 'rs6105735')
    single_snp_results <- data.frame(SNP=target_snps, Chr=NA, BP=NA,
                                     Lung_Overall=NA, Small=NA, Squam=NA, Adeno=NA, CAD=NA)

    # We first need to figure out where these SNPs are
    studies_list <- c('overall', 'cad')
    dir_list <- paste('/users/ryansun/desktop/LC/', studies_list, sep='')

    for (snp_it in 8:length(target_snps)) {
        found <- FALSE
        study_it <- 1
        cat('Searching for SNP ', snp_it, ' out of ', length(target_snps), '\n')

        while (!found) {
            # Set new study working directory, ideally just once per snp
            temp_success <- tryCatch(setwd(dir_list[study_it]), warning=function(w) w,
                                     error=function(e) e)
            if (class(temp_success)[1] %in% c('simpleError', 'simpleWarning')) {break}
            cat('Searching ', dir_list[study_it], '\n')

            # Search through all of chromosomes
            for (chr_it in 1:22) {
                temp_fname <- paste('LC_', studies_list[study_it], '_hg19_chr', chr_it,
                                    '.txt', sep='')
                # CAD different
                if (study_it == 2) {
                    temp_fname <- paste('CAD_hg19_chr', chr_it,
                                        '.txt', sep='')
                }

                # Open single chromosome SS file
                temp_SS_file <- fread(temp_fname, header=T)
                temp_row <- which(temp_SS_file$RS == target_snps[snp_it])

                # Found it, record
                if (length(temp_row) == 1) {
                    single_snp_results$Chr[snp_it] <- chr_it
                    single_snp_results$BP[snp_it] <- temp_SS_file$BP[temp_row]
                    cat('Found ', target_snps[snp_it], ' at ', chr_it, ':',
                        temp_SS_file$BP[temp_row], '\n')
                    found <- TRUE
                    break
                }

                # Not found, next chromosome
                cat('Not on chromosome ', chr_it, '\n')
            }

            # Try next study
            study_it <- study_it + 1
        }
    }


    ##################################################
    ##################################################
    # Now pull out the summary stats from each study
    studies_list <- c('overall', 'small', 'squam', 'adeno', 'cad')
    dir_list <- paste('/users/ryansun/desktop/LC/', studies_list, sep='')

    for (snp_it in 1:length(target_snps)) {
        study_it <- 1
        temp_chr <- single_snp_results$Chr[snp_it]
        if (is.na(temp_chr)) {next}
        cat('Searching for SNP ', snp_it, ' out of ', length(target_snps), '\n')

        for (study_it in 1:length(studies_list)) {
            setwd(dir_list[study_it])
            temp_fname <- paste('LC_', studies_list[study_it], '_hg19_chr', temp_chr,
                                    '.txt', sep='')
            # CAD different
            if (study_it == 5) {
                temp_fname <- paste('CAD_hg19_chr', temp_chr,
                                    '.txt', sep='')
            }

            # Open single chromosome SS file
            temp_SS_file <- fread(temp_fname, header=T)
            temp_row <- which(temp_SS_file$RS == target_snps[snp_it])

            # Not found?
            if (length(temp_row) == 0) {next}

            # Record
            if (study_it != 5) {
                single_snp_results[snp_it, 3+study_it] <- temp_SS_file$'P-value'[temp_row]
            } else {
                single_snp_results[snp_it, 3+study_it] <- temp_SS_file$p_dgc[temp_row]
            }

        }
    }

    setwd('/users/ryansun/documents/research/newsoftware/lung_cancer/lungCancerAssoc/data')
    write.table(single_snp_results, 'IL1B_pQTL_results.txt', append=F, quote=F, row.names=F, col.names=T)
}

