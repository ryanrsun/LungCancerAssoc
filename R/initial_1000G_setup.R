# We need to get the individual IDs of the people that we
# are using for covariance estimation / PCs.

initial_1000G_setup <- function() {
    # Italy, Spain, Great Britain, Utah European, Finland populations
    pops_to_use <- c('TSI', 'IBS', 'GBR', 'CEU', 'FIN')

    # Download a phenotype / sample info file
    system2(command='wget',
            args='ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped')

    # Open the file, get the individual IDs
    sample_info_tab <- fread('20130606_g1k.ped', header=T)

    # No fathers or children, just mothers and unrelated
    IDs_unrelated <- which(sample_info_tab$Population %in% pops_to_use &
                               sample_info_tab$Relationship == 'unrel')
    IDs_mother <- which(sample_info_tab$Population %in% pops_to_use &
                            sample_info_tab$Relationship == 'mother')
    IDs_to_keep <- c(IDs_unrelated, IDs_mother)
    IDs_to_keep <- sample_info_tab$"Individual ID"[IDs_to_keep]
    IDs_to_keep <- sort(IDs_to_keep)

    write.table(IDs_to_keep, 'Euro_keeplist.txt', sep='\n', append=F, quote=F, row.names=F, col.names=F)

    # GBR IDs
    GBR_IDs <-  which(sample_info_tab$Population == 'GBR')
    GBR_IDs <- sample_info_tab$"Individual ID"[GBR_IDs]
    write.table(GBR_IDs, 'GBR_keeplist.txt', sep='\n', append=F, quote=F, row.names=F, col.names=F)

    ##############################
    # Once we downloaded chr20 (wget from the 1000G FTP site, then use vcftools to
    # put it in plink format) and used eigenstrat to get the eigenvectors, plot by cohort
    evecs <- fread('chr20_1000G_euro.pca.evec', header=F, skip=1L)

    # Add the cohort information
    cohort_info <- data.frame(Individual_ID=paste(sample_info_tab$"Individual ID", sample_info_tab$"Individual ID", sep=':'),
                              Population=sample_info_tab$Population)
    cohort_info$Individual_ID <- as.character(cohort_info$Individual_ID)
    colnames(evecs)[1:3] <- c('Individual_ID', 'PC1', 'PC2')
    evecs_with_cohort <- merge(evecs, cohort_info, by="Individual_ID")

    # Check a random person
    temp_sub <- sample(x=1:nrow(evecs_with_cohort), size=5, replace=FALSE)
    temp_ID <- substr(evecs_with_cohort$Individual_ID[temp_sub], 1, gregexpr(':', evecs_with_cohort$Individual_ID[temp_sub])[[1]]-1)
    temp_rows <- which(sample_info_tab$`Individual ID` %in% temp_ID)
    evecs_with_cohort[temp_sub, ]
    sample_info_tab[temp_rows, ]

    # Plot
    ggplot(data=evecs_with_cohort, aes(x=PC1, y=PC2, color=Population)) +
        geom_point()




}
