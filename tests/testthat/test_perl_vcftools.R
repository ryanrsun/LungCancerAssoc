download_script <- function() {
    # Testing the LDAH gene
    pop_vec <- c('GBR', 'TSI', 'CEU', 'FIN', 'IBS')
    CHR <- 2
    start_bp <- 20784123
    end_bp <- 21122046
    chr_vcf_link <- paste('ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr', CHR,
                          '.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz', sep='')
    panel_link <- 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel'
    region_string <- paste(CHR, ':', start_bp, '-', end_bp, sep='')

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
}



test_that("1000G .pl script and vcftools produce the same genotypes", {

    # Not for CRAN
    skip_on_cran()

    # Testing the LDAH gene
    #download_script()

    CHR <- 2
    start_bp <- 20784123
    end_bp <- 21122046
    init_ped_name <- paste(CHR, '_', start_bp, '-', end_bp, '.ped', sep='')
    init_info_name <- paste(CHR, '_', start_bp, '-', end_bp, '.info', sep='')

    perl_map <- read.table(init_info_name, header=F)
    perl_ped <- read.table(init_ped_name, header=F)

    # Open the VCFtools file
    vcf_map <- read.table('LDAH.map', header=F)
    vcf_ped <- read.table('LDAH.ped', header=F)

    # Keep only the overlapping snps (some not in vcf because filters for MAF, missingness, etc.)
    dim(vcf_map)
    dim(perl_map)
    vcf_keep_rows <- which(vcf_map$V2 %in% perl_map$V1)
    perl_keep_rows <- which(perl_map$V1 %in% vcf_map$V2)
    length(vcf_keep_rows)
    length(perl_keep_rows)

    # Trim ped files
    perl_keep_cols1 <- c(1:6, perl_keep_rows*2+5)
    perl_keep_cols2 <- perl_keep_rows*2+6
    perl_keep_cols <- sort(c(perl_keep_cols1, perl_keep_cols2), decreasing=FALSE)
    vcf_keep_cols1 <- c(1:6, vcf_keep_rows*2+5)
    vcf_keep_cols2 <- vcf_keep_rows*2+6
    vcf_keep_cols <- sort(c(vcf_keep_cols1, vcf_keep_cols2), decreasing=FALSE)

    perl_ped <- perl_ped[, perl_keep_cols]
    vcf_ped <- vcf_ped[, vcf_keep_cols]

    # Order by individuals
    perl_keep_ind <- which(perl_ped$V2 %in% vcf_ped$V2)
    perl_ped <- perl_ped[perl_keep_ind, ]
    perl_ped <- perl_ped[order(perl_ped$V2, decreasing=FALSE), ]
    vcf_ped <- vcf_ped[order(vcf_ped$V2, decreasing=FALSE), ]

    # Loop and check
    for (i in 7:ncol(perl_ped)) {
        temp_col_perl <- perl_ped[, i]
        temp_col_vcf <- vcf_ped[, i]
        if (sum(temp_col_perl == temp_col_vcf) != length(temp_col_perl)) {
            cat(i)
        }
    }
})



