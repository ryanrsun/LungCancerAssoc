# Split files by chromosome
# Have squam, overall, adeno, and small
split_LC_SS <- function() {
    histology <- 'Squam'

    setwd(paste('/users/ryansun/desktop/LC/', histology, sep=''))

    all_SS <- fread(paste(histology, '.txt', sep=''), header=F, sep=',')
    all_SS <- all_SS[, 1:23, with=F]
    colnames(all_SS) <- c('chr:pos', 'RS', 'Chr', 'BP', 'reference_allele',
                          'effect_allele', 'EAF', 'OR_fixed', 'StdError_fixed', 'P-value',
                          'L95_fixed', 'U95_fixed', 'phete(Q)_Pvalue', 'I2', 'N_study', 'N_Cases',
                          'N_controls', 'Effects', 'OR_random', 'StdError_random', 'Pvalue_random',
                          'L95_random', 'U95_random')
    for (temp_chr in 1:22) {
        temp_data <- all_SS[which(all_SS$CHR == temp_chr), ]
        temp_data <- temp_data[order(temp_data$position, decreasing=F), ]
        temp_name <- paste('LC_', histology, '_hg19_chr', temp_chr, '.txt', sep='')
        write.table(temp_data, temp_name, append=F, quote=F, row.names=F, col.names=T, sep='\t')
        cat(dim(temp_data), '\n')
    }
}

# Split CAD
# Have squam, overall, adeno, and small
split_CAD_SS <- function() {
    setwd('/users/ryansun/desktop/LC/CAD')

    all_SS <- fread('cad.add.160614.website.txt', header=T)
    colnames(all_SS) <- c('RS', 'Chr', 'BP', 'effect_allele', 'noneffect_allele',
                          'effect_allele_freq', 'median_info', 'model', 'beta',
                          'se_dgc', 'p_dgc', 'het_pvalue', 'n_studies')

    for (temp_chr in 3:22) {
        temp_data <- all_SS[which(all_SS$Chr == temp_chr), ]
        temp_data <- temp_data[order(temp_data$BP, decreasing=F), ]
        temp_name <- paste('CAD_hg19_chr', temp_chr, '.txt', sep='')
        write.table(temp_data, temp_name, append=F, quote=F, row.names=F, col.names=T, sep='\t')
        cat(dim(temp_data), '\n')
    }
}


