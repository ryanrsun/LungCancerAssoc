# Use VCFtools to get the PLINK files (from 1000G vcf) for a list of genes

vcf_to_plink <- function() {
    args <- commandArgs(trailingOnly=TRUE)
    aID <- as.numeric(args[1])

    output_dir <- '/n/regal/xlin/ryansun/1000G_genes'
    gene_list <- read.table('LC_pathway_genes.txt')
    gene_info <- read.table('refgene_hg19_08302016.txt')
    colnames(gene_info) <- c('Gene', 'Chr', 'Start', 'End', 'other1', 'other2', 'other3')

    # Get gene info
    buffer <- 100000
    temp_gene <- as.character(gene_list[aID, 1])
    temp_row <- which(gene_info$Gene == temp_gene)
    if (length(temp_row) != 1) {stop(paste(temp_gene, ' not found', sep=''))}

    temp_chr <- as.numeric(as.character(gene_info$Chr[temp_row]))
    temp_start <- as.numeric(as.character(gene_info$Start[temp_row])) - buffer
    temp_end <- as.numeric(as.character(gene_info$End[temp_row])) + buffer

    input_name <- paste('ALL.chr', temp_chr, '.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz', sep='')
    temp_success <- system2(command='vcftools', args=c('--gzvcf', input_name, '--plink', '--out', temp_gene,
                                                       '--chr', temp_chr, '--from-bp', temp_start, '--to-bp', temp_end,
                                                       '--keep', 'Euro_keeplist.txt',
                                                       '--remove-indels', '--maf', '0.03', '--min-alleles', 2,
                                                       '--max-alleles', 2, '--max-missing', 1), wait=TRUE)

    cat(temp_success)

}
