list_all_genes_in_pathways <- function(x) {
    setwd('/users/ryansun/desktop')
    pathways_list <- read.table('IL1B_pathways.txt')

    all_genes <- apply(pathways_list[, 3:ncol(pathways_list)], 2, as.character)
    all_genes <- as.vector(all_genes)
    all_genes <- all_genes[-which(is.na(all_genes))]
    all_genes <- unique(all_genes)
    length(all_genes)
    head(all_genes)
    tail(all_genes)

    write.table(all_genes, 'LC_pathway_genes.txt', row.names=F, col.names=F, append=F, quote=F, sep='\t')
}
