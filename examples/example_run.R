start_example <- function()
{
    library(data.table)
    library(GBJ)

    setwd('/users/ryansun/documents/research/paper4/reference_data')
    pathways_tab <- fread('pathways_to_test_062817.txt', header=F)
    colnames(pathways_tab) <- c('Pathway_name', 'Pathway_description', paste('Gene', 1:(ncol(pathways_tab)-2), sep=''))
    pathways_tab_fname=NULL
    gene_tab <- fread('refgene_hg19_08302016.txt', header=F)
    colnames(gene_tab) <- c('Gene', 'CHR', 'Start', 'End', 'Other', 'Other2', 'Other3')
    gene_tab_fname=NULL
    SS_file=NULL
    SS_fname_root='summary_stats'
    evecs_tab=fread('gbr_eigenvectors.txt', header=F)
    evecs_tab_fname=NULL
    pathways_per_job=1
    gene_buffer=5000
    threshold_1000G=0.03
    prune_factor=0.5
    prune_limit=0.0625
    snp_limit=1000
    hard_snp_limit=1500
    prune_to_start=TRUE
    run_GHC=TRUE
    out_name_root='pathway_anal'
    refsnp_dir='/users/ryansun/documents/research/paper4/data/BCAC_hg19_summary_stats'
    input_dir='/users/ryansun/documents/research/paper4/data/BCAC_hg19_summary_stats'
    output_dir=NULL
    Snum=1
    aID=2
    checkpoint=TRUE


    libs='/n/home13/rsun/Rlibrary/3.3.1'
    library(GBJ, lib=libs)
    library(data.table, lib=libs)
    library(LungCancerAssoc, lib=libs)

    setwd('/n/home13/rsun/LungCancer/')
    pathways_tab <- fread('IL1B_pathways.txt', header=F)
    colnames(pathways_tab) <- c('Pathway_name', 'Pathway_description', paste('Gene', 1:(ncol(pathways_tab)-2), sep=''))
    pathways_tab_fname=NULL
    gene_tab <- fread('refgene_hg19_08302016.txt', header=F)
    colnames(gene_tab) <- c('Gene', 'CHR', 'Start', 'End', 'Other', 'Other2', 'Other3')
    gene_tab_fname=NULL
    SS_file=NULL
    SS_fname_root='LC_Squam_hg19_chr'
    evecs_tab=fread('chr20_1000G_euro.pca.evec', header=F, skip=1L)
    evecs_tab <- evecs_tab[, -22, with=F]
    evecs_tab_fname=NULL
    pathways_per_job=1
    gene_buffer=5000
    threshold_1000G=0.03
    prune_factor=0.5
    prune_limit=0.0625
    snp_limit=1000
    hard_snp_limit=1500
    prune_to_start=TRUE
    run_GHC=TRUE
    out_name_root='squam_anal'
    refsnp_dir='/n/regal/xlin/ryansun/compressed_1000G'
    input_dir='/n/home13/rsun/heart_disease_SS'
    output_dir='/n/home13/rsun/LungCancer'
    Snum=1
    aID=2
    checkpoint=TRUE



    run_pathwayanal(pathways_tab=pathways_tab, pathways_tab_fname=pathways_tab_fname,
         gene_tab=gene_tab, gene_tab_fname=gene_tab_fname,
         SS_file=SS_file, SS_fname_root=SS_fname_root,
         evecs_tab=evecs_tab, evecs_tab_fname=evecs_tab_fname,
         pathways_per_job=pathways_per_job, gene_buffer=gene_buffer, threshold_1000G=threshold_1000G,
         prune_factor=prune_factor, prune_limit=prune_limit, snp_limit=snp_limit, hard_snp_limit=hard_snp_limit,
         prune_to_start=prune_to_start, run_GHC=run_GHC, out_name_root=out_name_root,
         refsnp_dir=refsnp_dir, input_dir=input_dir, output_dir=output_dir,
         Snum=Snum, aID=aID, checkpoint=checkpoint)


}




