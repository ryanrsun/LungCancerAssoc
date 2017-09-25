# Read in single SNP results for one gene across a set of studies.
setwd("/Users/ryansun/Documents/Research/NewSoftware/Lung_cancer/LungCancerAssoc/data")
gene_tab <- fread('refgene_hg19_08302016.txt', header=F)
colnames(gene_tab) <- c('Gene', 'CHR', 'Start', 'End', 'Other1', 'Other2', 'Other3')

# Find IL1B
buffer <- 5000
IL1B_row <- which(gene_tab$Gene == 'IL1B')
IL1B_start <- as.numeric(gene_tab$Start[IL1B_row]) - buffer
IL1B_end <- as.numeric(gene_tab$End[IL1B_row]) + buffer
IL1B_chr <- as.numeric(gene_tab$CHR[IL1B_row])
cat(IL1B_chr, IL1B_start, IL1B_end, sep=':')
print(gene_tab[IL1B_row])

studies <- c('Squam', 'Small', 'Overall', 'Adeno', 'CAD')
IL1B_results <- NULL
for (study_it in 1:length(studies)) {

    # Set directory
    temp_wd <- paste('/users/ryansun/desktop/LC/', studies[study_it], sep='')
    setwd(temp_wd)

    # Open file
    if (study_it == 5) {
        temp_SS_fname <- paste('CAD_hg19_chr', IL1B_chr, '.txt', sep='')
    } else {
        temp_SS_fname <- paste('LC_', studies[study_it], '_hg19_chr', IL1B_chr, '.txt', sep='')
    }
    temp_SS_file <- fread(temp_SS_fname, header=T)

    # Get the IL1B snps for this study
    temp_snp_rows <- which(temp_SS_file$BP <= IL1B_end &
                               temp_SS_file$BP >= IL1B_start)

    # Not found
    if (length(temp_snp_rows) == 0) {next}

    # Record
    if (study_it == 5) {
        temp_study_results <- data.frame(RS=temp_SS_file$RS[temp_snp_rows],
                                         BP=temp_SS_file$BP[temp_snp_rows],
                                         P=temp_SS_file$"p_dgc"[temp_snp_rows],
                                         Study=studies[study_it])
    } else {
        temp_study_results <- data.frame(RS=temp_SS_file$RS[temp_snp_rows],
                                         BP=temp_SS_file$BP[temp_snp_rows],
                                         P=temp_SS_file$"P-value"[temp_snp_rows],
                                         Study=studies[study_it])
    }


    if (is.null(IL1B_results)) {
        IL1B_results <- temp_study_results
    } else {
        IL1B_results <- rbind(IL1B_results, temp_study_results)
    }

    # Checkpoint
    cat(study_it, '\n')
}

IL1B_results <- IL1B_results[order(IL1B_results$P, decreasing=FALSE), ]
IL1B_results
print(print(xtable(IL1B_results[1:10,], digits=3), include.rownames=F))
cat(length(which(IL1B_results$Study == 'CAD')),
    length(which(IL1B_results$Study == 'Small')),
    length(which(IL1B_results$Study == 'Squam')),
    length(which(IL1B_results$Study == 'Overall')),
    length(which(IL1B_results$Study == 'Adeno')))

setwd('/users/ryansun/dropbox/postdoc/2017-09-14')
write.table(IL1B_results, 'IL1B_single_snp.txt', append=F, quote=F,
            row.names=F, col.names=T, sep='\t')
#########################################################################



#########################################################################
#########################################################################
# General read results function for all pathways analyses
#########################################################################

read_results <- function(working_dir, fname_root, start_aID, end_aID,
                         save_dir=NULL, save_name='') {
    # Bind files
    setwd(working_dir)
    temp_results <- fread(paste(fname_root, start_aID, '.txt', sep=''))
    for (file_it in (start_aID+1):end_aID) {
        temp_fname <- paste(fname_root, file_it, '.txt', sep='')
        temp_file <- tryCatch(fread(temp_fname, header=T),
                              error=function(e) e)
        if (class(temp_file)[1] %in% c('simpleError')) {
            next
        }

        # Add to results
        mylist <- list(temp_results, temp_file)
        temp_results <- rbindlist(mylist)

        # Checkpoint
        if (file_it%%100 == 0) {
            closeAllConnections()
            cat(file_it)
        }
    }

    # Order, NA at the bottom
    temp_results <- temp_results[order(temp_results$GBJ_p, decreasing=FALSE), ]

    # Version with NA removed
    NA_rows <- which(is.na(temp_results$gbj))
    if (length(NA_rows) > 0) {
        NA_remove <- temp_results[-NA_rows, ]
    }

    # Save?
    if (!is.null(save_dir)) {
        setwd(save_dir)
        write.table(temp_results, save_name, append=F, quote=F, row.names=F,
                    col.names=T , sep='\t')
    }

    # Return full version, NA removed version
    return(list(full_results=temp_results,
                NA_remove=NA_remove,
                num_NA=length(NA_rows)))
}
#########################################################################


#################################################################
#################################################################
# Read pathway gene results for small cell
#################################################################
small_pathway_output <- read_results(working_dir='/users/ryansun/desktop/LC_results/small_pathway',
                                       fname_root='small_anal_S11_',
                                       start_aID=1, end_aID=258,
                                       save_dir='/users/ryansun/dropbox/postdoc/2017-09-14',
                                       save_name='small_pathway.txt')
small_pathway_output <- read_results(working_dir='/users/ryansun/desktop/LC_results/small_pathway',
                                     fname_root='small_anal_S11_',
                                     start_aID=1, end_aID=258,
                                     save_dir=NULL,
                                     save_name=NULL)
small_pathway <- small_pathway_output$full_results

# Print
threshold <- 1*10^(-3)
print_row <- max(which(small_pathway$GBJ_p <= threshold))
print(xtable(small_pathway[1:print_row, c(2,3,9,12), with=F],
             digits=c(0,0,0,-4,-4)), include.rownames=F)


#################################################################
#################################################################
# Read pathway gene results for squamous
#################################################################
squamous_pathway_output <- read_results(working_dir='/users/ryansun/desktop/LC_results/squam_pathway',
                                     fname_root='squam_anal_S21_',
                                     start_aID=1, end_aID=258,
                                     save_dir='/users/ryansun/dropbox/postdoc/2017-09-14',
                                     save_name='squam_pathway.txt')
squamous_pathway_output <- read_results(working_dir='/users/ryansun/desktop/LC_results/squam_pathway',
                                        fname_root='squam_anal_S21_',
                                        start_aID=1, end_aID=258,
                                        save_dir=NULL,
                                        save_name=NULL)
squamous_pathway <- squamous_pathway_output$full_results

# Print
threshold <- 1*10^(-4)
print_row <- max(which(squamous_pathway$GBJ_p <= threshold))
print(xtable(squamous_pathway[1:print_row, c(2,3,9,12), with=F],
             digits=c(0,0,0,-4,-4)), include.rownames=F)


#################################################################
#################################################################
# Read pathway gene results for CAD
#################################################################
CAD_pathway_output <- read_results(working_dir='/users/ryansun/desktop/LC_results/CAD_pathway',
                                        fname_root='CAD_anal_S31_',
                                        start_aID=1, end_aID=258,
                                        save_dir='/users/ryansun/dropbox/postdoc/2017-09-14',
                                        save_name='CAD_pathway.txt')
CAD_pathway_output <- read_results(working_dir='/users/ryansun/desktop/LC_results/CAD_pathway',
                                   fname_root='CAD_anal_S31_',
                                   start_aID=1, end_aID=258,
                                   save_dir=NULL,
                                   save_name=NULL)
CAD_pathway <- CAD_pathway_output$full_results

# Print
threshold <- 1*10^(-4)
print_row <- max(which(squamous_pathway$GBJ_p <= threshold))
print(xtable(squamous_pathway[1:print_row, c(2,3,9,12), with=F],
             digits=c(0,0,0,-4,-4)), include.rownames=F)





#################################################################
#################################################################
# Read single gene results for small cell
#################################################################
small_gene_output <- read_results(working_dir='/users/ryansun/desktop/LC_results/small_single',
                                        fname_root='small_singlegene_anal_S12_',
                                        start_aID=1, end_aID=265,
                                        save_dir='/users/ryansun/dropbox/postdoc/2017-09-14',
                                        save_name='small_single.txt')
small_gene_output <- read_results(working_dir='/users/ryansun/desktop/LC_results/small_single',
                                  fname_root='small_singlegene_anal_S12_',
                                  start_aID=1, end_aID=265,
                                  save_dir=NULL,
                                  save_name=NULL)
small_single <- small_gene_output$full_results

# Print
threshold <- 1*10^(-4)
print_row <- max(which(small_single$GBJ_p <= threshold))
print(xtable(small_single[1:print_row, c(1,3,9,12), with=F], digits=c(0,0,0,-4,-4)), include.rownames=F)

# Where is IL1B
temp_row <- which(small_single$Pathway_name == 'IL1B')
small_single[temp_row,]


#################################################################
#################################################################
# Read single gene results for squamous
#################################################################
squamous_gene_output <- read_results(working_dir='/users/ryansun/desktop/LC_results/squam_single',
                                  fname_root='squam_singlegene_anal_S22_',
                                  start_aID=1, end_aID=265,
                                  save_dir='/users/ryansun/dropbox/postdoc/2017-09-14',
                                  save_name='squam_single.txt')
squamous_gene_output <- read_results(working_dir='/users/ryansun/desktop/LC_results/squam_single',
                                     fname_root='squam_singlegene_anal_S22_',
                                     start_aID=1, end_aID=265,
                                     save_dir=NULL,
                                     save_name=NULL)
squamous_single <- squamous_gene_output$full_results

# Print
threshold <- 1*10^(-4)
print_row <- max(which(squamous_single$GBJ_p <= threshold))
print(xtable(squamous_single[1:print_row, c(1,3,9,12), with=F], digits=c(0,0,0,-4,-4)), include.rownames=F)

# Where is IL1B
temp_row <- which(squamous_single$Pathway_name == 'IL1B')
squamous_single[temp_row,]

#################################################################
#################################################################
# Read single gene results for CAD
#################################################################
CAD_gene_output <- read_results(working_dir='/users/ryansun/desktop/LC_results/CAD_single',
                                     fname_root='CAD_singlegene_anal_S32_',
                                     start_aID=1, end_aID=265,
                                     save_dir='/users/ryansun/dropbox/postdoc/2017-09-14',
                                     save_name='CAD_single.txt')
CAD_gene_output <- read_results(working_dir='/users/ryansun/desktop/LC_results/CAD_single',
                                fname_root='CAD_singlegene_anal_S32_',
                                start_aID=1, end_aID=265,
                                save_dir=NULL,
                                save_name=NULL)
CAD_single <- CAD_gene_output$full_results

# Print
threshold <- 1*10^(-4)
print_row <- max(which(CAD_single$GBJ_p <= threshold))
print(xtable(CAD_single[1:print_row, c(1,3,9,12), with=F], digits=c(0,0,0,-4,-4)), include.rownames=F)

# Where is IL1B
temp_row <- which(CAD_single$Pathway_name == 'IL1B')
CAD_single[temp_row,]


#################################################################
#################################################################
# Read single gene no prune results for Small
#################################################################
small_noprune_output <- read_results(working_dir='/users/ryansun/desktop/LC_results/Small_noprune',
                                fname_root='small_singlegene_noprune_S120_',
                                start_aID=1, end_aID=265,
                                save_dir='/users/ryansun/dropbox/postdoc/2017-09-14',
                                save_name='small_noprune.txt')
small_noprune_output <- read_results(working_dir='/users/ryansun/desktop/LC_results/Small_noprune',
                                     fname_root='small_singlegene_noprune_S120_',
                                     start_aID=1, end_aID=265,
                                     save_dir=NULL,
                                     save_name=NULL)
small_noprune <- small_noprune_output$full_results

# Print
threshold <- 1*10^(-4)
print_row <- max(which(small_noprune$GBJ_p <= threshold))
print(xtable(small_noprune[1:print_row, c(1,3,9,12), with=F], digits=c(0,0,0,-4,-4)), include.rownames=F)

# Where is IL1B
temp_row <- which(small_noprune$Pathway_name == 'IL1B')
small_noprune[temp_row,]



#################################################################
#################################################################
# Read single gene no prune results for squamous
#################################################################
squamous_noprune_output <- read_results(working_dir='/users/ryansun/desktop/LC_results/squam_noprune',
                                     fname_root='squam_singlegene_noprune_S220_',
                                     start_aID=1, end_aID=265,
                                     save_dir='/users/ryansun/dropbox/postdoc/2017-09-14',
                                     save_name='squamous_noprune.txt')
squamous_noprune_output <- read_results(working_dir='/users/ryansun/desktop/LC_results/squam_noprune',
                                        fname_root='squam_singlegene_noprune_S220_',
                                        start_aID=1, end_aID=265)
squamous_noprune <- squamous_noprune_output$full_results

# Print
threshold <- 1*10^(-4)
print_row <- max(which(squamous_noprune$GBJ_p <= threshold))
print(xtable(squamous_noprune[1:print_row, c(1,3,9,12), with=F], digits=c(0,0,0,-4,-4)), include.rownames=F)

# Where is IL1B
temp_row <- which(squamous_noprune$Pathway_name == 'IL1B')
squamous_noprune[temp_row,]


#################################################################
#################################################################
# Read single gene no prune results for CAD
#################################################################
CAD_noprune_output <- read_results(working_dir='/users/ryansun/desktop/LC_results/CAD_noprune',
                                        fname_root='CAD_singlegene_noprune_S320_',
                                        start_aID=1, end_aID=265,
                                        save_dir='/users/ryansun/dropbox/postdoc/2017-09-14',
                                        save_name='CAD_noprune.txt')
CAD_noprune_output <- read_results(working_dir='/users/ryansun/desktop/LC_results/CAD_noprune',
                                   fname_root='CAD_singlegene_noprune_S320_',
                                   start_aID=1, end_aID=265)
CAD_noprune <- CAD_noprune_output$full_results

# Print
threshold <- 1*10^(-4)
print_row <- max(which(CAD_noprune$GBJ_p <= threshold))
print(xtable(CAD_noprune[1:print_row, c(1,3,9,12), with=F], digits=c(0,0,0,-4,-4)), include.rownames=F)

# Where is IL1B
temp_row <- which(CAD_noprune$Pathway_name == 'IL1B')
CAD_noprune[temp_row,]


#################################################################
#################################################################
#################################################################
# Aggregate gene significance results across different phenotypes

setwd('/users/ryansun/dropbox/postdoc/2017-09-14')
CAD_noprune <- fread('CAD_noprune.txt', header=F)
small_noprune <- fread('small_noprune.txt', header=F)
squam_noprune <- fread('squamous_noprune.txt', header=F)

# Merge them
all_noprune <- merge(CAD_noprune[, c(1,9), with=F],
                     small_noprune[, c(1,9), with=F], by='Pathway_name', all=TRUE)
all_noprune <- merge(all_noprune,
                     squam_noprune[, c(1,9), with=F], by='Pathway_name', all=TRUE)

# Take a geometric average
avg_p <- apply(all_noprune[, 2:4, with=F], 1, function(x) {prod(x)/length(x)})
all_noprune <- cbind(all_noprune, avg_p)
all_noprune <- all_noprune[order(all_noprune$avg_p, decreasing=FALSE), ]


#################################################################
#################################################################
# Aggregate pathway significance results across different phenotypes

setwd('/users/ryansun/dropbox/postdoc/2017-09-14')
CAD_pathway <- fread('CAD_pathway.txt', header=T)
small_pathway <- fread('small_pathway.txt', header=T)
squam_pathway <- fread('squam_pathway.txt', header=T)

# Merge them
all_pathway <- merge(CAD_pathway[, c(1,9), with=F],
                     small_pathway[, c(1,9), with=F], by='Pathway_name', all=TRUE)
all_pathway <- merge(all_pathway,
                     squam_pathway[, c(1,9), with=F], by='Pathway_name', all=TRUE)

# Take a geometric average
avg_p <- apply(all_pathway[, 2:4, with=F], 1, function(x) {prod(x)/length(x)})
all_pathway <- cbind(all_pathway, avg_p)
all_pathway <- all_pathway[order(all_pathway$avg_p, decreasing=FALSE), ]

##########################################################################################
##########################################################################################

##########################################################################################
##########################################################################################
# Read in top gene p-values across one study - CAD
setwd("/Users/ryansun/Documents/Research/NewSoftware/Lung_cancer/LungCancerAssoc/data")
gene_tab <- fread('refgene_hg19_08302016.txt', header=F)
colnames(gene_tab) <- c('Gene', 'CHR', 'Start', 'End', 'Other1', 'Other2', 'Other3')

top_genes <- c('AIDA', 'SMARCA4', 'CDKN2B', 'APOC1', 'IL6R', 'VAMP8', 'CDKN2A',
               'ATP2B1', 'FURIN', 'APOE', 'GUCY1A3', 'FES', 'TGFB1', 'ATXN2',
               'COL4A2', 'ARNTL', 'SLC22A3', 'NOS3', 'SMAD3', 'AXL', 'IGF2BP1',
               'REST', 'APOB', 'SLC22A2', 'CHST6', 'DNM2', 'LPL', 'PTPN11')

top_gene_results <- NULL
buffer <- 5000
for (gene_it in 1:length(top_genes)) {
    # Find the gene
    temp_gene <- top_genes[gene_it]
    temp_row <- which(gene_tab$Gene == temp_gene)
    temp_start <- as.numeric(gene_tab$Start[temp_row]) - buffer
    temp_end <- as.numeric(gene_tab$End[temp_row]) + buffer
    temp_chr <- as.numeric(gene_tab$CHR[temp_row])
    cat(temp_chr, temp_start, temp_end, sep=':')
    print(gene_tab[temp_row])

    # Open squam cell
    temp_wd <- paste('/users/ryansun/desktop/LC/CAD', sep='')
    setwd(temp_wd)

    temp_SS_fname <- paste('CAD_hg19_chr', temp_chr, '.txt', sep='')
    temp_SS_file <- fread(temp_SS_fname, header=T)

    # Get the snps for this gene
    temp_snp_rows <- which(temp_SS_file$BP <= temp_end &
                               temp_SS_file$BP >= temp_start)

    # Record
    temp_study_results <- data.frame(RS=temp_SS_file$RS[temp_snp_rows],
                                     BP=temp_SS_file$BP[temp_snp_rows],
                                     Beta=temp_SS_file$beta[temp_snp_rows],
                                     SE=temp_SS_file$se_dgc[temp_snp_rows],
                                     P=temp_SS_file$p_dgc[temp_snp_rows],
                                     MAF=temp_SS_file$effect_allele_freq[temp_snp_rows],
                                     Gene=temp_gene)

    if (is.null(top_gene_results)) {
        top_gene_results <- temp_study_results
    } else {
        top_gene_results <- rbind(top_gene_results, temp_study_results)
    }
}



##########################################################################################
##########################################################################################
# Read in top gene p-values across one study - squamous cell LC
setwd("/Users/ryansun/Documents/Research/NewSoftware/Lung_cancer/LungCancerAssoc/data")
gene_tab <- fread('refgene_hg19_08302016.txt', header=F)
colnames(gene_tab) <- c('Gene', 'CHR', 'Start', 'End', 'Other1', 'Other2', 'Other3')

top_genes <- c('HLA-DQA1', 'TRIM38', 'HIST1H1B', 'ZKSCAN3', 'HIST1H2BB', 'SLC17A3', 'CYP2A6',
               'CHUK', 'EDN2', 'ABCG2', 'PKD2L1', 'XBP1', 'CHEK1', 'KDM4A',
               'ERLIN1', 'PEA15', 'TFAP2A')


top_gene_results <- NULL
buffer <- 5000
for (gene_it in 1:length(top_genes)) {
    # Find the gene
    temp_gene <- top_genes[gene_it]
    temp_row <- which(gene_tab$Gene == temp_gene)
    temp_start <- as.numeric(gene_tab$Start[temp_row]) - buffer
    temp_end <- as.numeric(gene_tab$End[temp_row]) + buffer
    temp_chr <- as.numeric(gene_tab$CHR[temp_row])
    cat(temp_chr, temp_start, temp_end, sep=':')
    print(gene_tab[temp_row])

    # Open squam cell
    temp_wd <- paste('/users/ryansun/desktop/LC/squam', sep='')
    setwd(temp_wd)

    temp_SS_fname <- paste('LC_squam_hg19_chr', temp_chr, '.txt', sep='')
    temp_SS_file <- fread(temp_SS_fname, header=T)

    # Get the snps for this gene
    temp_snp_rows <- which(temp_SS_file$BP <= temp_end &
                               temp_SS_file$BP >= temp_start)

    temp_beta <- log(temp_SS_file$OR_fixed[temp_snp_rows])
    temp_pvalue <- temp_SS_file$"P-value"[temp_snp_rows]
    temp_stat <- sqrt(qchisq(1-temp_pvalue, df=1))
    temp_se <- abs(temp_beta / temp_stat)


    # Record
    temp_study_results <- data.frame(RS=temp_SS_file$RS[temp_snp_rows],
                                     BP=temp_SS_file$BP[temp_snp_rows],
                                     OR=temp_SS_file$OR_fixed[temp_snp_rows],
                                     Beta=log(temp_SS_file$OR_fixed[temp_snp_rows]),
                                     SE=temp_se,
                                     P=temp_SS_file$"P-value"[temp_snp_rows],
                                     EAF=temp_SS_file$EAF[temp_snp_rows],
                                     Gene=temp_gene)

    if (is.null(top_gene_results)) {
        top_gene_results <- temp_study_results
    } else {
        top_gene_results <- rbind(top_gene_results, temp_study_results)
    }
}

top_rows <- which(top_gene_results$P < 1*10^(-5))
setwd('/users/ryansun/desktop')
write.table(top_gene_results[top_rows, ], 'top_squamous_snps_updated.txt',
            append=F, quote=F, row.names=F, col.names=T, sep='\t')

##################################################################
# Final squamous cells reports combining prune and no prune
squamous_noprune_output <- read_results(working_dir='/users/ryansun/desktop/LC_results/squam_noprune',
                                     fname_root='squam_singlegene_noprune_S220_',
                                     start_aID=1, end_aID=265,
                                     save_dir='/users/ryansun/dropbox/postdoc/2017-09-14',
                                     save_name='squam_noprune.txt')

squamous_gene_output <- read_results(working_dir='/users/ryansun/desktop/LC_results/squam_single',
                                     fname_root='squam_singlegene_anal_S22_',
                                     start_aID=1, end_aID=265,
                                     save_dir='/users/ryansun/dropbox/postdoc/2017-09-14',
                                     save_name='squam_single.txt')

setwd('/users/ryansun/dropbox/postdoc/2017-09-14')
squam_noprune <- fread('squam_noprune.txt', header=T)
squam_single <- fread('squam_single.txt', header=T)

# Only use pruned for those that didn't finish no prune
use_single_rows <- which(!(squam_single$Pathway_name %in% squam_noprune$Pathway_name))
squam_total <- rbind(squam_noprune, squam_single[use_single_rows, ])

# Write
setwd('/users/ryansun/desktop')
write.table(squam_total, 'all_squam_genes.txt', append=F, quote=F, row.names=F,
            col.names=T , sep='\t')


# What's missing?
setwd("/Users/ryansun/Documents/Research/NewSoftware/Lung_cancer/LungCancerAssoc/data")
single_gene_pathways <- fread('single_gene_pathways.txt', header=F)
missing_rows <- which(!(single_gene_pathways$V1 %in% squam_total$Pathway_name))
unique(ceiling(missing_rows / 100))



##################################################################
setwd('/users/ryansun/desktop/LC/CAD')
all_CAD <- fread('cad.add.160614.website.txt')
all_CAD <- all_CAD[, c(1,2,3,11), with=F]
colnames(all_CAD) <- c('RS', 'CHR', 'BP', 'p_CAD')
setwd('/users/ryansun/desktop/LC/squam')
all_squam <- fread('LC_Squam_hg19_chr1.txt')
all_squam <- all_squam[, c(2,3,4,10), with=F]
colnames(all_squam) <- c('RS', 'CHR_LC', 'BP_LC', 'p_LC')
for (i in 2:22) {
    temp_fname <- paste('LC_Squam_hg19_chr', i, '.txt', sep='')
    temp_squam <- fread(temp_fname)
    temp_squam <- temp_squam[, c(2,3,4,10), with=F]
    colnames(temp_squam) <- c('RS', 'CHR_LC', 'BP_LC', 'p_LC')
    mylist <- list(all_squam, temp_squam)
    all_squam <- rbindlist(mylist)
    cat(i)
}

# Merge, just leave RS and chr and BP and p-value
merged_data <- merge(all_squam, all_CAD, by='RS')
dim(merged_data)

# Just keep those with p-value < 0.01 in both cohorts
sig_rows <- which(merged_data$p_LC < 0.01 | merged_data$p_CAD < 0.01)
sig_data <- merged_data[sig_rows, ]
dim(sig_data)

# Now plot their p-values
ylabel <- expression(paste(-log[10], '(CAD_p)', sep=''))
xlabel <- expression(paste(-log[10], '(Squamous_p)', sep=''))
ggplot(sig_data, aes(x=-log10(p_LC), y=-log10(p_CAD))) +
    geom_point(color='blue') +
    xlab(xlabel) +
    ylab(ylabel) +
    ggtitle('All Common SNPs')
setwd('/users/ryansun/desktop')
ggsave('all_SNPs.png')

#########################################################################
# Now just those in LC risk genes
top_genes <- c('PSMA4', 'HFE', 'TERT', 'HLA-DQA1', 'TRIM38', 'HIST1H1B', 'ZKSCAN3', 'HIST1H2BB', 'SLC17A3', 'CYP2A6',
               'CHUK', 'EDN2', 'ABCG2', 'PKD2L1', 'XBP1', 'CHEK1', 'KDM4A',
               'ERLIN1', 'PEA15', 'TFAP2A')
setwd('/users/ryansun/documents/research/newsoftware/lung_cancer/lungcancerassoc/data')
gene_tab <- fread('refgene_hg19_08302016.txt')
colnames(gene_tab) <- c('Gene', 'CHR', 'Start', 'End', 'Other1', 'Other2', 'Other3')
gene_tab$CHR <- as.numeric(as.character(gene_tab$CHR))

top_snps <- NULL
for (gene_it in 1:length(top_genes)) {
    temp_gene <- top_genes[gene_it]
    temp_row <- which(gene_tab$Gene == temp_gene)
    temp_chr <- gene_tab$CHR[temp_row]
    temp_start <- gene_tab$Start[temp_row] - 5000
    temp_end <- gene_tab$End[temp_row] + 5000
    temp_snps <- which(merged_data$BP >= temp_start & merged_data$BP <= temp_end & merged_data$CHR == temp_chr)
    temp_data <- merged_data[temp_snps, ]
    temp_data <- cbind(temp_data, rep(temp_gene, length(temp_snps)))

    if (is.null(top_snps)) {
        top_snps <- temp_data
    } else {
        mylist <- list(top_snps, temp_data)
        top_snps <- rbindlist(mylist)
    }
    cat(gene_it, '\n')
}
dim(top_snps)

# Now plot their p-values
ylabel <- expression(paste(-log[10], '(CAD_p)', sep=''))
xlabel <- expression(paste(-log[10], '(Squamous_p)', sep=''))
ggplot(top_snps, aes(x=-log10(p_LC), y=-log10(p_CAD))) +
    geom_point(color='blue') +
    xlab(xlabel) +
    ylab(ylabel) +
    ggtitle('Top Squamous Cell Risk Genes')
setwd('/users/ryansun/desktop')
ggsave('top_squamous_risk_genes.png')


############################################################
# Now just those in the IL1B pathway list
setwd('/users/ryansun/documents/research/newsoftware/lung_cancer/lungcancerassoc/data')
IL1B_pathways <- fread('IL1B_pathways.txt')
all_genes <- unlist(IL1B_pathways[, 3:ncol(IL1B_pathways), with=F])
all_genes <- all_genes[-which(is.na(all_genes))]
top_genes <- unique(all_genes)
length(top_genes)
setwd('/users/ryansun/documents/research/newsoftware/lung_cancer/lungcancerassoc/data')
gene_tab <- fread('refgene_hg19_08302016.txt')
colnames(gene_tab) <- c('Gene', 'CHR', 'Start', 'End', 'Other1', 'Other2', 'Other3')
gene_tab$CHR <- as.numeric(as.character(gene_tab$CHR))
gene_tab <- gene_tab[which(gene_tab$Gene %in% top_genes), ]
gene_tab <- gene_tab[order(gene_tab$CHR), ]

IL1B_top_snps <- NULL
prev_chr <- -1
for (gene_it in 1:nrow(gene_tab)) {
    temp_gene <- gene_tab$Gene[gene_it]
    temp_chr <- gene_tab$CHR[gene_it]
    temp_start <- gene_tab$Start[gene_it] - 5000
    temp_end <- gene_tab$End[gene_it] + 5000

    # Go by chromosome to speed up
    if (temp_chr != prev_chr) {
        if (temp_chr > 22) {next}
        temp_chr_data <- merged_data[which(merged_data$CHR == temp_chr), ]
        prev_chr <- temp_chr
    }

    temp_snps <- which(temp_chr_data$BP >= temp_start & temp_chr_data$BP <= temp_end)
    if (length(temp_snps) == 0) {next}
    temp_data <- temp_chr_data[temp_snps, ]
    temp_data <- cbind(temp_data, rep(temp_gene, length(temp_snps)))

    if (is.null(IL1B_top_snps)) {
        IL1B_top_snps <- temp_data
    } else {
        mylist <- list(IL1B_top_snps, temp_data)
        IL1B_top_snps <- rbindlist(mylist)
    }
    if (gene_it%%200 == 0) {cat(gene_it, '\n')}
}
dim(IL1B_top_snps)

# Now plot their p-values
ylabel <- expression(paste(-log[10], '(CAD_p)', sep=''))
xlabel <- expression(paste(-log[10], '(Squamous_p)', sep=''))
ggplot(IL1B_top_snps, aes(x=-log10(p_LC), y=-log10(p_CAD))) +
    geom_point(color='blue') +
    xlab(xlabel) +
    ylab(ylabel) +
    ggtitle('All SNPs in IL1B Pathways')
setwd('/users/ryansun/desktop')
ggsave('IL1B_pathway_SNPs.png')

