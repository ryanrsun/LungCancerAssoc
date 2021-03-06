% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prune_snps.R
\name{prune_snps}
\alias{prune_snps}
\title{prune_snps.R}
\usage{
prune_snps(Snum, aID, part = NULL, fname_root, prune_R2, temp_Gmat,
  temp_Gmat_record, checkpoint)
}
\arguments{
\item{Snum}{The pruned file will be name S[Snum]_[aID].prune.out. The Snum and aID parameters
allow you to ensure that your pruned files do not overwrite each other if you have many jobs running at once.}

\item{aID}{The pruned file will be name S[Snum]_[aID].prune.out. The Snum and aID parameters
allow you to ensure that your pruned files do not overwrite each other if you have many jobs running at once.}

\item{part}{For run_pathwayanal_part.R the name will be S[Snum]_[aID]_[part].prune.out. Default is NULL.}

\item{fname_root}{The root of our .ped and .info files downloaded from 1000G, used to remove them at the end.}

\item{prune_R2}{Prune all SNPs that have pairwise correlation greater than this threshold.}

\item{temp_Gmat}{The genotypes of the SNPs, will remove columns in accordance with PLINK pruning.}

\item{temp_Gmat_record}{Info on the SNPs, will remove rows in accordance with PLINK pruning.}

\item{checkpoint}{A boolean, if TRUE, print out diagnostic/error messages.}
}
\value{
A list with the elements temp_Gmat (containing the genotypes at each qualifying SNP) and
temp_Gmat_record (containing the info on SNPs), or 1 if nothing to return.
}
\description{
We can only do set-based tests with SNPs that are in both the ped file (to estimate correlation)
and in the summary statistics file (obviously because we need the summary statistic). This function
tells us which SNPs are indeed in both. Use only with one region (contiguous length on one chromosome) at a time.
Make sure you have the PLINK binary in your working directory!
}

