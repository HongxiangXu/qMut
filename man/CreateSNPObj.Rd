% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Mutation_object.R
\name{CreateSNPObj}
\alias{CreateSNPObj}
\title{A SNP object}
\usage{
CreateSNPObj(
  snp_anno,
  snp_data,
  annotate_AA = T,
  metadata = NULL,
  remove_genes_dt = NA,
  verbose = F
)
}
\arguments{
\item{snp_anno}{Required. Generated from the Integrate_data() function.}

\item{snp_data}{Required. Generated from the Integrate_data() function.}

\item{metadata}{Optional.Samples in metadata without corresponding mutation identification will be automatically removed.}

\item{remove_genes_dt}{Optional. A dataframe containing removed genes. The column
of the dataframe should be matched with the example(Mtb_removed_gene_dt). Currently,
if you are working for mutation analysis of Mycobacterium tuberculosis, you can directly use presetted
remove_genes_dt through providing value "Mtb".}
}
\description{
This function create a basic SNP object designed for downstream
SNPs.analysis.
}
