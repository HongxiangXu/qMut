% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Mutation_object.R
\name{SNP}
\alias{SNP}
\title{The Mutation object for snp_data and snp_anno}
\description{
The Mutation object is created to apply all functions for mutation analysis
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{snp_list}}{A numeric field}

\item{\code{metadata}}{A numeric field}

\item{\code{modification_info}}{A non-numeric field}

\item{\code{nj_distance}}{NJ distance
Initialize a new MyClass object}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-Mutation Object-new}{\code{SNP$new()}}
\item \href{#method-Mutation Object-get_metadata}{\code{SNP$get_metadata()}}
\item \href{#method-Mutation Object-ML_data_dt}{\code{SNP$ML_data_dt()}}
\item \href{#method-Mutation Object-Evolution_dt}{\code{SNP$Evolution_dt()}}
\item \href{#method-Mutation Object-filter_snp}{\code{SNP$filter_snp()}}
\item \href{#method-Mutation Object-filter_sample}{\code{SNP$filter_sample()}}
\item \href{#method-Mutation Object-estimate_imp_all}{\code{SNP$estimate_imp_all()}}
\item \href{#method-Mutation Object-estimate_imp}{\code{SNP$estimate_imp()}}
\item \href{#method-Mutation Object-MutTree}{\code{SNP$MutTree()}}
\item \href{#method-Mutation Object-MutNumPlot}{\code{SNP$MutNumPlot()}}
\item \href{#method-Mutation Object-calculate_NJ}{\code{SNP$calculate_NJ()}}
\item \href{#method-Mutation Object-clone}{\code{SNP$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Mutation Object-new"></a>}}
\if{latex}{\out{\hypertarget{method-Mutation Object-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SNP$new(
  snp_anno,
  snp_data,
  metadata = NULL,
  remove_genes_dt = NA,
  annotate_AA = T,
  verbose = F
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{snp_anno}}{Initial value for the field snp_anno}

\item{\code{snp_data}}{Initial value for the field snp_data}

\item{\code{metadata}}{Initial value for the field metadata}

\item{\code{remove_genes_dt}}{Initial value for the field remove_genes_dt}

\item{\code{annotate_AA}}{Initial value for the field remove_genes_dt}

\item{\code{verbose}}{Initial value for the field remove_genes_dt
get_metadata}

\item{\code{modification_info}}{Initial value for the field remove_genes_dt}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Mutation Object-get_metadata"></a>}}
\if{latex}{\out{\hypertarget{method-Mutation Object-get_metadata}{}}}
\subsection{Method \code{get_metadata()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SNP$get_metadata()}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{Initialize}}{a new MyClass object}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Mutation Object-ML_data_dt"></a>}}
\if{latex}{\out{\hypertarget{method-Mutation Object-ML_data_dt}{}}}
\subsection{Method \code{ML_data_dt()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SNP$ML_data_dt(selected_pheno)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{selected_pheno}}{Initial value for the field selected_pheno}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
NULL
Initialize a new MyClass object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Mutation Object-Evolution_dt"></a>}}
\if{latex}{\out{\hypertarget{method-Mutation Object-Evolution_dt}{}}}
\subsection{Method \code{Evolution_dt()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SNP$Evolution_dt(x)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{Initial value for the field x
Getting SNPs in certain Genes, or certain SNP in all samples}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Mutation Object-filter_snp"></a>}}
\if{latex}{\out{\hypertarget{method-Mutation Object-filter_snp}{}}}
\subsection{Method \code{filter_snp()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SNP$filter_snp(gene, gene_type = "INDEX", upstream_length = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{gene}}{values for Index , Locus_tag or Gene}

\item{\code{gene_type}}{values for searching mutations by certain type, should be one of
"GENE","LOCUS_TAG","INDEX"}

\item{\code{upstream_length}}{whether to show the upstream mutation of genes}

\item{\code{length}}{Length of upstream and downstream regions of genes that should be included
Initialize a new MyClass object}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Mutation Object-filter_sample"></a>}}
\if{latex}{\out{\hypertarget{method-Mutation Object-filter_sample}{}}}
\subsection{Method \code{filter_sample()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SNP$filter_sample(sample_name = NULL, mutation_index = NULL, opposite = F)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{sample_name}}{sample names to filter}

\item{\code{mutation_index}}{mutation index.}

\item{\code{opposite}}{Used with mutation_index. When opposite is set to TRUE, samples
containing mutation_index will be removed; otherwise will be remained.
Getting SNPs counts in all samples}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Mutation Object-estimate_imp_all"></a>}}
\if{latex}{\out{\hypertarget{method-Mutation Object-estimate_imp_all}{}}}
\subsection{Method \code{estimate_imp_all()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SNP$estimate_imp_all()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Mutation Object-estimate_imp"></a>}}
\if{latex}{\out{\hypertarget{method-Mutation Object-estimate_imp}{}}}
\subsection{Method \code{estimate_imp()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SNP$estimate_imp(
  gene,
  gene_type = "GENE",
  metadata_col = NULL,
  upstream_length = NULL
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{gene}}{values for Index , Locus_tag or Gene}

\item{\code{gene_type}}{values for searching mutations by certain type, should be one of
"GENE","LOCUS_TAG","INDEX"}

\item{\code{metadata_col}}{Default is NULL. If the metadata_col was set, the estimate_imp function
will calculate significance difference of mutation counts among different metadata_col levels.
Should be colnames in metadata or the mutation INDEX.}

\item{\code{upstream_length}}{whether to show the upstream mutation of genes
MutTree
plot phylogenetic tree with heatmap of snps and lineage composition}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Mutation Object-MutTree"></a>}}
\if{latex}{\out{\hypertarget{method-Mutation Object-MutTree}{}}}
\subsection{Method \code{MutTree()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SNP$MutTree(locus_tag, limit = 5, size = 20, metadata_col = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{locus_tag}}{the locus_tag of wanted gene}

\item{\code{limit}}{the minimum frequence a sample need to have,defaults to 5}

\item{\code{size}}{the number of top occured mutation types to plot,defaults to 20}

\item{\code{metadata_col}}{Default is NULL. If the metadata_col was set, the estimate_imp function
will calculate significance difference of mutation counts among different metadata_col levels.
Should be colnames in metadata or the mutation INDEX.
MutNumPlot
plot mutation times and modification types in a single gene.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Mutation Object-MutNumPlot"></a>}}
\if{latex}{\out{\hypertarget{method-Mutation Object-MutNumPlot}{}}}
\subsection{Method \code{MutNumPlot()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SNP$MutNumPlot(
  locus_tag,
  limit = 1,
  label_num = 10,
  show_synomonous = T,
  Modification = NULL
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{locus_tag}}{the locus_tag of wanted gene}

\item{\code{limit}}{the minimum frequence a sample need to have,defaults to 1}

\item{\code{label_num}}{the number of top frequent mutationsto label, defaults to 10}

\item{\code{show_synomonous}}{Whether to show synomonous mutations.}

\item{\code{Modification}}{Modification types to show.
calculate_NJ
calculate sample distance for object based on Neighbor-Joining Algorithm.}

\item{\code{metadata_col}}{Default is NULL. If the metadata_col was set, the estimate_imp function}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Mutation Object-calculate_NJ"></a>}}
\if{latex}{\out{\hypertarget{method-Mutation Object-calculate_NJ}{}}}
\subsection{Method \code{calculate_NJ()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SNP$calculate_NJ(num_thread = 1, snp_range = c(0.05, 0.95))}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{num_thread}}{number of threads to use}

\item{\code{snp_range}}{SNPs within this frequency range were included in the NJ (Neighbor-Joining) distance calculation. This range was used to accelerate computational speed.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Mutation Object-clone"></a>}}
\if{latex}{\out{\hypertarget{method-Mutation Object-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{SNP$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
