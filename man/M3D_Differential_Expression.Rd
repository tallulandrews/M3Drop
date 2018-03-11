\name{M3DropFeatureSelection}
\alias{M3DropFeatureSelection}
\alias{bg__test_DE_K_equiv}
\title{Differentially Expressed Genes.}
\description{
	Use Michaelis-Menten curve to find differentially expressed (DE) genes. Functions tagged with "bg__" are not meant for direct usage and are not available in the Bioconductor release.
}
\usage{
  M3DropFeatureSelection(expr_mat, mt_method="bon", mt_threshold=0.05, suppress.plot=FALSE, xlim=NA)
  bg__test_DE_K_equiv(gene_info, fit=NA)
}
\arguments{
  \item{expr_mat}{a numeric matrix of normalized (not log-transformed) expression values, columns = samples, rows = genes.}
  \item{mt_method}{the multiple testing method used in p.adjust}
  \item{mt_threshold}{the threshold for identifying significantly DE genes.}
  \item{suppress.plot}{logical, whether to plot the fitted curve and highlight selected features.}
  \item{xlim}{specify the limits of the x-axis of the plot.}
  \item{fit}{Output from fitting the Michaelis-Menten equation (see: \code{bg__fit_MM}}
  \item{gene_info}{List of calculated gene-specific values output by \code{bg__calc_variables}}
}
\details{Fits a Michaelis-Menten function to the dropout-rate (if not provided) of the provided expression matrix. Identifies genes where the gene-specific K calculated as (S = mean expression, P = dropout rate): \deqn{K=\frac{S*P}{1-P}}{K = (S * P)/(1 - P)} is significantly larger than the K fitted to the entire dataset. Combines standard errors of the fitted K, the gene-specific dropout rate and the gene-specific average expression using error propagation rules. Determines the signficance of the gene-specific K using a Z-test of the log-transformed Ks with the propagated error then applies the specified multiple testing correction to identify DE genes. Plots the dropout rate vs gene expression with the fitted MM curve and highlights in purple the significantly DE genes.}
\value{ 
  \code{M3DropFeatureSelection} : a data.frame of significantly differentially expressed genes with columns: Gene, p.value, q.value
  \code{bg__test_DE_K_equiv} : a named list of containing: 
	pval : the significance of differential expression for each gene
	fold_change : ratio of gene-specific K to globally fit K
}
\examples{
	library(M3DExampleData)
	norm <- M3DropConvertData(Mmus_example_list$data, is.counts=TRUE)
	DE_genes <- M3DropFeatureSelection(norm, 
		mt_method="fdr", mt_threshold=0.01)
#  gene_info <- bg__calc_variables(Normalized_data$data[1:1000,])
#  DE_output <- bg__test_DE_K_equiv(gene_info)
}
\keyword{differential expression}
