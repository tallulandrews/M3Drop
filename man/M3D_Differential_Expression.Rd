\name{M3D_Differential_Expression}
\alias{M3D_Differential_Expression}
\alias{bg__test_DE_K_equiv}
\title{Differentially Expressed Genes.}
\description{
	Use Michaelis-Menten curve to find differentially expressed (DE) genes.
}
\usage{
  M3D_Differential_Expression(expr_mat, mt_method="bon", mt_threshold=0.05, suppress.plot=FALSE)
  bg__test_DE_K_equiv(expr_mat, fit=NA)
}
\arguments{
  \item{expr_mat}{a numeric matrix of expression values, columns = samples, rows = genes.}
  \item{mt_method}{the multiple testing method used in p.adjust}
  \item{mt_threshold}{the threshold for identifying significantly DE genes.}
  \item{suppress.plot}{logical, whether to plot the fitted curve and highlight DE genes.}
  \item{fit}{the output from fitting the Michaelis-Menten curve.}
}
\details{Fits a Michaelis-Menten function to the dropout-rate (if not provided) of the provided expression matrix. Identifies genes where the gene-specific K calculated as (S = mean expression, P = dropout rate): \deqn{K=\frac{S*P}{1-P}}{K = (S * P)/(1 - P)} is significantly larger than the K fitted to the entire dataset. Combines standard errors of the fitted K, the gene-specific dropout rate and the gene-specific average expression using error propagation rules. Determines the signficance of the gene-specific K using a Z-test with the propagated error then applies the specified multiple testing correction to identify DE genes. Plots the dropout rate vs gene expression with the fitted MM curve and highlights in purple the significantly DE genes.}
\value{ 
  \code{M3D_Differential Expression} : a data.frame of significantly differentially expressed genes with columns: Gene, p.value, q.value
  \code{bg__test_DE_K_equiv} : named list of raw p-values (pval) and relative K size (fold_change) for all genes.
}
\examples{
  library(M3DExampleData)
  DE_genes = M3D_Differential_Expression(Mmus_example_list$data, suppress.plot=TRUE)
  DE_genes = M3D_Differential_Expression(Mmus_example_list$data)
  DE_test_results = bg__test_DE_K_equiv(Mmus_example_list$data);
}
\seealso{
  \code{\link{bg__fit_MM}}
}
\keyword{differential expression}
