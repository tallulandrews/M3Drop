\name{M3Drop_Differential_Expression}
\alias{M3Drop_Differential_Expression}
\title{Differentially Expressed Genes.}
\description{
	Use Michaelis-Menten curve to find differentially expressed (DE) genes.
}
\usage{
  M3Drop_Differential_Expression(expr_mat, mt_method="bon", mt_threshold=0.05, suppress.plot=FALSE)
}
\arguments{
  \item{expr_mat}{a numeric matrix of expression values, columns = samples, rows = genes.}
  \item{mt_method}{the multiple testing method used in p.adjust}
  \item{mt_threshold}{the threshold for identifying significantly DE genes.}
  \item{suppress.plot}{logical, whether to plot the fitted curve and highlight DE genes.}
}
\details{Fits a Michaelis-Menten function to the dropout-rate (if not provided) of the provided expression matrix. Identifies genes where the gene-specific K calculated as (S = mean expression, P = dropout rate): \deqn{K=\frac{S*P}{1-P}}{K = (S * P)/(1 - P)} is significantly larger than the K fitted to the entire dataset. Combines standard errors of the fitted K, the gene-specific dropout rate and the gene-specific average expression using error propagation rules. Determines the signficance of the gene-specific K using a Z-test of the log-transformed Ks with the propagated error then applies the specified multiple testing correction to identify DE genes. Plots the dropout rate vs gene expression with the fitted MM curve and highlights in purple the significantly DE genes.}
\value{ 
  \code{M3Drop_Differential Expression} : a data.frame of significantly differentially expressed genes with columns: Gene, p.value, q.value
}
\examples{
  library(M3DExampleData)
  DE_genes = M3Drop_Differential_Expression(Mmus_example_list$data, suppress.plot=TRUE)
  DE_genes = M3Drop_Differential_Expression(Mmus_example_list$data)
}
\keyword{differential expression}
