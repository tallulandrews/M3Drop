\name{M3D_Differential_Expression}
\alias{M3D_Differential_Expression}
\title{Differentially Expressed Genes.}
\description{
	Use Michaelis-Menten curve to find differentially expressed (DE) genes.
}
\usage{M3D_Differential_Expression(expr_mat, weights=1, mt_method="bon", mt_threshold=0.05)}
\arguments{
  \item{expr_mat}{a numeric matrix of expression values, columns = samples, rows = genes.}
  \item{weights}{a numeric matrix of same dimension as expr_mat of weights for each expression value.}
  \item{mt_method}{the multiple testing method used in p.adjust}
  \item{mt_threshold}{the threshold for identifying significantly DE genes.}
}
\details{Fits a Michaelis-Menten function to the dropout-rate of the provided expression matrix. Identifies genes where the gene-specific K calculated as (S = mean expression): $$\frac{S}{S+K}$$ is significantly larger than the K fitted to the entire dataset. Combines standard errors of the fitted K, the gene-specific dropout rate and the gene-specific average expression using error propagation rules. Determines the signficance of the gene-specific K using a Z-test with the propagated error then applies the specified multiple testing correction to identify DE genes. Plots the dropout rate vs gene expression with the fitted MM curve and highlights in purple the significantly DE genes.}
\value{data.frame of significantly differentially expressed genes with columns: Gene, p.value, q.value}
\examples{
  M3D_Differential_Expression(Test_Cleaned$data)
}
\seealso{bg__fit_MM, bg__test_DE_K_equiv, p.adjust}
\keyword{differential expression}
