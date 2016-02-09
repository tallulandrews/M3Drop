\name{test_DE}
\alias{bg__test_DE_K_equiv}
\alias{bg__test_DE_P_equiv}
\alias{bg__test_DE_S_equiv}
\title{Identify Differentially Expressed Genes}
\description{Background function. Identifies differentially expressed genes by fitting the Michaelis-Menten model.}
\usage{bg__test_DE_K_equiv(expr_mat, weights = 1, fit = NA)}
\arguments{
  \item{expr_mat}{Numerical matrix of gene expression values. columns = cells, rows = genes.}
  \item{weights}{Numerical matrix of the same dimension as expr_mat with weights for each expression value.}
  \item{fit}{Output from fitting a Michaelis-Menten model.}
}
\details{Combines standard errors on various measures using error propatation rules. Identifies differentially expressed genes by fitting the Michaelis-Menten model then comparing the observed dropout-rate (P), mean expression (S) or Michalis constant (K) calculated for each gene. Calculates p-value using a Z-test to compare observed values to those predicted by the fitted model after combining standard errors on each measure (P, S, and K) using error propagation rules.}
\value{List containing:
  \item{pval}{Vector of p-values of test for each gene.}
  \item{fold_change}{Fold change of respective variable (P, S or K).}
}
\seealso{
  \code{\link{M3D_Differential Expression}}.
}
\examples{
  de_test = bg__test_DE_K_equiv(Test_Cleaned$data)
  de_test = bg__test_DE_P_equiv(Test_Cleaned$data)
  de_test = bg__test_DE_S_equiv(Test_Cleaned$data)
}
\keyword{background}
