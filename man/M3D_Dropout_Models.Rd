\name{M3D_Dropout_Models}
\alias{M3D_Dropout_Models}
\title{Get outliers from MM curve.}
\description{
	Identifies outliers left and right of a fitted Michaelis-Menten curve.
}
\usage{M3D_Dropout_Models(expr_mat, weights=1, xlim=NA)}
\arguments{
  \item{expr_mat}{a numeric matrix of expression values, columns = samples, rows = genes.}
  \item{weights}{a numeric matrix of same dimension as expr_mat of weights for each expression value.}
  \item{xlim}{ limits for x-axis of plot.}
}
\details{Plots the dropout-rate (P) vs average gene expression (S) for all genes. Fits three different models and adds the fitted curves to the plot. The three models are: the Michaelis-Menten function $$P=\frac{S}{S+K}$$[1], a logistic regression between P and log base 10 of S[2], and a double exponential $$ P = e^{\lambda S^2}$$[3]}
\value{List of fitted parameters, predicted values, and residuals for all three model fits (MM, Logistic, double exponential)}
\examples{
  M3D_Dropout_Models(Test_Cleaned$data)
}
\references{[1] Keener, J.; Sneyd, J. (2008). Mathematical Physiology: I: Cellular Physiology (2 ed.). Springer. ISBN 978-0-387-75846-6 [2] Kharchenko, PV; Silberstein, L; Scadden, DT. (2014) Bayesian approach to single-cell differential expression analysis. Nature Methods. 11:740-742 [3] Pierson, E; Yau, C. (2015) ZIFA: Dimensionality reduction for zero-inflated single-cell gene expression analysis. Genome Biology. 16:241 doi:10.1186/s13059-015-0805-z}
\keyword{Michaelis Menten, model fitting}
