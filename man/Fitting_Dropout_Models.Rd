\name{Fitting_Dropout_Models}
\alias{bg__fit_MM}
\alias{bg__fit_ZIFA}
\alias{bg__fit_logistic}
\title{Fit functions to the dropouts vs expression distribution.}
\description{
	Fits the modified Michaelis-Menten equation (MM), a logistic regession (logistic), or a double exponential (ZIFA) function to the relationship between mean expression and dropout-rate (proportion of zero values).
}
\usage{
  bg__fit_MM(p, s)
  bg__fit_logistic(p, s)
  bg__fit_ZIFA(p, s)
}
\arguments{
  \item{p}{a vector of dropout rates for each gene.}
  \item{s}{a vector of mean expression values for each gene. Must be the same order & length as p.}
}
\details{Fits one of different models to the relationship between dropout rate and mean expression. The three models are: 
		\code{bg__fit_MM} : the Michaelis-Menten function \deqn{P=1-\frac{S}{S+K}}{P = 1 - S/(K+S)} (see: [1]). Fit using \code{mle2} using normally distributed error.
		\code{bg__fit_logistic} : a logistic regression between P and log base 10 of S (used by [2]). Fit using \code{glm} (excludes genes where S == 0).
		\code{bg__fit_ZIFA} : a double exponential \deqn{P = e^{\lambda S^2}}{P = e^(-lambda*S^2)} (used by [3]). Fit using \code{lm} after log-transformation (genes were P == 0 are assigned a value of one tenth of the smallest P which is not 0).
}
\value{Named list including:
	K,fitted_err/B0,B1/lambda,fitted_err : the fitted parameters
	predictions : predicted values of p for each gene
	SSr/SAr :  sum of squared/absolute residuals
	model : vector of string descriptors of the fit
}
\examples{
#  library(M3DExampleData)
#  gene_info = bg__calc_variables(Mmus_example_list$data)
#  MM_fit = bg__fit_MM(gene_info$p, gene_info$s)
#  logistic_fit = bg__fit_logistic(gene_info$p, gene_info$s)
#  ZIFA_fit = bg__fit_ZIFA(gene_info$p, gene_info$s) 
}
\references{
  [1] Keener, J.; Sneyd, J. (2008). Mathematical Physiology: I: Cellular Physiology (2 ed.). Springer. ISBN 978-0-387-75846-6 
  [2] Kharchenko, PV; Silberstein, L; Scadden, DT. (2014) Bayesian approach to single-cell differential expression analysis. Nature Methods. 11:740-742 
  [3] Pierson, E; Yau, C. (2015) ZIFA: Dimensionality reduction for zero-inflated single-cell gene expression analysis. Genome Biology. 16:241 doi:10.1186/s13059-015-0805-z}
\keyword{Michaelis Menten, model fitting}
