\name{M3DropDropoutModels}
\alias{M3DropDropoutModels}
\title{Fit functions to the dropouts vs expression distribution.}
\description{
	Fits the modified Michaelis-Menten equation (MM), a logistic regession (logistic), or a double exponential (ZIFA) function to the relationship between mean expression and dropout-rate (proportion of zero values).
}
\usage{M3DropDropoutModels(expr_mat, xlim=NA, suppress.plot=FALSE)}
\arguments{
  \item{expr_mat}{a numeric matrix of normalized (not log-transformed) expression values, columns = samples, rows = genes.}
  \item{xlim}{limits for x-axis of plot.}
  \item{suppress.plot}{logical, whether to plot fit curves or not.}
}
\details{Plots the dropout-rate (P) vs average gene expression (S) for all genes. Fits three different models and adds the fitted curves to the plot. The three models are: 
	MMfit : the Michaelis-Menten function \deqn{P=1-\frac{S}{S+K}}{P = 1 - S/(K+S)} (see: [1]).
	LogiFit : a logistic regression between P and log base 10 of S (used by [2]).
	ExpoFit : a double exponential \deqn{P = e^{\lambda S^2}}{P = e^(-lambda*S^2)} (used by [3]).
}
\value{Invisibly, a list of output from each fit (MMfit, LogiFit, ExpoFit).}
\examples{
  library(M3DExampleData)
  norm <- M3DropConvertData(Mmus_example_list$data, is.counts=TRUE)
  M3DropDropoutModels(norm)
} 

\references{
  [1] Keener, J.; Sneyd, J. (2008). Mathematical Physiology: I: Cellular Physiology (2 ed.). Springer. ISBN 978-0-387-75846-6 
  [2] Kharchenko, PV; Silberstein, L; Scadden, DT. (2014) Bayesian approach to single-cell differential expression analysis. Nature Methods. 11:740-742 
  [3] Pierson, E; Yau, C. (2015) ZIFA: Dimensionality reduction for zero-inflated single-cell gene expression analysis. Genome Biology. 16:241 doi:10.1186/s13059-015-0805-z}
\keyword{Michaelis Menten, model fitting}
