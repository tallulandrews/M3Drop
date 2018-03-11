\name{M3DropGetExtremes}
\alias{M3DropGetExtremes}
\alias{bg__get_extreme_residuals}
\title{Get outliers from MM curve.}
\description{
	Identifies outliers left and right of a fitted Michaelis-Menten curve. Functions tagged with "bg__" are not meant for direct usage and are not available in the Bioconductor release.
}
\usage{
  M3DropGetExtremes(expr_mat, fdr_threshold=0.1, percent=NA, v_threshold=c(0.05,0.95), suppress.plot=FALSE)
  bg__get_extreme_residuals(expr_mat, fit=NA, fdr_threshold=0.1, percent=NA, v_threshold=c(0.05, 0.95), direction="right", suppress.plot=FALSE)
}
\arguments{
  \item{expr_mat}{a numeric matrix of normalized (not log-transformed) expression values, columns = samples, rows = genes.}
  \item{fit}{output from fitting the Michaelis-Menten equation, see \code{bg__fit_MM}}
  \item{fdr_threshold}{the threshold for identifying significant outliers after multiple testing correction.}
  \item{percent}{identify this percentage of data that is most extreme in each direction.}
  \item{v_threshold}{restrict to this range of dropout rates to avoid poorly fit regions of the data.}
  \item{direction}{"left" or "right", which tail of outliers to examine.}
  \item{suppress.plot}{logical, whether to plot the fitted Michaelis-Menten curve and highlight to identified most extreme outliers.}
}
\details{Fits a Michaelis-Menten function to the dropout-rate of the provided data, then identifies the most extreme left and/or right outliers from the curve. Horizontal residuals are calculated as : \deqn{\log_{10} S - \log_{10} \frac{K*(1-P)}{P}}{log_10(S) - log_10( (K * (1 - P)) / P )}. Extreme left[right] outliers are identified either as the \code{percent} smallest[largest] horizontal residuals. If \code{percent} is undefined (default) a normal distribution is fitted to the horizontal residuals and a Z-test is used to identify significant outliers after FDR multiple testing correction. 

Only genes with dropout rates within v_threshold will be considered to avoid the skewing of residuals due to the exponential parts of the MM curve near P = 0 & P = 1.

\code{M3DropGetExtremes} identifies both left and right residuals using the provided thresholds in each direction. Eg. will return the \code{percent} smallest and \code{percent} largest residuals. It also plots the fitted MM curve and highlights the left and right extreme outliers unless suppress.plot=TRUE .
\code{bg__get_extreme_residuals} identifies one-sided extreme outliers.
}
\value{
  \code{M3DropGetExtremes} List containing elements left and right, vectors of the names of the extreme genes to the left and right of the curve respectively.
  \code{bg__get_extreme_residuals} A vector of names of the extreme genes in the specified direction.
}
\examples{
  library(M3DExampleData)
  norm <- M3DropConvertData(Mmus_example_list$data, is.counts=TRUE)
  extreme_gene_lists <- M3DropGetExtremes(norm, fdr_threshold=0.1)
  extreme_gene_lists <- M3DropGetExtremes(norm, percent=0.01)
#  Lextremes <- bg__get_extreme_residuals(norm, fdr_threshold=0.1, direction="left")
}
\keyword{extremes, outliers, residuals}
