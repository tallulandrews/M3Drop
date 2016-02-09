\name{bg__get_extreme_residuals}
\alias{bg__get_extreme_residuals}
\title{Plot Dropouts vs Expression}
\description{
  Background function. Plots the dropout rate vs log(expression) for a gene expression matrix coloured by density.
}
\usage{bg__get_extreme_residuals(expr_mat, weights = 1, v_threshold=c(0.05,0.95), perc_most_extreme = NA, fdr_threshold=0.1, direction="right", suppress.plot=FALSE)}
\arguments{
  \item{expr_mat}{a numeric matrix of expression values, columns = samples, rows = genes.}
  \item{weights}{a numeric matrix of weights for each expression value.}
  \item{fit}{Output from fitting Michalis-Menten.}
  \item{v_threshold}{Range of dropout-rates to consider.}
  \item{perc_most_extreme}{Quantile most extreme genes to identify.}
  \item{fdr_threshold}{Use FDR to identify most extreme genes.}
  \item{direction}{Identify extreme genes to the "right" or "left" of the curve.}
  \item{suppress.plot}{Whether to make the plot or just calculate the requisite values.}
}
\details{Identifies extreme outliers from a fitted Michaelis-Menten curve. If perc_most_extreme is defined then returns that quantile of most extreme genes. Otherwise fits a normal distribution to horizontal residuals and uses fdr_threshold to identify extreme genes. Plots the fit of the normal distribution and the cut-off point unless suppress.plot = TRUE.}
\value{Vector of names of extreme genes.}
\seealso{
  \code{\link{M3D_Get_Extremes}}.
}
\examples{
  bg__get_extreme_residuals(Test_Cleaned$data)
}
\keyword{background}
