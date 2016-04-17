\name{M3Drop_Get_Extremes}
\alias{M3Drop_Get_Extremes}
\alias{bg__get_extreme_residuals}
\title{Get outliers from MM curve.}
\description{
	Identifies outliers left and right of a fitted Michaelis-Menten curve.
}
\usage{
  M3Drop_Get_Extremes(expr_mat, fdr_threshold=0.1, percent=NA, v_threshold=c(0.05,0.95), suppress.plot=FALSE)
  bg__get_extreme_residuals(expr_mat, fit=NA, v_threshold=c(0.05, 0.95), percent = NA, fdr_threshold = 0.1, direction="right", suppress.plot=FALSE)
}
\arguments{
  \item{expr_mat}{a numeric matrix of expression values, columns = samples, rows = genes.}
  \item{fdr_threshold}{the threshold for identifying significant outliers after multiple testing correction.}
  \item{percent}{identify this percentage of data that is most extreme in each direction.}
  \item{v_threshold}{restrict to this range of dropout rates to avoid poorly fit regions of the data.}
  \item{suppress.plot}{logical, whether to plot the fitted Michaelis-Menten curve and highlight to identified most extreme outliers.}
  \item{fit}{the output from fitting the Michaelis-Menten curve.}
  \item{direction}{identify outliers "left" or "right" of the curve.}
}
\details{Fits a Michaelis-Menten function to the dropout-rate of the provided data, then identifies the most extreme left and/or right outliers from the curve. Horizontal residuals are calculated as : \deqn{\log_{10} S - \log_{10} \frac{K*(1-P)}{P}}{log_10(S) - log_10( (K * (1 - P)) / P )}. Extreme left[right] outliers are identified either as the \code{percent} smallest[largest] horizontal residuals. If \code{percent} is undefined (default) a normal distribution is fitted to the horizontal residuals and a Z-test is used to identify significant outliers after FDR multiple testing correction. 

Only genes with dropout rates within v_threshold will be considered to avoid the skewing of residuals due to the exponential parts of the MM curve near P = 0 & P = 1.

The normal distribution fit can be visualized by running \code{bg__get_extreme_residuals} with suppress.plot=FALSE. 

\code{M3Drop_Get_Extremes} identifies both left and right residuals using the provided thresholds in each direction. Eg. will return the \code{percent} smallest and \code{percent} largest residuals. It also plots the fitted MM curve and highlights the left and right extreme outliers unless suppress.plot=TRUE .

}
\value{
  \code{M3Drop_Get_Extremes} List containing elements left and right, vectors of the names of the extreme genes to the left and right of the curve respectively.
  \code{bg__get_extreme_residuals} vector of names of the genes of the extreme genes.
}
\examples{
  library(M3DExampleData)
  extreme_gene_lists = M3Drop_Get_Extremes(Mmus_example_list$data, fdr_threshold=0.1, suppress.plot=TRUE)
  extreme_gene_lists = M3Drop_Get_Extremes(Mmus_example_list$data, percent=0.01) # plots fitted MM and highlights extreme genes.

  R_extreme_genes = bg__get_extreme_residuals(Mmus_example_list$data, direction="right", fdr_threshold=0.05) #plots the fitted normal distribution
}
\seealso{
  \code{\link{bg__fit_MM}}
  \code{\link{bg__horizontal_residuals_MM_log10}}
}
\keyword{extremes, outliers, residuals}
