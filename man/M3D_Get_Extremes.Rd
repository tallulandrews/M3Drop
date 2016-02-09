\name{M3D_Get_Extremes}
\alias{M3D_Get_Extremes}
\title{Get outliers from MM curve.}
\description{
	Identifies outliers left and right of a fitted Michaelis-Menten curve.
}
\usage{M3D_Get_Extremes(expr_mat, weights=1, fdr_threshold=0.1, percent=NA, v_threshold=c(0.05,0.95))}
\arguments{
  \item{expr_mat}{a numeric matrix of expression values, columns = samples, rows = genes.}
  \item{weights}{a numeric matrix of same dimension as expr_mat of weights for each expression value.}
  \item{fdr_threshold}{the threshold for identifying significant outliers after multiple testing correction.}
  \item{percent}{identify this percentage of data that is most extreme in each direction.}
  \item{v_threshold}{restrict to this range of dropout rates to avoid poorly fit regions of the data.}
}
\details{Fits a Michaelis-Menten function to the dropout-rate of the provided data, then identifies the most extreme outliers from the curve. Outliers are identified by fitting a normal distribution to the log residuals. A Z-test is used to identify significant outliers using an FDR multiple testing correction. If percent is provided this behaviour is overridden and the most extreme outliers are defined using quantiles.}
\value{List containing elements left and right, vectors of the extreme genes to the left and right of the curve respectively.}
\examples{
  M3D_Get_Extremes(Test_Cleaned$data, fdr_threshold=0.1)
  M3D_Get_Extremes(Test_Cleaned$data, percent=0.01)
}
\keyword{extremes}
