\name{bg__horizontal_residuals_MM_log10}
\alias{bg__horizontal_residuals_MM_log10}
\title{Calculate horizontal residuals}
\description{
	Calculates horizontal residuals from Michaelis-Menten curve.
}
\usage{
  bg__horizontal_residuals_MM_log10(K, p, s)
}
\arguments{
  \item{K}{fitted K parameter of the Michaelis-Menten equation}
  \item{p}{vector of dropout rates}
  \item{s}{vector of mean expression values}
}
\details{Horizontal residuals are calculated as : \deqn{\log_{10} S - \log_{10} \frac{K*(1-P)}{P}}{log_10(S) - log_10( (K * (1 - P)) / P )}. 
}
\value{
  A vector of horizontal residuals for each gene.
}
\examples{
  library(M3DExampleData)
  gene_info = bg__calc_variables(Mmus_example_list$data)
  residuals = bg__horizontal_residuals_MM_log10(K=9, p=gene_info$p, s=gene_info$s)
}
\seealso{
  \code{\link{M3D_Test_Shift}}
  \code{\link{M3D_Get_Extremes}}
}
\keyword{residuals}
