\name{bg__horizontal_residuals_MM_log10}
\alias{bg__horizontal_residuals_MM_log10}
\title{Calculate Horizontal Residuals}
\description{
	Calculates horizontal residuals from the Michaelis-Menten Function. Functions tagged with "bg__" are not meant for direct usage and are not available in the Bioconductor release.
}
\usage{
	bg__horizontal_residuals_MM_log10(K, p, s)
}
\arguments{
  \item{K}{fitted Michaelis-Menten constant.}
  \item{p}{Observed dropout rate.}
  \item{s}{Observed mean expression.}
}
\details{ Calculates the log-transformed horizontal residuals from the Michaelis-Menten function. Input values may be single values or vectors.

\deqn{r = log_{10}(S) - \frac{(1-P)K}{P}}{log10(S) - K*(1-P)/P}
	
}
\value{ Value of horizontal residual.}
\examples{
#  s <- c(10,100,100)
#  p <- c(0.9, 0.5, 0.1)
#  res <- bg__horizontal_residuals_MM_log10(10, p, s)
}
\keyword{background}
