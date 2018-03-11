\name{bg__default_mean2disp}
\alias{bg__default_mean2disp}
\title{Mean to Dispersion}
\description{
	Default function to calculate the negative binomial dispersion parameter given the mean expression of the gene.
}
\usage{
	bg__default_mean2disp(mu, coeffs=c(3.967816,-1.855054))
}
\arguments{
  \item{mu}{mean expression of the gene(s)}
  \item{coeffs}{coefficients of linear regression between log(cv2) and log10(mean expression) for observed data.}
}
\details{
Given the coefficients of a linear regression between the natural logarithm of the squared coefficient of variation and log base 10 of mean expression, calculates the dispersion parameter (1/size) of the negative binomial distribution for a given expression level.

Default coefficients were fit to the Buettner at al.[1] single cell mESC data.
}
\value{the dispersion parameter (1/size) for the (mu,size) paramterization of the \code{rnbinom} function.
}
\examples{
# mu <- 100
# r <- bg__default_mean2disp(mu)
# n_cells <- 50
# sim <- rnbinom(n_cells, size=1/r, mu=mu)
}
\keyword{simulations}
\references{
[1] Buettner et al. (2015) Computational analysis of cell-to-cell heterogeneity in single-cell RNA-sequencing data reveals hidden subpopulations of cells. Nature Biotechnology 33 : 155-160.
}
