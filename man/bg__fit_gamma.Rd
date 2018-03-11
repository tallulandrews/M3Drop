\name{bg__fit_gamma}
\alias{bg__fit_gamma}
\title{Fit Gamma Distribution}
\description{
	Fits the parameters of a gamma distribution to a set of observations
}
\usage{
	bg__fit_gamma(x)
}
\arguments{
  \item{x}{vector of observations to fit.}
}
\details{
	scale = var(x)/mean(x)
	shape = mean(x)/scale
}
\value{list with two entries: shape and size.}
\examples{
# dat <- rgamma(100, 0.1, 1)
# params <- bg__fit_gamma(dat)
}
\keyword{simulations}
