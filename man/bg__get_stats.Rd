\name{bg__get_stats}
\alias{bg__get_stats}
\title{Calculate Simulation Statistics}
\description{
	Given significant genes and a ground truth calculates the false discovery rate and false negative rate of the significance calls.
}
\usage{	bg__get_stats(sig, TP, ngenes)}
\arguments{
  \item{sig}{vector of gene IDs/names which were called as significant by some method.}
  \item{TP}{vector of gene IDs/names which known to be positives (ground truth).}
  \item{ngenes}{number of total genes under consideration.}
}
\details{False discovery rate (FDR) is the proportion of all genes considered significant which are not in the ground truth. False negative rate (FNR) is the proportion of all ground truth genes which are not considered significant.
}
\value{vector of two elements: (FDR, FNR).}
\examples{
#	Calls <- c(1,3,5,8,10)
#	Truth <- c(2,3,5,9)
#	total=10
#	stats <- bg__get_stats(Calls, Truth, total)
}
\keyword{simulations, statistics}
