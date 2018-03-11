\name{bg__shift_size}
\alias{bg__shift_size}
\title{Shift Size Parameter}
\description{
	Shifts a fitted size paramter according to the a power-law relationship.
}
\usage{
	bg__shift_size(mu_all, size_all, mu_group, coeffs) 
}
\arguments{
	\item{mu_all}{vector of original mus}
	\item{size_all}{vector of original sizes}
	\item{mu_group}{vector of new mus}
	\item{coeffs}{coefficients of the power-law fit}
}
\details{
	Wrapper of short function which converts values to their logs then shifts them using the parameters of a linear fit before un-logging them.
}
\value{vector of new size parameters}
\examples{
#require("M3DExampleData")
#counts <- NBumiConvertData(Mmus_example_list$data, is.counts=TRUE);
#fit <- NBumiFitModel(counts);
#ceoffs <- NBumiFitDispVsMean(fit, suppress.plot=TRUE)
#new_size <- bg__shift_size(rowMeans(counts), fit$size, rowMeans(counts)*2, coeffs)
}
\keyword{simulations}
