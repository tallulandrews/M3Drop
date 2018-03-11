\name{NBumiFeatureSelectionHighVar}
\alias{NBumiFeatureSelectionHighVar}
\title{Variance-based Feature Selection}
\description{
	Ranks genes by residual dispersion from mean-dispersion power-law relationship.
}
\usage{
	NBumiFeatureSelectionHighVar(fit)
	
}
\arguments{
  \item{fit}{output from NBumiFitModel or NBumiFitBasicModel.}
}
\details{
	Uses linear regression on log-transformed mean expression and fitted dispersions. Ranks genes by the residuals of this fit, negative = high variance, positive = low variance.
}
\value{
	Sorted vector of residuals.
}
\examples{
	library(M3DExampleData)
	counts <- NBumiConvertData(Mmus_example_list$data)
	fit <- NBumiFitModel(counts);
	HVGs <- names(NBumiFeatureSelectionHighVar(fit)[1:2000]);
}
\keyword{single cell}
\keyword{feature selection}
\keyword{depth-adjusted negative binomial}
