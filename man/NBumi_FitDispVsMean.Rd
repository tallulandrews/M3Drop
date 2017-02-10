\name{NBumiFitDispVsMean}
\alias{NBumiFitDispVsMean}
\title{Fit function between mean and dispersion}
\description{
	Fits a power-law relationship between mean expression and dispersion.
}
\usage{
	NBumiFitDispVsMean(fit, suppress.plot=TRUE)
}
\arguments{
  \item{fit}{output from NBumiFitModel or NBumiFitBasicModel.}
  \item{suppress.plot}{Whether to plot the calculated fit.}
}
\details{
	Fits a power-law relationship between mean expression and fitted gene-specific dispersions, using linear regression on ln transformed values. 
	Lowly expressed genes, mean expression < 2^4, are excluded as long as at least 2000 genes remain for fitting.
}
\value{
	Coefficients of linear regression.
}
\examples{
	library(M3DExampleData)
	counts <- as.matrix(Mmus_example_list$data);
	counts <- counts[rowSums(counts) > 0,];
	fit <- NBumiFitModel(counts);
	coeffs <- NBumiFitDispVsMean(fit, suppress.plot=TRUE);
}
