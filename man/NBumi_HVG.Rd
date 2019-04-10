\name{NBumiHVG}
\alias{NBumiHVG}
\title{Variance-based Feature Selection Accounting for Library Size and Sample Variance}
\description{
	Tests for significantly high variability in droplet-based datasets.
}
\usage{
	NBumiHVG(counts, fit, fdr_thresh=0.05, suppress.plot=FALSE, method=c("DANB", "basic"))
}
\arguments{
  \item{counts}{raw count matrix (e.g. from NBumiConvertData).}
  \item{fit}{output from NBumiFitModel or NBumiFitBasicModel.}
  \item{fdr_thresh}{multiple testing correction threshold to apply to filter output.}
  \item{suppress.plot}{whether to plot mean vs variance & selected features.}
  \item{method}{whether to use DANB dispersions or raw sample variances.}
}
\details{
	Assumes a constant dispersion parameter due to technical noise (see: [1]), which is
	estimated using linear regression. Gene-specific observed sample variances are 
	compared to their respective expected variances. Significance is evaluated using
	a Z-test with the expected variance of the sample variance for a Negative Binomial
	(see: [2]). 

	The method argument controls whether the expected and observed variances are adjusted
	to account for uneven library sizes between cells. The default "DANB" option does
	the adjustment.
}
\value{
	a data.frame with columns: Gene, effect.size, p.value, q.value
}
\examples{
	library(M3DExampleData)
	counts <- NBumiConvertData(Mmus_example_list$data)
	fit <- NBumiFitModel(counts);
	HVGs <- NBumiHVG(counts, fit, suppress.plot=TRUE);
	HVGs_uncorrected <- NBumiHVG(counts, fit, suppress.plot=TRUE, method="basic");
}
\keyword{single cell}
\keyword{feature selection}
\keyword{depth-adjusted negative binomial}
\references{
	[1] Svensson, V. (2019) Droplet scRNA-seq is not zero-inflated. bioRxiv. doi: https://doi.org/10.1101/582064
	[2] Rose, C. and Smith, M. D. Mathematical Statistics with Mathematica. New York: Springer-Verlag, 2002. p264
}
