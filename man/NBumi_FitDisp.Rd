\name{bg__fit_size_to_var}
\alias{bg__fit_size_to_var}
\title{Fit gene-specific dispersion}
\description{
	Fits dispersion for a specific gene for the depth-adjusted negative binomial. Functions tagged with "bg__" are not meant for direct usage and are not available in the Bioconductor release.
}
\usage{
	bg__fit_size_to_var(obs, mu_vec, max_size, min_size=10^-10, convergence=0.001)
}
\arguments{
  \item{obs}{observed variance corrected for library size.}
  \item{mu_vec}{expectation of depth-adjusted negative-binomial for each observation}
  \item{max_size}{maximum dispersion for genes that exhibit Poissonian expression.}
  \item{min_size}{initial dispersion estimate.}
  \item{convergence}{acceptable error on dispersion estimate.}
}
\details{
	Minimized the difference between the mean of expected variances for each observation based on the depth-adjust negative binomial(DANB) and 
	the observed variance. The observed variance is corrected for library size by subtracting the observation specific mean as per the DANB from 
	the observed UMI counts prior to calculating the variance. The optimization is done using an iterative method similar to Newton's method, which ends when the difference between estimates is less than "convergence" or if the estimated dispersion exceeds "max_size" or after 1000 iterations. 
}
\value{
	Estimated dispersion.
}
\examples{
#	library(M3DExampleData)
#	counts <- as.matrix(Mmus_example_list$data);
#	counts <- counts[rowSums(counts) > 0,];
#	fit <- NBumiFitModel(counts);
#	coeffs <- NBumiFitDispVsMean(fit, suppress.plot=TRUE);
}
