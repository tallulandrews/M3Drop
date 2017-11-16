\name{bg__var_vs_drop}
\alias{bg__var_vs_drop}
\title{Variance vs Dropout Rate}
\description{
	Fixes the mean expression level and simulates a single gene with different fold changes across two equally sized subpopulations. 
}
\usage{
	bg__var_vs_drop(pop_size, fixed_mean, K=10.3, dispersion_from_mean=bg__default_mean2disp, suppress.plot=TRUE)
}
\arguments{
  \item{pop_size}{size of each subpopulation, i.e. 1/2 the total number of cells.}
  \item{fixed_mean}{expected mean expression across all cells.}
  \item{K}{K of the Michaelis-Menten function.}
  \item{dispersion_from_mean}{function which calculates the dispersion parameter of the negative binomial distribution given the mean expression of the gene.}
  \item{suppress.plot}{Whether to plot variance vs fold change and dropout rate vs fold change plots.}
}
\details{
	Simulates a single gene that is differentially expressed across two equally sized cell populations (\code{pop_size}) with a constant expected mean expression across all cells. Simulations are performed for fold changes of every integer value from 1 (no DE) to 100. Simulations use a negative binomial distribution inflated with zeros according to a Michaelis-Menten function. 

	Default mean2disp relationship and K were fit to the Buettner et al. [1] mESC data.

	Calculates pearson correlations between observed sample variance and dropout rate and fold change of the differential expression. Sample variance is decomposed into between subpopulation and within subpopulation variance ANOVA. Optionally plots the relationship between observed sample variance and dropout rate and the fold change. 
}
\value{
Name list of values:
	var_r = pearson correlation between observed sample variance and fold change.
	drop_r = pearson correlation between observed dropout rate and fold change.
	vars = vector of observed sample variances
	drops = vector of observed dropout rates
	fc = vector of expected fold changes
	Vbtw = vector of between subpopulation variances
	Vwithin = vector of within subpopulation variances
}
\examples{
#	mu100_N200 <- bg__var_vs_drop(100,100, suppress.plot=FALSE)
#	mu100_N2000 <- bg__var_vs_drop(1000,100)
#	mu1000_N200 <- bg__var_vs_drop(1000,100)
#	c(mu100_N200$var_r,mu100_N2000$var_r,mu1000_N200$var_r)
}
\keyword{simulations, variance, dropouts, differential expression}
\references{
[1] Buettner et al. (2015) Computational analysis of cell-to-cell heterogeneity in single-cell RNA-sequencing data reveals hidden subpopulations of cells. Nature Biotechnology 33 : 155-160.
}
