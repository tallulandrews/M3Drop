\name{NBumiCheckFit}
\alias{NBumiCheckFit}
\alias{NBumiCheckFitFS}
\title{Check Fit Quality}
\description{
	Checks the quality of the fit of the depth-adjusted negative binomial model.
}
\usage{
	NBumiCheckFit(counts, fit, suppress.plot=FALSE)
	NBumiCheckFitFS(counts, fit, suppress.plot=FALSE)
}
\arguments{
  \item{counts}{a numeric matrix of raw UMI counts, columns = samples, rows = genes.}
  \item{fit}{output from NBumiFitModel or NBumiFitBasicModel.}
  \item{suppress.plot}{Whether to plot the observed vs expected frequency of dropouts for genes & cells.}
}
\details{
	Calculates expected dropouts for genes and cells and compares to observed values. Optionally plots observed vs expected dropouts for both genes and cells.
	\code{NBumiCheckFit} uses depth-adjusted negative binomial with gene-specific dispersions.
	\code{NBumiCheckFitFS} uses depth-adjusted negative binomial with dispersions calculated from the power-law between gene-specific dispersion and mean expression (as is used for feature selection).
}
\value{
	Invisibly, named list of output: 
		gene_error = sum of squared error between observed and expected gene-specific total dropouts
		cell_error = sum of squared error between observed and expected cell-specific total dropouts
		exp_ps = gene by sample matrix of probability of a dropout for the negative binomial model fitted to each observation.
}
\examples{
	library(M3DExampleData)
	counts <- NBumiConvertData(Mmus_example_list$data)
	fit <- NBumiFitModel(counts);
	fit_quality <- NBumiCheckFitFS(counts, fit);
}
