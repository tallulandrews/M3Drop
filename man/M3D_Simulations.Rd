\name{bg__MakeSimData}
\alias{bg__MakeSimData}
\alias{bg__MakeSimDE}
\alias{bg__MakeSimDVar}
\alias{bg__MakeSimHVar}
\title{Make Simulated Data}
\description{
	Makes simulated data based on a negative binomial distribution inflated with zeros based on the Michaelis-Menten equation.
}
\usage{
	bg__MakeSimData(dispersion_fun=bg__default_mean2disp, n_cells=300, dispersion_factor=1, base_means=10^rnorm(25000, 1, 1), K=10.3)
	bg__MakeSimDE(dispersion_fun=bg__default_mean2disp, fold_change=10, frac_change=0.1, n_cells=300, sub_pop=0.5, dispersion_factor=1, base_means=10^rnorm(25000,1,1), K=10.3)
	bg__MakeSimDVar(dispersion_fun=bg__default_mean2disp, fold_change=10, frac_change=0.1, n_cells=300, sub_pop=0.5, dispersion_factor=1, base_means=10^rnorm(25000,1,1), K=10.3)
	bg__MakeSimHVar(dispersion_fun=bg__default_mean2disp, fold_change=10, frac_change=0.1, n_cells=300, dispersion_factor=1, base_means=10^rnorm(25000,1,1), K=10.3)
}
\arguments{
  \item{dispersion_fun}{a function which takes mean experssion and returns the dispersion parameter of the negative binomial distribution.}
  \item{n_cells}{total number of cells (columns) in the simulated dataset.}
  \item{sub_pop}{proportion of cells with changed expression.}
  \item{frac_change}{proportion of genes with changed expression.}
  \item{fold_change}{fold change in dispersion or mean expression.}
  \item{dispersion_factor}{a factor that multiplies the calculated mean-specific dispersion for all genes.}
  \item{base_means}{a vector of background mean expression values.}
  \item{K}{K of the Michaelis-Menten function}
}
\details{Generates simulated single-cell gene expression data using a zero-inflated negative binomial distribution. A user-supplied function relates the dispersion parameter (1/size of the R parameterization of the negative binomial distribution). Zeros are added based on a Michaelis-Menten function. 

Default values of base_means, K, and dispersion_fun were fit to the Buettner et al. 2015 data [1].


\code{bg__MakeSimData} generates simulated single-cell data for a single homogeneous population.

\code{bg__MakeSimDE} generates simulated single-cell data for two different populations where a proportion of genes have a \code{fold_change} difference in the mean for population "2".

\code{bg__MakeSimDVar} generates simulated single-cell data for two different populations where a proportion of genes have a \code{fold_change} difference in the dispersion for population "2".

\code{bg__MakeSimHVar} generates simulated single-cell data for a single homogeneous population where a proportion of genes have a \code{fold_change} increase in dispersion over the expectation given the mean expression of the gene.
}
\value{
	\code{bg__MakeSimData} : a gene expression matrix where rows are genes, columns are cells.
	\code{bg__MakeSimDE}, \code{bg__MakeSimDVar}, \code{bg__MakeSimHVar} : a list of three named items:
		data : the gene expression matrix where rows are genes, columns are cells
		cell_labels : a vector of 1 or 2 indicating which cells are the unchanged ("1") or changed ("2") population.
		TP : a vector of row IDs of those genes that change (true positives).
}
\examples{
#  means = c(1,2,5,10,20,50,100,200,500,1000,2000,5000)
#  population1 <- bg__MakeSimData(n_cells=10, base_means=means)
#  population2 <- bg__MakeSimData(n_cells=10, base_means=means*2, dispersion_factor=0.5)
#  sim_DE <- bg__MakeSimDE(n_cells=100, base_means=means)
#  sim_DVar <- bg__MakeSimDVar(n_cells=100, sub_pop=0.25, base_means=means)
#  sim_HVar <- bg__MakeSimHVar(base_means=means, fold_change=3)
}
\keyword{simulations, differential expression, differential variance}
\references{
[1] Buettner et al. (2015) Computational analysis of cell-to-cell heterogeneity in single-cell RNA-sequencing data reveals hidden subpopulations of cells. Nature Biotechnology 33 : 155-160.
}
