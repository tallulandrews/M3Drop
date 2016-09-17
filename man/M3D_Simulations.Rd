\name{M3DropMakeSimData}
\alias{M3DropMakeSimData}
\alias{M3DropMakeSimDE}
\alias{M3DropMakeSimDVar}
\alias{M3DropMakeSimHVar}
\title{Make Simulated Data}
\description{
	Makes simulated data based on a negative binomial distribution inflated with zeros based on the Michaelis-Menten equation.
}
\usage{
	M3DropMakeSimData(dispersion_fun=bg_mean2disp, n_cells=300, dispersion_factor=1, base_means=10^rnorm(25000, 1, 1), K=10.3)
	M3DropMakeSimDE(dispersion_fun=bg_mean2disp, fold_change=10, frac_change=0.1, n_cells=300, sub_pop=0.5, dispersion_factor=1, base_means=10^rnorm(25000,1,1), K=10.3)
	M3DropMakeSimDVar(dispersion_fun=bg_mean2disp, fold_change=10, frac_change=0.1, n_cells=300, sub_pop=0.5, dispersion_factor=1, base_means=10^rnorm(25000,1,1), K=10.3)
	M3DropMakeSimHVar(dispersion_fun=bg_mean2disp, fold_change=10, frac_change=0.1, n_cells=300, dispersion_factor=1, base_means=10^rnorm(25000,1,1), K=10.3)
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


\code{M3DropMakeSimData} generates simulated single-cell data for a single homogeneous population.

\code{M3DropMakeSimDE} generates simulated single-cell data for two different populations where a proportion of genes have a \code{fold_change} difference in the mean for population "2".

\code{M3DropMakeSimDVar} generates simulated single-cell data for two different populations where a proportion of genes have a \code{fold_change} difference in the dispersion for population "2".

\code{M3DropMakeSimHVar} generates simulated single-cell data for a single homogeneous population where a proportion of genes have a \code{fold_change} increase in dispersion over the expectation given the mean expression of the gene.
}
\value{
	\code{M3DropMakeSimData} : a gene expression matrix where rows are genes, columns are cells.
	\code{M3DropMakeSimDE}, \code{M3DropMakeSimDVar}, \code{M3DropMakeSimHVar} : a list of three named items:
		data : the gene expression matrix where rows are genes, columns are cells
		cell_labels : a vector of 1 or 2 indicating which cells are the unchanged ("1") or changed ("2") population.
		TP : a vector of row IDs of those genes that change (true positives).
}
\examples{
  means = c(1,2,5,10,20,50,100,200,500,1000,2000,5000)
  population1 <- M3DropMakeSimData(n_cells=10, base_means=means)
  population2 <- M3DropMakeSimData(n_cells=10, base_means=means*2, dispersion_factor=0.5)
  sim_DE <- M3DropMakeSimDE(n_cells=100, base_means=means)
  sim_DVar <- M3DropMakeSimDVar(n_cells=100, sub_pop=0.25, base_means=means)
  sim_HVar <- M3DropMakeSimHVar(base_means=means, fold_change=3)
}
\keyword{simulations, differential expression, differential variance}
\references{
[1] Buettner et al. (2015) Computational analysis of cell-to-cell heterogeneity in single-cell RNA-sequencing data reveals hidden subpopulations of cells. Nature Biotechnology 33 : 155-160.
}
