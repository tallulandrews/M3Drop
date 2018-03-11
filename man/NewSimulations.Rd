\name{Simulation Trifecta}
\alias{M3DropSimulationTrifecta}
\alias{NBumiSimulationTrifecta}
\title{Make Simulated Data from a provided scRNASeq dataset.}
\description{
	Fits either a zero-inflated negative binomial model (M3Drop) or the depth-adjusted negative binomial model (NBumi) to a provided dataset then simulates data from that model including differentially expressed (DE), differentially variable (DV), or globally unusually variable (HV) genes. 
}
\usage{
	M3DropSimulationTrifecta(original_data, n_genes=25000, n_cells=250, sub_pop_prop=0.5)
	NBumiSimulationTrifecta(original_data, n_genes=25000, n_cells=250, sub_pop_prop=0.5)
}
\arguments{
  \item{original_data}{the expression matrix of a scRNASeq dataset to base the simulations on. Should be normalized (not log-transformed) values for M3Drop or raw counts for NBumi.}
  \item{n_genes}{number of genes to simulated for each expression matrix.}
  \item{n_cells}{number of cells to simulate for each condition (total columns of final matrices = 2*n_cells).}
  \item{sub_pop_prop}{proportion of cells in one of the sub-populations.}
}
\details{Generates simulated single-cell gene expression data based on an existing dataset. Three expression matrices are produced each with the same cell and gene-specific parameters but where the log2-fold change is applied to either the mean (DE), or variance (DV) in one half of the cells or is applied to the variance across all the cells (HV). 

Mean expression for each simulated gene is drawn from a log-normal distribution fit to the original dataset. These means are then bottom thresholded to ensure all genes have a mean expression >= 10^-5, and top thresholded to ensure no gene has a mean expression greater than the largest mean expression in the original dataset.

\code{M3DropSimulationTrifecta}

\code{NBumiSimulationTrifecta} : Cell-specific library sizes are drawn from a gamma distribution fit to the original data.

}
\value{
	a named list of output including:
		truth - the true log (base 2) fold changes in expression level or variability for each gene.
		groups - a vector specifying the group ID for each cell for the DE and DV genes (1 = control, 2 = different).
		de - the count matrix containing genes differentially expressed across the two groups.
		dv - the count matrix containing genes with differential variability across the two groups.
		hv - the count matrix containing genes with globally unusual variability.
		
}
\examples{
	library(M3DExampleData)
	counts <- NBumiConvertData(Mmus_example_list$data)
	norm <- M3DropConvertData(Mmus_example_list$data, is.log=FALSE, is.counts=TRUE)
	ZINB_sim <- M3DropSimulationTrifecta(norm)
	DANB_sim <- NBumiSimulationTrifecta(counts)
}
\keyword{simulations, differential expression, differential variance}
