\name{NBumiFitModel}
\alias{NBumiFitModel}
\alias{NBumiFitBasicModel}
\title{Fit Depth-Adjusted Negative Binomial Model}
\description{
	Fits means and dispersions for depth-adjusted or basic negative binomial models to a read/UMI count matrix.
}
\usage{
	NBumiFitModel(counts)
	NBumiFitBasicModel(counts)
	
}
\arguments{
  \item{counts}{a numeric matrix of raw UMI counts, columns = samples, rows = genes.}
}
\details{

\code{NBumiFitModel} Fits a depth-adjusted negative binomial model for each expression observation. Each expression is modelled using a negative binomial distribution with mean equal to t_i*t_j/T, where t_i is the total counts for sample i, t_j is the total counts for gene j and T is the total counts. Dispersions (R size parameter) are fit such that: var_j(counts_ij-mu_ij) = sum(mu_ij+mu_ij^2/size_j). Cases where genes exhibit poissonian behavior (size->infinity) are assigned a size of 10^10. 
\code{NBumiFitBasicModel} Fits a basic negative binomial model for each expression observation. Each expression is modelled using a negative binomial distribution with mean equal to t_j/n, where t_j is the total counts for gene j and n is the number of cells. Dispersions (R size parameter) are fit such that: var_j(counts_ij) = sum(mu_j+mu_j^2/size_j). Cases where genes exhibit poissonian behavior (size->infinity) are assigned a size of 10^10. 

}
\value{
	A named list of: 
		mus = a genes by samples matrix of means for the observation-specific negative binomial
		sizes = a vector of dispersions for each gene.
		vals = named list of summary statistics of the original count matrix: 
			tis = total molecules per cell, 
			tjs = total molecules per gene, 
			dis = total dropouts per cell, 
			djs = total dropouts per genes, 
			total = total molecules in dataset, 
			nc = number of cells, 
			ng = number of genes
}
\examples{
	library(M3DExampleData)
	counts <- NBumiConvertData(Mmus_example_list$data)
	fit <- NBumiFitModel(counts);
}
\keyword{single cell}
