\name{Brennecke_getVariableGenes}
\alias{Brennecke_getVariableGenes}
\title{Identify Highly Variable Genes}
\description{
	Implements the method of Brennecke et al. (2013) to identify highly variable genes.
}
\usage{Brennecke_getVariableGenes(expr_mat, spikes=NA, suppress.plot=FALSE, fdr=0.1, minBiolDisp=0.5)}
\arguments{
  \item{expr_mat}{a numeric matrix of expression values, columns = samples, rows = genes.}
  \item{spikes}{a vector of gene names of row numbers of spike-in genes which are subject to only technical variance.}
  \item{suppress.plot}{Whether to make the plot or just calculate the requisite values.}
  \item{fdr}{Use FDR to identify significantly highly variable genes.}
  \item{minBiolDisp}{Minimum percentage of variance due to biological factors.}
}
\details{Identifies significantly highly variable genes as detailed in Brennecked et al [1]. If spike-ins are provided they are used fit a function to the relationship between gene expression and variance due to technical factors. If spike-ins are not provided then all genes are used in the fitting. }
\value{Vector of names of highly variable genes.}
\references{Brennecke et al. (2013) Accounting for technical noise in single-cell RNA-seq experiments. Nature Methods 10, 1093-1095. doi:10.1038/nmeth.2645}
\examples{
  library(M3DExampleData)
  Brennecke_getVariableGenes(Mmus_example_list$data)
  Brennecke_getVariableGenes(Mmus_example_list$data, spikes=5550:5600)
}
\keyword{highly variable genes}
