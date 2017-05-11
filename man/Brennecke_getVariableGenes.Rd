\name{BrenneckeGetVariableGenes}
\alias{BrenneckeGetVariableGenes}
\title{Identify Highly Variable Genes}
\description{
	Implements the method of Brennecke et al. (2013) to identify highly variable genes.
}
\usage{BrenneckeGetVariableGenes(expr_mat, spikes=NA, suppress.plot=FALSE, fdr=0.1, minBiolDisp=0.5, fitMeanQuantile=0.8)}
\arguments{
  \item{expr_mat}{a numeric matrix of normalized or raw (not log-transformed) expression values, columns = samples, rows = genes.}
  \item{spikes}{a vector of gene names of row numbers of spike-in genes which are subject to only technical variance.}
  \item{suppress.plot}{Whether to make the plot or just calculate the requisite values.}
  \item{fdr}{Use FDR to identify significantly highly variable genes.}
  \item{minBiolDisp}{Minimum percentage of variance due to biological factors.}
  \item{fitMeanQuantile}{Threshold for genes to be used in fitting. May need to be decreased in low-depth/umi-tagged datasets to achieve good fit.}
}
\details{Identifies significantly highly variable genes as detailed in Brennecked et al [1]. If spike-ins are provided they are used fit a function to the relationship between gene expression and variance due to technical factors. If spike-ins are not provided then all genes are used in the fitting. }
\value{Vector of names of highly variable genes.}
\references{Brennecke et al. (2013) Accounting for technical noise in single-cell RNA-seq experiments. Nature Methods 10, 1093-1095. doi:10.1038/nmeth.2645}
\examples{
  library(M3DExampleData)
  HVG <- BrenneckeGetVariableGenes(Mmus_example_list$data)
  HVG_spike <- BrenneckeGetVariableGenes(Mmus_example_list$data, spikes=5550:5600)
}
\keyword{highly variable genes}
