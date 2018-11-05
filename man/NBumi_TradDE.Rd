\name{bg__nbumiGroupDE}
\alias{bg__nbumiGroupDE}
\title{Perform Traditional Differential Expression}
\description{
	Performs traditional (i.e. compare defined groups) differential expression using a negative binomial model with MM zero-inflation. Functions tagged with "bg__" are not meant for direct usage and are not available in the Bioconductor release.
}
\usage{
}
\arguments{
  \item{counts}{a numeric matrix of raw UMI counts, columns = samples, rows = genes.}
  \item{fit}{fit NB UMI model from NBumiFitModel}
  \item{groups}{a vector of biological group IDs for each cell(columns).}
}
\details{
THIS FUNCTION SHOULD NOT BE USED.

\code{unfinished__nbumiGroupDE} Uses a log-likelihood ratio test to perform model selection between a model of constant mean expression vs a model of different mean expression across the biological groups. Probabilities of observing the data given the model are calculated using a negative binomial distribution with means adjusted for the total molecules detected per cell and dispersion fit to observed variance and adjusted to the mean of each group based on a globally fit power-law relationship. Significance is evaluated using the chi-square distribution.

}
\value{
	A table of fold mean expression differences for each biological group relative to the global mean expression level with raw p-values and FDR corrected p-values for each gene.
}
\examples{
#library(M3DExampleData)
#counts <- as.matrix(Mmus_example_list$data);
#counts <- counts[rowSums(counts) > 0,];
#counts <- counts[1:1000,]
#fit <- NBumiFitModel(counts);
#DE_output <- bg__nbumiGroupDE(counts, fit, Mmus_example_list$labels)
}
\keyword{differential expression}
