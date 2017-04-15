\name{bg__NBumiGroupDE}
\alias{bg__NBumiGroupDE}
\title{Perform Traditional Differential Expression}
\description{
	Performs traditional (i.e. compare defined groups) differential expression using a negative binomial model with MM zero-inflation.
}
\usage{
	bg__NBumiGroupDE(counts, fit, groups)
}
\arguments{
  \item{counts}{a numeric matrix of raw UMI counts, columns = samples, rows = genes.}
  \item{fit}{fit NB UMI model from NBumiFitModel}
  \item{groups}{a vector of biological group IDs for each cell(columns).}
}
\details{
THIS FUNCTION SHOULD BE CONSIDERED UNSTABLE.

\code{bg__NBGroupDE} Uses a log-likelihood ratio test to perform model selection between a model of constant mean expression vs a model of different mean expression across the biological groups. Probabilities of observing the data given the model are calculated using a negative binomial distribution with means adjusted for the total molecules detected per cell and dispersion fit to observed variance and adjusted to the mean of each group based on a globally fit power-law relationship. Significance is evaluated using the chi-square distribution.

}
\value{
	A table of fold mean expression differences for each biological group relative to the global mean expression level with raw p-values and FDR corrected p-values for each gene.
}
\examples{
#library(M3DExampleData)
#counts <- as.matrix(Mmus_example_list$data);
#counts <- counts[rowSums(counts) > 0,];
#fit <- NBumiFitModel(counts);
#DE_output <- NBumiGroupDE(counts, fit, Mmus_example_list$labels)
#DE_output <- NBumiCGroupDE(counts, fit, Mmus_example_list$labels)
}
\keyword{differential expression}
