\name{NBumiCompareModels}
\alias{NBumiCompareModels}
\title{Compare negative binomial models}
\description{
	Compares the fit of the depth-adjusted negative binomial model and basic negative binomial model.
}
\usage{
	NBumiCompareModels(counts, size_factor=colSums(counts)/median(colSums(counts)))
}
\arguments{
  \item{counts}{a numeric matrix of raw UMI counts, columns = samples, rows = genes.}
  \item{size_factor}{a calculated size factor for library size normalization.}
}
\details{
	Compares the fit of the depth-adjusted negative binomial model and basic negative binomial model. Depth
	-adjusted negative binomial is fit to raw molecule counts. Basic negative binomial is fit to library-size 
	normalized counts. The absolute error between observed gene-specific dropouts and expectations given each 
	model is calculated. And a plot of fitted and observed mean-expression vs dropouts is created.
}
\value{
	Vector of errors for each model.
}
\examples{
	library(M3DExampleData)
	counts <- as.matrix(Mmus_example_list$data);
	counts <- counts[rowSums(counts) > 0,];
	coeffs <- NBumiCompareModels(counts);
}
