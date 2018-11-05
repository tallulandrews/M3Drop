\name{NBumiCompareModels}
\alias{NBumiCompareModels}
\title{Compare negative binomial models}
\description{
	Compares the fit of the depth-adjusted negative binomial model and basic negative binomial model.
}
\usage{
	NBumiCompareModels(counts, size_factor=(Matrix::colSums(counts)/median(Matrix::colSums(counts))))
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
	A named list containing: 
		errors : Vector of errors for each model.
		basic_fit : object for the basic negative binomial.
		adjusted_fit : object for the depth-adjusted negative binomial.
}
\examples{
	library(M3DExampleData)
	counts <- NBumiConvertData(Mmus_example_list$data)
	out <- NBumiCompareModels(counts);
	out$errors
}
