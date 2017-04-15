\name{PoissonUMIFeatureSelectionDropouts}
\alias{PoissonUMIFeatureSelectionDropouts}
\title{Dropout-based Feature Selection}
\description{
	Ranks genes by significance of increase in dropouts compared to expectation.
}
\usage{
	PoissonUMIFeatureSelectionDropouts(fit)
	
}
\arguments{
  \item{fit}{output from NBumiFitModel or NBumiFitBasicModel.}
}
\details{
	Calculates dropout probability for each observation using depth-adjusted negative binomial means and dispersions equal to the mean (Poisson). Total dropouts per gene are modelled using the normal 
	approximation of the sum of bernoulli variables. And significance is evaluated using a Z-test.
}
\value{
	Sorted vector of p-values
}
\examples{
library(M3DExampleData)
counts <- as.matrix(Mmus_example_list$data);
counts <- counts[rowSums(counts) > 0,];
fit <- NBumiFitModel(counts);
Dropout_features <- names(PoissonUMIFeatureSelectionDropouts(fit)[1:2000]);
}
\keyword{single cell}
\keyword{feature selection}
\keyword{depth-adjusted negative binomial}
