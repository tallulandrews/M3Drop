\name{NBumiFeatureSelectionOther}
\alias{obsolete__nbumiFeatureSelectionDropouts}
\alias{obsolete__nbumiFeatureSelectionHighVarDist2Med}
\title{Other Feature Selection Methods}
\description{
	Ineffective alternative feature selection methods based on the depth-adjusted negative binomial model. Functions tagged with "bg__" are not meant for direct usage and are not available in the Bioconductor release.
}
\usage{
	
}
\arguments{
  \item{fit}{output from NBumiFitModel or NBumiFitBasicModel.}
  \item{window_size}{window for calculating the moving median.}
}
\details{

	Calculates dropout probability for each observation using depth-adjusted negative binomial means and dispersions calculated 
	from a fitted power-law relationship between mean and dispersion. Total dropouts per gene are modelled using the normal 
	approximation of the sum of bernoulli variables. And significance is evaluated using a Z-test.

	\code{obsolete__nbumiFeatureSelectionDropouts} Ranks genes by significance of increase in dropouts compared to expectation allowing for gene-specific dispersions.
	\code{obsolete__nbumiFeatureSelectionHighVarDist2Med} Ranks genes by the distance to median of log-transformed estimated dispersions.
}
\value{
	Sorted vector of p-values/distances.
}
\examples{
library(M3DExampleData)
#counts <- as.matrix(Mmus_example_list$data);
#counts <- counts[rowSums(counts) > 0,];
#fit <- NBumiFitModel(counts);
#Dropout_features <- names(obsolete__nbumiFeatureSelectionDropouts(fit)[1:2000]);
#dist2med_features <- names(obsolete__nbumiFeatureSelectionHighVarDist2Med(fit)[1:2000]);
}
\keyword{single cell}
\keyword{feature selection}
\keyword{depth-adjusted negative binomial}
