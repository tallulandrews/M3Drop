\name{NBumiFeatureSelectionCombinedDrop}
\alias{NBumiFeatureSelectionCombinedDrop}
\title{Dropout-based Feature Selection}
\description{
	Ranks genes by significance of increase in dropouts compared to expectation.
}
\usage{
	NBumiFeatureSelectionCombinedDrop(fit, ntop=NULL, method="fdr", qval.thresh=2, suppress.plot=TRUE)
	
}
\arguments{
  \item{fit}{output from NBumiFitModel or NBumiFitBasicModel.}
  \item{ntop}{number of top ranked genes to return}
  \item{method}{correction method for multiple comparisons (check ?p.adjust.methods for more details)}
  \item{qval.thresh}{significant threshold}
  \item{suppress.plot}{logical, whether to plot the fitted curve and highlight selected features}
}
\details{
	Calculates dropout probability for each observation using depth-adjusted negative binomial means and dispersions calculated 
	from a fitted power-law relationship between mean and dispersion. Total dropouts per gene are modelled using the normal 
	approximation of the sum of bernoulli variables. And significance is evaluated using a Z-test.

	If provided, ntop will overrule the significance threshold.
}
\value{
	dataframe with columns: 
		Gene 
		effect_size (difference between observed and expected dropout rate)
		p.value
		q.value (corrected by adjustment method specifed by the method argument)
}
\examples{
	library(M3DExampleData)
	counts <- NBumiConvertData(Mmus_example_list$data)
	fit <- NBumiFitModel(counts);
	Drop_features <- names(NBumiFeatureSelectionCombinedDrop(fit, qval.thresh=0.05));
}
\keyword{single cell}
\keyword{feature selection}
\keyword{depth-adjusted negative binomial}
