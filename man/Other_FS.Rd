\name{Other Feature Selection Methods}
\alias{irlbaPcaFS}
\alias{giniFS}
\alias{corFS}
\title{Other Feature Selection Methods}
\description{
	Other feature selection methods which only rank genes.
}
\usage{
	irlbaPcaFS(expr_mat, pcs=c(2,3))
	giniFS(expr_mat, suppress.plot=TRUE)
	corFS(expr_mat, dir=c("both", "pos", "neg"), fdr=NULL)
}
\arguments{
  \item{expr_mat}{normalized but not log-transformed gene expression matrix, rows=genes, cols=cells.}
  \item{pcs}{which principle components to use to score genes.}
  \item{suppress.plot}{whether to plot the gene expression vs Gini score.}
  \item{dir}{direction of correlation to consider.}
  \item{fdr}{apply an fdr threshold for significant features.}
}
\details{

	\code{irlbaPcaFS} Features are ranked by the sum of the magnitude of the loadings for the specified principle components. PCA is performed using the irlba package using sparse matricies for speed.
	\code{giniFS} Fits a loess curve between the maximum expression value and gini-index of each gene. Genes are ranked by p-value from a normal distribution fit to the residuals of the curve. As proposed by GiniClust [1].
	\code{corFS} Calculates all gene-gene correlations then ranks genes by the magnitude of the most positive or negative correlation. "both" will rank by the average of the magnitudes of the most positive & negative correlation.
}
\value{
	Sorted vector of scores for each gene from best features to worst features.
}
\examples{
library(M3DExampleData)
norm <- M3DropConvertData(Mmus_example_list$data[1:500,]);
Features_gini <- giniFS(norm, suppress.plot=FALSE);
Features_cor <- corFS(norm);
Features_pca <- irlbaPcaFS(norm);
}
\keyword{single cell}
\keyword{feature selection}
