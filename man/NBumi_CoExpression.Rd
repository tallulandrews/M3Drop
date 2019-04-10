\name{NBumiCoexpression}
\alias{NBumiCoexpression}
\title{Variance-based Feature Selection}
\description{
	Ranks genes by residual dispersion from mean-dispersion power-law relationship.
}
\usage{
	NBumiCoexpression(counts, fit, gene_list=NULL, method=c("both", "on", "off"))
}
\arguments{
  \item{counts}{raw count matrix (e.g. from NBumiConvertData).}
  \item{fit}{output from NBumiFitModel or NBumiFitBasicModel.}
  \item{gene_list}{set of gene names to test coexpression of.}
  \item{method}{type of coexpression to test (see: Details).}
}
\details{
	Tests for co-expression using the normal approximation of a binomial test.
	Co-expression is defined according to the method argument as follows:
	\describe{
		\item{on}{two genes are both >0 in more cells than expected.}
		\item{off}{two genes are both 0 in more cells than expected.}
		\item{both}{two genes are either both >0 or both 0 in more cells than expected.}
	}
	In all cases the null expectation is calculated using the fit library-size adjusted negative binomial model.
	This remove spurious co-expression due to cells with differing numbers of detected genes.
}
\value{
	a matrix of Z-scores for each pair of genes in the provided gene list.
}
\examples{
	library(M3DExampleData)
	counts <- NBumiConvertData(Mmus_example_list$data)
	fit <- NBumiFitModel(counts);
	genes <- c("Sox2", "Eomes", "Zscan4d", "Obox1", "Obox3")
	co <- NBumiCoexpression(counts, fit, genes, method="both");
	on <- NBumiCoexpression(counts, fit, genes, method="on");
	off <- NBumiCoexpression(counts, fit, genes, method="off");
}
\keyword{single cell}
\keyword{feature selection}
\keyword{depth-adjusted negative binomial}
