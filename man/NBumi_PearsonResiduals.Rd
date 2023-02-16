\name{NBumiPearsonResiduals}
\alias{NBumiPearsonResiduals}
\alias{NBumiPearsonResidualsApprox}
\title{Calculate Pearson Residuals}
\description{
	Uses the NBumi depth-adjusted negative binomial model to calculate Pearson Residuals as an approach to normalizing 10X/UMI tagged data.
}
\usage{
	NBumiPearsonResiduals(counts, fits=NULL)
	NBumiPearsonResidualsApprox(counts, fits=NULL)
}
\arguments{
  \item{counts}{a numeric matrix of raw UMI counts, columns = samples/cells, rows = genes.}
  \item{fits}{the output from NBumiFitModel.}
}
\details{
Calculates a unique expectation for each gene expression observation based on the depth-adjusted negative binomial model (see: NBumiFitModel for details). This expection (mu) is equal to t_i*t_j/T, where t_i is the total counts for sample i, t_j is the total counts for gene j and T is the total counts.

\code{NBumiPearsonResidualApprox} Pearson residuals are approximated as: (counts - expectation) / sqrt(expectation). This assumes a Poisson-distribution of counts.

\code{NBumiPearsonResiduals} Pearson residuals are approximated as: (counts - mu) / sqrt(mu + mu^2/size). This uses a negative-binomial distribution of counts.

}
\value{
	a matrix of pearson residuals of equal dimension as your original matrix.
}
\examples{
	library(M3DExampleData)
	fit <- NBumiFitModel(counts);
	pearson1 <- NBumiPearsonResiduals(counts,fit)
	pearson2 <- NBumiPearsonResidualsApprox(counts,fit)
}
\keyword{single cell}
\seealso{NBumiFitModel}
