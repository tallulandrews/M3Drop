\name{Imputation}
\title{Normalized Data using the DANB model}
\alias{NBumiImputeNorm}
\description{
	Normalizes data to a common library size, imputing zeros as needed.
}
\usage{
	NBumiImputeNorm(counts, fit, total_counts_per_cell=median(fit$vals$tis))
}
\arguments{
  \item{counts}{raw count matrix, rows=genes, cols=cells}
  \item{fit}{Output from \code{NBumiFitModel}}
  \item{total_counts_per_cell}{reference library size to normalize all cells to.}
}
\details{
	Converts raw counts into positions in the CDF for the depth-adjusted negative binomial model fit to each observation. Adjusts the DANB parameters (mean and size) for the new library size. Then calculate the normalized counts for the equivalent position in the CDF for the NB using the new parameters. 
}
\value{
	Normalized count matrix.
}

\examples{
library(M3DExampleData)
counts <- as.matrix(Mmus_example_list$data);
counts <- counts[rowSums(counts) > 0,];
fit <- NBumiFitModel(counts);
normed_counts <- NBumiImputeNorm(counts, fit, 1000000)
}
\keyword{single cell}
\keyword{normalization}
