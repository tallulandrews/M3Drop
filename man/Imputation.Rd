\name{NBumiImputeZeros}
\alias{NBumiImputeZeros}
\title{Impute Zero values using the DANB model}
\description{
	Imputes zero counts using the expectation from the depth-adjusted negative binomial model. 
}
\usage{
	NBumiImputeZeros(counts, fit)
}
\arguments{
  \item{counts}{a numeric matrix of raw UMI counts, columns = samples, rows = genes.}
  \item{fit}{fit NB UMI model from NBumiFitModel}
}
\details{
	Replaces all zero counts with the expectation from the depth-adjusted negative binomial model, accounting for differences in sequencing depth and gene-specific expression levels. 
}
\value{
	Gene expression matrix with imputed zeros.
}
\examples{
library(M3DExampleData)
counts <- as.matrix(Mmus_example_list$data);
counts <- counts[rowSums(counts) > 0,];
counts <- counts[1:1000,]
fit <- NBumiFitModel(counts);
imputed <- NBumiImputeZeros(counts, fit)
}
\keyword{differential expression}
