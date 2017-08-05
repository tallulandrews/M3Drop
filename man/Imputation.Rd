\name{NBumiImputeZeros}
\alias{NBumiImputeZeros}
\alias{M3DropImputZeros}
\title{Impute Zero values using the DANB model or M3Drop model}
\description{
	Imputes zero counts using the expectation from the depth-adjusted negative binomial model or Michaelis-Menten function. 
}
\usage{
	NBumiImputeZeros(counts, fit)
	M3DropImputeZeros(expr_mat)
}
\arguments{
  \item{counts}{a numeric matrix of raw UMI counts, columns = samples, rows = genes.}
  \item{fit}{fit NB UMI model from NBumiFitModel}
  \item{expr_mat}{a numeric matrix of normalized (not log-transformed) expression values, columns = samples, rows = genes.}
}
\details{
\code{NBumiImputeZeros} : Replaces all zero counts with the expectation from the depth-adjusted negative binomial model, accounting for differences in sequencing depth and gene-specific expression levels. 

\code{M3DropImputZeros} : Fits and calculates expected mean expression level of each gene based on the observed dropout rate using the Michaelis-Menten equation:
\deqn{S=\frac{K}{P}+K}{S = K/P + K} 
Then for all expression values which are less than the respective expected mean expression are replaced with the average of the observed expression and expected mean expression. Thus all zeros remain below all other observations at one half of the expected mean expression of the gene. Genes with no dropouts are assigned a dropout frequency of 1/(n_cells*2).
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

Normalized_data <- M3DropCleanData(Mmus_example_list$data,
                labels = Mmus_example_list$labels,
                is.counts=TRUE, min_detected_genes=2000)

imputed <- M3DropImputeZeros(Normalized_data$data);
}
\keyword{imputation}
