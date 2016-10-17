\name{bg__calc_variables}
\alias{bg__calc_variables}
\title{Calculate Gene Variables}
\description{
	Calculates a suite of gene-specific variables including: mean, dropout rate, and their standard errors.
}
\usage{
	bg__calc_variables(expr_mat)
}
\arguments{
  \item{expr_mat}{a numeric matrix of normalized (not log-transformed) expression values, columns = samples, rows = genes.}
}
\details{ Calculates mean expression and its standard error ( sd/sqrt(n) ), and the dropout rate, proportion of observations that are zero, and its standard error ( sqrt(p*(1-p)/n) ) for each gene in the dataset. Performs various checks to ensure expression matrix is of suitable format for M3Drop.
}
\value{
	Named list of calculated values:
	s = vector of mean expression for each gene
	s_stderr = vector of mean expression standard error for each gene
	p = vector of dropout rate for each gene
	p_stderr = vector of dropout rate standard error for each gene
}
\examples{
  library(M3DExampleData)
  Normalized_data <- M3DropCleanData(Mmus_example_list$data,
                  labels = Mmus_example_list$labels,
                  is.counts=TRUE, min_detected_genes=2000)
  vals <- bg__calc_variables(Normalized_data$data[1:10,])
}
\keyword{background}
