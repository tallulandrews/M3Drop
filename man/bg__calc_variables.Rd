\name{bg__calc_variables}
\alias{bg__calc_variables}
\title{Calculate Dropout Rate and Mean Expression}
\description{Background function. Calculates dropout rate and mean expression for each gene in a given gene expression matrix.}
\usage{bg__calc_variables(expr_mat, weights = 1)}
\arguments{
  \item{expr_mat}{Numerical matrix of gene expression values. columns = samples/cells, rows = genes.}
  \item{weights}{Numerical matrix of the same dimension as expr_mat containing weights for each expression value.}
}
\details{Calculates weighted or unweighted dropout rate (frequency of zeros) and mean expression, as well as the standard error of the mean expression for each gene.}
\value{List containing:
  \item{s}{mean expression of each gene.}
  \item{p}{dropout rate (proportion of zeros) of each gene.}
  \item{s_stderr}{standard error of mean expression of each gene.}
  \item{s_stderr_nozero}{standard error of mean expression of each gene excluding zeros.}
  \item{p_stderr}{standard error of dropout rate of each gene.}
}
\seealso{
  \code{\link{bg__dropout_plot_base}}.
}
\examples{
  gene_info = bg__calc_variables(Test_Cleaned$data)
}
\keyword{plot, background}
