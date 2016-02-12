\name{bg__dropout_plot_base}
\alias{bg__dropout_plot_base}
\alias{bg__calc_variables}
\title{Plot Dropouts vs Expression}
\description{
  Background Functions. Plots the dropout rate vs log(expression) for a gene expression matrix coloured by density.
}
\usage{
  bg__dropout_plot_base(expr_mat, xlim = NA, suppress.plot = FALSE)
  bg__calc_variables(expr_mat)
}
\arguments{
  \item{expr_mat}{a numeric matrix of expression values, columns = samples, rows = genes.}
  \item{xlim}{limits for xaxis of plot.}
  \item{suppress.plot}{Whether to make the plot or just calculate the requisite values.}
}
\details{
  \code{bg__calc_variables} calculates the mean expression and dropout rate for each gene.

  \code{bg__dropout_plot_base} Plots dropout rate vs log average expression for each gene. Points coloured by density. 
}
\value{
  \code{bg__calc_variables} returns a list of the vector of dropout rates (p) and the vector of average gene expression (s)

  \code{bg__dropout_plot_base} returns a named list including:
  \item{p}{vector of dropout rate of each gene.}
  \item{s}{vector of mean expression values of each gene.}
  \item{xes}{vector of log_10 mean expression values.}
  \item{data}{the original gene expression matrix.}
  \item{order}{vector of the reordering of genes from lowest expressed to highest expressed.}
}
\examples{
  library(M3DExampleData)
  bg__calc_variables(Mmus_example_list$data)
  bg__dropout_plot_base(Mmus_example_list$data)
}
\keyword{background}
