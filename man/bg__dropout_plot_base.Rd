\name{bg__dropout_plot_base}
\alias{bg__dropout_plot_base}
\title{Plot Dropouts vs Expression}
\description{
  Background function. Plots the dropout rate vs log(expression) for a gene expression matrix coloured by density.
}
\usage{bg__dropout_plot_base(expr_mat, weights = 1, xlim = NA, suppress.plot = FALSE)}
\arguments{
  \item{expr_mat}{a numeric matrix of expression values, columns = samples, rows = genes.}
  \item{weights}{a numeric matrix of weights for each expression value.}
  \item{xlim}{limits for xaxis of plot.}
  \item{suppress.plot}{Whether to make the plot or just calculate the requisite values.}
}
\details{
  Plots dropout rate vs log average expression for each gene. Points coloured by density. If weights are provided then weighted dropout rate and weighted average expression is calculated instead. 
}
\value{Returns object of type M3Dbaseplot.}
\seealso{
  \code{\link{M3D_Dropout_Models}}.
  \code{\link{M3D_Differential_Expression}}.
  \code{\link{M3D_Get_Extremes}}.
}
\examples{
  bg__dropout_plot_base(Test_Cleaned$data)
}
\keyword{plot, background}
