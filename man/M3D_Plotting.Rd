\name{M3DropPlottingFunctions}
\alias{bg__dropout_plot_base}
\alias{bg__add_model_to_plot}
\alias{bg__highlight_genes}
\title{Make M3Drop Plots}
\description{
	Background functions used for making M3Drop plots. Functions tagged with "bg__" are not meant for direct usage and are not available in the Bioconductor release.
}
\usage{
	bg__dropout_plot_base(expr_mat, xlim=NA, suppress.plot=FALSE)
	bg__add_model_to_plot(fitted_model, base_plot, lty=1, lwd=1, col="dodgerblue", legend_loc="topright")
	bg__highlight_genes(base_plot, expr_mat, genes, col="darkorange", pch=16)
}
\arguments{
  \item{expr_mat}{a numeric matrix of normalized (not log-transformed) expression values, columns = samples, rows = genes.}
  \item{xlim}{limits of x-axis of plot, for comparing distributions across many datasets.}
  \item{suppress.plot}{Whether to make the plot or just return values.}
  \item{fitted_model}{output from fitting a model to dropout rate vs mean expression, see \code{bg__fit_MM},\code{bg__fit_ZIFA},\code{bg__fit_logistic}}
  \item{base_plot}{output from \code{bg__dropout_plot_base}}
  \item{legend_loc}{coordinates of top-right corner of the legend.}
  \item{genes}{list of genes to be highlighted.}
  \item{col,lty,lwd,pch}{Graphical parameters passed to plotted points & lines respectively.}
}
\details{ 
\code{bg__dropout_plot_base} plots grey-scale distribution of dropout-rate vs log10(mean expression) for each each. Colour indicates density of points light-grey = low, black = high.

\code{bg__add_model_to_plot} adds a line based on the predicted dropout rate for each gene from a particular model.

\code{bg__highlight_genes} highlights specified genes in the given colour.

}
\value{
\code{bg__dropout_plot_base} A named list of output:
	p = dropout rate for each gene
	s = mean expression of each gene
	xes = log10 transformed mean expression of each gene
	data = original expression matrix
	order = ordering of xes from smallest to largest
\code{bg__add_model_to_plot} Invisibly location of the plotted legend.
\code{bg__highlight_genes} Invisibly the coordinates of the highlighted genes.
}
\keyword{plotting, background}
