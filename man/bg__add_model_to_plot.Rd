\name{bg__add_model_to_plot}
\alias{bg__add_model_to_plot}
\title{Add MM Curve to Plot}
\description{
  Background Function. Add a line for the fit MM curve to an existing plot of dropouts vs gene expression.
}
\usage{
  bg__add_model_to_plot(fitted_model, base_plot, lty=1, lwd=1, col="black", legend_loc="topright")
}
\arguments{
  \item{fitted_model}{output from fitting the MM function.}
  \item{base_plot}{output from \code{bg__dropout_plot_base}.}
  \item{lty}{line type.}
  \item{lwd}{line width.}
  \item{col}{colour of line.}
  \item{legend_loc}{location to put the legend.}
}
\details{Adds a line to the existing plot with an appropriate legend.}
\value{Invisibly, the coordinates of the plotted legend.}
\examples{
  library(M3DExampleData)
  BasePlot = bg__dropout_plot_base(Mmus_example_list$data)
  fit = bg__fit_MM(BasePlot$p, BasePlot$s)
  position = bg__add_model_to_plot(fit, BasePlot)
}
\keyword{background}
