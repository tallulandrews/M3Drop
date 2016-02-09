\name{bg__highlight_genes}
\alias{bg__highlight_genes}
\title{Highlight Genes on Dropout plot}
\description{Background function. Highlights a set of genes using a specified colour.}
\usage{bg__highlight_genes(base_plot, genes, colour = "purple", pch = 16)}
\arguments{
  \item{base_plot}{An object of class M3Dbaseplot.}
  \item{genes}{A vector of genes to highlight. May be numeric, logical or character.}
  \item{colour}{Colour to highlight the genes with.}
  \item{pch}{Point type to highlight the genes with.}
}
\details{If the genes vector is numeric or logical, highlights genes corresponding to the specified rows of the expression matrix used to build base_plot. If the genes vector contains characters or a type that can be coerced to characters, then will highlight the matching genes from the expression matrix used to build base_plot.}
\value{None}
\seealso{
  \code{\link{M3D_Differential_Expression}}.
  \code{\link{M3D_Get_Extremes}}.
}
\examples{
  Base_Plot = bg__dropout_plot_base(Test_Cleaned$data)
  bg__highlight_genes(Base_Plot, c("MRPL18","SOD2","PSMB1","FPR2","TCEB2","ATP6V0C"))
}
\keyword{plot, background}
