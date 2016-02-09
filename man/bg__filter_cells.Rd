\name{bg__filter_cells}
\alias{bg__filter_cells}
\title{Remove Low Quality Cells}
\description{Background function. Removes cells with few detected genes.}
\usage{bg__filter_cells(data, labels = NA, suppress.plot = FALSE, threshold = NA)}
\arguments{
  \item{expr_mat}{Numerical matrix of gene expression values. columns = cells, rows = genes.}
  \item{labels}{Known population IDs of cells.}
  \item{suppress.plot}{Don't create plot.}
  \item{threshold}{minimum number of detected genes to retain cell.}
}
\details{If threshold is provided then removes all cells with fewer than that number of detected genes. Otherwise, a normal distribution is fit to the distribution of number of detected genes per cell and those with significantly few detected genes (5\% FDR) are removed. Creates a plot of the fit normal distribution unless suppress.plot = TRUE.}
\value{List containing:
  \item{expr_mat}{Expression matrix after removing low quality cells.}
  \item{labels}{Known population IDs of cells after removing low quality cells.}
}
\seealso{
  \code{\link{M3D_Clean_Data}}.
}
\examples{
  data = bg__filter_cells(Test_Cleaned$data, Test_Cleaned$labels)
}
\keyword{normalization, background}
