\name{M3DropCleanData}
\alias{M3DropCleanData}
\title{Filter Expression Data}
\description{
	Filters and normalizes a given expression matrix. Removes low quality cells and undetected genes, and normalizes counts to counts per million. 
}
\usage{
  M3DropCleanData(expr_mat, labels = NA, is.counts = TRUE, suppress.plot = FALSE, pseudo_genes = NA, min_detected_genes = NA)
}
\arguments{
  \item{expr_mat}{a numeric matrix of expression values, columns = samples/cells, rows = genes.}
  \item{labels}{a vector of length equal to the number of columns of expr_mat with names or group IDs for each cell.}
  \item{is.counts}{logical, whether the provided data is unnormalized read/fragment counts.}
  \item{suppress.plot}{logical, whether to plot the distribution of number of detected genes per cell.}
  \item{pseudo_genes}{a vector of gene names of known pseudogenes which will be removed from the cleaned data.}
  \item{min_detected_genes}{minimum number of genes/cell for a cell to be included in the cleaned data.}
}
\details{Retains genes detected (expression>0) in more than 3 cells and with mean normalized expression >= 10^-5. If min_detected_genes is defined all cells not reaching the threshold are removed. Otherwise, fits a normal distribution to the distribution of detected genes/cell and removes those cells with significantly few detected genes (FDR 5\%). This fit is plotted for visual inspection. If is.counts==TRUE then each column is converted to counts per million (ignoring ERCC spike-ins if present).}
\value{list with elements: data, the normalized filtered expression matrix; and labels, labels of the remaining cells.}
\examples{
  library(M3DExampleData)
  cpm = M3DropCleanData(Mmus_example_list$data, Mmus_example_list$labels, is.counts=TRUE, min_detected_genes=2000) # Remove all cells with < 2000 detected genes and convert to cpm 
  filtered_only = M3DropCleanData(Mmus_example_list$data, Mmus_example_list$labels, is.counts=FALSE) # Removes cells with significantly few detected genes (FDR=5%)
}
\keyword{normalization, quality control}
