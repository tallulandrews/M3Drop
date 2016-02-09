\name{M3D_Expression_Heatmap}
\alias{M3D_Expression_Heatmap}
\title{Plot Expression Heatmap}
\description{
	Plots a heatmap of the expression of a set of genes using gplots' heatmap.2
}
\usage{M3D_Expression_Heatmap(Genes, expr_mat, cell_labels=NA, interesting_genes=NA, marker_genes=NA, outlier_cells=NA)}
\arguments{
  \item{Genes}{a character vector of gene names to be plot.}
  \item{expr_mat}{a numeric matrix of expression values, columns = samples, rows = genes.}
  \item{cell_labels}{factor of labels for each cell in the expression matrix that will to coloured in a top bar of the heatmap.}
  \item{interesting_genes}{list of vectors of gene names that will be coloured on a side bar of the heatmap.}
  \item{marker_genes}{a character vector of gene names to be labelled on the heatmap.}
  \item{outlier_cells}{a character vector of cells to be labelled on the heatmap.}
}
\details{Modifies the gplots function heatmap.2 to replace the row dendrogram with a legend of the colours used in the columns colour bar (cell_labels) and use a custom colour scalling. Expression is displayed as Z-scores coloured blue-white-red centered at 0 and binned in the range [-2,2].}
\value{output from heatmap.2 call.}
\examples{
  M3D_Expression_Heatmap(head(rownames(Test_Cleaned$data)),Test_Cleaned$data, cell_labels = Test_Cleaned$labels)
}
\seealso{bg__expression_heatmap, heatmap.2}
\keyword{heatmap}
