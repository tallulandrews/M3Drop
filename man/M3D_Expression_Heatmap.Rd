\name{M3Drop_Expression_Heatmap}
\alias{M3Drop_Expression_Heatmap}
\alias{bg__expression_heatmap}
\title{Plot Heatmap of Gene Expression}
\description{
	Plots a customized heatmap of scaled log expression values.
}
\usage{
  M3Drop_Expression_Heatmap(genes, expr_mat, cell_labels=NA, interesting_genes=NA, key_genes=genes, key_cells=NA)
   bg__expression_heatmap(genes, expr_mat, cell_labels=NA, gene_labels=NA, key_genes=genes, key_cells=NA)
}
\arguments{
  \item{genes}{a character vector of gene names to be plot.}
  \item{expr_mat}{a numeric matrix of expression values, columns = samples, rows = genes.}
  \item{cell_labels}{factor of labels for each cell in the expression matrix that will to coloured in a top bar of the heatmap.}
  \item{interesting_genes}{list of vectors of gene names that will be used to colour the bar at the top of the heatmap.}
  \item{gene_labels}{factor of labels for each gene in the expression matrix t will be used to colour the bar at the side of the heatmap.}
  \item{key_genes}{a character vector of gene names to be labelled on the heatmap.}
  \item{key_cells}{a character vector of cells to be labelled on the heatmap.}
}
\details{Modifies the gplots function heatmap.2 to replace the row dendrogram with a legend of the colours used in the columns colour bar (cell_labels) and use a custom colour scalling. Expression is displayed as Z-scores of log transformed expression (adding a pseudocount of 1) coloured blue-white-red centered at 0 and binned in the range [-2,2].}
\value{Invisibly, output from heatmap.2 call.}
\examples{
  library(M3DExampleData)
  M3Drop_Expression_Heatmap(head(rownames(Mmus_example_list$data),20),Mmus_example_list$data, cell_labels = Mmus_example_list$labels)
}
\keyword{heatmap}
