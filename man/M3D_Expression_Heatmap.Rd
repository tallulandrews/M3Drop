\name{M3DropExpressionHeatmap}
\alias{M3DropExpressionHeatmap}
\alias{bg__expression_heatmap}
\title{Plot Heatmap of Gene Expression}
\description{
	Plots a customized heatmap of scaled log expression values. Functions tagged with "bg__" are not meant for direct usage and are not available in the Bioconductor release.
}
\usage{
  M3DropExpressionHeatmap(genes, expr_mat, cell_labels=NA, interesting_genes=NA, key_genes=genes, key_cells=NA)
  bg__expression_heatmap(genes, expr_mat, cell_labels=NA, gene_labels=NA, key_genes=genes, key_cells=NA)
}
\arguments{
  \item{genes}{a character vector of gene names to be plot.}
  \item{expr_mat}{a numeric matrix of normalized (not log-transformed) expression values, columns = samples, rows = genes.}
  \item{cell_labels}{factor of labels for each cell in the expression matrix that will be used to colour the bar on the top of the heatmap.}
  \item{interesting_genes}{list of vectors of gene names that will be used to colour the bar to the left of the heatmap.}
  \item{key_genes}{a character vector of gene names to be labelled on the heatmap.}
  \item{key_cells}{a character vector of cells to be labelled on the heatmap. Unlabelled cells will be assigned a numerical index}
  \item{gene_labels}{factor of labels for each gene that will be used to colour the bar on the left of the heatmap.}
}
\details{Modifies the gplots function heatmap.2 to replace the row dendrogram with a legend of the colours used in the columns colour bar (cell_labels) and use a custom colour scalling. Expression is displayed as Z-scores of log transformed expression (adding a pseudocount of 1) coloured blue-white-red centered at 0 and binned in the range [-2,2].

\code{M3DropExpressionHeatmap} is a wrapper for \code{bg__expression_heatmap} that checks and reformats provided arguments. 
}
\value{Invisibly, output from heatmap.2 call.}
\examples{
  library(M3DExampleData)
  M3DropExpressionHeatmap(head(rownames(Mmus_example_list$data),20),Mmus_example_list$data, cell_labels = Mmus_example_list$labels)
#  bg__expression_heatmap(head(rownames(Mmus_example_list$data),20),Mmus_example_list$data, cell_labels = Mmus_example_list$labels)
}
\keyword{heatmap}
