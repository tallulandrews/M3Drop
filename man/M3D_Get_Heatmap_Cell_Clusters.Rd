\name{M3DropGetHeatmapClusters}
\alias{M3DropGetHeatmapClusters}
\alias{M3DropGetHeatmapNames}
\alias{M3DropGetHeatmapClusters}
\title{Extracts clusters/ordered names from heatmap output}
\description{Extracts the groupings correponding to the given number of clusters from heatmap output.}
\usage{M3DropGetHeatmapClusters(heatout, k, type="cell")
	M3DropGetHeatmapNames(heatout, type="cell")}
\arguments{
  \item{heatout}{Output from a gene-expression heatmap.}
  \item{k}{Number of clusters.}
  \item{type}{One of "cell" or "gene" indicating whether to get clusters or names of cells(columns) or genes(rows).}
}
\details{\code{M3DropGetHeatmapClusters} : Traverses down the row or column dendrogram and cuts at the first point where there are at least k clusters.
	 \code{M3DropGetHeatmapNames}: gets the names of the cells/genes in the order that they appear on the dendrogram.
}
\value{A vector of cluster labels for each cell.}
\examples{
  library(M3DExampleData)
  genes <- rownames(Mmus_example_list$data)[1:20]
  heatmap_out <- M3DropExpressionHeatmap(genes, Mmus_example_list$data)
  clusters <- M3DropGetHeatmapClusters(heatmap_out, k=5)
  heatmap_gene_labels <- M3DropGetHeatmapNames(heatmap_out, type="gene")
}
\keyword{heatmap, cluster}
