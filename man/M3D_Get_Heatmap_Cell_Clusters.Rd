\name{M3D_Get_Heatmap_Cell_Clusters}
\alias{M3D_Get_Heatmap_Cell_Clusters}
\title{Extracts cell clusters from heatmap output}
\description{Extracts the clustering correponding to the given number of clusters from heatmap output.}
\usage{M3D_Get_Heatmap_Cell_Clusters(heatmap_output, k)}
\arguments{
  \item{heatmap_output}{Output from a gene-expression heatmap.}
  \item{k}{Number of clusters.}
}
\details{Uses \code{cutree} to extract clusters from the column dendrogram of a heatmap.}
\value{A vector of cluster labels for each cell.}
\examples{
  library(M3DExampleData)
  heatmap_out = M3D_Expression_Heatmap(rownames(Mmus_example_list$data)[1:20], Mmus_example_list$data)
  clusters = M3D_Get_Heatmap_Cell_Clusters(heatmap_out, k=5)
}
\keyword{heatmap}
