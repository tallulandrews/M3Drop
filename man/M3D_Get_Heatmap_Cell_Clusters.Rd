\name{M3Drop_Get_Heatmap_Cell_Clusters}
\alias{M3Drop_Get_Heatmap_Cell_Clusters}
\title{Extracts cell clusters from heatmap output}
\description{Extracts the clustering correponding to the given number of clusters from heatmap output.}
\usage{M3Drop_Get_Heatmap_Cell_Clusters(heatout, k)}
\arguments{
  \item{heatout}{Output from a gene-expression heatmap.}
  \item{k}{Number of clusters.}
}
\details{Traverses down the dendrogram and cuts at the first point where there are at least k clusters.}
\value{A vector of cluster labels for each cell.}
\examples{
  library(M3DExampleData)
  heatmap_out = M3Drop_Expression_Heatmap(rownames(Mmus_example_list$data)[1:20], Mmus_example_list$data)
  clusters = M3Drop_Get_Heatmap_Cell_Clusters(heatmap_out, k=5)
}
\keyword{heatmap}
