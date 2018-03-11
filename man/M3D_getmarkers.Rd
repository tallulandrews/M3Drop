\name{M3DropGetMarkers}
\alias{M3DropGetMarkers}
\title{Identify marker genes}
\description{
	Calculates area under the ROC curve for each gene to predict the best group of cells from all other cells.
}
\usage{
  M3DropGetMarkers(expr_mat, labels)
}
\arguments{
  \item{expr_mat}{a numeric matrix of normalized expression values, columns = samples, rows = genes.}
  \item{labels}{a vector of group ids for each cell/sample.}
}
\details{Uses the ROCR package to calculate the AUC for each gene for the group with the highest average rank. Significant is calculated using a Wilcox rank-sum test. 
}
\value{
  A dataframe with a row for each gene and columns: AUC, Group (which label this gene had the highest average rank for), and pval (uncorrected p-value of prediction).
}
\examples{
  library(M3DExampleData)
  norm <- M3DropConvertData(Mmus_example_list$data, is.counts=TRUE)
  marker_gene_table <- M3DropGetMarkers(norm, Mmus_example_list$labels)
}
\keyword{markers}
