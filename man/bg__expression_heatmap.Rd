\name{bg__expression_heatmap}
\alias{bg__expression_heatmap}
\title{Gene Expression Heatmap}
\description{Background function. Makes a heatmap of gene expression for a given set of genes..}
\usage{bg__expression_heatmap(genes, data, cell_labels = NA, gene_labels = NA, key_genes = NA, key_cells = NA)}
\arguments{
  \item{genes}{A vector of genes to highlight. May be numeric, logical or character.}
  \item{data}{A numeric matrix of expression values. columns = samples/cells, rows = genes}
  \item{cell_labels}{A vector of known cell populations.}
  \item{gene_labels}{A vector of known gene groups.}
  \item{key_genes}{A vector of important genes that will be labelled on the plot}
  \item{key_cells}{A vector of important cells that will be labelled on the plot}
}
\details{ Customizes the gplots \link{heatmap.2} function for displaying single-cell RNASeq data. Includes a legend for the colour bar used to label cell populations. }
\value{list of \link{heatmap.2} output.}
\seealso{
  \code{\link{M3D_Expression_Heatmap}}.
  \code{\link{heatmap.2}}.
}
\examples{
  bg__expression_heatmap(c("MRPL18","SOD2","PSMB1","FPR2","TCEB2","ATP6V0C"), Test_Cleaned$data, cell_labels = Test_Cleaned$labels)
}
\keyword{plot, background, heatmap}
