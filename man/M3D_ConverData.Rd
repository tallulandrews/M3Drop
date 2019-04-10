\name{M3DropConvertData}
\alias{M3DropConvertData}
\title{Convert Data to be suitable for M3Drop}
\description{
	Recognizes a variety of R objects/classes and extracts expression 
	matrices from them then converts that to a normalized but 
	not-log transformed matrix appropriate for input into M3Drop functions.
}
\usage{
  M3DropConvertData(input, is.log=FALSE, is.counts=FALSE, pseudocount=1)
}
\arguments{
  \item{input}{ a matrix, data.frame or object }
  \item{is.log}{ has the data been log-transformed? (assumes log-base 2 with pseudocount of 1) }
  \item{is.counts}{ is the data raw unnormalized counts? (raw counts will be CPM normalized) }
  \item{pseudocount}{ pseudocount added before log-transformation }
}
\details{
	You must have loaded the respective packages (in parentheses) into your namespace before running this function on the respective objects. Note that to maintain scalability sparse matrices will remain as such.

	Supported classes/objects:
		\describe{

		\item{SCESet (scater <= 1.4.0)}{uses "exprs" or if unavailable then "counts"}

		\item{SingleCellExperiment (scater >= 1.6.0)}{uses "normcounts" if available, then "logcounts", which is assumed to be log-normalized, then "counts"}

		\item{CellDataSet (monocle)}{uses "exprs", specify log/counts using arguments}

		\item{ExpressionSet (Biobase)}{uses "exprs", specify log/counts using arguments}

		\item{seurat (Seurat)}{uses "raw.data" as counts.}
		}

		Matrix/Dataframe classes :

		\describe{
		\item{dgCMatrix (Matrix)}{specify log/counts using arguments}

		\item{data.table (data.table)}{specify log/counts using arguments}

		\item{DataTable (S4Vectors)}{specify log/counts using arguments}

		\item{DataFrame (S4Vectors)}{specify log/counts using arguments}

		\item{AnnotatedDataFrame (Biobase)}{specify log/counts using arguments}

		\item{matrix (base-r)}{specify log/counts using arguments}

		\item{data.frame (base-r)}{specify log/counts using arguments}
		}

}
\value{A normalized but not log-transformed matrix appropriate for input into M3Drop functions.}
\examples{
	# Simulated raw count matrix:
	set.seed(42)
	counts <- matrix(rpois(200, lambda=3), ncol=10)
	expr_mat <- M3DropConvertData(counts, is.counts=TRUE)

	# log normalized data frame
	lognorm <-log2( t(t(counts)/colSums(counts)*100)+1 )
	lognorm <- as.data.frame(lognorm)
	expr_mat <- M3DropConvertData(lognorm)

	# Sparse matrix
	require("Matrix")
	counts <- Matrix(counts, sparse=TRUE)
	expr_mat <- M3DropConvertData(counts, is.counts=TRUE)

	# SingleCellExperiment Object
	require("SingleCellExperiment")
	SCE <- SingleCellExperiment(assays=list(counts=counts))
	expr_mat <- M3DropConvertData(SCE)

	# monocle Object
	require("monocle")
	obj <- suppressWarnings(newCellDataSet(as.matrix(lognorm)))
	expr_mat <- M3DropConvertData(obj, is.log=TRUE)
}
\keyword{normalization}
