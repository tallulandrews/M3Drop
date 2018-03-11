\name{NBumiConvertData}
\alias{NBumiConvertData}
\title{Convert Data to be suitable for NBumi}
\description{
	Recognizes a variety of R objects/classes and extracts expression 
	matrices from them then converts that to a count matrix for input 
	into NBumi functions.
}
\usage{
  NBumiConvertData(input, is.log=FALSE, is.counts=FALSE, pseudocount=1)
}
\arguments{
  \item{input}{a matrix, data.frame or object }
  \item{is.log}{ has the data been log-transformed? (assumes log-base 2 with pseudocount of 1) }
  \item{is.counts}{ is the data raw unnormalized counts? (raw counts will be CPM normalized) }
  \item{pseudocount}{ pseudocount added before log-transformation }
}
\details{
	You must have loaded the respective packages (in parentheses) into your namespace before running this function on the respective objects. Note that to maintain scalability sparse matrices will remain as such.

	Supported classes/objects:
		SCESet (scater <= 1.4.0) : uses "counts" or if unavailable then "exprs" 
		SingleCellExperiment (scater >= 1.6.0) : uses "counts", if unavailable then "logcounts", which is assumed to be log-normalized.
		CellDataSet (monocle) : uses "exprs", specify log/counts using arguments
		ExpressionSet (Biobase) : uses "exprs", specify log/counts using arguments
		seurat (Seurat) : uses "raw.data" as counts.

		Matrix/Dataframe classes :
		dgCMatrix (Matrix) : specify log/counts using arguments
		data.table (data.table) : specify log/counts using arguments
		DataTable (S4Vectors) : specify log/counts using arguments
		DataFrame (S4Vectors): specify log/counts using arguments
		AnnotatedDataFrame (Biobase) : specify log/counts using arguments
		matrix (base-r) : specify log/counts using arguments
		data.frame (base-r) : specify log/counts using arguments

	Counts are rounded up to integers if necessary, if counts are 
	unavailable then this will attempt to convert log2 normalized expression
	to counts by de-logging, subtracting the pseudocount, and then 
	un-normalizing by rescaling cells based on their relative number of 
	detected genes, finally expression is rounded up to integers for use 
	as counts. 

}
\value{A count matrix appropriate for input into NBumi functions.}
\examples{
	# Simulated raw count matrix:
	set.seed(42)
	counts <- matrix(rpois(200, lambda=3), ncol=10)
	input_counts <- NBumiConvertData(counts, is.counts=TRUE)

	# log normalized data frame
	lognorm <-log2( t(t(counts)/colSums(counts)*100)+1 )
	lognorm <- as.data.frame(lognorm)
	input_counts <- NBumiConvertData(lognorm)

	# Sparse matrix
	require("Matrix")
	counts <- Matrix(counts, sparse=TRUE)
	input_counts <- NBumiConvertData(counts, is.counts=TRUE)

	# SingleCellExperiment Object
	require("SingleCellExperiment")
	SCE <- SingleCellExperiment(assays=list(counts=counts))
	input_counts <- NBumiConvertData(SCE)

	# monocle Object
	require("monocle")
	obj <- suppressWarnings(newCellDataSet(as.matrix(lognorm)))
	input_counts <- NBumiConvertData(obj, is.log=TRUE)
}
\keyword{normalization}
