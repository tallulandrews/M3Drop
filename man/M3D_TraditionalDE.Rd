\name{M3DropTraditionalDE}
\alias{bg__m3dropTraditionalDE}
\alias{bg__m3dropTraditionalDEShiftDisp}
\alias{bg__fitdispersion}
\alias{bg__get_mean2disp}
\title{Perform Traditional Differential Expression}
\description{
	Performs traditional (i.e. compare defined groups) differential expression using a negative binomial model with MM zero-inflation. Functions tagged with "bg__" are not meant for direct usage and are not available in the Bioconductor release.
}
\usage{
	bg__get_mean2disp(expr_mat)
	bg__fitdispersion(expr_mat)
}
\arguments{
  \item{expr_mat}{a numeric matrix of library-size normalized expression values, columns = samples, rows = genes.}
}
\details{
THESE FUNCTIONS SHOULD NOT BE USED.

\code{unfinished__m3dropTraditionalDE} : Uses a log-likelihood ratio test to perform model selection between a model of constant mean expression vs a model of different mean expression across the biological groups. 
Probabilities of observing the data given the model are calculated using a zero-inflated negative binomial distribution. Global relationships between mean and dispersion (power-law) as well as mean and 
dropouts (Michaelis-Menten) for genes are fit independently for each batch. Dispersions are fixed for each batch and calculated using the fitted power-law using the global mean expression of each gene.
Significance is evaluated using the chi-square distribution.

\code{unfinished__m3dropTraditionalDEShiftDisp}: Uses a log-likelihood ratio test to perform model selection between a model of constant mean expression vs a model of different mean expression across the biological groups. 
Probabilities of observing the data given the model are calculated using a zero-inflated negative binomial distribution. Global relationships between mean and dispersion (power-law) as well as mean and 
dropouts (Michaelis-Menten) for genes are fit independently for each batch. Dispersions are shifted from the global variance according to the mean expression for each biological group, using batch-specific power-laws.
Significance is evaluated using the chi-square distribution.


\code{broken__m3dropCTraditionalDE} : Uses a log-likelihood ratio test to perform model selection between a model of constant mean expression vs a model of different mean expression across the biological groups. 
Probabilities of observing the data given the model are calculated using a zero-inflated negative binomial distribution. Global relationships between mean and dispersion (power-law) as well as mean and 
dropouts (Michaelis-Menten) for genes are fit to the full dataset. Significance is evaluated using the chi-square distribution.

\code{bg__get_mean2disp} fits a power-law relationship between the squared coefficient of variation and mean expression of each gene, which is used to calculate the expected dispersion parameter for the negative binomial distribution from a given mean expression value.

\code{bg__fitdispersion} estimates gene-specific dispersions from the mean and variance of gene expression values \deqn{r = \frac{\mu^2}{\sigma^2-\mu}}{r = mu^2/(var-mu)}. Then fits a power-law relationship between the estimated dispersion and mean exprssion.
}

\value{
	\code{bg__m3dropTraditionalDE} : a table of observed mean expression levels for each biological group and each batch with raw p-values and FDR corrected p-values for each gene.
	\code{bg__m3dropTraditionalDEShiftDisp} : a table of observed mean expression levels for each biological group and each batch with raw p-values and FDR corrected p-values for each gene.
	\code{broken__m3dropCTraditionalDE} : a table of observed mean expression levels for each biological group with raw p-values and FDR corrected p-values for each gene.
	\code{bg__get_mean2disp} : a function which calculates the expected dispersion given the mean expression.
	\code{bg__fitdispersion} : exponent of the power-law relationship between dispersion and mean expression.
}
\examples{
library(M3DExampleData)
#Normalized_data <- M3DropCleanData(Mmus_example_list$data,
#                  labels = Mmus_example_list$labels,
#                  is.counts=TRUE, min_detected_genes=2000)
#DE_output <- bg__m3dropTraditionalDE(Normalized_data$data[1:100,], Normalized_data$labels)
#DE_shifted_output <- bg__m3dropTraditionalDEShiftDisp(Normalized_data$data[1:100,], Normalized_data$labels)
#DE_output_batches <- bg__m3dropTraditionalDE(Normalized_data$data[1:100,], Normalized_data$labels, batches=Normalized_data$labels) # each biological condition was performed separately
}
\keyword{differential expression}
