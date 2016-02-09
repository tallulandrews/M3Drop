\name{bg__invert_MM}
\alias{bg__invert_MM}
\title{Inverse Michaelis-Menten}
\description{Background function. Calculates expression level given dropout rate and K.}
\usage{bg__invert_MM(K, p)}
\arguments{
  \item{K}{Michaelis constant}
  \item{p}{Scalar or vector of dropout rates.}
}
\details{Inverts the Michaelis-Menten equation. Calculates expression level as K*(1-p)/p.}
\value{expression level for each value of p.}
\seealso{
  \code{\link{bg__horizontal_residuals_MM_log10}}.
  \code{\link{bg__test_DE_S_equiv}}.
}
\examples{
  Base_Plot = bg__dropout_plot_base(Test_Cleaned$data)
  bg__invert_MM(Base_Plot$K, Base_Plot$P)
}
\keyword{MMmodel, background}
