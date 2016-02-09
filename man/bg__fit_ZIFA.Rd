\name{bg__fit_ZIFA}
\alias{bg__fit_ZIFA}
\title{Fit Double Exponential}
\description{Background function. Fits double exponential model to dropout rate vs expression data.}
\usage{bg__fit_ZIFA(p, s)}
\arguments{
  \item{p}{Vector of dropout rates.}
  \item{s}{Vector of mean expression.}
}
\details{Fits double exponential model[1]: p = e^{-\lambda s^2} using maximum likelihood estimation, using mle2 from bbmle package. }
\value{List containing:
  \item{lambda}{Fitted double exponential coefficient}
  \item{Lerr}{Standard error of lambda}
  \item{fitted_err}{Fitted standard deviation around the curve.}
  \item{predictions}{predicted dropout rate.}
  \item{model}{text description of the fit.}
  \item{SSr}{Sum of squared residuals.}
  \item{SAr}{Sum of absolute residuals.}
}
\seealso{
  \code{\link{M3D_Dropout_Models}}.
}
\examples{
  Base_Plot = bg__dropout_plot_base(Test_Cleaned$data)
  bg__fit_ZIFA(Base_Plot$P, Base_Plot$S)
}
\keyword{fit, ZIFAmodel, background}
\references{Pierson, E and Yau, C. (2015) ZIFA: Dimensionality reduction for zero-inflated single-cell gene expression analysis. Genome Biology 16:241. doi:10.1186/s13059-015-0805-z}
