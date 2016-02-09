\name{bg__fit_MM}
\alias{bg__fit_MM}
\title{Fit Michaelis-Menten}
\description{Background function. Fits Michaelis constant to dropout rate vs expression data.}
\usage{bg__fit_MM(p, s)}
\arguments{
  \item{p}{Vector of dropout rates.}
  \item{s}{Vector of mean expression.}
}
\details{Fits Michaelis-Menten constant using maximum likelihood estimation, using mle2 from bbmle package.}
\value{List containing:
  \item{K}{Fitted Michaelis constant}
  \item{Kerr}{Standard Error of Michaelis constant}
  \item{fitted_err}{Fitted standard deviation around the curve.}
  \item{predictions}{predicted dropout rate.}
  \item{model}{text description of the fit.}
  \item{SSr}{Sum of squared residuals.}
  \item{SAr}{Sum of absolute residuals.}
}
\seealso{
  \code{\link{M3D_Dropout_Models}}.
  \code{\link{M3D_Differential Expression}}.
  \code{\link{M3D_Get_Extremes}}.
}
\examples{
  Base_Plot = bg__dropout_plot_base(Test_Cleaned$data)
  bg__fit_MM(Base_Plot$P, Base_Plot$S)
}
\keyword{fit, MMmodel, background}
