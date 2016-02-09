\name{bg__fit_logistic}
\alias{bg__fit_logistic}
\title{Fit Logistic Regression}
\description{Background function. Fits logistic regression to dropout rate vs log(expression).}
\usage{bg__fit_logistic(p, s)}
\arguments{
  \item{p}{Vector of dropout rates.}
  \item{s}{Vector of mean expression.}
}
\details{Fits logistic regression p~log(s) using base glm function using family="binomial".}
\value{List containing:
  \item{predictions}{predicted dropout rate.}
  \item{B0}{Fitted intercept.}
  \item{B1}{Fitted slope.}
  \item{model}{text description of the fit.}
  \item{SSr}{Sum of squared residuals.}
  \item{SAr}{Sum of absolute residuals.}
}
\seealso{
  \code{\link{M3D_Dropout_Models}}.
}
\examples{
  Base_Plot = bg__dropout_plot_base(Test_Cleaned$data)
  bg__fit_logistic(Base_Plot$P, Base_Plot$S)
}
\keyword{fit,  background}
