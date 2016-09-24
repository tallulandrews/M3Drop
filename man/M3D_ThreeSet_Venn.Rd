\name{M3DropThreeSetVenn}
\alias{M3DropThreeSetVenn}
\title{Three-way Venn Diagram}
\description{
	Plot an area-proportional three-set Venn Diagram with labels.
}
\usage{
  M3DropThreeSetVenn(set1, set2, set3, names)
}
\arguments{
  \item{set1}{a vector of items in the first set.}
  \item{set2}{a vector of items in the second set.}
  \item{set3}{a vector of items in the third set.}
  \item{names}{a vector of names of each set}
}
\details{Approximates area-proportional three-set Venn Diagram with code by David J. States (available at: http://tolstoy.newcastle.edu.au/R/help/03a/1115.html). Then places labels within each circle and overlap-region using code by Tallulah Andrews. 
}
\value{None}
\examples{
  SetA <- c(1:20)
  SetB <- c(15:30)
  SetC <- c(5,10,15,20,25,30,35,40,45,50,55,60)
  M3DropThreeSetVenn(SetA, SetB, SetC, names=c("A","B","C"))
}
\keyword{Venn Diagram}
