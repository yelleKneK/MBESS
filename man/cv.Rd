\name{cv}
\alias{cv}
\title{Function to calculate the regular (which is also biased) estimate of the coefficient of variation or the unbiased estimate of the coefficient of variation.}
\description{Returns the estimated coefficient of variation or the unbiased estimate of the coefficient of variation.}
\usage{
cv(C.of.V=NULL, mean=NULL, sd=NULL, N=NULL, unbiased=FALSE)
}
\arguments{
  \item{C.of.V}{Usual estimate of the coefficient of variation (\code{C.of.V=sd/mean})}
  \item{mean}{observed mean}
  \item{sd}{observed standard deviation (based on \code{N}-1 in the denominator of the variance)}
  \item{N}{sample size}

  \item{unbiased}{return the unbiased estimate of the coefficient of variation}
}
\details{
A function to calculate the usual estimate of the coefficient of variation or its unbiased estimate.
}
\value{
Returns the estimated coefficient of variation (regular but biased estimate or unbiased estimate.
}
%\references{}

\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})}
%\note{}
\seealso{\code{\link{ci.cv}}}
\examples{
cv(mean=100, sd=15)
cv(mean=100, sd=15, N=50, unbiased=TRUE)
cv(C.of.V=.15, N=2, unbiased=TRUE)}
\keyword{design}
\keyword{htest}
