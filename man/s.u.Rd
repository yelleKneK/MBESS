\name{s.u}
\alias{s.u}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Unbiased estimate of the population standard deviation}
\description{
Transforms the usual (and biased) estimate of the standard deviation into an unbiased estimator.}
\usage{
s.u(s=NULL, N=NULL, X=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{the usual estimate of the standard deviation (i.e., the square root of the unibased estimate of the variance)}
  \item{N}{sample size \code{s} is based}
  \item{X}{vector of scores in which the unbiased estimate of the standard deviation should be calculated}
}
\details{
Returns the unbiased estimate for the standard deviation.
}
\value{
The unbiased estimate for the standard deviation.
}
\references{
Holtzman, W. H. (1950). The unbiased estimate of the population variance and standard deviation. \emph{American 
Journal of Psychology}, \emph{63}, 615--617.
}
\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})}
%\note{}
%\seealso{}
\examples{
set.seed(113)
X <- rnorm(10, 100, 15)

# Square root of the unbiased estimate of the variance (not unbiased)
var(X)^.5

# One way to implement the function.
s.u(s=var(X)^.5, N=length(X))

# Another way to implement the function.
s.u(X=X)
}
\keyword{design}
\keyword{htest}
