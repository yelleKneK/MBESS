\name{ancova.random.data}
\alias{ancova.random.data}
\title{ Generate random data for an ANCOVA model }
\description{
Generate random data for a simple (one-response-one-covariate) ANCOVA model considering the covariate as random. Data can be generated in the
contexts of both randomized design (same population covariate mean across groups) and non-randomized 
design (different population covariate means across groups).
}
\usage{
ancova.random.data(mu.y, mu.x, sigma.y, sigma.x, rho, J, n, randomized = TRUE)
}
\arguments{
  \item{mu.y}{ a vector of the population group means of the response variable }
  \item{mu.x}{ the population mean of the covariate (in the randomized design context), or a vector of the population group
  means of the covariate (in the non-randomized design context) }
  \item{sigma.y}{ the population standard deviation of the response (outcome) variable }
  \item{sigma.x}{ the population standard deviation of the covariate}
  \item{rho}{ the population correlation coefficient between the response and the covariate }
  \item{J}{ the number of groups }
  \item{n}{ the number of sample size \emph{per group} }
  \item{randomized}{ a logical statement of whether randomized design is used}
}
\details{
This function uses a multivariate normal distribution to generate the random data; the covariate is considered 
as a random variable in the model. This function uses \code{mvrnorm} in the \code{MASS} package in an internal function, and
thus it requires the \code{MASS} package be installed.

This function assumes homogeneous covariance matrix among groups, in both the randomized design and
non-randomized design contexts. 
}
\value{
This function returns an \eqn{n} by \eqn{J2} matrix, where \eqn{n} and \eqn{J} are as defined 
in the argument. The first \eqn{J} columns of the matrix contains the random data for the response, and
the second \eqn{J} columns of the matrix contains the random data for the covariate.
}
%\references{  }
\author{ Keke Lai (University of California-Merced) and Ken Kelley (University of Notre Dame) <kkelley@nd.edu>}
%\note{}
\seealso{ \code{mvrnorm} in the \code{MASS} package }
\examples{
random.data <- ancova.random.data(mu.y=c(3,5), mu.x=10, sigma.y=1, 
sigma.x=2, rho=.8, J=2, n=20)
}

\keyword{ design }
