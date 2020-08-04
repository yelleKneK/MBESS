\name{ss.aipe.cv}
\alias{ss.aipe.cv}
\title{Sample size planning for the coefficient of variation given the goal of Accuracy in Parameter Estimation approach to sample 
size planning}

\description{Determines the necessary sample size so that the expected confidence interval width
for the coefficient of variation will be sufficiently narrow, optionally with a desired degree of 
certainty that the interval will not be wider than desired. The value of \code{C.of.V} should be positive.}
\usage{
ss.aipe.cv(C.of.V = NULL, width = NULL, conf.level = 0.95, 
degree.of.certainty = NULL, assurance=NULL, certainty=NULL, 
mu = NULL, sigma = NULL, alpha.lower = NULL, alpha.upper = NULL, 
Suppress.Statement = TRUE, sup.int.warns = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{C.of.V}{population coefficient of variation on which the sample size procedure is based}
  \item{width}{desired (full) width of the confidence interval}
  \item{conf.level}{confidence interval coverage; 1-Type I error rate}
  \item{degree.of.certainty}{value with which confidence can be placed that describes the likelihood of obtaining a confidence interval less than the value specified (e.g., .80, .90, .95)}
  \item{assurance}{ an alias for \code{degree.of.certainty} }
  \item{certainty}{an alias for \code{degree.of.certainty}}
  \item{mu}{population mean (specified with \code{sigma} when \code{C.of.V} is not specified)}
  \item{sigma}{population standard deviation (specified with \code{mu} when \code{C.of.V}) is not specified)}
  \item{alpha.lower}{Type I error for the lower confidence limit}
  \item{alpha.upper}{Type I error for the upper confidence limit}
  \item{Suppress.Statement}{Suppress a message restating the input specifications}
  \item{sup.int.warns}{suppress internal function warnings (e.g., warnings associated with \code{qt})}
  \item{\dots}{for modifying parameters of functions this function calls }
}
%\details{}
\value{Returns the necessary sample size given the input specifications.}
%\references{}
\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})}
%\note{}
\seealso{\code{ss.aipe.cv.sensitivity}, \code{cv}}
\examples{
# Suppose one wishes to have a confidence interval with an expected width of .10 
# for a 99\% confidence interval when the population coefficient of variation is .25.
ss.aipe.cv(C.of.V=.1, width=.1, conf.level=.99)

# Ensuring that the confidence interval will be sufficiently narrow with a 99\% 
# certainty for the situation above.
ss.aipe.cv(C.of.V=.1, width=.1, conf.level=.99, degree.of.certainty=.99)
}
\keyword{design}
\keyword{htest}
