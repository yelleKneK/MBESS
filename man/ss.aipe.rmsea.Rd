\name{ss.aipe.rmsea}
\alias{ss.aipe.rmsea}
\title{ Sample size planning for RMSEA in SEM}
\description{
Sample size planning for the population root mean square error of approximation (RMSEA) from the accuracy in parameter estimation (AIPE) perspective. The sample size is planned so that the expected width of a confidence interval for the population RMSEA is no larger than desired.
}
\usage{
ss.aipe.rmsea(RMSEA, df, width, conf.level = 0.95)
}
\arguments{
  \item{RMSEA}{the input RMSEA value }
  \item{df}{degrees of freedom of the model}
  \item{width}{desired confidence interval width}
  \item{conf.level}{desired confidence level (e.g., .90, .95, .99, etc.)}
}
%\details{}
\value{
Returns the necessary total sample size in order to achieve the desired 
degree of accuracy (i.e., the sufficiently narrow confidence interval). 
}
%\references{  }
\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu}) and Keke Lai}
%\note{}
\seealso{ \code{\link{ci.rmsea}}}
\examples{
\dontrun{
# ss.aipe.rmsea(RMSEA=.035, df=50, width=.05, conf.level=.95)
}
}
\keyword{ design }
