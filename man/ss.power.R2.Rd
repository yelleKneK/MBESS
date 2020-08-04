\name{ss.power.R2}
\alias{ss.power.R2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to plan sample size so that the test of the squared multiple correlation coefficient is sufficiently powerful.}
\description{
Function for determining the necessary sample size for the test of the squared multiple correlation 
coefficient or for determining the statistical power given a specified sample size for 
the squared multiple correlation coefficient in models where the regressors are regarded as fixed.
}
\usage{
ss.power.R2(Population.R2 = NULL, alpha.level = 0.05, desired.power = 0.85,
p, Specified.N = NULL, Cohen.f2 = NULL, Null.R2 = 0, 
Print.Progress = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Population.R2}{Population squared multiple correlation coefficient}
  \item{alpha.level}{Type I error rate}
  \item{desired.power}{desired degree of statistical power}
  \item{p}{the number of predictor variables}
  \item{Specified.N}{the sample size used to calculate power (rather than determine necessary sample size)}
  \item{Cohen.f2}{Cohen's (1988) effect size for multiple regression: \code{Population.R2}/(1-\code{Population.R2})}
  \item{Null.R2}{value of the null hypothesis that the squared multiple correlation will be evaluated against (this will typically be zero)}
  \item{Print.Progress}{if the progress of the iterative procedure is printed to the screen as the iterations are occuring}
  \item{\dots}{possible additional parameters for internal functions}
}
\details{
Determine the necessary sample size given a particular \code{Population.R2}, \code{alpha.level}, \code{p}, and \code{desired.power}. Alternatively, given \code{Population.R2}, \code{alpha.level}, \code{p}, and \code{Specified.N}, the function can be used to determine the statistical power.
}
\value{
  \item{Sample.Size}{returns either \code{Necessary.Sample.Size} or \code{Specified.Sample.Size}, depending on if sample size is being determined for a desired degree of statistical power analysis or if statistical power is being determined given a specified sample size, respectively}
  \item{Actual.Power}{Actual power of the situation described}
}
\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})}
\note{
When determining sample size for a desired degree of power, there will always be a 
slightly larger degree of actual power. This is the case because the algorithm employed 
determines sample size until the actual power is no less than the 
desired power (given sample size is a whole number power will almost certainly 
not be exactly the specified value). This is the same as other statistical power 
procedures that return whole numbers for necessary sample size.
}
\seealso{\code{ss.aipe.R2}, \code{ss.power.reg.coef}, \code{conf.limits.ncf}}
\examples{
# ss.power.R2(Population.R2=.5, alpha.level=.05, desired.power=.85, p=5)
# ss.power.R2(Cohen.f2=1, alpha.level=.05, desired.power=.85, p=5)
# ss.power.R2(Population.R2=.5, Specified.N=15, alpha.level=.05, 
# desired.power=.85, p=5)
# ss.power.R2(Cohen.f2=1, Specified.N=15, alpha.level=.05, desired.power=.85, p=5)
}
\keyword{design}
