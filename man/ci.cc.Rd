\name{ci.cc}
\alias{ci.cc}
\title{
Confidence interval for the population correlation coefficient
}
\description{
This function is used to form a confidence interval for the population correlation coefficient. Note that this appraoch assumes that the variables the sample correlation coefficient are based are assumed to be bivariate normally distributed (e.g., Hays, 1994, Chapter 14). 
}
\usage{ci.cc(r, n, conf.level = 0.95, alpha.lower = NULL, alpha.upper = NULL)}
\arguments{
  \item{r}{observed value of the correlation coefficient (specifically the zero-order Pearson product-moment correlation coefficient)}
  \item{n}{sample size}
  \item{conf.level}{desired confidence level, where the error rate is the same on each side}
\item{alpha.lower}{the Type I error rate for the lower confidence interval limit}
\item{alpha.upper}{the Type I error rate for the upper confidence interval limit}
}
\details{
Note that this appraoch to confidence intervals does will not generally lead to a symmetric confidence interval. The function first transforms $r$ into emph{Z}\' , forms a confidence interval for the population value (i.e., $zeta$), and then transforms the confidence limits for $zeta$ into the scale of the correlation coefficient. 
}
\value{
\item{Lower.Limit }{lower limit of the confidence interval}
\item{Estimated.Correlation }{observed value of the correlation coefficient}
\item{Upper.Limit }{upper limit of the confidence interval}
}
\references{
Kelley, K. (2007). Confidence intervals for standardized effect sizes: Theory, application, 
and implementation.\emph{ Journal of Statistical Software, 20}(8), 1--24.

Hays, W. L. (1994). \emph{Statistics} (5th ed). Fort Worth, TX: Harcourt Brace College Publishers)
}
\author{ Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu}) }

\note{
This confidence interval assumes that the two variables the correlation is based are bivariate normal. See Hays (2004, Chapter 14) for details. 
}

\seealso{
\code{\link{transform_Z.r}}, \code{\link{transform_r.Z}}
}
\examples{
# Example, from Hayes. Suppose n=100 and r=.35. 
ci.cc(r=.35, n=100, conf.level=.95)

# Here is another way to enter the above example. 
ci.cc(r=.35, n=100, conf.level=NULL, alpha.lower=.025, alpha.upper=.025)

# Here are examples of one-sided confidence intervals. 
ci.cc(r=.35, n=100, conf.level=NULL, alpha.lower=0, alpha.upper=.05)
ci.cc(r=.35, n=100, conf.level=NULL, alpha.lower=.05, alpha.upper=0)
}
\keyword{ Design }
