\name{ci.cv}
\alias{ci.cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Confidence interval for the coefficient of variation}
\description{
Function to calculate the confidence interval for the population coefficient of variation using the noncentral \code{t}-distribution.}
\usage{
ci.cv(cv=NULL, mean = NULL, sd = NULL, n = NULL, data = NULL, 
conf.level = 0.95, alpha.lower = NULL, alpha.upper = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cv}{coefficient of variation}
  \item{mean}{sample mean}
  \item{sd}{sample standard deviation (square root of the unbiased estimate of the variance)}
  \item{n}{sample size}
  \item{data}{vector of data for which the confidence interval for the coefficient of variation is to be calculated}
  \item{conf.level}{desired confidence level (1-Type I error rate)}
  \item{alpha.lower}{the proportion of values beyond the lower limit of the confidence interval (cannot be used with \code{conf.level}).}
  \item{alpha.upper}{the proportion of values beyond the upper limit of the confidence interval (cannot be used with \code{conf.level}).}
  \item{...}{ allows one to potentially include parameter values for inner functions}
}
\details{
Uses the noncentral \emph{t}-distribution to calculate the confidence interval for the population coefficient of variation.
}
\value{
\item{Lower.Limit.CofV}{Lower confidence interval limit}
\item{Prob.Less.Lower}{Proportion of the distribution beyond \code{Lower.Limit.CofV}}
\item{Upper.Limit.CofV}{Upper confidence interval limit}
\item{Prob.Greater.Upper}{Proportion of the distribution beyond \code{Upper.Limit.CofV}} 
\item{C.of.V}{Observed coefficient of variation}

}
\references{
Johnson, B. L., & Welch, B. L. (1940). Applications of the non-central \emph{t}-distribution. \emph{Biometrika}, 31, 362--389.

Kelley, K. (2007). Sample size planning for the coefficient of variation from the accuracy in parameter estimation approach. \emph{Behavior Research Methods, 39} (4), 755--766.
  
Kelley, K. (2007). Constructing confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20} (8), 1--24.

McKay, A. T. (1932). Distribution of the coefficient of variation and the extended \emph{t} distribution, \emph{Journal of the Royal Statistical Society}, \emph{95}, 695--698.
}

\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{\code{\link{cv}}}
\examples{
set.seed(113)
N <- 15
X <- rnorm(N, 5, 1)
mean.X <- mean(X)
sd.X <- var(X)^.5

ci.cv(mean=mean.X, sd=sd.X, n=N, alpha.lower=.025, alpha.upper=.025,
conf.level=NULL)
ci.cv(data=X, conf.level=.95)
ci.cv(cv=sd.X/mean.X, n=N, conf.level=.95)

}
\keyword{models}% at least one, from doc/KEYWORDS
\keyword{htest}% __ONLY ONE__ keyword per line
