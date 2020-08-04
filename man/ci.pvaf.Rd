\name{ci.pvaf}
\alias{ci.pvaf}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Confidence Interval for the Proportion of Variance Accounted for (in the dependent variable by knowing the levels of the factor) }
\description{
Function to obtain the exact confidence limits for the proportion of variance of the dependent variable accounted for by knowing the levels of the factor (or the grouping factor in a single factor design) group status in a fixed factor analysis of variance.
}
\usage{
ci.pvaf(F.value = NULL, df.1 = NULL, df.2 = NULL, N = NULL, 
conf.level = 0.95, alpha.lower = NULL, alpha.upper = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{F.value}{ observed \emph{F}-value from fixed effects analysis of variance}
  \item{df.1}{ numerator degrees of freedom }
  \item{df.2}{ denominator degrees of freedom }
  \item{N}{ sample size }
  \item{conf.level}{ confidence interval coverage (i.e., 1-Type I error rate); default is .95 }
  \item{alpha.lower}{ Type I error for the lower confidence limit }
  \item{alpha.upper}{ Type I error for the upper confidence limit }
  \item{\dots}{ allows one to potentially include parameter values for inner functions }
}
\details{
The confidence level must be specified in one of following two ways: using confidence interval coverage (\code{conf.level}), or lower and upper confidence 
limits (\code{alpha.lower} and \code{alpha.upper}).

This function uses the confidence interval transformation principle (Steiger, 2004) to transform the confidence limits for the noncentrality parameter to the confidence limits for the population proportion of variance accounted for by knowing the group status. The confidence interval for the noncentral \emph{F} parameter can be obtained from the
 function \code{conf.limits.ncf} in MBESS, which is used within this function.
  
}
\value{
Returns the confidence interval for the proportion of variance of the dependent variable accounted for by knowing group status in a fixed factor analysis of variance (using a noncentral \emph{F}-distribution).
\item{Lower.Limit.Proportion.of.Variance.Accounted.for}{The lower confidence limit for the proportion of variance accounted for in the deviation by group status.}
\item{Upper.Limit.Proportion.of.Variance.Accounted.for}{The upper confidence limit for the proportion of variance accounted for in the deviation by group status.}
}
\references{ 
Fleishman, A. I. (1980). Confidence intervals for correlation ratios. \emph{Educational and Psychological Measurement, 40}, 659--670.

Kelley, K. (2007). Constructing confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20} (8), 1--24.

Steiger, J. H. (2004). Beyond the \emph{F} Test: Effect size confidence intervals and tests of close fit in the 
Analysis of Variance and Contrast Analysis. \emph{Psychological Methods, 9}, 164--182.}


\author{ Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu}) }
\note{ 
This function can be used for single or factorial ANOVA designs.
}
\seealso{\code{conf.limits.ncf} }
\examples{
\dontrun{
## Bargman (1970) gave an example in which a 5-group ANOVA with 11 subjects in each 
## group is conducted and the observed F value is 11.2213. This example was used 
## in Venables (1975),  Fleishman (1980), and Steiger (2004). If one wants to calculate the 
## exact confidence interval for the proportion of variance accounted for in that example, 
## this function can be used.

ci.pvaf(F.value=11.221, df.1=4, df.2=50, N=55)

ci.pvaf(F.value=11.221, df.1=4, df.2=50, N=55, conf.level=.90)

ci.pvaf(F.value=11.221, df.1=4, df.2=50, N=55,  alpha.lower=0, alpha.upper=.05)
}
}
\keyword{ design }% at least one, from doc/KEYWORDS
