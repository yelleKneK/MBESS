\name{ci.snr}
\alias{ci.snr}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Confidence Interval for the Signal-To-Noise Ratio }
\description{
  Function to obtain the exact confidence interval for the signal-to-noise ratio (i.e., the variance of the specific factor over the error variance).
}
\usage{
ci.snr(F.value = NULL, df.1 = NULL, df.2 = NULL, N = NULL, conf.level = 0.95,
 alpha.lower = NULL, alpha.upper = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{F.value}{ observed \emph{F}-value from the analysis of variance}
  \item{df.1}{ numerator degrees of freedom }
  \item{df.2}{ denominator degrees of freedom }
  \item{N}{ sample size }
  \item{conf.level}{ confidence interval coverage (i.e., 1 - Type I error rate), default is .95 }
  \item{alpha.lower}{ Type I error for the lower confidence limit }
  \item{alpha.upper}{ Type I error for the upper confidence limit }
  \item{\dots}{ allows one to potentially include parameter values for inner functions }
}
\details{
The confidence level must be specified in one of following two ways: using 
confidence interval coverage (\code{conf.level}), or lower and upper confidence 
limits (\code{alpha.lower} and \code{alpha.upper}).

This function uses the confidence interval transformation principle (Steiger, 2004) to transform the confidence limits for the noncentality parameter to the confidence limits for the population's signal-to-noise ratio. The confidence interval for noncentral \emph{F} parameter can be obtained 
from the \code{conf.limits.ncf} function in MBESS, which is used internally within this function.
  
}
\value{
  
Returns the confidence limits for the signal-to-noise ratio.
\item{Lower.Limit.Signal.to.Noise.Ratio}{lower limit for signal to noise ratio}
\item{Upper.Limit.Signal.to.Noise.Ratio}{upper limit for signal to noise ratio}

}

\references{ 
Kelley, K. (2007). Constructing confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20} (8), 1--24.

Fleishman, A. I. (1980). Confidence intervals for correlation ratios. \emph{Educational and Psychological Measurement, 40}, 659--670.

Steiger, J. H. (2004). Beyond the \emph{F} Test: Effect size confidence intervals and tests of close fit in the Analysis of Variance and Contrast Analysis.  \emph{Psychological Methods, 9}, 164--182. }

\author{ Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu}) }
\note{ 
The signal to noise ratio is defined as the variance due to the particular factor over the error variance (i.e., the mean square error).
}
\seealso{ \code{ci.srsnr}, \code{conf.limits.ncf} }
\examples{

## Bargman (1970) gave an example in which a 5-group ANOVA with 11 subjects in each 
## group is conducted and the observed F value is 11.2213. This example was 
## used in Venables (1975),  Fleishman (1980), and Steiger (2004). If one wants to calculate 
## the exact confidence interval for the signal-to-noise ratio of that example, this
## function can be used.

ci.snr(F.value=11.221, df.1=4, df.2=50, N=55)

ci.snr(F.value=11.221, df.1=4, df.2=50, N=55, conf.level=.90)

ci.snr(F.value=11.221, df.1=4, df.2=50, N=55,  alpha.lower=.02, alpha.upper=.03)
}

\keyword{ design }% at least one, from doc/KEYWORDS
