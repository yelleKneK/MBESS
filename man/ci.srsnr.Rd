\name{ci.srsnr}
\alias{ci.srsnr}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Confidence Interval for the Square Root of the Signal-To-Noise Ratio }
\description{
  Function to calculate the exact confidence interval for the square root of the signal-to-noise ratio.
}
\usage{
ci.srsnr(F.value = NULL, df.1 = NULL, df.2 = NULL, N = NULL, 
conf.level = 0.95, alpha.lower = NULL, alpha.upper = NULL, ...)
}

\arguments{
  \item{F.value}{ observed \emph{F}-value from the analysis of variance }
  \item{df.1}{ numerator degrees of freedom }
  \item{df.2}{ denominator degrees of freedom }
  \item{N}{ sample size }
  \item{conf.level}{ confidence interval coverage (i.e., 1 - Type I error rate); default is .95 }
  \item{alpha.lower}{ Type I error for the lower confidence limit }
  \item{alpha.upper}{ Type I error for the upper confidence limit }
  \item{\dots}{ allows one to potentially include parameter values for inner functions }
}
\details{
The confidence level must be specified in one of following two ways: using 
confidence interval coverage (\code{conf.level}), or lower and upper confidence 
limits (\code{alpha.lower} and \code{alpha.upper}).

The square root of the signal-to-noise ratio is defined as the standard deviation due to the particular factor over the 
standard deviation of the error (i.e., the square root of the mean square error). This function uses the confidence 
interval transformation principle (Steiger, 2004) to transform the confidence limits for the noncentality 
parameter to the confidence limits for square root of signal-to-noise ratio. The confidence interval 
for noncentral \emph{F} parameter can be abtained from function \code{conf.limits.ncf} in MBESS. 
}
\value{
 
Returns the square root of the confidence limits for the signal to noise ratio.

\item{Lower.Limit.of.the.Square.Root.of.the.Signal.to.Noise.Ratio}{lower limit of the square root of the signal to noise ratio}
\item{Upper.Limit.of.the.Square.Root.of.the.Signal.to.Noise.Ratio}{upper limit of the square root of the signal to noise ratio}
}
\references{
Fleishman, A. I. (1980). Confidence intervals for correlation ratios. \emph{Educational and Psychological Measurement, 40}, 659--670.

Kelley, K. (2007). Constructing confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20} (8), 1--24.


Steiger, J. H. (2004). Beyond the \emph{F} Test: Effect size confidence intervals and tests of close fit in the Analysis of Variance and Contrast Analysis. \emph{Psychological Methods, 9}, 164--182.  }

\author{ Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu}) }

% \note{ }
\seealso{ \code{ci.snr}, \code{conf.limits.ncf} }
\examples{
## To illustrate the calculation of the confidence interval for noncentral 
## F parameter,Bargman (1970) gave an example in which a 5-group ANOVA with 
## 11 subjects in each group is conducted and the observed F value is 11.2213. 
## This exmaple continued to be used in Venables (1975),  Fleishman (1980), 
## and Steiger (2004). If one wants to calculate the exact confidence interval 
## for square root of the signal-to-noise ratio of that example, this 
## function can be used.

ci.srsnr(F.value=11.221, df.1=4, df.2=50, N=55)

ci.srsnr(F.value=11.221, df.1=4, df.2=50, N=55, conf.level=.90)

ci.srsnr(F.value=11.221, df.1=4, df.2=50, N=55,  alpha.lower=.02, alpha.upper=.03)
}
\keyword{ design }% at least one, from doc/KEYWORDS
