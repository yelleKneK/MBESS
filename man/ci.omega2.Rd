\name{ci.omega2}
\alias{ci.omega2}

\title{
Confidence Interval for omega-squared (\eqn{\omega^2}) for fixed-effects ANOVA and ANCOVA designs (and partial omega-squared \eqn{\omega^2_p} for multifactor ANOVA and ANCOVA designs)}
\description{
Function to obtain the exact confidence interval using the non-central $F$ distribution for omega-squared or partial omega-squared.}
\usage{
ci.omega2(F.value = NULL, df.1 = NULL, df.2 = NULL, N = NULL, conf.level = 0.95, 
alpha.lower = NULL, alpha.upper = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{F.value}{
The value of the $F$-statistic for the analysis of (co)variace model (ANOVA) or, in the case of a multifactor ANOVA, the $F$-statistic for the particular factor.)
}
  \item{df.1}{
numerator degrees of freedom
}
  \item{df.2}{
denominator degrees of freedom
}
  \item{N}{
total sample size (i.e., the number of individual entities in the data)
}
  \item{conf.level}{ confidence interval coverage (i.e., 1-Type I error rate), default is .95 }
  \item{alpha.lower}{ Type I error for the lower confidence limit }
  \item{alpha.upper}{ Type I error for the upper confidence limit }
  \item{\dots}{ allows one to potentially include parameter values for inner functions }
}
\details{
The confidence level must be specified in one of following two ways: using 
confidence interval coverage (\code{conf.level}), or lower and upper confidence 
limits (\code{alpha.lower} and \code{alpha.upper}). The value returned is the confidence 
interval limits for the population \eqn{\omega^2} (or partial \eqn{\omega^2}).

This function uses the confidence interval transformation principle (Steiger, 2004) to transform the confidence limits for the noncentality parameter to the confidence limits for the population's (partial) omega-squared (\eqn{\omega^2}). The confidence interval for the noncentral \emph{F}-parameter can be obtained
from the \code{conf.limits.ncf} function in MBESS, which is used internally within this function.
}
\value{
Returns the confidence limits for (partial) omega-sqaured.
\item{lower_Limit_omega2}{lower limit for omega-squared}
\item{lower_Limit_omega2}{upper limit for omega-squared}
}
\references{
Fleishman, A. I. (1980). Confidence intervals for correlation ratios. \emph{Educational and Psychological Measurement, 40}, 659--670.

Kelley, K. (2007). Constructing confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20} (8), 1--24.

Steiger, J. H. (2004). Beyond the \emph{F} Test: Effect size confidence intervals and tests of close fit in the Analysis of Variance and Contrast Analysis. \emph{Psychological Methods, 9}, 164--182.
}
\author{ Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu}) }


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{ci.srsnr}, \code{ci.snr}, \code{conf.limits.ncf}
}
\examples{
## To illustrate the calculation of the confidence interval for noncentral 
## F parameter,Bargman (1970) gave an example in which a 5-group ANOVA with 
## 11 subjects in each group is conducted and the observed F value is 11.2213. 
## This exmaple continued to be used in Venables (1975),  Fleishman (1980), 
## and Steiger (2004). If one wants to calculate the exact confidence interval 
## for omega-squared of that example, this function can be used.

ci.omega2(F.value=11.221, df.1=4, df.2=50, N=55)

ci.omega2(F.value=11.221, df.1=4, df.2=50, N=55, conf.level=.90)

  }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{design}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
