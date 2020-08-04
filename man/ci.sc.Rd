\name{ci.sc}
\alias{ci.sc}

\title{ Confidence Interval for a Standardized Contrast in a Fixed Effects ANOVA}
\description{
Function to obtain the confidence interval for a standardized contrast in a fixed effects analysis of variance context.
}
\usage{
ci.sc(means = NULL, s.anova = NULL, c.weights = NULL, n = NULL, 
N = NULL, Psi = NULL, ncp = NULL, conf.level = 0.95, 
alpha.lower = NULL, alpha.upper = NULL, df.error = NULL, ...)
}

\arguments{
  \item{means}{ a vector of the group means or the means of the particular
   level of the effect (for fixed effect designs) }
  
  \item{s.anova}{the standard deviation of the errors from the ANOVA model (i.e., the square root of the mean square error)}
  \item{c.weights}{the contrast weights (chose weights so that the positive \emph{c}-weights sum to 1 and the negative \emph{c}-weights sum to -1; i.e., use fractional values not integers).}
  \item{n}{ sample sizes \emph{per group} or sample sizes for the level of the particular factor (if length 1 it is 
  assumed that the sample size \emph{per group} or for the level of the particular factor are are equal) }
  
  \item{N}{ total sample size }
  \item{Psi}{ the (unstandardized) contrast effect, obtained by multiplying the \emph{j}th mean 
  by the \emph{j}th contrast weight (this is the unstandardized effect) }
  
  \item{ncp}{ the noncentrality parameter from the \emph{t}-distribution }
  \item{conf.level}{ desired level of confidence for the computed interval (i.e., 1 - the Type I error rate) }
  \item{alpha.lower}{ the Type I error rate for the lower confidence interval limit }
  \item{alpha.upper}{ the Type I error rate for the upper confidence interval limit }
  \item{df.error}{ the degrees of freedom for the error. In one-way designs, this is simply 
  \emph{N}-length (means) and need not be specified; it must be specified if the design has multiple factors. }
  
  \item{\dots}{ optional additional specifications for nested functions }
}
%\details{}
\value{
  \item{Lower.Conf.Limit.Standardized.Contrast }{the lower confidence limit for the standardized contrast}
  \item{Standardized.contrast }{standardized contrast}
  \item{Upper.Conf.Limit.Standardized.Contrast}{the upper confidence limit for the standardized contrast}
}

\references{ 
Kelley, K. (2007). Constructing confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20} (8), 1--24.

Lai, K., & Kelley, K. (2007). Sample size planning for standardized ANCOVA and ANOVA 
contrasts: Obtaining narrow confidence intervals. \emph{Manuscript submitted for publication}.

Steiger, J. H. (2004). Beyond the \emph{F} Test: Effect size confidence intervals and tests of close fit in the 
Analysis of Variance and Contrast Analysis. \emph{Psychological Methods, 9}, 164--182.  

}

\author{ Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu}) }

\note{
Be sure to use the standard deviation and not the error variance for \code{s.anova}, not the square of this value (the error variance) which would come from the source table (i.e., do not use the variance of the error but rather use its square root, the standard deviation).

Be sure to use the error variance and not its square root (i.e., use the variance of the  standard deviation of the errors).
Be sure to use the standard deviations of errors for \code{s.anova} and \code{s.ancova}, not the square of these values (i.e., do not use the variance of the errors).

Be sure to use fractional \emph{c}-weights when doing complex contrasts (not integers) to specify \code{c.weights}. For exmaple, in an ANCOVA of four groups, if the user wants to compare the mean of group 1 and 2 with the mean of group 3 and 4, \code{c.weights} should be specified as c(0.5, 0.5, -0.5, -0.5) rather than c(1, 1, -1, -1). Make sure the sum of the contrast weights are zero.


}
\seealso{ \code{conf.limits.nct}, \code{ci.src}, \code{ci.smd}, \code{ci.smd.c}, \code{ci.sm}, \code{ci.c} }
\examples{
# Here is a four group example. Suppose that the means of groups 1--4 are 2, 4, 9, 
# and 13, respectively. Further, let the error variance be .64 and thus the standard
# deviation would be .80 (note we use the standard deviation in the function, not the 
# variance). The standardized contrast of interest here is the average of groups 1 and 4
# versus the average of groups 2 and 3. 

ci.sc(means=c(2, 4, 9, 13), s.anova=.80, c.weights=c(.5, -.5, -.5, .5), 
n=c(3, 3, 3, 3), N=12, conf.level=.95)


# Here is an example with two groups. 
ci.sc(means=c(1.6, 0), s.anova=.80, c.weights=c(1, -1), n=c(10, 10), N=20, conf.level=.95)
}
\keyword{ design}
