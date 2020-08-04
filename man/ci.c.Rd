\name{ci.c}
\alias{ci.c}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Confidence interval for a contrast in a fixed effects ANOVA }
\description{
Function to calculate the exact confidence interval for a contrast in a fixed effects analysis of variance context. This function assumes homogeneity of variance (as does the ANOVA upon which 's.anova' is based).
}
\usage{
ci.c(means = NULL, s.anova = NULL, c.weights = NULL, n = NULL, 
N = NULL, Psi = NULL, conf.level = 0.95, alpha.lower = NULL, 
alpha.upper = NULL, df.error = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{means}{a vector of the group means or the means of the particular level of the effect (for fixed effect designs) }
  \item{s.anova}{the standard deviation of the errors from the ANOVA model (i.e., the square root of the mean square error)}
  \item{c.weights}{the contrast weights (choose weights so that the positive \emph{c}-weights sum to 1 and the negative \emph{c}-weights sum to -1; i.e., use fractional values not integers).}
  \item{n}{ sample sizes \emph{per group} or level of the particular factor (if length 1 it is assumed that the per group/level sample sizes are equal) }
  \item{N}{ total sample size }
  \item{Psi}{ the (unstandardized) contrast effect, obtained by multiplying the \emph{j}th mean by the jth contrast weight (this is the unstandardized effect) }
  \item{conf.level}{ confidence interval coverage (i.e., 1- Type I error rate); default is .95 }
  \item{alpha.lower}{ Type I error for the lower confidence limit }
  \item{alpha.upper}{ Type I error for the upper confidence limit}
  \item{df.error}{ the degrees of freedom for the error. In one-way designs, this is simply \emph{N}-length (means) and need not be specified; it must be specified if the design has multiple factors. }
  \item{\dots}{ allows one to potentially include parameter values for inner functions }
}
%\details{}
\value{
Returns the confidence limits for the contrast: 
  
  \item{Lower.Conf.Limit.Contrast }{the lower confidence limit for the contrast effect}
  \item{Contrast }{the value of the estimated unstandardized contrast effect}
  \item{Upper.Conf.Limit.Contrast}{the upper confidence limit for the contrast effect}
}
\references{
Kelley, K. (2007). Constructing confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20} (8), 1--24.

Steiger, J. H. (2004). Beyond the \emph{F} Test: Effect size confidence intervals and tests of close fit in the 
analysis of variance and contrast analysis. \emph{Psychological Methods, 9}, 164--182.
}

\author{ Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu}) }

\note{Be sure to use the standard deviation and not the error variance for \code{s.anova}, not the square of this value (the error variance) which would come from the source table (i.e., do not use the variance of the error but rather use its square root, the standard deviation).

Be sure to use fractional \emph{c}-weights when doing complex contrasts (not integers) to specify \code{c.weights}. For example, in an ANCOVA of four groups, if the user wants to compare the mean of group 1 and 2 with the mean of group 3 and 4, \code{c.weights} should be specified as c(0.5, 0.5, -0.5, -0.5) rather than c(1, 1, -1, -1). Make sure the sum of the contrast weights is zero.

}

\seealso{ \code{conf.limits.nct}, \code{ci.sc}, \code{ci.src}, \code{ci.smd}, \code{ci.smd.c}, \code{ci.sm}}

\examples{
# Here is a four group example. Suppose that the means of groups 1--4 are 2, 4, 9, 
# and 13, respectively. Further, let the error variance be .64 and thus the standard
# deviation would be .80 (note we use the standard deviation in the function, not the 
# variance). The contrast of interest here is the average of groups 1 and 4 versus the 
# average of groups 2 and 3. 

ci.c(means=c(2, 4, 9, 13), s.anova=.80, c.weights=c(.5, -.5, -.5, .5), 
n=c(3, 3, 3, 3), N=12, conf.level=.95)


# Here is an example with two groups. 
ci.c(means=c(1.6, 0), s.anova=.80, c.weights=c(1, -1), n=c(10, 10), N=20, conf.level=.95)


# An example given by Maxwell and Delaney (2004, pp. 155--171) :
# 24 subjects of mild hypertensives are assigned to one of four treatments: drug 
# therapy, biofeedback, dietary modification, and a treatment combining all the 
# three previous treatments. Subjects' blood pressure is measured two weeks
# after the termination of treatment. Now we want to form a 95% level
# confidence interval for the difference in blood pressure between subjects
# who received drug treatment and those who received biofeedback treatment 

## Drug group's mean = 94; group size=4
## Biofeedback group's mean = 91; group size=6 
## Diet group's mean = 92; group size=5
## Combination group's mean = 83; group size=5
## Mean Square Within (i.e., 'error.variance') = 67.375

ci.c(means=c(94, 91, 92, 83), s.anova=sqrt(67.375), c.weights=c(1, -1, 0, 0), 
n=c(4, 6, 5, 5), N=20, conf.level=.95)

}
\keyword{design}% at least one, from doc/KEYWORDS

