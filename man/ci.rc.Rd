\name{ci.rc}
\alias{ci.rc}
\title{ Confidence Interval for a Regression Coefficient }
\description{
A function to calculate a confidence interval for the population regression coefficient of interest
using the standard approach and the noncentral approach when the regression coefficients are
standardized.
}

\usage{
ci.rc(b.k, SE.b.k = NULL, s.Y = NULL, s.X = NULL, N, K, R2.Y_X = NULL, 
R2.k_X.without.k = NULL, conf.level = 0.95, R2.Y_X.without.k = NULL, 
t.value = NULL, alpha.lower = NULL, alpha.upper = NULL, 
Noncentral = FALSE, Suppress.Statement = FALSE, ...)
}

\arguments{
  \item{b.k}{ value of the regression coefficient for the \emph{k}th predictor variable }
  \item{SE.b.k}{ standard error for the \emph{k}th predictor variable }
  \item{s.Y}{ standard deviation of \emph{Y}, the dependent variable }
  \item{s.X}{ standard deviation of \emph{X}, the predictor variable of interest }
  \item{N}{ sample size }
   \item{K}{the number of predictors }
  \item{R2.Y_X}{ the squared multiple correlation coefficient predicting \emph{Y} from the \emph{k} predictor variables }
  \item{R2.k_X.without.k}{ the squared multiple correlation coefficient predicting the \emph{k}th predictor variable
(i.e., the predictor of interest) from the remaining \emph{K}-1 predictor variables }
  \item{conf.level}{ desired level of confidence for the computed interval (i.e., 1 - the Type I error rate)}
  \item{R2.Y_X.without.k}{ the squared multiple correlation coefficient predicting \emph{Y} from the \emph{K}-1 predictor
variable with the \emph{k}th predictor of interest excluded }
  \item{t.value}{ the \emph{t}-value evaluating the null hypothesis that the population regression coefficient
for the \emph{k}th predictor equals zero }
  \item{alpha.lower}{ the Type I error rate for the lower confidence interval limit }
  \item{alpha.upper}{ the Type I error rate for the upper confidence interval limit }
  \item{Noncentral}{ \code{TRUE} or \code{FALSE} statement specifying whether or not the noncentral approach to
confidence intervals should be used }
  \item{Suppress.Statement}{\code{TRUE} or \code{FALSE} statement specifying whether or not a statement should be printed
that identifies the type of confidence interval formed}
  \item{\dots}{ optional additional specifications for nested functions }
}

\details{
This function calls upon \code{ci.reg.coef} in MBESS, but has a different naming system. See \code{ci.reg.coef} for more details.

For standardized variables, do not specify the standard deviation of the variables and input the
standardized regression coefficient for \code{b.k}.
}

\value{
Returns the confidence limits for the standardized regression coefficients of interest from the standard
approach to confidence interval formation or from the noncentral approach to confidence interval
formation using the noncentral \emph{t}-distribution.
}

\references{ 
Kelley, K. (2007). Confidence intervals for standardized effect sizes: Theory, application, 
and implementation.\emph{ Journal of Statistical Software, 20}(8), 1--24.

Kelley, K. & Maxwell, S. E. (2003). Sample size for Multiple Regression: Obtaining regression
coefficients that are accurate, not simply significant.\emph{ Psychological Methods, 8}, 305--321.

Kelley, K. & Maxwell, S. E. (2008). Power and accuracy for omnibus and targeted effects:
Issues of sample size planning with applications to Multiple Regression. \emph{Handbook of Social Research
Methods}, J. Brannon, P. Alasuutari, and L. Bickman (Eds.). New York, NY: Sage Publications.

Smithson, M. (2003). \emph{Confidence intervals}. New York, NY: Sage Publications. 

Steiger, J. H. (2004). Beyond the \emph{F} Test: Effect size confidence intervals and tests of close fit in the 
Analysis of Variance and Contrast Analysis. \emph{Psychological Methods, 9}, 164--182.}

\author{ Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu}) }

\note{ Not all of the values need to be specified, only those that contain all of the necessary information
in order to compute the confidence interval (options are thus given for the values that need to be
specified).
}

\seealso{ \code{ss.aipe.reg.coef}, \code{conf.limits.nct}, \code{ci.reg.coef}, \code{ci.src} }

%\examples{}

\keyword{ design }% at least one, from doc/KEYWORDS
