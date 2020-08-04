\name{ci.src}
\alias{ci.src}

\title{ Confidence Interval for a Standardized Regression Coefficient }
\description{
Function to obtain the confidence interval for a standardized regression coefficient.
}
\usage{
ci.src(beta.k = NULL, SE.beta.k = NULL, N = NULL, K = NULL, R2.Y_X = NULL, 
R2.k_X.without.k = NULL, conf.level = 0.95, R2.Y_X.without.k = NULL, 
t.value = NULL, b.k = NULL, SE.b.k = NULL, s.Y = NULL, s.X = NULL, 
alpha.lower = NULL, alpha.upper = NULL, Suppress.Statement = FALSE, ...)
}

\arguments{
  \item{beta.k}{ the standardized regression coefficient }
  \item{SE.beta.k}{the standard error of the standarized regression coefficient }
  \item{N}{ sample size }
  \item{K}{ the number of predictors }
  \item{R2.Y_X}{ the squared multiple correlation coefficient predicting \emph{Y} from the \emph{k} predictor variables }
  \item{R2.k_X.without.k}{ the squared multiple correlation coefficient predicting the \emph{k}th predictor variable
(i.e., the predictor of interest) from the remaining \emph{p}-1 predictor variables }
  \item{conf.level}{ desired level of confidence for the computed interval (i.e., 1 - the Type I error rate)}
  \item{R2.Y_X.without.k}{ the squared multiple correlation coefficient predicting \emph{Y} from the \emph{p}-1 predictor
variable with the \emph{k}th predictor of interest excluded }
  \item{t.value}{ the \emph{t}-value evaluating the null hypothesis that the population regression coefficient
for the \emph{k}th predictor equals zero }
  \item{b.k}{ the unstandardized regression coefficient }
  \item{SE.b.k}{ the standard error of the unstandardized regression coefficient }
  \item{s.Y}{ standard deviation of \emph{Y}, the dependent variable }
  \item{s.X}{ standard deviation of \emph{X}, the predictor variable of interest }
  \item{alpha.lower}{  the Type I error rate for the lower confidence interval limit }
  \item{alpha.upper}{  the Type I error rate for the upper confidence interval limit }
  \item{Suppress.Statement}{ \code{TRUE} or \code{FALSE} statement specifying whether or not a statement should be printed
that identifies the type of confidence interval formed }
  \item{\dots}{ optional additional specifications for nested functions }
}
\details{
For standardized variables, do not specify the standard deviation of the variables and input the
standardized regression coefficient for \code{b.k}.
}
\value{
Returns the confidence limits specified for the regression coefficient of interest from the standard
approach to confidence interval formation or from the noncentral approach to confidence interval
formation using the noncentral \emph{t}-distribution.

}
\references{ 
Kelley, K. (2007). Constructing confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20} (8), 1--24.

Kelley, K., & Maxwel, S. E. (2003). Sample size for Multiple Regression: Obtaining regression
coefficients that are accurate, not simply significant. \emph{Psychological Methods, 8}, 305--321.

Kelley, K., & Maxwell, S. E. (2008). Sample Size Planning with applications to multiple regression: Power and accuracy for omnibus and targeted effects. In P. Alasuuta, J. Brannen, & L. Bickman (Eds.), \emph{The Sage handbook of social research methods} (pp. 166--192). Newbury Park, CA: Sage.

Smithson, M. (2003). \emph{Confidence intervals}. New York, NY: Sage Publications. 

Steiger, J. H. (2004). Beyond the \emph{F} Test: Effect size confidence intervals and tests of close fit in the 
Analysis of Variance and Contrast Analysis. \emph{Psychological Methods, 9}, 164--182.
}

\author{ Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu}) }

\note{ 
This function calls upon \code{ci.reg.coef} in MBESS, but has a different naming scheme. See  \code{ci.reg.coef} for more details.

To form a confidence interval for the unstandardized regression coefficient, use \code{ci.rc}. This function is used to form a confidence interval for the 
standardized regression coefficient.

Not all of the values need to be specified, only those that contain all of the necessary information
in order to compute the confidence interval (options are thus given for the values that need to be
specified).

}

\seealso{ \code{ss.aipe.reg.coef}, \code{conf.limits.nct}, \code{ci.reg.coef}, \code{ci.rc} }

%\examples{}

\keyword{ design }
