\name{ci.reg.coef}
\alias{ci.reg.coef}
\title{Confidence interval for a regression coefficient}
\description{A function to calculate a confidence interval around the population 
regression coefficient of interest using the standard approach and the noncentral
approach when the regression coefficients are standardized.}
\usage{
ci.reg.coef(b.j, SE.b.j=NULL, s.Y=NULL, s.X=NULL, N, p, R2.Y_X=NULL,
R2.j_X.without.j=NULL, conf.level=0.95, R2.Y_X.without.j=NULL, 
t.value=NULL, alpha.lower=NULL, alpha.upper=NULL, Noncentral=FALSE, 
Suppress.Statement=FALSE, ...)}
\arguments{
  \item{b.j}{value of the regression coefficient for the \emph{j}th predictor variable}
  \item{SE.b.j}{standard error for the \emph{j}th predictor variable}
  \item{s.Y}{standard deviation of \emph{Y}, the dependent variable}
  \item{s.X}{standard deviation of \eqn{X_j}, the predictor variable of interest}
  \item{N}{sample size}
  \item{p}{the number of predictors}
  \item{R2.Y_X}{the squared multiple correlation coefficient predicting \code{Y} from the \code{p} predictor variables}
  \item{R2.j_X.without.j}{the squared multiple correlation coefficient predicting the \emph{\code{j}}th predictor variable (i.e., the predictor of interest) from the remaining \code{p}-1 predictor variables}
  \item{conf.level}{desired level of confidence for the computed interval (i.e., 1 - the Type I error rate)}
  \item{R2.Y_X.without.j}{the squared multiple correlation coefficient predicting \code{Y} from the \emph{\code{p}}-1 predictor variable with the \code{j}th predictor of interest excluded}
  \item{t.value}{the \emph{t}-value evaluating the null hypothesis that the population regression coefficient for the \code{j}th predictor equals zero}
  \item{alpha.lower}{the Type I error rate for the lower confidence interval limit}
  \item{alpha.upper}{the Type I error rate for the upper confidence interval limit}
  \item{Noncentral}{\code{TRUE} or \code{FALSE}, specifying whether or not the noncentral approach to confidence intervals should be used}
  \item{Suppress.Statement}{\code{TRUE}/\code{FALSE} statement specifying whether or not a statement should be printed that identifies the type of confidence interval formed}
  \item{...}{optional additional specifications for nested functions}
}
\details{
For standardized variables, do not specify the standard deviation of the variables and input the standardized
regression coefficient for \code{b.j}.}

\value{Returns the confidence limits specified for the regression coefficient of interest from the standard approach to confidence interval formation or from the noncentral approach to confidence interval formation using the noncentral \emph{t}-distribution.}

\references{

Kelley, K. & Maxwell, S. E. (2003). Sample size for Multiple Regression: Obtaining regression coefficients that are accurate, not simply significant. 
\emph{Psychological Methods, 8}, 305--321.

Kelley, K. & Maxwell, S. E. (2008). Sample Size Planning with applications to multiple regression: Power and accuracy for omnibus and targeted effects. In P. Alasuuta, J. Brannen, & L. Bickman (Eds.), \emph{The Sage handbook of social research methods} (pp. 166--192). Newbury Park, CA: Sage.

Smithson, M. (2003). \emph{Confidence intervals}. New York, NY: Sage Publications. 
}

\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu}) }

\note{Not all of the values need to be specified, only those that
contain all of the necessary information in order to compute the
confidence interval (options are thus given for the values that need
to be specified).

The function \code{ci.rc} in MBESS also calculates the confidence interval
for the population (unstandardized) regression coefficient. The 
function \code{ci.src} also calculates the confidence interval
for the population (standardized) regression coefficient. These two 
functions perform the same tasks as \code{ci.reg.coef} does and 
are preferred to it because of simpler arguments.
}

\seealso{
\code{ss.aipe.reg.coef}, \code{conf.limits.nct}, \code{ci.rc}, \code{ci.src}
}

%\examples{}
\keyword{htest}
