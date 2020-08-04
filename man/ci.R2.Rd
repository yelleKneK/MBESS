\name{ci.R2}
\alias{ci.R2}
\title{Confidence interval for the population squared multiple correlation coefficient}
\description{A function to calculate the confidence interval for the population squared multiple correlation coefficient.}
\usage{
ci.R2(R2 = NULL, df.1 = NULL, df.2 = NULL, conf.level = .95, 
Random.Predictors=TRUE, Random.Regressors, F.value = NULL, N = NULL, 
p = NULL, K, alpha.lower = NULL, alpha.upper = NULL, tol = 1e-09)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R2}{squared multiple correlation coefficient}
  \item{df.1}{numerator degrees of freedom}
  \item{df.2}{denominator degrees of freedom}
  \item{conf.level}{confidence interval coverage; 1-Type I error rate}
  \item{Random.Predictors}{whether or not the predictor variables are random or fixed (random is default)}
  \item{Random.Regressors}{an alias for \code{Random.Predictors}; \code{Random.Regressors} overrides 
  \code{Random.Predictors}}
  \item{F.value}{obtained \emph{F}-value}
  \item{N}{sample size}
  \item{p}{number of predictors}
  \item{K}{alias for \code{p}, the number of predictors}
  \item{alpha.lower}{Type I error for the lower confidence limit}
  \item{alpha.upper}{Type I error for the upper confidence limit}
  \item{tol}{tolerance for iterative convergence}}

\details{
This function can be used with random predictor variables (\code{Random.Predictors=TRUE}) or when predictor
variables are fixed (\code{Random.Predictors=FALSE}). In many applications of multiple regression, 
predictor variables are random, which is the default in this function. 

For random predictors, the function implements the procedure of Lee (1971), which was implemented by 
Algina and Olejnik (2000; specifically in their \emph{ci.smcc.bisec.sas} SAS script). When \code{Random.Predictors=TRUE},
the function implements code that is in part based on the Alginia and Olejnik (2000) SAS script. 

When \code{Random.Predictors=FALSE}, and thus the predictors are planned and thus fixed in 
hypothetical replications of the study, the confidence limits are based on a
noncentral \eqn{F}-distribution (see \code{conf.limits.ncf}).
}

\value{
  \item{Lower.Conf.Limit.R2 }{upper limit of the confidence interval around the population multiple correlation coefficient}
  \item{Prob.Less.Lower }{proportion of the distribution less than \code{Lower.Conf.Limit.R2}}
  \item{Upper.Conf.Limit.R2 }{upper limit of the confidence interval around the population multiple correlation coefficient}
  \item{Prob.Greater.Upper }{proportion of the distribution greater than \code{Upper.Conf.Limit.R2}}
}
\references{
Algina, J. & Olejnik, S. (2000). Determining Sample Size for Accurate Estimation of 
the Squared Multiple Correlation Coefficient. \emph{Multivariate Behavioral Research, 35}, 
119--136.

Kelley, K. (2007). Constructing confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20} (8), 1--24.

Lee, Y. S. (1971). Some results on the sampling distribution of the multiple correlation coefficient.
\emph{Journal of the Royal Statistical Society, B, 33}, 117--130.

Smithson, M. (2003). \emph{Confidence intervals}. New York, NY: Sage Publications. 

Steiger, J. H. & Fouladi, R. T. (1992) R2: A computer program for interval estimation, power calculation,
and hypothesis testing for the squared multiple correlation. \emph{Behavior research methods, instruments and computers, 4}, 581--582.

}
\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})}

\seealso{\code{ss.aipe.R2}, \code{conf.limits.ncf}}
\examples{
# For random predictor variables.
# ci.R2(R2=.25, N=100, K=5, conf.level=.95, Random.Predictors=TRUE)

# ci.R2(F.value=6.266667, N=100, K=5, conf.level=.95, Random.Predictors=TRUE)

# For fixed predictor variables.
# ci.R2(R2=.25, N=100, K=5, conf.level=.95, Random.Predictors=TRUE)

# ci.R2(F.value=6.266667, N=100, K=5, conf.level=.95, Random.Predictors=TRUE)

# One sided confidence intervals when predictors are random.
# ci.R2(R2=.25, N=100, K=5, alpha.lower=.05, alpha.upper=0, conf.level=NULL,
# Random.Predictors=TRUE)

# ci.R2(R2=.25, N=100, K=5, alpha.lower=0, alpha.upper=.05, conf.level=NULL,
# Random.Predictors=TRUE)

# One sided confidence intervals when predictors are fixed.
# ci.R2(R2=.25, N=100, K=5, alpha.lower=.05, alpha.upper=0, conf.level=NULL,
# Random.Predictors=FALSE)

# ci.R2(R2=.25, N=100, K=5, alpha.lower=0, alpha.upper=.05, conf.level=NULL,
# Random.Predictors=FALSE)
}
\keyword{multivariate}
\keyword{htest}
\keyword{regression}
