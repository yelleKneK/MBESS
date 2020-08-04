\name{ci.R}
\alias{ci.R}

\title{ Confidence interval for the multiple correlation coefficient }

\description{
A function to obtain the confidence interval for the population multiple correlation coefficient when predictors are random (the default) or fixed.
}

\usage{
ci.R(R = NULL, df.1 = NULL, df.2 = NULL, conf.level = 0.95, 
Random.Predictors = TRUE, Random.Regressors, F.value = NULL, 
N = NULL, K=NULL, alpha.lower = NULL, alpha.upper = NULL, ...)
}

\arguments{
  \item{R}{multiple correlation coefficient }
  \item{df.1}{ numerator degrees of freedom }
  \item{df.2}{denominator degrees of freedom }
  \item{conf.level}{confidence interval coverage (i.e., 1- Type I error rate); default is .95 }
  \item{Random.Predictors}{whether or not the predictor variables are random or fixed (random is default) }
  \item{Random.Regressors}{an alias for \code{Random.Predictors}; \code{Random.Regressors} overrides
   \code{Random.Predictors} }
  \item{F.value}{ obtained \emph{F}-value }
  \item{N}{sample size }
  \item{K}{number of predictors }
  \item{alpha.lower}{Type I error for the lower confidence limit }
  \item{alpha.upper}{Type I error for the upper confidence limit }
  \item{\dots}{allows one to potentially include parameter values for inner functions }
}

\details{
This function is based on the function \code{ci.R2} in MBESS package.

This function can be used with random predictor variables (\code{Random.Predictors=TRUE}) or
when predictor variables are fixed (\code{Random.Predictors=FALSE}). In many applications in the behavioral,
educational, and social sciences, predictor variables are random, which is the default for this
function.

For random predictors, the function implements the procedure of Lee (1971), which was implemented
by Algina and Olejnik (2000; specifically in their \emph{ci.smcc.bisec.sas} SAS script). When 
\code{Random.Predictors=TRUE}, the function implements code that is in part based on the Alginia and
Olejnik (2000) SAS script.

When \code{Random.Predictors=FALSE}, and thus the predictors are planned and thus fixed in hypothetical
replications of the study, the confidence limits are based on a noncentral \emph{F}-distribution (see
\code{conf.limits.ncf}).
}
\value{
  \item{Lower.Conf.Limit.R }{lower limit of the confidence interval around the population multiple correlation coefficient}
  \item{Prob.Less.Lower }{proportion of the distribution less than \code{Lower.Conf.Limit.R}}
  \item{Upper.Conf.Limit.R}{upper limit of the confidence interval around the population multiple correlation coefficient}
  \item{Prob.Greater.Upper}{proportion of the distribution greater than \code{Upper.Conf.Limit.R}}
}
\references{
Algina, J. & Olejnik, S. (2000). Determining sample size for accurate estimation of the squared
multiple correlation coefficient. \emph{Multivariate Behavioral Research, 35}, 119--136.

Lee, Y. S. (1971). Some results on the sampling distribution of the multiple correlation coefficient.
\emph{Journal of the Royal Statistical Society, B, 33}, 117--130.

Smithson, M. (2003). \emph{Confidence intervals}. New York, NY: Sage Publications.

Steiger, J. H. (2004). Beyond the \emph{F} Test: Effect size confidence intervals and tests of close fit in the 
Analysis of Variance and Contrast Analysis. \emph{Psychological Methods, 9}, 164--182.

Steiger, J. H. & Fouladi, R. T. (1992). R2: A computer program for interval estimation, power
calculation, and hypothesis testing for the squared multiple correlation. \emph{Behavior research methods,
instruments and computers, 4}, 581--582.

}

\author{ Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})}
%\note{ }
\seealso{ \code{ci.R2}, \code{ss.aipe.R2}, \code{conf.limits.nct} }
%\examples{}
\keyword{ design }% at least one, from doc/KEYWORDS
\keyword{ regression }% __ONLY ONE__ keyword per line
