\name{conf.limits.nct}
\alias{conf.limits.nct}

\title{Confidence limits for a noncentrality parameter from a t-distribution}

\description{
Function to determine the noncentrality parameters necessary to form a confidence interval around the population noncentrality parameter and related parameters.
}

\usage{
conf.limits.nct(ncp, df, conf.level = 0.95, alpha.lower = NULL, 
alpha.upper = NULL, t.value, tol = 1e-09, sup.int.warns = TRUE, 
...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{ncp}{the noncentrality parameter (e.g., observed \emph{t}-value) of interest.}
  \item{df}{the degrees of freedom.}
  \item{conf.level}{the level of confidence for a symmetric confidence interval.}
  \item{alpha.lower}{the proportion of values beyond the lower limit of the confidence interval (cannot be used with \code{conf.level}).}
  \item{alpha.upper}{the proportion of values beyond the upper limit of the confidence interval (cannot be used with \code{conf.level}).}
  \item{t.value}{alias for \code{ncp}}
  \item{tol}{is the tolerance of the iterative method for determining the critical values.}
  \item{sup.int.warns}{Suppress internal warnings (from internal functions): \code{TRUE} or \code{FALSE}}
  \item{\dots}{ allows one to potentially include parameter values for inner functions}
}
\details{
Function for finding the upper and lower confidence limits for a noncentral parameter from a noncentral \emph{t}-distribution with \code{df} degrees of freedom. 
This function is especially helpful when forming confidence intervals around standardized mean differences (i.e., Cohen's \emph{d}; Glass's \emph{g}; Hedges' \emph{g}), standardized regression coefficients, and 
coefficients of variations. The \code{Lower.Limit} and the \code{Upper.Limit} values correspond to the noncentral parameters of a \emph{t}-distribution with \code{df} degrees of
freedom whose upper and lower tails contain the desired proportion of the respective noncentral \emph{t}-distribution.
When \code{ncp} is zero, the \code{Lower.Limit} and \code{Upper.Limit} are simply the desired quantiles of the
central \emph{t}-distribution with \code{df} degrees of freedom.

Note that the confidence interval limit(s) are found twice, using two different methods. The first method uses the \code{optimize} function, whereas the second method uses the \code{nlm} function. The best of the two methods, if not equal and numerically exact, is taken. This does not concern the user. 
}

\value{ 
  \item{Lower.Limit}{Value of the distribution with \code{Lower.Limit} noncentral value that has at its specified quantile \code{F.value}}
  \item{Prob.Less.Lower}{Proportion of the distribution beyond (i.e., less than) \code{Lower.Limit}}
  \item{Upper.Limit}{Value of the distribution with \code{Upper.Limit} noncentral value that has at its specified quantile \code{F.value}}
  \item{Prob.Greater.Upper}{Proportion of the distribution beyond (i.e., larger than) \code{Upper.Limit}}
  }
\references{
Cumming, G. & Finch, S. (2001). A primer on the understanding, use, and calculation of confidence intervals that are
based on central and noncentral distributions, \emph{Educational and Psychological Measurement, 61}, 532--574.

Kelley, K. (2005). The effects of nonnormal distributions on confidence intervals around the standardized mean
difference: Bootstrap and parametric confidence intervals, \emph{Educational and Psychological Measurement, 65}, 51--69.

Kelley, K. (2007). Constructing confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20} (8), 1--24.

Steiger, J. & Fouladi, T. (1997). Noncentrality interval estimation and the evaluation of statistical models. In L. Harlow, 
S. Muliak, & J. Steiger (Eds.), \emph{What if there were no significance tests?}. Mahwah, NJ: Lawrence Erlbaum.

}
\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})}

\section{Warning}{At the present time, the largest \code{ncp} that R can accurately handle is 37.62.
}

\seealso{\code{pt}, \code{qt}, \code{ci.smd}, \code{ci.smd.c}, \code{ss.aipe}, \code{conf.limits.ncf}, \code{conf.limits.nc.chisq}}
\examples{

# Suppose observed t-value based on 'df'=126 is 2.83. Finding the lower 
# and upper critical values for the population noncentrality parameter 
# with a symmetric confidence interval with 95\% confidence is given as:
conf.limits.nct(ncp=2.83, df=126, conf.level=.95)

# Modifying the above example so that a nonsymmetric 95% confidence interval
# can be formed:
conf.limits.nct(ncp=2.83, df=126, alpha.lower=.01, alpha.upper=.04,
conf.level=NULL)

# Modifying the above example so that a single-sided 95% confidence interval
# can be formed:
conf.limits.nct(ncp=2.83, df=126, alpha.lower=0, alpha.upper=.05,
conf.level=NULL)

}
\keyword{models}
\keyword{htest}
