\name{ci.smd.c}
\alias{ci.smd.c}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Confidence limits for the standardized mean difference using the control
group standard deviation as the divisor.}
\description{
Function to calculate the confidence limits for the standardized mean difference using the control group standard deviation 
as the divisor (Glass's \emph{g}).}

\usage{
ci.smd.c(ncp = NULL, smd.c = NULL, n.C = NULL, n.E = NULL, 
conf.level = 0.95, alpha.lower = NULL, alpha.upper = NULL, 
tol = 1e-09, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ncp}{ is the estimated noncentrality parameter, this is generally the observed \emph{t}-statistic from comparing the control and experimental group (assuming homogeneity of variance)}
  \item{smd.c}{ is the standardized mean difference (using the control group standard deviation in the denominator)}
  \item{n.C}{ is the sample size for the control group}
  \item{n.E}{ is the sample size for experimental group}
  \item{conf.level}{ is the confidence level (1-Type I error rate)}
  \item{alpha.lower}{ is the Type I error rate for the lower tail}
  \item{alpha.upper}{ is the Type I error rate for the upper tail}
  \item{tol}{ is the tolerance of the iterative method for determining the critical values}
  \item{\dots}{ Potentially include parameter for inner functions}
}

\value{
  \item{Lower.Conf.Limit.smd.c}{The lower bound of the computed confidence interval}
  \item{smd.c}{The standardized mean difference based on the control group standard deviation}
  \item{Upper.Conf.Limit.smd.c}{The upper bound of the computed confidence interval}
}

\references{
Cohen, J. (1988). \emph{Statistical power analysis for the behavioral sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.

Cumming, G. & Finch, S. (2001). A primer on the understanding, use, and calculation of confidence intervals that are
based on central and noncentral distributions, \emph{Educational and Psychological Measurement, 61}, 532--574.

Glass, G. V. (1976). Primary, secondary, and meta-analysis of research. \emph{Educational Researcher, 5}, 3--8.

Hedges, L. V. (1981). Distribution theory for Glass's Estimator of effect size and related estimators. \emph{Journal of Educational Statistics, 2}, 107--128.

Kelley, K. (2007). Constructing confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20} (8), 1--24.

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval estimation and the evaluation of
statistical methods. In L. L. Harlow, S. A. Mulaik, & J. H. Steiger (Eds.), \emph{What if there were
no significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.
}

\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})}

\section{Warning}{
This function uses \code{conf.limits.nct}, which has as one of its arguments \code{tol} 
(and can be modified with \code{tol} of the present function). 
If the present function fails to converge (i.e., if it runs but does not report a solution), 
it is likely that the \code{tol} value is too restrictive and should be increased by a factor of 10, but probably by no more than 100. 
Running the function \code{conf.limits.nct} directly will report the actual probability values of the limits found. This should be 
done if any modification to \code{tol} is necessary in order to ensure acceptable confidence limits for the noncentral-\emph{t} 
parameter have been achieved.
}

\seealso{\code{smd.c}, \code{smd}, \code{ci.smd}, \code{conf.limits.nct}}
\examples{
ci.smd.c(smd.c=.5, n.C=100, n.E=100, conf.level=.95)
}
\keyword{models}
\keyword{htest}
