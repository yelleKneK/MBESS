\name{ci.smd}
\alias{ci.smd}
\title{Confidence limits for the standardized mean difference.}
\description{
Function to calculate the confidence limits for the population standardized 
mean difference using the square root of the pooled variance as the divisor.
This function is thus used to determine the confidence bounds for the population 
quantity of what is generally referred to as Cohen's \emph{d} (delta being 
that population quantity).}
\usage{
ci.smd(ncp=NULL, smd=NULL, n.1=NULL, n.2=NULL, conf.level=.95, 
alpha.lower=NULL, alpha.upper=NULL, tol=1e-9, ...)

}
%- maybe also 'usage' for other objects documented here. 
\arguments{
  \item{ncp}{ is the estimated noncentrality parameter, this is generally the observed \emph{t}-statistic from comparing the two groups and assumes homogeneity of variance}
  \item{smd}{ is the standardized mean difference (using the pooled standard deviation in the denominator)}
  \item{n.1}{ is the sample size for Group 1}
  \item{n.2}{ is the sample size for Group 2}
  \item{conf.level}{ is the confidence level (1-Type I error rate)}
  \item{alpha.lower}{ is the Type I error rate for the lower tail}
  \item{alpha.upper}{ is the Type I error rate for the upper tail}
  \item{tol}{ is the tolerance of the iterative method for determining the critical values}
  \item{\dots}{ allows one to potentially include parameter values for inner functions}
}

\value{
  \item{Lower.Conf.Limit.smd}{The lower bound of the computed confidence interval}
  \item{smd}{The standardized mean difference}
  \item{Upper.Conf.Limit.smd}{The upper bound of the computed confidence interval}
}

\references{
Cohen, J. (1988) Statistical power analysis for the behavioral sciences (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.

Cumming, G. & Finch, S. (2001). A primer on the understanding, use, and calculation of confidence intervals that are
based on central and noncentral distributions, \emph{Educational and Psychological Measurement, 61}, 532--574.

Hedges, L. V. (1981). Distribution theory for Glass's Estimator of effect size and related estimators. \emph{Journal of Educational Statistics, 2}, 107--128.

Kelley, K. (2007). Constructing confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20} (8), 1--24.

Kelley, K., Maxwell, S. E., & Rausch, J. R. (2003). Obtaining Power or Obtaining Precision: Delineating Methods
of Sample-Size Planning, \emph{Evaluation and the Health Professions, 26}, 258--287.

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval estimation and the evaluation of
statistical methods. In L. L. Harlow, S. A. Mulaik,&J.H. Steiger (Eds.), \emph{What if there were
no significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.
}

\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})}

\section{Warning}{
This function uses \code{conf.limits.nct}, which has as one of its arguments \code{tol} 
(and can be modified with \code{tol} of the present function). 
If the present function fails to converge (i.e., if it runs but does not report a solution), 
it is likely that the \code{tol} value is too restrictive and should be increased by a factor of 10, but probably by no more than 100. 
Running the function \code{conf.limits.nct} directly will report the actual probability values of the limits found. This should be 
done if any modification to \code{tol} is necessary in order to ensure acceptable confidence limits for the noncentral-\emph{t} parameter have been achieved.
}

\seealso{\code{\link{smd}}, \code{\link{smd.c}}, \code{\link{ci.smd.c}}, \code{\link{conf.limits.nct}}}

\examples{
# Steiger and Fouladi (1997) example values.
ci.smd(ncp=2.6, n.1=10, n.2=10, conf.level=1-.05)
ci.smd(ncp=2.4, n.1=300, n.2=300, conf.level=1-.05)
}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{htest}% __ONLY ONE__ keyword per line
