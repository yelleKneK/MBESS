\name{ss.aipe.sm}
\alias{ss.aipe.sm}
\title{Sample size planning for Accuracy in Parameter Estimation (AIPE) of the standardized mean}
\description{
A function to calculate the appropriate sample size for the standardized mean such that the 
width of the confidence interval is sufficiently narrow. 
}
\usage{
ss.aipe.sm(sm, width, conf.level = 0.95, assurance = NULL, certainty=NULL, ...)
}
\arguments{
  \item{sm}{the population standardized mean }
  \item{width}{ the desired full width of the obtained confidence interval }
  \item{conf.level}{the desired confidence interval coverage, (i.e., 1 - Type I error rate)}
  \item{assurance}{ parameter to ensure that the obtained confidence interval width is 
  narrower than the desired width with a specified degree of certainty (must be \code{NULL} or between zero and unity) }
  \item{certainty}{an alias for \code{assurance}}
  \item{\dots}{ allows one to potentially include parameter values for inner functions }
}

%\details{}

\value{
  \item{n }{the necessary sample size in order to achieve the desired degree of accuracy (i.e., the sufficiently narrow confidence interval)}
}

\references{ 
Cumming, G. & Finch, S. (2001). A primer on the understanding, use, and calculation of confidence intervals that are
based on central and noncentral distributions, \emph{Educational and Psychological Measurement, 61}, 532--574.

Hedges, L. V. (1981). Distribution theory for Glass's Estimator of effect size and related estimators. \emph{Journal of Educational Statistics, 2}, 107--128.

Kelley, K. (2005). The effects of nonnormal distributions on confidence intervals around the standardized mean
difference: Bootstrap and parametric confidence intervals, \emph{Educational and Psychological Measurement, 65}, 51--69.

Kelley, K. (2007). Constructing confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20} (8), 1--24.

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the standardized mean difference: 
Accuracy in Parameter Estimation via narrow confidence intervals. \emph{Psychological Methods, 11(4)}, 363--385.

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval estimation and the evaluation of
statistical methods. In L. L. Harlow, S. A. Mulaik,& J.H. Steiger (Eds.), \emph{What if there were
no significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum. }


\author{ Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu}); Keke Lai}

%\note{}

\seealso{ \code{conf.limit.nct}, \code{ci.sm} }

\examples{
# Suppose the population mean is believed to be 20, and the population
# standard deviation is believed to be 2; thus the population standardized
# mean is believed to be 10. To determine the necessary sample size for a 
# study so that the full width of the 95 percent confidence interval 
# obtained in the study will be, with 90% assurance, no wider than 2.5, 
# the function should be specified as follows. 

# ss.aipe.sm(sm=10, width=2.5, conf.level=.95, assurance=.90)
}

\keyword{ design }
