\name{ss.aipe.sc}
\alias{ss.aipe.sc}
\title{Sample size planning for Accuracy in Parameter Estimation (AIPE) of the standardized contrast in ANOVA}
\description{
A function to calculate the appropriate sample size per group for the standardized contrast in ANOVA such that the 
width of the confidence interval is sufficiently narrow. 
}
\usage{
ss.aipe.sc(psi, c.weights, width, conf.level = 0.95, 
assurance = NULL, certainty = NULL, ...)
}
\arguments{
  \item{psi}{ population standardized contrast }
  \item{c.weights}{ the contrast weights }
  \item{width}{ the desired full width of the obtained confidence interval }
  \item{conf.level}{the desired confidence interval coverage, (i.e., 1 - Type I error rate)}
  \item{assurance}{ parameter to ensure that the obtained confidence interval width is narrower than the desired width with a specified degree of certainty (must be NULL or between zero and unity) }
  \item{certainty}{an alias for \code{assurance}}
  \item{\dots}{ allows one to potentially include parameter values for inner functions }
}
%\details{}
\value{
\item{n}{necessary sample size \emph{per group}}
}

\references{
Cumming, G. & Finch, S. (2001). A primer on the understanding, use, and calculation of confidence intervals that are
based on central and noncentral distributions, \emph{Educational and Psychological Measurement, 61}, 532--574.

Hedges, L. V. (1981). Distribution theory for Glass's Estimator of effect size and related estimators. \emph{Journal of Educational Statistics, 2}, 107--128.

Kelley, K. (2005). The effects of nonnormal distributions on confidence intervals around the standardized mean
difference: Bootstrap and parametric confidence intervals, \emph{Educational and Psychological Measurement, 65}, 51--69.

Kelley, K. (2007). Constructing confidence intervals for standardized effect sizes: Theory, application, 
and implementation. \emph{Journal of Statistical Software, 20} (8), 1--24.

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the standardized mean difference: 
Accuracy in Parameter Estimation via narrow confidence intervals. \emph{Psychological Methods, 11(4)}, 363--385.

Lai, K., & Kelley, K. (2007). Sample size planning for standardized ANCOVA and ANOVA 
contrasts: Obtaining narrow confidence intervals. \emph{Manuscript submitted for publication}.

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval estimation and the evaluation of
statistical methods. In L. L. Harlow, S. A. Mulaik, & J.H. Steiger (Eds.), \emph{What if there were
no significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum. }


\author{ Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu}); Keke Lai }
%\note{ }

\seealso{ \code{ci.sc}, \code{conf.limits.nct}, \code{ss.aipe.c} }
\examples{
# Suppose the population standardized contrast is believed to be .6 
# in some 5-group ANOVA model. The researcher is interested in comparing 
# the average of means of group 1 and 2 with the average of group 3 and 4.

# To calculate the necessary sample size per gorup such that the width 
# of 95 percent confidence interval of the standardized 
# contrast is, with 90 percent assurance, no wider than .4:

# ss.aipe.sc(psi=.6, c.weights=c(.5, .5, -.5, -.5, 0), width=.4, assurance=.90) 
}

\keyword{design}

