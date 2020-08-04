\name{ss.aipe.sm.sensitivity}
\alias{ss.aipe.sm.sensitivity}
\title{Sensitivity analysis for sample size planning for the standardized mean from the Accuracy in Parameter Estimation (AIPE) 
Perspective }
\description{Performs a sensitivity analysis when planning sample size from the Accuracy in Parameter Estimation (AIPE) 
Perspective for the standardized mean.
}
\usage{
ss.aipe.sm.sensitivity(true.sm = NULL, estimated.sm = NULL, 
desired.width = NULL, selected.n = NULL, assurance = NULL, 
certainty=NULL, conf.level = 0.95, G = 10000, print.iter = TRUE, 
detail = TRUE, ...)
}
\arguments{
  \item{true.sm}{ population standardized mean }
  \item{estimated.sm}{ estimated standardized mean }
  \item{desired.width}{ desired full width of the confidence interval for the population standardized mean }
  \item{selected.n}{ selected sample size to use in order to determine distributional properties of a given value of sample size }
  \item{assurance}{ parameter to ensure that the obtained confidence interval width is narrower 
  than the desired width with a specified degree of certainty (must be \code{NULL} or between zero and unity) }
  \item{certainty}{an alias for \code{assurance}}
  \item{conf.level}{the desired confidence interval coverage, (i.e., 1 - Type I error rate)}
  \item{G}{ number of generations (i.e., replications) of the simulation}
  \item{print.iter}{ to print the current value of the iterations}
  \item{detail}{whether the user needs a detailed (\code{TRUE}) or brief (\code{FALSE}) report of the simulation results; the 
    detailed report includes all the raw data in the simulations}
  \item{\dots}{ allows one to potentially include parameter values for inner functions }
}

%\details{}

\value{
\item{sm.obs}{vector of the observed standardized mean}
\item{Full.Width}{vector of the full confidence interval width}
\item{Width.from.sm.obs.Lower}{vector of the lower confidence interval width}
\item{Width.from.sm.obs.Upper}{vector of the upper confidence interval width}
\item{Type.I.Error.Upper}{iterations where a Type I error occurred on the upper end of the confidence interval}
\item{Type.I.Error.Lower}{iterations where a Type I error occurred on the lower end of the confidence interval}
\item{Type.I.Error}{iterations where a Type I error happens}
\item{Lower.Limit}{the lower limit of the obtained confidence interval}
\item{Upper.Limit}{the upper limit of the obtained confidence interval}
\item{replications}{number of replications of the simulation}
\item{True.sm}{the population standardized mean}
\item{Estimated.sm}{the estimated standardized mean}
\item{Desired.Width}{the desired full confidence interval width}
\item{assurance}{parameter to ensure that the obtained confidence interval width is narrower than the desired width with a specified degree of certainty}
\item{Sample.Size}{the sample size used in the simulation}
\item{mean.full.width}{mean width of the obtained full confidence intervals}
\item{median.full.width}{median width of the obtained full confidence intervals}
\item{sd.full.width}{standard deviation of the widths of the obtained full confidence intervals}
\item{Pct.Width.obs.NARROWER.than.desired}{percentage of the obtained full confidence interval widths that are narrower than the desired width }
\item{mean.Width.from.sm.obs.Lower}{mean lower width of the obtained confidence intervals}
\item{mean.Width.from.sm.obs.Upper}{mean upper width of the obtained confidence intervals}
\item{Type.I.Error.Upper}{Type I error rate from the upper side}
\item{Type.I.Error.Lower}{Type I error rate from the lower side}
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

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval estimation and the evaluation of
statistical methods. In L. L. Harlow, S. A. Mulaik, & J.H. Steiger (Eds.), \emph{What if there were
no significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum. }

\author{ Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu}); Keke Lai}
%\note{ }

\seealso{ \code{ss.aipe.sm}}

\examples{
# Since 'true.sm' equals 'estimated.sm', this usage
# returns the results of a correctly specified situation.
# Note that 'G' should be large (10 is used to make the 
# example run easily)
# Res.1 <- ss.aipe.sm.sensitivity(true.sm=10, estimated.sm=10, 
# desired.width=.5, assurance=.95, conf.level=.95, G=10,
# print.iter=FALSE)

# Lists contained in Res.1.
# names(Res.1) 

#Objects contained in the 'Results' lists.
# names(Res.1$Results) 

#How many obtained full widths are narrower than the desired one?
# Res.1$Summary$Pct.Width.obs.NARROWER.than.desired

# True standardized mean difference is 10, but specified at 12.
# Change 'G' to some large number (e.g., G=20)
# Res.2 <- ss.aipe.sm.sensitivity(true.sm=10, estimated.sm=12, 
# desired.width=.5, assurance=NULL, conf.level=.95, G=20)

# The effect of the misspecification on mean confidence intervals is:
# Res.2$Summary$mean.full.width

}

\keyword{design}
\keyword{htest}
