\name{ss.aipe.sc.sensitivity}
\alias{ss.aipe.sc.sensitivity}
\title{Sensitivity analysis for sample size planning for the standardized ANOVA contrast from 
the Accuracy in Parameter Estimation (AIPE) Perspective}
\description{
Performs a sensitivity analysis when planning sample size from the Accuracy in Parameter Estimation (AIPE) 
Perspective for the standardized ANOVA contrast.
}
\usage{
ss.aipe.sc.sensitivity(true.psi = NULL, estimated.psi = NULL, c.weights, 
desired.width = NULL, selected.n = NULL, assurance = NULL, certainty=NULL, 
conf.level = 0.95, G = 10000, print.iter = TRUE, detail = TRUE, ...)
}
\arguments{
  \item{true.psi}{population standardized contrast }
  \item{estimated.psi}{ estimated standardized contrast}
  \item{c.weights}{ the contrast weights }
  \item{desired.width}{ the desired full width of the obtained confidence interval }
  \item{selected.n}{ selected sample size to use in order to determine distributional properties of at a given value of sample size }
  \item{assurance}{ parameter to ensure that the obtained confidence interval width is narrower than the desired width with a specified degree of certainty (must be NULL or between zero and unity) }
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
\item{psi.obs}{observed standardized contrast in each iteration}
\item{Full.Width}{vector of the full confidence interval width}
\item{Width.from.psi.obs.Lower}{vector of the lower confidence interval width}
\item{Width.from.psi.obs.Upper}{vector of the upper confidence interval width}
\item{Type.I.Error.Upper}{iterations where a Type I error occurred on the upper end of the confidence interval}
\item{Type.I.Error.Lower}{iterations where a Type I error occurred on the lower end of the confidence interval}
\item{Type.I.Error}{iterations where a Type I error happens}
\item{Lower.Limit}{the lower limit of the obtained confidence interval}
\item{Upper.Limit}{the upper limit of the obtained confidence interval}
\item{replications}{number of replications of the simulation}
\item{True.psi}{population standardized contrast}
\item{Estimated.psi}{estimated standardized contrast}
\item{Desired.Width}{the desired full width of the obtained confidence interval}
\item{assurance}{the value assigned to the argument \code{assurance}}
\item{Sample.Size.per.Group}{sample size per group}
\item{Number.of.Groups}{number of groups}
\item{mean.full.width}{mean width of the obtained full conficence intervals}
\item{median.full.width}{median width of the obtained full confidence intervals}
\item{sd.full.width}{standard deviation of the widths of the obtained full confidence intervals}
\item{Pct.Width.obs.NARROWER.than.desired}{percentage of the obtained full confidence interval widths that are narrower than the desired width}
\item{mean.Width.from.psi.obs.Lower}{mean lower width of the obtained confidence intervals}
\item{mean.Width.from.psi.obs.Upper}{mean upper width of the obtained confidence intervals}
\item{Type.I.Error.Upper}{Type I error rate from the upper side}
\item{Type.I.Error.Lower}{Type I error rate from the lower side}

}

\references{ 
Cumming, G. & Finch, S. (2001). A primer on the understanding, use, and calculation of confidence intervals that are
based on central and noncentral distributions, \emph{Educational and Psychological Measurement, 61}, 532--574.

Hedges, L. V. (1981). Distribution theory for Glass's Estimator of effect size and related estimators. \emph{Journal of Educational Statistics, 2}, 107--128.

Kelley, K. (2007). Constructing confidence intervals for standardized effect sizes: Theory, application, 
and implementation. \emph{Journal of Statistical Software, 20} (8), 1--24.

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the standardized mean difference: 
Accuracy in Parameter Estimation via narrow confidence intervals. P\emph{sychological Methods, 11} (4), 363--385.

Lai, K., & Kelley, K. (2007). Sample size planning for standardized ANCOVA and ANOVA 
contrasts: Obtaining narrow confidence intervals. \emph{Manuscript submitted for publication}.

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval estimation and the evaluation of
statistical methods. In L. L. Harlow, S. A. Mulaik, & J.H. Steiger (Eds.), \emph{What if there where
no significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum. }


\author{ Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu}); Keke Lai (University of California -- Merced) }

%\note{}

\seealso{\code{ss.aipe.sc}, \code{ss.aipe.c}, \code{conf.limits.nct}}

%\examples{}

\keyword{design}
