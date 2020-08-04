\name{ss.aipe.sc.ancova.sensitivity}
\alias{ss.aipe.sc.ancova.sensitivity}
\title{Sensitivity analysis for the sample size planning method for standardized ANCOVA contrast }
\description{
Sensitivity analysis for the sample size planning method with the goal to obtain sufficiently narrow confidence intervals for standardized 
ANCOVA complex contrasts.
}
\usage{
ss.aipe.sc.ancova.sensitivity(true.psi = NULL, estimated.psi = NULL, 
c.weights, desired.width = NULL, selected.n = NULL, mu.x = 0, 
sigma.x = 1, rho, divisor = "s.ancova", assurance = NULL, 
conf.level = 0.95, G = 10000, print.iter = TRUE, detail = TRUE, ...)
}
\arguments{
  \item{true.psi}{ the population standardized ANCOVA contrast }
  \item{estimated.psi}{ the estimated standardized ANCOVA contrast }
  \item{c.weights}{ the contrast weights}
  \item{desired.width}{ the desired full width of the obtained confidence interval }
  \item{selected.n}{ selected sample size to use in order to determine distributional properties of a given value of sample size}
  \item{mu.x}{ the population mean for the covariate }
  \item{sigma.x}{ the population standard deviation of the covariate }
  \item{rho}{ the population correlation coefficient between the response and the covariate }
  \item{divisor}{ which error standard deviation to be used in standardizing the contrast; the value can be
  either \code{"s.ancova"} or \code{"s.anova"} }
  \item{assurance}{ parameter to ensure that the obtained confidence interval width is narrower than the 
  desired width with a specified degree of certainty (must be \code{NULL} or between zero and unity) }
  \item{conf.level}{ the desired confidence interval coverage, (i.e., 1 - Type I error rate) }
  \item{G}{ number of generations (i.e., replications) of the simulation }
  \item{print.iter}{ to print the current value of the iterations }
  \item{detail}{ whether the user needs a detailed (\code{TRUE}) or brief (\code{FALSE}) report of the simulation results; the 
    detail report includes all the raw data in the simulations }
  \item{\dots}{ allows one to potentially include parameter values for inner functions }
}
\details{
The sample size planning method this function is based on is developed in the context of simple (i.e., one-response-one-covariate) 
ANCOVA model and randomized design (i.e., same population covariate mean across groups).

An ANCOVA contrast can be standardized in at least two ways: (a) divided by the error standard deviation of the
ANOVA model, (b) divided by the error standard deviation of the ANCOVA model. This function can be used to analyze
both types of standardized ANCOVA contrasts.

The population mean and standard deviation of the covariate does not affect the sample size planning 
procedure; they can be specified as any values that are considered as reasonable by the user.
}
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
\item{Sample.Size.per.Group}{sample size \emph{per group}}
\item{Number.of.Groups}{number of groups}
\item{mean.full.width}{mean width of the obtained full confidence intervals}
\item{median.full.width}{median width of the obtained full confidence intervals}
\item{sd.full.width}{standard deviation of the widths of the obtained full confidence intervals}
\item{Pct.Width.obs.NARROWER.than.desired}{percentage of the obtained full confidence interval widths that are narrower than the desired width}
\item{mean.Width.from.psi.obs.Lower}{mean lower width of the obtained confidence intervals}
\item{mean.Width.from.psi.obs.Upper}{mean upper width of the obtained confidence intervals}
\item{Type.I.Error.Upper}{Type I error rate from the upper side}
\item{Type.I.Error.Lower}{Type I error rate from the lower side}
\item{Type.I.Error}{Type I error rate}
}
\references{ 
Kelley, K. (2007). Constructing confidence intervals for standardized effect sizes: Theory, application, 
and implementation. \emph{Journal of Statistical Software, 20} (8), 1--24.

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the standardized mean difference: 
Accuracy in Parameter Estimation via narrow confidence intervals. \emph{Psychological Methods, 11} (4), 363--385.

Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for ANCOVA and ANOVA contrasts: Sample size planning via narrow confidence intervals. 
\emph{British Journal of Mathematical and Statistical Psychology, 65}, 350--370.

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval estimation and the evaluation of
statistical methods. In L. L. Harlow, S. A. Mulaik, & J.H. Steiger (Eds.), \emph{What if there were
no significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum. }


\author{ Keke Lai }

%\note{ }

\seealso{ \code{\link{ss.aipe.sc.ancova}}; \code{\link{ss.aipe.sc.sensitivity}} }
%\examples{}

\keyword{ design }
