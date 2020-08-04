\name{ss.aipe.sc.ancova}
\alias{ss.aipe.sc.ancova}
\title{Sample size planning from the AIPE perspective for standardized ANCOVA contrasts }
\description{
Sample size planning from the accuracy in parameter estimation (AIPE) perspective for standardized ANCOVA contrasts.
}
\usage{
ss.aipe.sc.ancova(Psi = NULL, sigma.anova = NULL, sigma.ancova = NULL,
psi = NULL, ratio = NULL, rho = NULL, divisor = "s.ancova", 
c.weights, width, conf.level = 0.95, assurance = NULL, ...)
}
\arguments{
  \item{Psi}{ the population unstandardized ANCOVA (adjusted) contrast }
  \item{sigma.anova}{ the population error standard deviation of the ANOVA model }
  \item{sigma.ancova}{ the population error standard deviation of the ANCOVA model  }
  \item{psi}{ the population standardized ANCOVA (adjusted) contrast }
  \item{ratio}{ the ratio of \code{sigma.ancova} over \code{sigma.anova} }
  \item{rho}{ the population correlation coefficient between the response and the covariate }
  \item{divisor}{ which error standard deviation to be used in standardizing the contrast; the value can be
  either \code{"s.ancova"} or \code{"s.anova"} }
  \item{c.weights}{ contrast weights }
  \item{width}{ the desired full width of the obtained confidence interval}
  \item{conf.level}{ the desired confidence interval coverage, (i.e., 1 - Type I error rate) }
  \item{assurance}{ parameter to ensure that the obtained confidence interval width is narrower 
  than the desired width with a specified degree of certainty (must be \code{NULL} or between zero and unity)}
  \item{\dots}{ allows one to potentially include parameter values for inner functions }
}
\details{
The sample size planning method this function is based on is developed in the context of simple (i.e., one-response-one-covariate) 
ANCOVA model and randomized design (i.e., same population covariate mean across groups).

An ANCOVA contrast can be standardized in at least two ways: (a) divided by the error standard deviation of the
ANOVA model, (b) divided by the error standard deviation of the ANCOVA model. This function can be used to analyze
both types of standardized ANCOVA contrasts.

Not all of the arguments about the effect sizes need to be specified. If \code{divisor="s.ancova"} is 
used in the argument, then input either (a) \code{psi}, or (b) \code{Psi} and \code{s.ancova}. 
If \code{divisor="s.anova"} is used in the argument, possible specifications 
are (a) \code{Psi}, \code{s.ancova}, and \code{s.anova}; (b) \code{psi}, and \code{ratio}; 
(c) \code{psi}, and \code{rho}. 
}
\value{
This function returns the sample size \emph{per group}.
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
no significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.}


\author{ Keke Lai (University of California--Merced) }

\note{ 
When \code{divisor="s.anova"} and the argument \code{assurance} is specified, the necessary
sample size \emph{per group} returned by the function with \code{assurance} specified is slightly underestimated. 
The method to obtain exact sample size in the above situation has not been developed yet. A practical solution is
to use the sample size returned as the starting value to conduct a priori Montre Carlo simulations with 
function \code{\link{ss.aipe.sc.ancova.sensitivity}}, as discussed in Lai & Kelley (under review).
 }

\seealso{ \code{ss.aipe.sc}, \code{ss.aipe.sc.ancova.sensitivity} }
\examples{
\dontrun{
ss.aipe.sc.ancova(psi=.8, width=.5, c.weights=c(.5, .5, 0, -1))

ss.aipe.sc.ancova(psi=.8, ratio=.6, width=.5, 
c.weights=c(.5, .5, 0, -1), divisor="s.anova")

ss.aipe.sc.ancova(psi=.5, rho=.4, width=.3, 
c.weights=c(.5, .5, 0, -1), divisor="s.anova")
}

}

\keyword{ design }
