\name{ss.aipe.c.ancova.sensitivity}
\alias{ss.aipe.c.ancova.sensitivity}
\title{Sensitivity analysis for sample size planning for the (unstandardized) contrast in randomized ANCOVA
 from the Accuracy in Parameter Estimation (AIPE) Perspective }
\description{Performs a sensitivity analysis when planning sample size from 
the Accuracy in Parameter Estimation (AIPE) Perspective for the (unstandardized) contrast in randomized ANCOVA design.
}
\usage{
ss.aipe.c.ancova.sensitivity(true.error.var.ancova = NULL, 
est.error.var.ancova = NULL, true.error.var.anova = NULL, 
est.error.var.anova = NULL, rho, est.rho = NULL, G = 10000, 
mu.y, sigma.y, mu.x, sigma.x, c.weights, width,
conf.level = 0.95, assurance = NULL, certainty=NULL)
}
\arguments{
  \item{true.error.var.ancova}{population error variance of the ANCOVA model }
  \item{est.error.var.ancova}{estimated error variance of the ANCOVA model }
  \item{true.error.var.anova}{population error variance of the ANOVA model (i.e., excluding the covariate) }
  \item{est.error.var.anova}{estimated error variance of the ANOVA model (i.e., excluding the covariate) }
  \item{rho}{population correlation coefficient of the response and the covariate }
  \item{est.rho}{ estimated correlation coefficient of the response and the covariate }
  \item{G}{ number of generations (i.e., replications) of the simulation}
  \item{mu.y}{vector that contains the response's population mean of each group}
  \item{sigma.y}{the population standard deviation of the response}
  \item{mu.x}{the population mean of the covariate}
  \item{sigma.x}{the population standard deviation of the covariate }
  \item{c.weights}{ the contrast weights }
  \item{width}{ the desired full width of the obtained confidence interval }
  \item{conf.level}{the desired confidence interval coverage, (i.e., 1 - Type I error rate)}
  \item{assurance}{ parameter to ensure that the obtained confidence interval width is narrower than the desired width with a specified degree of certainty (must be NULL or between zero and unity) }
  \item{certainty}{an alias for \code{assurance}}
}

\details{
The arguments \code{mu.y}, \code{mu.x}, \code{sigma.y}, and \code{sigma.x} are used to generate random data in the simulations
for the sensitivity analysis. The value of \code{mu.y} should be the same as the square root of \code{true.error.var.anova}

So far this function is based on one-covariate randomized ANCOVA design only. The argument \code{mu.x} should be
a single number, because it is assumed that the population mean of the covariate is equal across groups in randomized
ANCOVA. 
}
\value{
\item{Psi.obs}{the observed (unstandardized) contrast}
\item{se.Psi}{the standard error of the observed (unstandardized) contrast}
\item{se.Psi.restricted}{the standard error of the observed (unstandardized) contrast calculated by ignoring the covariate}
\item{se.res.over.se.full}{the ratio of contrast's full standard error over the restricted one in each iteration}
\item{width.obs}{full confidence interval width}
\item{Type.I.Error}{Type I error happens in each iteration}
\item{Type.I.Error.Upper}{Type I error happens in the upper end in each iteration}
\item{Type.I.Error.Lower}{Type I error happens in the lower end in each iteration}
\item{Type.I.Error}{percentage of Type I error happened in the entire simulation}
\item{Type.I.Error.Upper}{percentage of Type I error happened in the upper end in the entire simulation}
\item{Type.I.Error.Lower}{percentage of Type I error happened in the lower end in the entire simulation}
\item{width.NARROWER.than.desired}{percentage of obtained widths that are narrower than the desired width}
\item{Mean.width.obs}{mean width of the obtained full confidence intervals}
\item{Median.width.obs}{median width of the obtained full confidence intervals}
\item{Mean.se.res.vs.se.full}{the mean of the ratios of contrast's full standard error over the restricted one}
\item{Psi.pop}{population (unstandardized) contrast}
\item{Contrast.Weights}{contrast weights}
\item{mu.y}{the response's population mean of each group}
\item{mu.x}{the population mean of the covariate}
\item{sigma.x}{the population standard deviation of the covariate}
\item{Sample.Size.per.Group}{sample size per group}
\item{conf.level}{the desired confidence interval coverage, (i.e., 1 - Type I error rate)}
\item{assurance}{specified \code{assurance}}
\item{rho}{population correlation coefficient of the response and the covariate}
\item{est.rho}{estimated correlation coefficient of the response and the covariate}
\item{true.error.var.ANOVA}{population error variance of the ANOVA model}
\item{est.error.var.ANOVA}{estimated error variance of the ANOVA model}
}

%\references{}

\author{Keke Lai (University of Notre Dame; \email{Lai.15@ND.Edu})}

%\note{ }

%\seealso{ }

\examples{
\dontrun{
ss.aipe.c.ancova.sensitivity(true.error.var.ancova=30, 
est.error.var.ancova=30, rho=.2, mu.y=c(10,12,15,13), mu.x=2, 
G=1000, sigma.x=1.3, sigma.y=2, c.weights=c(1,0,-1,0), width=3)

ss.aipe.c.ancova.sensitivity(true.error.var.anova=36, 
est.error.var.anova=36, rho=.2, est.rho=.2, G=1000, 
mu.y=c(10,12,15,13), mu.x=2, sigma.x=1.3, sigma.y=6, 
c.weights=c(1,0,-1,0), width=3, assurance=NULL)
}
}

\keyword{design}
