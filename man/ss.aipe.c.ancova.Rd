\name{ss.aipe.c.ancova}
\alias{ss.aipe.c.ancova}
\title{ Sample size planning for a contrast in randomized ANCOVA from the Accuracy in Parameter Estimation (AIPE) perspective}
\description{
A function to calculate the appropriate sample size per group for the (unstandardized) contrast, in one-covariate randomized
ANCOVA, so that the width of the confidence interval is sufficiently narrow. 
}
\usage{
ss.aipe.c.ancova(error.var.ancova = NULL, error.var.anova = NULL, 
rho = NULL, c.weights, width, conf.level = 0.95, 
assurance = NULL, certainty = NULL)
}
\arguments{
  \item{error.var.ancova}{the population error variance of the ANCOVA model (i.e., the mean square within of the ANCOVA model)}
  \item{error.var.anova}{the population error variance of the ANOVA model (i.e., the mean square within of the ANOVA model)}
  \item{rho}{the population correlation coefficient of the response and the covariate}
  \item{c.weights}{ the contrast weights }
  \item{width}{ the desired full width of the obtained confidence interval }
  \item{conf.level}{the desired confidence interval coverage, (i.e., 1 - Type I error rate)}
  \item{assurance}{ parameter to ensure that the obtained confidence interval width is narrower than the desired width with a specified degree of certainty (must be NULL or between zero and unity) }
  \item{certainty}{an alias for \code{assurance}}
}

\details{
Either the error variance of the ANCOVA model or of the ANOVA model can be used to plan the appropriate
sample size per group. When using the error variance of the ANOVA model to plan sample size, the correlation coefficient
of the response and the covariate is also needed.}

\value{
  \item{n}{the necessary sample size \emph{per group}}
}

\references{ 
Kelley, K., Maxwell, S. E., & Rausch, J. R. (2003). Obtaining power or obtaining 
precision: Delineating methods of sample size planning. \emph{Evaluation and the Health Professions, 26}, 258-287.

Maxwell, S. E., & Delaney, H. D. (2004). \emph{Designing experiments and analyzing data: A model 
comparison perspective}. Mahwah, NJ: Erlbaum.
}

\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu}); Keke Lai <Lai.15@ND.Edu> }

%\note{}

\seealso{\code{ci.c.ancova}, \code{ci.sc.ancova}, \code{ss.aipe.c} }

\examples{
# Suppose the population error variance of some three-group ANOVA model
# is believed to be 40, and the population correlation coefficient 
# of the response and the covariate is 0.22. The researcher is 
# interested in the difference between the mean of group 1 and 
# the average of means of group 2 and 3. To plan the sample size so 
# that, with 90 percent certainty, the obtained 95 percent full 
# confidence interval width is no wider than 3:

ss.aipe.c.ancova(error.var.anova=40, rho=.22, 
c.weights=c(1, -0.5, -0.5), width=3, assurance=.90)
}

\keyword{design }
