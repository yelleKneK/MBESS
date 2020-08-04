\name{ss.aipe.c}
\alias{ss.aipe.c}
\title{Sample size planning for an ANOVA contrast from the Accuracy in Parameter Estimation (AIPE) perspective}

\description{
A function to calculate the appropriate sample size \emph{per group} for the (unstandardized) ANOVA contrast 
so that the width of the confidence interval is sufficiently narrow. }

\usage{
ss.aipe.c(error.variance = NULL, c.weights, width, conf.level = 0.95, 
assurance = NULL, certainty = NULL, MSwithin = NULL, SD = NULL, ...)
}
\arguments{
  \item{error.variance}{the common error variance; i.e., the mean square error}
  \item{c.weights}{ the contrast weights }
  \item{width}{ the desired full width of the obtained confidence interval }
  \item{conf.level}{the desired confidence interval coverage, (i.e., 1 - Type I error rate)}
  \item{assurance}{ parameter to ensure that the obtained confidence interval width is narrower than the desired width with a specified degree of certainty (must be NULL or between zero and unity) }
  \item{certainty}{an alias for \code{assurance}}
  \item{MSwithin}{an alias for \code{error.variance} }
  \item{SD}{the standard deviation of the common error in ANOVA model}
  \item{\dots}{ allows one to potentially include parameter values for inner functions }
}

%\details{}

\value{
    \item{n}{the necessary sample size \emph{per group}}
}

\references{
Kelley, K., Maxwell, S. E., & Rausch, J. R. (2003). Obtaining power or obtaining 
precesion: Delineating methods of sample size planning. \emph{Evaluation and the Health Professions, 26}, 258--287.

Maxwell, S. E., & Delaney, H. D. (2004). \emph{Designing experiments and analyzing data: A model 
comparison perspective}. Mahwah, NJ: Erlbaum.
}

\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu}), Keke Lai }

\note{
Be sure to use the error varaince and not its square root (i.e., the standard deviation of the errors).
}

\seealso{\code{ss.aipe.sc}, \code{ss.aipe.c.ancova}, \code{ci.c}}
\examples{
# Suppose the population error variance of some three-group ANOVA model
# is believed to be 40. The researcher is interested in the difference 
# between the mean of group 1 and the average of means of group 2 and 3. 
# To plan the sample size so that, with 90 percent certainty, the 
# obtained 95 percent full confidence interval width is no wider than 3:

ss.aipe.c(error.variance=40, c.weights=c(1, -0.5, -0.5), width=3, assurance=.90)
}

\keyword{design}
