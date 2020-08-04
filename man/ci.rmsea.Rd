\name{ci.rmsea}
\alias{ci.rmsea}
\title{Confidence interval for the population root mean square error of approximation}
\description{
Confidence interval for the population root mean square error of approximation (RMSEA).
}

\usage{
ci.rmsea(rmsea, df, N, conf.level = 0.95, alpha.lower = NULL, 
alpha.upper = NULL)
}
\arguments{
  \item{rmsea}{observed root mean square error of approximation}  
  \item{df}{degrees of freedom of the model}
  \item{N}{sample size}
  \item{conf.level}{desired confidence level (e.g., .90, .95, .99)}
  \item{alpha.lower}{ the Type I error rate for the lower tail }
  \item{alpha.upper}{ the Type I error rate for the upper tail }
}
\details{
Provides a confidence interval for the population root mean square error of approximation (RMSEA) using the noncentral chi-square distribution (e.g., Steiger & Lind, 1980).
}
\value{
returns the upper and lower limit as well as the observed value of the RMSEA.
}
\references{ 
Steiger, J. H., & Lind, J. C. (1980). \emph{Statistically-based tests for the number of common 
factors}. Paper presented at the annual Spring meeting of the Psychometric Society, Iowa City, IA.
}
\author{ Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu}) }
%\note{ }
%\seealso{  }
%\examples{}

\keyword{ design }
