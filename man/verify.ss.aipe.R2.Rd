\name{verify.ss.aipe.R2}
\alias{verify.ss.aipe.R2}
\title{Internal MBESS function for verifying the sample size in ss.aipe.R2}
\description{
Internal function called upon by \code{ss.aipe.R2} when \code{verify.ss=TRUE}. This function then calls 
upon \code{ss.aipe.R2.sensitivity} for the simulation study.
}
\usage{
verify.ss.aipe.R2(Population.R2 = NULL, conf.level = 0.95, width = NULL, 
Random.Predictors = TRUE, which.width = "Full", p = NULL, n = NULL, 
degree.of.certainty = NULL, g = 500, G = 10000, print.iter=FALSE, ...)
}
\arguments{
  \item{Population.R2}{value of the population multiple correlation coefficient}
  \item{conf.level}{confidence interval level (e.g., .95, .99, .90); 1-Type I error rate}
  \item{width}{width of the confidence interval (see \code{which.width})}
  \item{Random.Predictors}{whether or not the predictor variables are random (set to \code{TRUE}) or are fixed (set to \code{FALSE})}
  \item{which.width}{defines the width that \code{width} refers to}
  \item{p}{the number of predictor variables}
  \item{n}{starting sample size (i.e., from \code{ss.aipe.R2})}
  \item{degree.of.certainty}{value with which confidence can be placed that describes the likelihood of obtaining a confidence interval less than the value specified (e.g., .80, .90, .95)}
  \item{g}{simulations for the preliminary sample size (much smaller than \code{G})}
  \item{G}{number of replications for the actual Monte Carlo simulation (should be large)}
  \item{print.iter}{specify whether or not the internal iterations should be printed}
  \item{\dots}{additional arguments passed to internal functions}
}
\details{
This function is internal to MBESS and is called upon when \code{verify.ss=TRUE} in
the \code{ss.aipe.R2} function. Although users can use \code{verify.ss.aipe.R2} directly, it is not 
recommended.
}
\value{
Returns the exact (provided \code{G} is large enough) sample size necessary to satisfy the 
conditions specified.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})}
\keyword{design}
