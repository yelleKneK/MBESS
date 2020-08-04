\name{ss.aipe.pcm}
\alias{ss.aipe.pcm}

\title{Sample size planning for polynomial change models in longitudinal study

}
\description{
This function plans sample size with respect to the group-by-time interaction in the context of a longitudinal design with two groups. It plans sample size from the accuracy in parameter estimation (AIPE) perspective, where the goal is to obtain a sufficiently narrow confidence interval for the fixed effect polynomial
change coefficient parameter (e.g., linear, quadratic, etc.). The sample size returned can be one such that (a) the expected confidence interval width is sufficiently narrow, or (b) the observed confidence interval will be sufficiently narrow with a specified high degree of assurance (e.g., .99, .95, .90, etc.). This function accompanies Kelley and Rausch (2011).
}
\usage{
ss.aipe.pcm(true.variance.trend, error.variance, 
variance.true.minus.estimated.trend = NULL, duration, frequency, 
width, conf.level = 0.95, trend = "linear", assurance = NULL)
}

\arguments{
  \item{true.variance.trend}{The variance of the individuals' true change coefficients (i.e., \eqn{\sigma^2_{\upsilon_m}} in Kelley & Rausch, 2011) for the polynomial trend (e.g., linear, quadratic, etc.) of interest.}
  \item{error.variance}{The true error variance (i.e., \eqn{\sigma^2_{\epsilon}} in Kelley & Rausch, 2011).}
  \item{variance.true.minus.estimated.trend}{The variance of the difference between the \eqn{m}th true change coefficient minus the \eqn{m}th estimated change coefficient (i.e., \eqn{\sigma^2_{\hat{\pi}_{m} - \pi_{m}}} from Equation 19 in Kelley & Rausch, 2011).}
  \item{duration}{The duration of the study.}
  \item{frequency}{The number of times measurement occurs within each unit of time. }
  \item{width}{width of the confidence interval}
  \item{conf.level}{The desired level of confidence for the confidence interval that will be computed at the completion of the study. }
  \item{trend}{The polynomial trend (1st-3rd) of interest specified as "linear", "quadratic", or "cubic".}
  \item{assurance}{Value with which confidence can be placed that describes the likelihood of obtaining a confidence interval less than the value specified (e.g, .80, .90, .95)}
}

%\details{ }

\value{
Returns the necessary sample size for the combination of the desired goals and values of the population parameters for a specific design. 
}
\references{Kelley, K., & Rausch, J. R. (2011). Accuracy in parameter estimation for polynomial change models. \emph{Psychological Methods}.}

\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})}

\note{Like in all formal sample size planning methods that require the value of one or more population parameter(s), if the population parameters are incorrectly specified, there is no guarantee that the sample size this function returns will be accurate. Of course, the further away from the true values, the further away the true sample size will tend to be.

The number of timepoints in a study (say \eqn{M}) is defined by \eqn{f \times D + 1}, where \eqn{f} is the frequency and \eqn{D} is the duration.}

%\seealso{}

\examples{ 
\dontrun{
# An example used in Kelley and Rausch for the expected confidence interval 
# width (returns 278). Thus, a necessary sample size of 278 is required when 
# the duration of the study will be 4 units and the frequency of measurement 
# occasions is 1 year in order for the expected confidence interval 
# width to be 0.025 units.  

ss.aipe.pcm(true.variance.trend=0.003, error.variance=0.0262, duration=4, 
frequency=1, width=0.025, conf.level=.95)

# Now, when incorporating an assurance parameter (returns 316). 
# Thus, a necessary sample size of 316 will ensure that the 95\% confidence 
# interval will be sufficiently narrow (i.e., have a width less than .025 units) 
# at least 99\% of the time.

ss.aipe.pcm(true.variance.trend=.003, error.variance=.0262, duration=4, 
frequency=1, width=.025, conf.level=.95, assurance=.99)
}
}

\keyword{design}
\keyword{multivariate}
