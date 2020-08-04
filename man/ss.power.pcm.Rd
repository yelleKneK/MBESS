\name{ss.power.pcm}
\alias{ss.power.pcm}

\title{Sample size planning for power for polynomial change models}
\description{ Returns power given the sample size, or sample size given the desired power, for polynomial change models  (currently only linear, that is, straight-line, change models)}
\usage{ss.power.pcm(beta, tau, level.1.variance, frequency, duration, desired.power = NULL, 
N = NULL, alpha.level = 0.05, standardized = TRUE, directional = FALSE)
}

\arguments{
  \item{beta}{ the level two regression coefficient for the group by time (linear) interaction; 
  where "X" is coded -.5 and .5 for the two groups. }
  \item{tau}{ the true variance of the individuals' slopes }
  \item{level.1.variance}{ level one variance }
  \item{frequency}{ frequency of measurements per unit of time duration of the study in the
  particular units (e.g., age, hours, grade level, years, etc.)  }
  \item{duration}{ time in some number of units (e.g., years) }
  \item{desired.power}{ desired power }
  \item{N}{ total sample size (one-half in each of the two groups) }
  \item{alpha.level}{ Type I error rate }
  \item{standardized}{ the standardized slope is the unstandardized slope divided by the square root of tau, the variance of the unique effects for beta.  }
  \item{directional}{ should a one (\code{TRUE}) or two (\code{FALSE}) tailed test be performed.}
}
%\details{}
%\value{}

\references{ Raudenbush, S. W., & X-F., Liu. (2001). Effects of study duration, frequency of observation, 
and sample size on power in studies of group differences in polynomial change. \emph{Psychological Methods, 6}, 387--401.  }

\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})}

%\note{}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }

\examples{
# Example from Raudenbush and Liu (2001)
ss.power.pcm(beta=-.4, tau=.003, level.1.variance=.0262, frequency=2, duration=2, 
desired.power=.80, alpha.level=.05, standardized=TRUE, directional=FALSE)
ss.power.pcm(beta=-.4, tau=.003, level.1.variance=.0262, frequency=2, duration=2,
N=238, alpha.level=.05, standardized=TRUE, directional=FALSE)


# The standardized effect size is obtained as beta/sqrt(tau): -.4/sqrt(.003) = -.0219.
# ss.power.pcm(beta=-.0219, tau=.003, level.1.variance=.0262, frequency=2, duration=2, 
# desired.power=.80, alpha.level=.05, standardized=FALSE, directional=FALSE)
ss.power.pcm(beta=-.0219, tau=.003, level.1.variance=.0262, frequency=2, duration=2, 
N=238, alpha.level=.05, standardized=FALSE, directional=FALSE)

}

\keyword{ design }

