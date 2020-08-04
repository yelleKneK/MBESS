\name{conf.limits.ncf}
\alias{conf.limits.ncf}
\title{Confidence limits for noncentral F parameters}
\description{
Function to determine the noncentral parameter that leads to the observed \emph{F}-value, 
so that a confidence interval around the population \emph{F}-value can be conducted. Used for forming confidence intervals around noncentral parameters (given the monotonic relationship between the \emph{F}-value and the noncentral value).
}
\usage{
conf.limits.ncf(F.value = NULL, conf.level = .95, df.1 = NULL, 
df.2 = NULL, alpha.lower = NULL, alpha.upper = NULL, tol = 1e-09,
Jumping.Prop = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{F.value}{the observed \emph{F}-value}
  \item{conf.level}{the desired degree of confidence for the interval}
  \item{df.1}{the numerator degrees of freedom}
  \item{df.2}{the denominator degrees of freedom}
  \item{alpha.lower}{Type I error for the lower confidence limit}
  \item{alpha.upper}{Type I error for the upper confidence limit}
  \item{tol}{tolerance for iterative convergence}
  \item{Jumping.Prop}{Value used in the iterative scheme to determine the noncentral 
  parameters necessary for confidence interval construction using noncentral 
  \emph{F}-distributions (\code{0 < Jumping.Prop < 1}) (users should not need to change this value)}
}
\details{
This function is the relied upon by the \code{ci.R2} and \code{ss.aipe.R2}. If the function fails 
(or if a function relying upon this function fails), adjust the \code{Jumping.Prop}
(to a smaller value).
}
\value{
  \item{Lower.Limit}{Value of the distribution with \code{Lower.Limit} noncentral value that has at its specified quantile \code{F.value}}
  \item{Prob.Less.Lower}{Proportion of cases falling below \code{Lower.Limit}}
  \item{Upper.Limit}{Value of the distribution with \code{Upper.Limit} noncentral value that has at its specified quantile \code{F.value}}
  \item{Prob.Greater.Upper}{Proportion of cases falling above \code{Upper.Limit}}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu}); Keke Lai (University of Califonia-Merced)}

\seealso{\code{ss.aipe.R2}, \code{ci.R2}, \code{conf.limits.nct}}
\examples{
conf.limits.ncf(F.value = 5, conf.level = .95, df.1 = 5, 
df.2 = 100)

# A one sided confidence interval.
conf.limits.ncf(F.value = 5, conf.level = NULL, df.1 = 5, 
df.2 = 100, alpha.lower = .05, alpha.upper = 0, tol = 1e-09,
Jumping.Prop = 0.1)
}
\keyword{design}
\keyword{multivariate}
\keyword{regression}
