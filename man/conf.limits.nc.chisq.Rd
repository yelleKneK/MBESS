\name{conf.limits.nc.chisq}
\alias{conf.limits.nc.chisq}
\title{Confidence limits for noncentral chi square parameters}
\description{
Function to determine the noncentral parameter that leads to the observed \code{Chi.Square}-value, 
so that a confidence interval for the population noncentral chi-squrae value can be formed.
}
\usage{conf.limits.nc.chisq(Chi.Square=NULL, conf.level=.95, df=NULL, 
alpha.lower=NULL, alpha.upper=NULL, tol=1e-9, Jumping.Prop=.10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Chi.Square}{the observed chi-square value}
  \item{conf.level}{the desired degree of confidence for the interval}
  \item{df}{the degrees of freedom}
  \item{alpha.lower}{Type I error for the lower confidence limit}
  \item{alpha.upper}{Type I error for the upper confidence limit}
  \item{tol}{tolerance for iterative convergence}
  \item{Jumping.Prop}{Value used in the iterative scheme to determine the noncentral 
  parameters necessary for confidence interval construction using noncentral 
  chi square-distributions (\code{0 < Jumping.Prop < 1})}
}
\details{
If the function fails (or if a function relying upon this function fails), adjust the \code{Jumping.Prop}
(to a smaller value).
}
\value{
  \item{Lower.Limit}{Value of the distribution with \code{Lower.Limit} noncentral value that has at its specified quantile \code{Chi.Square}}
  \item{Prob.Less.Lower}{Proportion of cases falling below \code{Lower.Limit}}
  \item{Upper.Limit}{Value of the distribution with \code{Upper.Limit} noncentral value that has at its specified quantile \code{Chi.Square}}
  \item{Prob.Greater.Upper}{Proportion of cases falling above \code{Upper.Limit}}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu});
Keke Lai (University of California--Merced)
}

\seealso{\code{conf.limits.nct}, \code{conf.limits.ncf}}
\examples{
# A typical call to the function.
conf.limits.nc.chisq(Chi.Square=30, conf.level=.95, df=15)

# A one sided (upper) confidence interval.
 conf.limits.nc.chisq(Chi.Square=30, alpha.lower=0, alpha.upper=.05, 
 conf.level=NULL, df=15)
}
\keyword{design}
\keyword{multivariate}
\keyword{regression}
