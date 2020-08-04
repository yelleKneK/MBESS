\name{t.and.smd.conversion}
\alias{delta2lambda}
\alias{lambda2delta}
\title{Conversion functions for noncentral t-distribution}
\description{
Functions useful for converting a standardized mean difference to a noncentrality parameter, and vice versa.
}
\usage{
lambda2delta(lambda, n.1, n.2)
delta2lambda(delta, n.1, n.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda}{noncentral value from a \emph{t}-distribution}
  \item{delta}{population value of the standardized mean difference}
  \item{n.1}{sample size in group 1}
  \item{n.2}{sample size in group 2}
}
\details{
Although \code{lambda} is the population noncentral value, an estimate of it is the observed value of a
\emph{t}-statistic. Likewise, delta can be estimated as the observed standardized mean difference. Thus, the observed 
standardized mean difference can be converted to the observed \emph{t}-value. These functions are especially helpful in the
context of forming confidence intervals for the population standardized mean difference. 
}
\value{
Either the value of \code{delta} given \code{lambda} or \code{lambda} given \code{delta} (and the \emph{per group} sample sizes).
}

\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})}

\seealso{\code{smd}, \code{ci.smd}, \code{ss.aipe.smd}}
\examples{
lambda2delta(lambda=2, n.1=113, n.2=113)
delta2lambda(delta=.266076, n.1=113, n.2=113)
}
\keyword{design}
