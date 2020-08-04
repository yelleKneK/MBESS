\name{smd.c}
\alias{smd.c}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Standardized mean difference using the control group as the basis of 
standardization}
\description{
Function to calculate the standardized mean difference (regular or unbiased) using the control group standard deviation 
as the basis of standardization (for either raw data or summary measures).
}
\usage{
smd.c(Group.T = NULL, Group.C = NULL, Mean.T = NULL, Mean.C = NULL, 
s.C = NULL, n.C = NULL, Unbiased=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Group.T}{Raw data for the treatment group.}
  \item{Group.C}{Raw data for the control group.}
  \item{Mean.T}{The mean of the treatment group.}
  \item{Mean.C}{The mean of the control group.}
  \item{s.C}{The standard deviation of the control group (i.e., the square root of the unbiased estimator of the population variance).}
  \item{n.C}{The sample size of the control group.}
  \item{Unbiased}{Returns the unbiased estimate of the standardized mean difference using the standard deviation of the control group.}
}
\details{

When \code{Unbiased=TRUE}, the unbiased estimate of the standardized mean difference (using the control group as the
basis of standardization) is returned (Hedges, 1981). Although
the unbiased estimate of the standardized mean difference is not often reported, at least at the present time, it is
nevertheless made available to those who are interested in calculating this quantity. 

}
\value{
  Returns the estimated standardized mean difference using the control group standard deviation as the basis of standardization.
}
\references{

Hedges, L. V. (1981). Distribution theory for Glass's Estimator of effect size and related estimators. \emph{Journal of Educational Statistics, 2}, 107--128.

Glass, G. (1976). Primary, secondary, and meta-analysis of research. \emph{Educational Researcher, 5}, 3--8.  

}
\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})}

\seealso{
\code{smd}, \code{conf.limits.nct}}
\examples{
# Generate sample data.
set.seed(113)
g.T <- rnorm(n=25, mean=.5, sd=1)
g.C <- rnorm(n=25, mean=0, sd=1)
smd.c(Group.T=g.T, Group.C=g.C)

M.T <- .66745
M.C <- .24878
sd.c <- 1.1311
n.c <- 25
smd.c(Mean.T=M.T, Mean.C=M.C, s=sd.c)
smd.c(Mean.T=M.T, Mean.C=M.C, s=sd.c, n.C=n.c, Unbiased=TRUE)
}
\keyword{models}
\keyword{htest}
