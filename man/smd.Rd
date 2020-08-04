\name{smd}
\alias{smd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Standardized mean difference}
\description{
Function to calculate the standardized mean difference (regular or unbiased) using either raw data or summary measures.
}
\usage{
smd(Group.1 = NULL, Group.2 = NULL, Mean.1 = NULL, Mean.2 = NULL, 
s.1 = NULL, s.2 = NULL, s = NULL, n.1 = NULL, n.2 = NULL,
Unbiased=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Group.1}{Raw data for group 1.}
  \item{Group.2}{Raw data for group 2.}
  \item{Mean.1}{The mean of group 1.}
  \item{Mean.2}{The mean of group 2.}
  \item{s.1}{The standard deviation of group 1 (i.e., the square root of the unbiased estimator of the population variance).}
  \item{s.2}{The standard deviation of group 2 (i.e., the square root of the unbiased estimator of the population variance).}
  \item{s}{The pooled group standard deviation (i.e., the square root of the unbiased estimator of the population variance).}
  \item{n.1}{The sample size within group 1.}
  \item{n.2}{The sample size within group 2.}
  \item{Unbiased}{Returns the unbiased estimate of the standardized mean difference.}
}
\details{

When \code{Unbiased=TRUE}, the unbiased estimate of the standardized mean difference is returned (Hedges, 1981).

}
\value{
  Returns the estimated standardized mean difference.
}
\references{

Cohen, J. (1988). \emph{Statistical power analysis for the behavioral sciences} (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.

Cumming, G. & Finch, S. (2001). A primer on the understanding, use, and calculation of confidence intervals that are
based on central and noncentral distributions, \emph{Educational and Psychological Measurement, 61}, 532--574.

Hedges, L. V. (1981). Distribution theory for Glass's Estimator of effect size and related estimators. \emph{Journal of Educational Statistics, 2}, 107--128.

Kelley, K. (2005) The effects of nonnormal distributions on confidence intervals around the standardized mean
difference: Bootstrap and parametric confidence intervals, \emph{Educational and Psychological Measurement, 65}, 51--69.

Kelley, K. (2007). Constructing confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20} (8), 1--24.

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval estimation and the evaluation of
statistical methods. In L. L. Harlow, S. A. Mulaik, & J. H. Steiger (Eds.), \emph{What if there were
no significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.

}
\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})}

\seealso{
\code{smd.c}, \code{conf.limits.nct}, \code{ss.aipe}
}
\examples{
# Generate sample data.
set.seed(113)
g.1 <- rnorm(n=25, mean=.5, sd=1)
g.2 <- rnorm(n=25, mean=0, sd=1)
smd(Group.1=g.1, Group.2=g.2)

M.x <- .66745
M.y <- .24878
sd <- 1.048
smd(Mean.1=M.x, Mean.2=M.y, s=sd)

M.x <- .66745
M.y <- .24878
n1 <- 25
n2 <- 25
sd.1 <- .95817
sd.2 <- 1.1311
smd(Mean.1=M.x, Mean.2=M.y, s.1=sd.1, s.2=sd.2, n.1=n1, n.2=n2)

smd(Mean.1=M.x, Mean.2=M.y, s.1=sd.1, s.2=sd.2, n.1=n1, n.2=n2, 
Unbiased=TRUE)

}
\keyword{models}
\keyword{htest}
