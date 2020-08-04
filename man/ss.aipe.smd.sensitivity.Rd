\name{ss.aipe.smd.sensitivity}
\alias{ss.aipe.smd.sensitivity}
\title{Sensitivity analysis for sample size given the Accuracy in Parameter Estimation approach for the standardized mean difference.}
\description{
Performs sensitivity analysis for sample size determination for the standardized mean difference 
given a population and a standardized mean difference. Allows one to determine the effect of being
wrong when estimating the population standardized mean difference in terms of the 
width of the obtained (two-sided) confidence intervals. 
}
\usage{
ss.aipe.smd.sensitivity(true.delta = NULL, estimated.delta = NULL, 
desired.width = NULL, selected.n=NULL, assurance=NULL, certainty = NULL, 
conf.level = 0.95, G = 10000, print.iter = TRUE, ...)
}
\arguments{
  \item{true.delta}{population standardized mean difference}
  \item{estimated.delta}{estimated standardized mean difference; can be \code{true.delta} to perform standard simulations}
  \item{desired.width}{describe full width for the confidence interval around the population standardized mean difference}
  \item{selected.n}{selected sample size to use in order to determine distributional properties of at a given value of sample size}
\item{assurance}{parameter to ensure confidence interval width with a specified degree of certainty (must 
be \code{NULL} or between zero and unity)}
  \item{certainty}{an alias for \code{assurance}}
  \item{conf.level}{the desired degree of confidence (i.e., 1-Type I error rate).}
  \item{G}{number of generations (i.e., replications) of the simulation}
  \item{print.iter}{to print the current value of the iterations}
  \item{...}{for modifying parameters of functions this function calls}
}
\details{
For sensitivity analysis when planning sample size given the desire to obtain narrow confidence intervals
for the population standardized mean difference. Given a population value and an estimated value, one can determine
the effects of incorrectly specifying the population standardized mean difference (\code{true.delta}) on the 
obtained widths of the confidence intervals. Also, one can evaluate the percent of the confidence intervals
that are less than the desired width (especially when modifying the \code{certainty} parameter); see \code{ss.aipe.smd})

Alternatively, one can specify \code{selected.n} to determine the results at a particular sample size (when doing this \code{estimated.delta} cannot be specified).
}
\value{
  \item{Results}{list of the results in \code{G}-length vector form}
  \item{Specifications}{specification of the function}
  \item{Summary}{summary measures of some important descriptive statistics}
\item{d}{contained in \code{Results} list: vector of the observed \code{d} values}
\item{Full.Width}{contained in \code{Results} list: vector of }
\item{Width.from.d.Upper}{contained in \code{Results} list: vector of the observed upper widths of the confidence interval (upper limit minus observed standardized mean difference)}
\item{Width.from.d.Lower}{contained in \code{Results} list: vector of the observed lower widths of the confidence interval (standardized mean difference minus lower limit)}
\item{Type.I.Error.Upper}{contained in \code{Results} list: iterations where a Type I error occurred on the upper end of the confidence interval}
\item{Type.I.Error.Lower}{contained in \code{Results} list: iterations where a Type I error occurred on the lower end of the confidence interval}
\item{Type.I.Error}{contained in \code{Results} list: iterations where a Type I error occurred}
\item{Upper.Limit}{contained in \code{Results} list: vector of the obtained upper limits from the simulation}
\item{Low.Limit}{contained in \code{Results} list: vector of the obtained lower limits from the simulation}
\item{replications}{contained in \code{Specifications} list: number of generations (i.e., replication) of the simulation}
\item{true.delta}{contained in \code{Specifications} list: population value of the standardized mean difference}
\item{estimated.delta}{contained in \code{Specifications} list: value of the population (mis)specified for purposes of sample size planning}
\item{desired.width}{contained in \code{Specifications} list: desired full width of the confidence interval around the population standardized mean difference}
\item{certainty}{contained in \code{Specifications} list: desired degree of certainty that the obtained confidence interval width is less than the value specified}
\item{n.j}{contained in \code{Specifications} list: sample size per group given the specifications}
\item{mean.full.width}{contained in \code{Summary} list: mean width of the obtained confidence intervals}
\item{median.full.width}{contained in \code{Summary} list: median width of the obtained confidence intervals}
\item{sd.full.width}{contained in \code{Summary} list: standard deviation of the obtained confidence intervals}
\item{Pct.Less.Desired}{contained in \code{Summary} list: Percent of the confidence widths less than the width specified.}
\item{mean.Width.from.d.Lower}{contained in \code{Summary} list:mean width of the lower portion of the confidence interval (from d)}
\item{mean.Width.from.d.Upper}{contained in \code{Summary} list:mean width of the upper portion of the confidence interval (from d)}
\item{Type.I.Error.Upper}{contained in \code{Summary} list: Type I error rate from the upper side}
\item{Type.I.Error.Lower}{contained in \code{Summary} list: Type I error rate from the lower side}

}
\references{
Cumming, G. & Finch, S. (2001). A primer on the understanding, use, and calculation of confidence intervals that are
based on central and noncentral distributions, \emph{Educational and Psychological Measurement, 61}, 532--574.

Hedges, L. V. (1981). Distribution theory for Glass's Estimator of effect size and related estimators. \emph{Journal of Educational Statistics, 2}, 107--128.

Kelley, K. (2005). The effects of nonnormal distributions on confidence intervals around the standardized mean
difference: Bootstrap and parametric confidence intervals, \emph{Educational and Psychological Measurement, 65}, 51--69.

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval estimation and the evaluation of
statistical methods. In L. L. Harlow, S. A. Mulaik, & J.H. Steiger (Eds.), \emph{What if there were
no significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.
}
\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})}
\note{Returns three lists, where each list has multiple components.}
\seealso{\code{ss.aipe.smd}}
\examples{
# Since 'true.delta' equals 'estimated.delta', this usage 
# returns the results of a correctly specified situation.
# Note that 'G' should be large (50 is used to make the example run easily)
# Res.1 <- ss.aipe.smd.sensitivity(true.delta=.5, estimated.delta=.5, 
# desired.width=.30, certainty=NULL, conf.level=.95, G=50,
# print.iter=FALSE)

# Lists contained in Res.1.
# names(Res.1) 

#Objects contained in the 'Results' lists.
# names(Res.1$Results) 

#Extract d from the Results list of Res.1.
# d <- Res.1$Results$d 

# hist(d)

# Pull out summary measures
# Res.1$Summary

# True standardized mean difference is .4, but specified at .5.
# Change 'G' to some large number (e.g., G=5,000)
# Res.2 <- ss.aipe.smd.sensitivity(true.delta=.4, estimated.delta=.5, 
# desired.width=.30, certainty=NULL, conf.level=.95, G=50, 
# print.iter=FALSE)

# The effect of the misspecification on mean confidence intervals is:
# Res.2$Summary$mean.full.width

# True standardized mean difference is .5, but specified at .4.
# Res.3 <- ss.aipe.smd.sensitivity(true.delta=.5, estimated.delta=.4, 
# desired.width=.30, certainty=NULL, conf.level=.95, G=50, 
# print.iter=FALSE)

# The effect of the misspecification on mean confidence intervals is:
# Res.3$Summary$mean.full.width
}
\keyword{design}
\keyword{htest}
