\name{aipe.smd}
\alias{ss.aipe.smd.lower}
\alias{ss.aipe.smd.upper}
\alias{ss.aipe.smd.full}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sample size planning for the standardized mean different from the accuracy 
in parameter estimation approach}
\description{
A set of functions that \code{ss.aipe.smd} calls upon to calculate the appropriate sample size 
for the standardized mean difference such that the expected value of the confidence interval 
is sufficiently narrow.
}


\usage{
ss.aipe.smd.full(delta, conf.level, width, ...)
ss.aipe.smd.lower(delta, conf.level, width, ...)
ss.aipe.smd.upper(delta, conf.level, width, ...)
}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{the population value of the standardized mean difference}
  \item{conf.level}{the desired degree of confidence (i.e., 1-Type I error rate)}
  \item{width}{desired width of the specified (i.e., \code{Lower}, \code{Upper}, \code{Full}) region of the confidence interval}
  \item{...}{specify additional parameters in functions these functions call upon}
  }

\value{
  \item{n}{The necessary sample size \emph{per group} in order to satisfy the specified goals.}
}

\references{
Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.

Cumming, G. & Finch, S. (2001). A primer on the understanding, use, and calculation of confidence intervals that are
based on central and noncentral distributions, \emph{Educational and Psychological Measurement, 61}, 532--574.

Hedges, L. V. (1981). Distribution theory for Glass's Estimator of effect size and related estimators. \emph{Journal of Educational Statistics, 2}, 107--128.

Kelley, K. (2005). The effects of nonnormal distributions on confidence intervals around the standardized mean
difference: Bootstrap and parametric confidence intervals, \emph{Educational and Psychological Measurement, 65}, 51--69.

Kelley, K., Maxwell, S. E., & Rausch, J. R. (2003). Obtaining Power or Obtaining Precision: Delineating Methods
of Sample-Size Planning, \emph{Evaluation and the Health Professions, 26}, 258--287.

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the standardized mean difference: 
Accuracy in Parameter Estimation via narrow confidence intervals. \emph{Psychological Methods, 11(4)}, 363--385.

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval estimation and the evaluation of
statistical methods. In L. L. Harlow, S. A. Mulaik, & J.H. Steiger (Eds.), \emph{What if there were
no significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum.
}

\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})}

\note{The function \code{ss.aipe.smd} is the function users should generally use. The function
 \code{ss.aipe.smd} calls upon these functions as needed. They can be thought of loosely 
 as internal MBESS functions.}

\section{Warning}{The returned value is the sample size \emph{per group}. Currently only
\code{ss.aipe.smd.full} returns the exact value. However, \code{ss.aipe.smd.lower} and \code{ss.aipe.smd.upper}
provide approximate sample size values.} 

\seealso{\code{ss.aipe.smd}}

\keyword{design}
\keyword{htest}
