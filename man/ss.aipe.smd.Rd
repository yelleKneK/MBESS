\name{ss.aipe.smd}
\alias{ss.aipe.smd}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sample size planning for the standardized mean difference from the 
Accuracy in Parameter Estimation (AIPE) perspective
}

\description{
A function to calculate the appropriate sample size for the standardized mean difference such that 
the expected value of the confidence interval is sufficiently narrow, optionally with a 
\code{degree.of.certainty}.}

\usage{ss.aipe.smd(delta, conf.level, width, which.width="Full", 
degree.of.certainty=NULL, assurance=NULL, certainty=NULL, ...)

}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{the population value of the standardized mean difference}
  \item{conf.level}{the desired degree of confidence (i.e., 1-Type I error rate)}
  \item{width}{desired width of the specified (i.e., \code{Full}, \code{Lower}, and \code{Upper} widths) region of the confidence interval}
  \item{which.width}{the width that the \code{width} argument refers identifies the width of interest (i.e., \code{Full}, \code{Lower}, and \code{Upper} widths)}
  \item{degree.of.certainty}{parameter to ensure confidence interval width with a specified degree of certainty}
  \item{assurance}{ an alias for \code{degree.of.certainty} }
\item{certainty}{an alias for \code{degree.of.certainty}}
  \item{...}{for modifying parameters of functions this function calls upon}
}

\value{
Returns the necessary sample size \emph{per group} in order to achieve the desired 
degree of accuracy (i.e., the sufficiently narrow confidence interval).
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

Steiger, J. H., & Fouladi, R. T. (1997) Noncentrality interval estimation and the evaluation of
statistical methods. In L. L. Harlow, S. A. Mulaik, & J. H. Steiger (Eds.), \emph{What if there where
no significance tests?} (pp. 221-257). Mahwah, NJ: Lawrence Erlbaum.
}

\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})}

\note{The function \code{ss.aipe.smd} is the preferred function, and is the one that is recommended for widespread use.
The functions \code{ss.aipe.smd.lower}, \code{ss.aipe.smd.upper} and 
\code{ss.aipe.smd.full} are called from the \code{ss.aipe.smd} function.}

\section{Warning}{Finding sample size for lower and uppper confidence limits is approximate, but very close to being exact. The \code{pt()} function is limited to accurate values 
when the the noncentral parameter is less than 37.62.} 

\seealso{\code{smd}, \code{smd.c}, \code{ci.smd}, \code{ci.smd.c}, 
\code{conf.limits.nct}, \code{power.t.test}, \code{ss.aipe.smd.lower}, 
\code{ss.aipe.smd.upper}, \code{ss.aipe.smd.full}}

\examples{
# ss.aipe.smd(delta=.5, conf.level=.95, width=.30)
# ss.aipe.smd(delta=.5, conf.level=.95, width=.30, degree.of.certainty=.8)
# ss.aipe.smd(delta=.5, conf.level=.95, width=.30, degree.of.certainty=.95)
}
\keyword{design}
\keyword{htest}
