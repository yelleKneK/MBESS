\name{Expected.R2}
\alias{Expected.R2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Expected value of the squared multiple correlation coefficient}
\description{
Returns the expected value of the squared multiple correlation coefficient given the population squared multiple correlation coefficient, sample size, and the number of predictors}
\usage{
Expected.R2(Population.R2, N, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Population.R2}{population squared multiple correlation coefficient}
  \item{N}{sample size}
  \item{p}{the number of predictor variables}
}
\details{
Uses the hypergeometric function as discussed in section 28 of Stuart, Ord, and Arnold (1999) in order to obtain the \emph{correct}
value for the squared multiple correlation coefficient. Many times an exact value is given that ignores the hypergeometric function. 
This function yields the correct value. 
}
\value{
Returns the expected value of the squared multiple correlation coefficient.
}
\references{
Olkin, I. & Pratt, J. W. (1958). Unbiased estimation of certain correlation coefficients. \emph{Annals of Mathematical statistics, 29}, 201--211.

Stuart, A., Ord, J. K., & Arnold, S. (1999). \emph{Kendall's advanced theory of statistics: {C}lassical inference and the linear model} (Volume 2A, 2nd Edition).
New York, NY: Oxford University Press. 
}
\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})}

\note{
Uses package \code{gsl} and its \code{hyperg_2F1} function.
}
\seealso{\code{ss.aipe.R2}, \code{ci.R2}, \code{Variance.R2}}
\examples{
# library(gsl)
# Expected.R2(.5, 10, 5)
# Expected.R2(.5, 25, 5)
# Expected.R2(.5, 50, 5)
# Expected.R2(.5, 100, 5)
# Expected.R2(.5, 1000, 5)
# Expected.R2(.5, 10000, 5)
}
\keyword{design}% at least one, from doc/KEYWORDS
