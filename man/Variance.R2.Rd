\name{Variance.R2}
\alias{Variance.R2}
\title{Variance of squared multiple correlation coefficient}
\description{
Function to determine the variance of the squared multiple correlation coefficient given the population squared multiple correlation coefficient,
sample size, and the number of predictors.
}
\usage{
Variance.R2(Population.R2, N, p)
}
\arguments{
  \item{Population.R2}{population squared multiple correlation coefficient}
  \item{N}{sample size}
  \item{p}{the number of predictor variables}
}
\details{
Uses the hypergeometric function as discussed in and section 28 of Stuart, Ord, and Arnold (1999) in order to obtain the \emph{correct}
value for the variance of the squared multiple correlation coefficient.}
\value{
Returns the variance of the of the squared multiple correlation coefficient.
}
\references{
Stuart, A., Ord, J. K., & Arnold, S. (1999). \emph{Kendall's advanced theory of statistics: {C}lassical inference and the linear model} (Volume 2A, 2nd Edition).
New York, NY: Oxford University Press.}

\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})}
\note{
Uses package \code{gsl} and its \code{hyperg_2F1} function.
}
\seealso{\code{Expected.R2}, \code{ci.R2}, \code{ss.aipe.R2}}
\examples{
# library(gsl)
# Variance.R2(.5, 10, 5)
# Variance.R2(.5, 25, 5)
# Variance.R2(.5, 50, 5)
# Variance.R2(.5, 100, 5)
}
\keyword{design}
