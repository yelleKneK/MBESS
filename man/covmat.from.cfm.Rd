\name{covmat.from.cfm}
\alias{covmat.from.cfm}

\title{ Covariance matrix from confirmatory (single) factor model.}
\description{
  Function calculates a covariance matrix using the specified \code{Lambda} and \code{Psi.Square} values from a confirmatory
factor model approach (McDonald, 1999).
}
\usage{
covmat.from.cfm(Lambda, Psi.Square, tol.det = 1e-05)
}

\arguments{
  \item{Lambda}{the vector of population factor loadings}
  \item{Psi.Square}{the vector of population error variances }
  \item{tol.det}{the specified tolerance for the determinant }
}

%\details{}

\value{
  \item{Population.Covariance }{the population covariance matrix}
  \item{True.Covariance}{the true covariance matrix}
  \item{True.Covariance}{the error covariance matrix}

}
\references{McDonald, R. P. (1999). \emph{Test theory: A unified approach}. Mahwah, NJ: Erlbaum.}

\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu}); 
Leann Terry (Indiana University; \email{ljterry@Indiana.Edu})
}

%\note{ }

\seealso{\code{\link{CFA.1}};\code{\link[sem]{sem}} }

\examples{
# General Congeneric
# covmat.from.cfm(Lambda=c(.8, .9, .6, .8), Psi.Square=c(.6, .2, .1, .3), tol.det=.00001)


# True-score equivalent
# covmat.from.cfm(Lambda=c(.8, .8, .8, .8), Psi.Square=c(.6, .2, .1, .3), tol.det=.00001)


# Parallel 
# covmat.from.cfm(Lambda=c(.8, .8, .8, .8), Psi.Square=c(.2, .2, .2, .2), tol.det=.00001)
}

\keyword{multivariate}
