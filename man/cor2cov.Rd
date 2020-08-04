\name{cor2cov}
\alias{cor2cov}
\title{ Correlation Matrix to Covariance Matrix Conversion }
\description{
Function to convert a correlation matrix to a covariance matrix.
}
\usage{
cor2cov(cor.mat, sd, discrepancy=1e-5)
}
\arguments{
  \item{cor.mat}{ the correlation matrix to be converted }
  \item{sd}{ a vector that contains the standard deviations of the variables in the correlation matrix }
  \item{discrepancy}{a neighborhood of 1, such that numbers on the main diagonal of the correlation matrix
  will be considered as equal to 1 if they fall in this neighborhood}
}

\details{The correlation matrix to convert can be either symmetric or triangular. The covariance matrix
returned is always a symmetric matrix.}

%\value{}

%\references{  }

\author{ Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu}), Keke Lai }

\note{ The correlation matrix input should be a square matrix, and the length of \code{sd} should be equal to
the number of variables in the correlation matrix (i.e., the number of rows/columns). Sometimes the correlation
matrix input may not have exactly 1's on the main diagonal, due to, eg, rounding; \code{discrepancy} specifies
the allowable discrepancy so that the function still considers the input as a correlation matrix and can 
proceed (but the function does not change the numbers on the main diagonal).    
}

%\seealso{  }

%\examples{}

\keyword{ design }
