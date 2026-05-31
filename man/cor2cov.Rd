\name{cor2cov}
\alias{cor2cov}
\title{ Correlation Matrix to Covariance Matrix Conversion }
\description{
Function to convert a correlation matrix to a covariance matrix.
}
\usage{
cor2cov(cor.mat, sd)
}
\arguments{
  \item{cor.mat}{ the correlation matrix to be converted }
  \item{sd}{ a vector that contains the standard deviations of the variables in the correlation matrix }
}

\details{The correlation matrix to convert can be either symmetric or triangular. The covariance matrix
returned is always a symmetric matrix.}

%\value{}

%\references{  }

\author{ Ken Kelley (University of Notre Dame; \email{kkelley@nd.edu}), Keke Lai }

\note{ The correlation matrix input should be a square matrix, and the length of \code{sd} should be equal to
the number of variables in the correlation matrix (i.e., the number of rows/columns).  
}

%\seealso{  }

%\examples{}

\keyword{ design }
