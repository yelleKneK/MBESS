\name{CFA.1}
\alias{CFA.1}

\title{ One-factor confirmatory factor analysis model}

\description{Returns the MLE estimates and the estimated asymptotic covariance matrix of parameter estimates for one-factor confirmatory factor analysis model
}

\usage{
CFA.1(S, N, equal.loading = FALSE, equal.error = FALSE, package="lavaan", 
	se="standard", ...)
}

\arguments{
  \item{S}{ covariance matrix of the indicators }
  \item{N}{ total sample size }
  \item{equal.loading}{logical statement indicating whether the path coefficients are the same}
  \item{equal.error}{ logical statement indicating whether the manifest variables have the same error variances}
  \item{package}{the package used in confirmatory factor analysis (\code{sem} or \code{lavaan}}
  \item{se}{See the \code{\link[lavaan]{cfa}} and check the \code{se} argument}
  \item{\dots}{Additional arguments for the \code{\link[lavaan]{cfa}} function}
}

\value{
  \item{Model}{the factor analysis model specified by the user}
  \item{Factor.Loadings}{factor loadings}
  \item{Indicator.var}{the error variances of the indicator variables}
  \item{Parameter.cov}{the covariance matrix of the parameters}
  \item{converged}{\code{TRUE} or \code{FALSE} statement on if the model converged}
  \item{package}{notes the package used to get the output}
}
\note{The output will differ slightly, both in form and potentially values, based on which package \pkg{lavaan} or \pkg{sem} is used.}

\author{ Keke Lai (University of California-Merced) and Ken Kelley (University of Notre Dame)}

\seealso{ \code{\link[sem]{sem}}, \code{covmat.from.cfm} }
\examples{
\dontrun{
cov.mat<- matrix(
c(1.384, 1.484, 1.988, 2.429, 3.031,
1.484, 2.756, 2.874, 3.588, 4.390,
1.988, 2.874, 4.845, 4.894, 6.080,
2.429, 3.588, 4.894, 6.951, 7.476,
3.031, 4.390, 6.080, 7.476, 10.313), nrow=5)


CFA.1(N=300, S=cov.mat, package="lavaan")

CFA.1(N=300, S=cov.mat, package="sem")

}
}

\keyword{ multivariate}
