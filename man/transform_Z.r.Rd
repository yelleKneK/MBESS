\name{transform_Z.r}
\alias{transform_Z.r}
\title{
Transform Fischer's \emph{Z} into the scale of a correlation coefficient 
}
\description{
A function to transform Fischer's \emph{Z} into the scale of a correlation coefficient. 
}
\usage{
transform_Z.r(Z)
}
\arguments{
  \item{Z}{Fisher's Z or Fisher's \emph{Z}\' value.}
}
\details{
This function is typically used in the context of forming a confidence interval for a population correlation coefficient. Note that, in that situation, the two variables are assumed to follow a bivariate normal distribution (e.g., Hays, 1994).
}
\value{
returns a value on the scale of a correlation coefficieint from a value of Fisher's \emph{Z}. 
}
\references{
Kelley, K. (2007). Confidence intervals for standardized effect sizes: Theory, application, 
and implementation.\emph{ Journal of Statistical Software, 20}(8), 1--24.

Hays, W. L. (1994). \emph{Statistics} (5th ed). Fort Worth, TX: Harcourt Brace College Publishers)
}
\author{ Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu}) }

\seealso{
\code{\link{transform_r.Z}}, \code{\link{ci.cc}}
}

\examples{
# From Hays (1994, pp. 649--650)
transform_Z.r(0.3654438)
}
\keyword{ design }
