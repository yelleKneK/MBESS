\name{transform_r.Z}
\alias{transform_r.Z}
\title{Transform a correlation coefficient (r) into the scale of Fischer's \emph{Z}
}
\description{
This function transform a correlation coefficient into the scale of Fischer's \emph{Z}.}
\usage{
transform_r.Z(r)
}
\arguments{
  \item{r}{correlation coefficient (between two variables)}
}
\details{
This function is typically used in the context of forming a confidence interval for a population correlation coefficient. Note that, in that situation, the two variables are assumed to follow a bivariate normal distribution (e.g., Hays, 1994). 
}
\value{
returns a value on the scale of Fisher's \emph{Z}, also called Fisher's emph{Z}\', from a given correlation value. 
}
\references{
Kelley, K. (2007). Confidence intervals for standardized effect sizes: Theory, application, 
and implementation.\emph{ Journal of Statistical Software, 20}(8), 1--24.

Hays, W. L. (1994). \emph{Statistics} (5th ed). Fort Worth, TX: Harcourt Brace College Publishers)
}
\author{ Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu}) }

\seealso{
\code{\link{transform_Z.r}}, \code{\link{ci.cc}}
}
\examples{
# From Hays (1994, pp. 649--650)
transform_r.Z(.35)
}

\keyword{ design }
