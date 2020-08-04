\name{signal.to.noise.R2}
\alias{signal.to.noise.R2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Signal to noise using squared multiple correlation coefficient}
\description{
Function that calculates five different signal-to-noise ratios using the squared multiple correlation coefficient.
}
\usage{
signal.to.noise.R2(R.Square, N, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R.Square}{usual estimate of the squared multiple correlation coefficient (with no adjustments)}
  \item{N}{sample size}
  \item{p}{number of predictors}
}
\details{
The method of choice is \code{phi2.UMVUE.NL}, but it requires \code{p} of 5 or more. In situations where \code{p} < 5, it is suggested that \code{phi2.UMVUE.L} be used. 
}
\value{
\item{phi2.hat}{Basic estimate of the signal-to-noise ratio using the usual estimate of the squared multiple correlation coefficient: \code{phi2.hat}=\code{R.Square}/(1-\code{R.Square})}
\item{phi2.adj.hat}{Estimate of the signal-to-noise ratio using the usual adjusted R Square in place of \emph{R}-Square: \code{phi2.hat}=\emph{Adj.R2}/(1-\emph{Adj.R2})}
\item{phi2.UMVUE}{Muirhead's (1985) unique minimum variance unbiased estimate of the signal-to-noise ratio (Muirhead uses theta-U): see reference or code for formula}
\item{phi2.UMVUE.L}{Muirhead's (1985) unique minimum variance unbiased linear estimate of the signal-to-noise ratio (Muirhead uses theta-L): see reference or code for formula}
\item{phi2.UMVUE.NL}{Muirhead's (1985) unique minimum variance unbiased nonlinear estimate of the signal-to-noise ratio (Muirhead uses theta-NL); requires the number of predictors to be greater than five: see reference or code for formula}
}
\references{

Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum.

Muirhead, R. J. (1985). Estimating a particular function of the multiple correlation coefficient. \emph{Journal of the American Statistical Association,  80}, 923--925.
}
\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})}

\seealso{
\code{ci.R2}, \code{ss.aipe.R2}
}
\examples{
signal.to.noise.R2(R.Square=.5, N=50, p=2)
signal.to.noise.R2(R.Square=.5, N=50, p=5)
signal.to.noise.R2(R.Square=.5, N=100, p=2)
signal.to.noise.R2(R.Square=.5, N=100, p=5)
}
\keyword{models}
\keyword{htest}
