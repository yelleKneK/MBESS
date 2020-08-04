\name{ss.aipe.cv.sensitivity}
\alias{ss.aipe.cv.sensitivity}

\title{Sensitivity analysis for sample size planning given the Accuracy in Parameter Estimation approach for the coefficient of variation.}
\description{
Performs sensitivity analysis for sample size determination for the coefficient of variation
given a population coefficient of variation (or population mean and standard deviation) and goals for the 
sample size procedure. Allows one to determine the effect of being wrong when estimating the 
population coefficient of variation in terms of the width of the obtained (two-sided) confidence intervals. 
The values of \code{True.C.of.V} and \code{Estimated.C.of.V} should be positive.}
\usage{
ss.aipe.cv.sensitivity(True.C.of.V = NULL, Estimated.C.of.V = NULL,
width = NULL, degree.of.certainty = NULL, assurance=NULL, certainty=NULL, 
mean = 100, Specified.N = NULL, conf.level = 0.95, 
G = 1000, print.iter = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{True.C.of.V}{population coefficient of variation}
  \item{Estimated.C.of.V}{estimated coefficient of variation}
  \item{width}{desired confidence interval width}
  \item{degree.of.certainty}{parameter to ensure confidence interval width with a specified degree of certainty (must be NULL or between zero and unity)}
  \item{assurance}{ the alias for \code{degree.of.certainty} }
\item{certainty}{an alias for \code{degree.of.certainty}}
  \item{mean}{Some arbitrary value that the simulation uses to generate data (the variance of the data is determined by the mean and the coefficient of variation)}
  \item{Specified.N}{selected sample size to use in order to determine distributional properties of at a given value of sample size (not used with \code{Estimated.C.of.V})}
  \item{conf.level}{the desired degree of confidence (i.e., 1-Type I error rate).}
  \item{G}{number of generations (i.e., replications) of the simulation}
  \item{print.iter}{to print the current value of the iterations}
}
\details{

For sensitivity analysis when planning sample size given the desire to obtain narrow confidence intervals
for the population coefficient of variation. Given a population value and an estimated value, one can determine
the effects of incorrectly specifying the population coefficient of variation (\code{True.C.of.V}) on the 
obtained widths of the confidence intervals. Also, one can evaluate the percent of the confidence intervals
that are less than the desired width (especially when modifying the \code{degree.of.certainty} parameter); see \code{ss.aipe.cv})

Alternatively, one can specify \code{Specified.N} to determine the results at a particular sample size (when doing this \code{Estimated.C.of.V} cannot be specified).
}
\value{
  \item{Data.from.Simulation}{list of the results in matrix form}
  \item{Specifications}{specification of the function}
  \item{Summary.of.Results}{summary measures of some important descriptive statistics}  
}

%\references{ ~put references to the literature/web site here ~ }

\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})}
\note{Returns three lists, where each list has multiple components.}

\seealso{\code{cv}, \code{ss.aipe.cv}}
%\examples{}

\keyword{design}
\keyword{htest}
