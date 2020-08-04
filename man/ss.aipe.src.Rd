\name{ss.aipe.src}
\alias{ss.aipe.src}

\title{ sample size necessary for the accuracy in parameter estimation approach for a standardized regression coefficient of interest }
\description{
A function used to plan sample size from the accuracy in parameter estimation approach for a standardized
regression coefficient of interest given the input specification.
}
\usage{
ss.aipe.src(Rho2.Y_X = NULL, Rho2.k_X.without.k = NULL, K = NULL, 
beta.k = NULL, width, which.width = "Full", sigma.Y = 1, sigma.X.k = 1, 
RHO.XX = NULL, Rho.YX = NULL, which.predictor = NULL, 
alpha.lower = NULL, alpha.upper = NULL, conf.level = .95, 
degree.of.certainty = NULL, assurance=NULL, certainty=NULL, 
Suppress.Statement = FALSE)
}

\arguments{
  \item{Rho2.Y_X}{ Population value of the squared multiple correlation coefficient }
  \item{Rho2.k_X.without.k}{ Population value of the squared multiple correlation coefficient predicting the
\emph{k}th predictor variable from the remaining \emph{p}-1 predictor variables }
  \item{K}{ the number of predictor variables }
  \item{beta.k}{ the regression coefficient for the \emph{k}th predictor variable (i.e., the predictor of
interest) }
  \item{width}{ the desired width of the confidence interval }
  \item{which.width}{ which width (\code{"Full"}, \code{"Lower"}, or \code{"Upper"}) the width refers 
  to (at present, only \code{"Full"} can be specified) }
  \item{sigma.Y}{ the population standard deviation of \emph{Y} (i.e., the dependent variables) }
  \item{sigma.X.k}{ the population standard deviation of the \emph{k}th \emph{X} variable (i.e., the 
  predictor variable of interest) }
  \item{RHO.XX}{ Population correlation matrix for the \emph{p} predictor variables }
  \item{Rho.YX}{ Population \emph{p} length vector of correlation between the dependent variable (\emph{Y})
and the \emph{p} independent variables }
  \item{which.predictor}{ identifies which of the \emph{p} predictors is of interest }
  \item{alpha.lower}{ Type I error rate for the lower confidence interval limit }
  \item{alpha.upper}{ Type I error rate for the upper confidence interval limit }
  \item{conf.level}{ desired level of confidence for the computed interval (i.e., 1 - the Type I error
rate) }
  \item{degree.of.certainty}{ degree of certainty that the obtained confidence interval will be sufficiently narrow, which 
  yields an approximate sample size to be verified with function \code{ss.aipe.reg.coef.sensitivity}
to determine if it is appropriate. }
 \item{assurance}{ an alias for \code{degree.of.certainty} }
\item{certainty}{an alias for \code{degree.of.certainty}}
  \item{Suppress.Statement}{ \code{TRUE/FALSE} statement whether or not a sentence describing the situation defined
is printed with the necessary sample size }
}
\details{
Not all of the arguments need to be specified, only those that provide all of the necessary information
so that the sample size can be determined for the conditions specified.
}
\value{
Returns the necessary sample size in order for the goals of accuracy in parameter estimation to be
satisfied for the confidence interval for a particular regression coefficient given the input specifications.
}
\references{ 
Kelley, K. & Maxwell, S. E. (2003). Sample size for Multiple Regression: Obtaining regression
coefficients that are accurate, not simply significant. \emph{Psychological Methods, 8}, 305--321.

Kelley, K. & Maxwell, S. E. (2008). Sample Size Planning with applications to multiple regression: Power and accuracy for omnibus and targeted effects. In P. Alasuuta, J. Brannen, & L. Bickman (Eds.), \emph{The Sage handbook of social research methods} (pp. 166--192). Newbury Park, CA: Sage.
}
\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})}

\note{This function calls upon \code{\link{ss.aipe.reg.coef}} in MBESS but has a different naming 
scheme. See \code{\link{ss.aipe.reg.coef}} for more details.}

\section{Warning }{As discussed in Kelley and Maxwell (2008), the sample size planning approach from the AIPE perspective used in this function is only an approximation.}

\seealso{ \code{ss.aipe.reg.coef.sensitivity}, \code{conf.limits.nct}, 
\code{ss.aipe.reg.coef}, \code{ss.aipe.rc} }
\examples{
# Exchangable correlation structure
# Rho.YX <- c(.3, .3, .3, .3, .3)
# RHO.XX <- rbind(c(1, .5, .5, .5, .5), c(.5, 1, .5, .5, .5), c(.5, .5, 1, .5, .5),
# c(.5, .5, .5, 1, .5), c(.5, .5, .5, .5, 1))

# ss.aipe.src(width=.1, which.width="Full", sigma.Y=1, sigma.X=1, RHO.XX=RHO.XX,
# Rho.YX=Rho.YX, which.predictor=1, conf.level=1-.05)

# ss.aipe.src(width=.1, which.width="Full", sigma.Y=1, sigma.X=1, RHO.XX=RHO.XX,
# Rho.YX=Rho.YX, which.predictor=1, conf.level=1-.05, degree.of.certainty=.85)
}
\keyword{ design }
