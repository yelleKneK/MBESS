\name{ss.aipe.R2}
\alias{ss.aipe.R2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Sample Size Planning for Accuracy in Parameter Estimation
for the multiple correlation coefficient. }
\description{
Determines necessary sample size for the multiple correlation coefficient so that the 
confidence interval for the population multiple correlation coefficient is 
sufficiently narrow. Optionally, there is a certainty parameter that allows one to be 
a specified percent certain that the observed interval will be no wider than desired.
}
\usage{
ss.aipe.R2(Population.R2 = NULL, conf.level = 0.95, width = NULL, 
Random.Predictors = TRUE, Random.Regressors, which.width = "Full", p = NULL, 
K, degree.of.certainty = NULL, assurance=NULL, certainty=NULL, 
verify.ss = FALSE, Tol = 1e-09, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Population.R2}{ value of the population multiple correlation coefficient }
  \item{conf.level}{ confidence interval level (e.g., .95, .99, .90); 1-Type I error rate }
  \item{width}{width of the confidence interval (see \code{which.width}) }
  \item{Random.Predictors}{ whether or not the predictor variables are random (set to \code{TRUE}) or are fixed (set to \code{FALSE}) }
  \item{Random.Regressors}{ an alias for \code{Random.Predictors}; \code{Random.Regressors} 
  overrides \code{Random.Predictors} }
  \item{which.width}{ defines the width that \code{width} refers to }
  \item{p}{ the number of predictor variables }
  \item{K}{ an alias for \code{p}; \code{K} overrides \code{p} }
  \item{degree.of.certainty}{value with which confidence can be placed that describes the likelihood of obtaining a confidence interval less than the value specified (e.e.g, .80, .90, .95)}
  \item{assurance}{ an alias for \code{degree.of.certainty} }
\item{certainty}{an alias for \code{degree.of.certainty}}
  \item{verify.ss}{evaluates numerically via an internal Monte Carlo simulation the exact sample size given the specifications}
  \item{Tol}{the tolerance of the iterative function \code{conf.limits.nct} for convergence}  
  \item{\dots}{for modifying the parameters of functions this function calls upon}
}
\details{
This function determines a necessary sample size so that the expected confidence 
interval width for the squared multiple correlation coefficient is sufficiently narrow (when \code{degree.of.certainty=NULL}) so that the obtained confidence 
interval is no larger than the value specified with some desired degree of 
certainty (i.e., a probability that the obtained width is less than the specified 
width). The method depends on whether or not the regressors are regarded as fixed 
or random. This is the case because the distribution theory for the two cases is 
different and thus the confidence interval procedure is conditional on the type of 
regressors. The default methods are approximate but can be made exact with the 
specification of \code{verify.ss=TRUE}, which performs an a priori Monte Carlo simulation study. Kelley (2007) and Kelley & Maxwell (2008) detail the methods used in the 
function, with the former focusing on random regressors and the latter on fixed regressors.

It is recommended that the option \code{verify.ss} should always be used! Doing so uses the method implied sample size as an estimate and then evaluates with an internal Monte Carlo simulation (i.e., via "brute-force" methods) the exact sample size given the goals specified. When \code{verify.ss=TRUE}, the default number of iterations is 10,000 but this can be changed by specifying G=5000 (or some other value; 10000 is the recommended) When \code{verify.ss=TRUE} is specified, an internal function \code{verify.ss.aipe.R2} calls upon the \code{ss.aipe.R2.sensitivity} function for purposes of the 
internal Monte Carlo simulation study. See the \code{verify.ss.aipe.R2} function for arguments that can be passed from \code{ss.aipe.R2} to \code{verify.ss.aipe.R2}.
}

\value{
 \item{Required.Sample.Size}{sample size that should be used given the conditions specified.}
%\item{Expected.Width}{given the sample size, the expected (half)-width of interest.}
}

\references{ 
Algina, J. & Olejnik, S. (2000). Determining sample size for accurate estimation of the squared
multiple correlation coefficient. \emph{Multivariate Behavioral Research, 35}, 119--136.

Steiger, J. H. & Fouladi, R. T. (1992). R2: A computer program for interval estimation, power calculation,
and hypothesis testing for the squared multiple correlation. \emph{Behavior research methods, instruments and computers, 4}, 581--582.

Kelley, K. (2007). Sample size planning for the squared multiple correlation coefficient:
Accuracy in parameter estimation via narrow confidence intervals, 
\emph{manuscripted submitted for publication}.

Kelley, K. & Maxwell, S. E. (2008). Power and accuracy for omnibus and targeted effects: 
Issues of sample size planning with applications to multiple regression. In P. Alasuuta, J. Brannen, & L. Bickman (Eds.),
\emph{Handbook of Social Research Methods} (pp. 166--192). Newbury Park, CA: Sage.
}

\author{ Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu}) }

\note{ This function without \code{verify.SS=FALSE} can be slow to converge when
 \code{verify.SS=TRUE}, the function can take some time to converge (e.g., 15 minutes). 
 Most times this will not be the case, but it is possible in some situations.
}


\seealso{ \code{ci.R2}, \code{conf.limits.nct}, \code{ss.aipe.R2.sensitivity} }


\examples{
\dontrun{
# Returned sample size should be considered approximate; exact sample
# size is obtained by specifying the argument 'verify.ss=TRUE' (see below).
# ss.aipe.R2(Population.R2=.50, conf.level=.95, width=.10, which.width="Full",
# p=5, Random.Predictors=TRUE)
# Uncomment to run in order to get exact sample size.
# ss.aipe.R2(Population.R2=.50, conf.level=.95, width=.10, which.width="Full",
# p=5, Random.Predictors=TRUE, verify.ss=TRUE)


# Same as above, except the predictor variables are considered fixed.
# Returned sample size should be considered approximate; exact sample
# size is obtained by specifying the argument 'verify.ss=TRUE'.
# ss.aipe.R2(Population.R2=.50, conf.level=.95, width=.10, which.width="Full",
# p=5, Random.Predictors=FALSE)
# Uncomment to run in order to get exact sample size.
#ss.aipe.R2(Population.R2=.50, conf.level=.95, width=.10, which.width="Full",
#p=5, Random.Predictors=FALSE, verify.ss=TRUE)


# Returned sample size should be considered approximate; exact sample
# size is obtained by specifying the argument 'verify.ss=TRUE'.
# ss.aipe.R2(Population.R2=.50, conf.level=.95, width=.10, which.width="Full",
# p=5, degree.of.certainty=.85, Random.Predictors=TRUE)
# Uncomment to run in order to get exact sample size.
#ss.aipe.R2(Population.R2=.50, conf.level=.95, width=.10, which.width="Full",
#p=5, degree.of.certainty=.85, Random.Predictors=TRUE, verify.ss=TRUE)


# Same as above, except the predictor variables are considered fixed.
# Returned sample size should be considered approximate; exact sample
# size is obtained by specifying the argument 'verify.ss=TRUE'.
# ss.aipe.R2(Population.R2=.50, conf.level=.95, width=.10, which.width="Full",
# p=5, degree.of.certainty=.85, Random.Predictors=FALSE)
# Uncomment to run in order to get exact sample size.
#ss.aipe.R2(Population.R2=.50, conf.level=.95, width=.10, which.width="Full",
#p=5, degree.of.certainty=.85, Random.Predictors=FALSE, verify.ss=TRUE)
}
}


\keyword{design}
