\name{power.density.equivalence.md}
\alias{power.density.equivalence.md}

\title{ Density for power of two one-sided tests procedure (TOST) for equivalence }
\description{
A function to calculate density for the power of the two one-sided tests prodedure (TOST). (See package \code{equivalence}, function \code{tost}.)
}
\usage{
power.density.equivalence.md(power_sigma, alpha = alpha, theta1 = theta1, 
theta2 = theta2, diff = diff, sigma = sigma, n = n, nu = nu)
}

\arguments{
  \item{power_sigma}{x-value for integration}
  \item{alpha}{ \code{alpha} level for the 2 \emph{t}-tests (usually \code{alpha}=0.05). Confidence interval for full test is at level (1-2*\eqn{\alpha}) }
  \item{theta1}{ lower limit of equivalence interval on appropriate scale (regular or log) }
  \item{theta2}{ upper limit of equivalence interval on appropriate scale }
  \item{diff}{ true difference (ratio on log scale) in treatment means on appropriate scale }
  \item{sigma}{ sqrt(error variance) as fraction (root MSE from ANOVA, or coefficient of variation) }
  \item{n}{ number of subjects per treatment (number of total subjects for crossover design) }
  \item{nu}{ degrees of freedom for sigma }
}

%\details{}

\value{
\item{power_density}{density at diff for power of TOST: the probability that the confidence interval will lie within ['theta1', 'theta2']}
}
\references{ 
Diletti, E., Hauschke D. & Steinijans, V.W. (1991). Sample size determination of bioequivalence assessment by means of confidence intervals, \emph{International Journal of Clinical Pharmacology, Therapy and Toxicology, 29}, No. 1, 1--8.

Phillips, K.F. (1990). Power of the Two One-Sided Tests Procedure in Bioquivalence, \emph{Journal of Pharmacokinetics and Biopharmaceutics, 18}, No. 2, 139--144.

Schuirmann, D.J. (1987). A comparison of the two one-sided tests procedure and the power approach for assessing the equivalence of average bioavailability,
\emph{Journal of Pharmacokinetics and Biopharmaceutics, 15}. 657--680.
}

\author{ Kem Phillips; \email{kemphillips@comcast.net} }

%\note{ }


\seealso{ \code{\link{power.equivalence.md.plot}}, \code{\link{power.density.equivalence.md}}}

\examples{
\dontrun{
# This function is called by power.equivalence.md within 
# the integrate function. It is integrated over 
# appropriate limits to compute the power. Use

power.density.equivalence.md(.1, alpha=.05, theta1=-.2, theta2=.2, diff=.05, 
	sigma= .20, n=24, nu=22)

# The usage for the logarithmic scale is the same, except that 
# theta1, theta2, and diff must be on that scale. That is, use log(.8), etc.
}
}
\keyword{ design }
