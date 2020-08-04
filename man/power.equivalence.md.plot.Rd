\name{power.equivalence.md.plot}
\alias{power.equivalence.md.plot}

\title{ Plot power of Two One-Sided Tests Procedure (TOST) for Equivalence}
\description{
A function to plot the power of the two one-sided tests prodedure (TOST) for various alternatives. (See also package \code{equivalence}, function \code{tost}.)
}
\usage{
power.equivalence.md.plot(alpha, logscale, theta1, theta2, sigma, n, nu, title2)
}

\arguments{
  \item{alpha}{ \emph{alpha} level for the 2 \emph{t}-tests (usually \emph{alpha}=0.05). 
  Confidence interval for full test is at level 1- 2*\code{alpha} }
  \item{logscale}{ whether to use logarithmic scale \code{TRUE} or not \code{FALSE} }
  \item{theta1}{ lower limit of equivalence interval}
  \item{theta2}{ upper limit of equivalence interval }
  \item{sigma}{ \code{sqrt}(error variance) as fraction (root MSE from ANOVA, or coefficient of variation) }
  \item{n}{ number of subjects per treatment (number of total subjects for crossover design) }
  \item{nu}{ degrees of freedom for \code{sigma} }
  \item{title2}{ Title appearing at bottom of plot }
}

%\details{}

\value{
  \item{power}{Plot of power of TOST (probability that (1-2*\code{alpha}) confidence 
  interval will lie within (\code{theta1}, \code{theta2}) given \code{sigma}, \code{n}, 
  and \code{nu}.  Also returns matrix of 201 differences between \code{theta1} 
  and \code{theta2} as first column, and power values corresponding to \code{n} for other columns.}
}


\references{ 
Diletti, E., Hauschke D. & Steinijans, V.W. (1991) Sample size determination of bioequivalence assessment by means of confidence intervals, \emph{International Journal of Clinical Pharmacology, Therapy and Toxicology, 29}, No. 1, 1-8.

Phillips, K.F. (1990) Power of the Two One-Sided Tests Procedure in Bioquivalence,
\emph{Journal of Pharmacokinetics and Biopharmaceutics, 18}, No. 2, 139-144.

Schuirmann, D.J. (1987) A comparison of the two one-sided tests procedure and the power approach for assessing the equivalence of average bioavailability,
\emph{Journal of Pharmacokinetics and Biopharmaceutics, 15}. 657-680.
}

\author{
Kem Phillips; \email{kemphillips@comcast.net}
 }
 
%\note{ }

%\seealso{ \code{\link{power.equivalence.md.plot}}, \code{\link{power.density.equivalence.md}}}

\examples{
\dontrun{
# Suppose that two formulations of a drug are to be compared 
# on the regular scale using a two-period crossover design, 
# with theta1 = -0.20, theta2 = 0.20, rm(CV) = 0.20, and 
# we choose 
n<-c(9,12,18,24,30,40,60)

# corresponding to 
nu<-c(7,10,16,22,28,38,58)

# degrees of freedom.  We need to test bioequivalence at the 
# .05 significance level, which corresponds to having a .90 confidence
# interval lying within (-0.20, 0.20). This corresponds to 
# Phillips (1990),  Figure 3.  Use

power.equivalence.md.plot(.05, FALSE, -.2, .2, .20, n, nu, 'Phillips Figure 3')

# If the formulations are compared on the logarithmic scale with 
# theta1 = 0.80, theta2 = 1.25, and 

n<-c(8,12,18,24,30,40,60)

# corresponding  to 
nu<-c(6,10,16,22,28,38,58)

# degrees of freedom. This corresponds to Diletti, Figure 1c. Use

power.equivalence.md.plot(.05, TRUE, .8, 1.25, .20, n, nu, 'Diletti, Figure 1c')
}
}

\keyword{ design}
