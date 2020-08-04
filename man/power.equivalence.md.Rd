\name{power.equivalence.md}
\alias{power.equivalence.md}
\title{ Power of Two One-Sided Tests Procedure (TOST) for Equivalence }
\description{
A function to calculate the power of the two one-sided tests prodedure (TOST). This is 
the probability that a confidence interval lies within a specified equivalence 
interval. (See also package \code{equivalence}, function \code{tost}.)

}
\usage{
power.equivalence.md(alpha, logscale, ltheta1, ltheta2, ldiff, sigma, n, nu)
}

\arguments{
  \item{alpha}{ \emph{alpha} level for the 2 emph{t}-tests (usually \emph{alpha}=0.05). 
  Confidence interval for full test is at level 1- 2*\code{alpha} }
  \item{logscale}{whether to use logarithmic scale (\code{TRUE}) or not (\code{FALSE}) }
  \item{ltheta1}{ lower limit of equivalence interval}
  \item{ltheta2}{ upper limit of equivalence interval }
  \item{ldiff}{ true difference (ratio on log scale) in treatment means }
  \item{sigma}{ \code{sqrt}(error variance) as fraction (root MSE from ANOVA, or coefficient of variation) }
  \item{n}{ number of subjects per treatment (number of total subjects for crossover design) }
  \item{nu}{ degrees of freedom for \code{sigma} }
}

%\details{}

\value{
  \item{power}{Power of TOST; the probability that the confidence interval will lie within ['theta1', 'theta2'] given \code{sigma}, \code{n}, and \code{nu}}
  }
 
 
\references{ 
Diletti, E., Hauschke D. & Steinijans, V.W. (1991). Sample size determination of 
bioequivalence assessment by means of confidence intervals, 
\emph{International Journal of Clinical Pharmacology, Therapy and Toxicology, 29}, No. 1, 1--8.

Phillips, K.F. (1990). Power of the Two One-Sided Tests Procedure in Bioquivalence,
\emph{Journal of Pharmacokinetics and Biopharmaceutics, 18}, No. 2, 139--144.

Schuirmann, D.J. (1987). A comparison of the two one-sided tests procedure and the power approach for assessing the equivalence of average bioavailability,
\emph{Journal of Pharmacokinetics and Biopharmaceutics, 15}. 657--680.
}


\author{
Kem Phillips; \email{kemphillips@comcast.net}
 }


%\note{}

%\seealso{ \code{\link{power.equivalence.md.plot}}, \code{\link{power.density.equivalence.md}}}

\examples{
 # Suppose that two formulations of a drug are to be compared on 
 # the regular scale using a two-period crossover design, with 
 # theta1 = -0.20, theta2 = 0.20, rm{CV} = 0.20, the 
 # difference in the mean bioavailability is 0.05 (5 percent), and we choose 
 # n=24, corresponding to 22 degrees of freedom.  We need to test 
 # bioequivalence at the 5 percent significance level, which corresponds to 
 # having a 90 percent confidence interval lying within (-0.20, 0.20). Then 
 # the power will be 0.8029678.  This corresponds to Phillips (1990), 
 # Table 1, 5th row, 5th column, and Figure 3.  Use
 
power.equivalence.md(.05, FALSE, -.2, .2, .05, .20, 24, 22)


# If the formulations are compared on the logarithmic scale with 
# theta1 = 0.80, theta2 = 1.25, n=18 (16 degrees of freedom), and 
# a ratio of test to reference of 1.05. Then the power will be 0.7922796.
# This corresponds to Diletti, Table 1, power=.80, CV=.20, ratio=1.05, and Figure 1c. Use

power.equivalence.md(.05, TRUE, .8, 1.25, 1.05, .20, 18, 16)
}


\keyword{design}
