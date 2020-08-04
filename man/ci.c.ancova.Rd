\name{ci.c.ancova}
\alias{ci.c.ancova}
\title{ Confidence interval for an (unstandardized) contrast in ANCOVA with one covariate}
\description{
To calculate the confidence interval for an unstandardized contrast in the one-covariate ANCOVA.
}
\usage{
ci.c.ancova(Psi, adj.means, s.ancova = NULL, c.weights, n,
cov.means, SSwithin.x, conf.level = 0.95, ...)
}
\arguments{
  \item{Psi}{ the unstandardized contrast of adjusted means}
  \item{adj.means}{the vector that contains the adjusted mean of each group on the dependent variable}
  \item{s.ancova}{the standard deviation of the errors from the ANCOVA model (i.e., the square root of the mean square error from ANCOVA)}
  \item{c.weights}{ the contrast weights }
  \item{n}{ either a single number that indicates the sample size \emph{per group} or a vector that contains the
    sample size of each group}
  \item{cov.means}{ a vector that contains the group means of the covariate}
  \item{SSwithin.x}{ the sum of squares within groups obtained from the summary table for ANOVA on the covariate }
  \item{conf.level}{the desired confidence interval coverage, (i.e., 1 - Type I error rate)}
  \item{\dots}{ allows one to potentially include parameter values for inner functions }
}

%\details{}

\value{
\item{lower.limit}{the lower confidence limit of the (unstandardized) ANCOVA contrast}
\item{upper.limit}{the upper confidence limit of the (unstandardized) ANCOVA contrast}
}

\references{
Kelley, K. (2007). Constructing confidence intervals for standardized effect sizes: Theory, application, and implementation. \emph{Journal of Statistical Software, 20} (8), 1--24.

Maxwell, S. E., & Delaney, H. D. (2004). \emph{Designing experiments and analyzing data: A model comparison perspective}. Mahwah, NJ: Erlbaum. 
}

\author{Keke Lai (University of California--Merced) and Ken Kelley (University of Notre Dame; \email{kkelley@nd.edu})}

\note{ 
Be sure to use the standard deviation and not the error variance for \code{s.ancova}, not the square of this value which would come from the source table (i.e., do not use the variance of the error but rather use the square root).

If \code{n} receives a single number, that number is considered as the sample size \emph{per group}. If \code{n}
receives a vector, the vector is considered as the sample size of each group.

Be sure to use fractions not the integers to specify \code{c.weights}. For example, in an ANCOVA of four groups, 
if the user wants to compare the mean of group 1 and 2 with the mean of group 3 and 4, \code{c.weights} should
be specified as c(0.5, 0.5, -0.5, -0.5) rather than c(1, 1, -1, -1). Make sure the sum of the contrast weights 
are zero.
}

\seealso{\code{ci.c}, \code{ci.sc.ancova} }

\examples{
# Maxwell & Delaney (2004, pp. 428-468) offer an example that 30 depressive 
# individuals are randomly assigned to three groups, 10 in each, and ANCOVA 
# is performed on the posttest scores using the participants' pretest 
# scores as the covariate. The means of pretest scores of group 1 to 3 are 
# 17, 17.7, and 17.4, respectively, and the adjusted means of groups 1 to 3 
# are 7.5, 12, and 14, respectively. The error variance in ANCOVA is 29, 
# and the sum of squares within groups from ANOVA on the covariate is 
# 313.37. 

# To obtain the confidence interval for adjusted mean of group 1 versus 
# group 2:
ci.c.ancova(adj.means=c(7.5, 12, 14), s.ancova=sqrt(29), c.weights=c(1, -1, 0), 
n=10, cov.means=c(17, 17.7, 17.4), SSwithin.x=313.37)}

\keyword{ design }
