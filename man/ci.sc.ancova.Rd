\name{ci.sc.ancova}
\alias{ci.sc.ancova}
\title{ Confidence interval for a standardized contrast in ANCOVA with one covariate}
\description{
Calculate the confidence interval for a standardized contrast in ANCOVA with one covariate. The standardizer
(i.e., the divisor) can be either the error standard deviation of the ANOVA model (i.e., the model excluding the covariate)
or of the ANCOVA model.
}
\usage{
ci.sc.ancova(Psi=NULL, adj.means=NULL, s.anova = NULL, s.ancova, 
standardizer = "s.ancova", c.weights, n, cov.means, SSwithin.x, 
conf.level = 0.95)
}
\arguments{
  \item{Psi}{unstandardized contrast of adjusted means }
  \item{adj.means}{the vector that contains the adjusted mean of each group on the dependent variable}
  \item{s.anova}{the standard deviation of the errors from the ANOVA model (i.e., the square root of the mean square error from ANOVA)}
  \item{s.ancova}{the standard deviation of the errors from the ANCOVA model (i.e., the square root of the mean square error from ANCOVA)}
  \item{standardizer}{which error standard deviation the user wants to use, the value of which can be
  either \code{"s.ancova"} or \code{"s.anova"}}
  \item{c.weights}{the contrast weights (chose weights so that the positive \emph{c}-weights sum to 1 and the negative \emph{c}-weights sum to -1; i.e., use fractional values not integers).}
  \item{n}{either a single number that indicates the sample size per group, or a vector that contains the
    sample size of each group}
  \item{cov.means}{ a vector that contains the group means of the covariate}
  \item{SSwithin.x}{the sum of squares within groups obtained from the summary table for ANOVA on the covariate}
  \item{conf.level}{the desired confidence interval coverage, (i.e., 1 - Type I error rate)}
}

%\details{}

\value{
\item{standardizer}{the divisor used in the standardization}
\item{psi.limit.lower}{the lower confidence limit of the standardized contrast}
\item{psi}{the estimated contrast}
\item{psi.limit.upper}{the upper confidence limit of the standardized contrast}
}

\references{ 
Kelley, K. (2007). Constructing confidence intervals for standardized effect sizes: Theory, application, 
and implementation. \emph{Journal of Statistical Software, 20} (8), 1--24.

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the standardized mean difference: 
Accuracy in Parameter Estimation via narrow confidence intervals. \emph{Psychological Methods, 11}, 363--385.

Lai, K., & Kelley, K. (2012). Accuracy in parameter estimation for ANCOVA and ANOVA contrasts: Sample size planning via narrow confidence intervals. \emph{British Journal of Mathematical and Statistical Psychology, 65}, 350--370.

Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval estimation and the evaluation of
statistical methods. In L. L. Harlow, S. A. Mulaik, & J.H. Steiger (Eds.), \emph{What if there were
no significance tests?} (pp. 221--257). Mahwah, NJ: Lawrence Erlbaum. }

\author{Keke Lai (University of California--Merced) and Ken Kelley \email{kkelley@nd.edu}}

\note{
Be sure to use the standard deviations and not the error variances for \code{s.anova} and \code{s.ancova}, not the squares of these values which would come from the source tables (i.e., do not use the variance of the errors but rather use its square root, the standard deviation).

If \code{n} receives a single number, that number is considered as the sample size per group. If \code{n}
is assigned to a vector, the vector is considered as the sample size of each group.

Be sure to use fractional \emph{c}-weights when doing complex contrasts (not integers) to specify \code{c.weights}. For example, in an ANCOVA of four groups, if the user wants to compare the mean of group 1 and 2 with the mean of group 3 and 4, \code{c.weights} should be specified as c(0.5, 0.5, -0.5, -0.5) rather than c(1, 1, -1, -1). Make sure the sum of the contrast weights are zero.

The argument to be assigned to \code{standardizer} must be either \code{"s.ancova"} or \code{"s.anova"}. 
}

\seealso{\code{ci.c.ancova}, \code{ci.sc} }

\examples{
# Maxwell & Delaney (2004, pp. 428--468) offer an example that 30 depressive 
# individuals are randomly assigned to three groups, 10 in each, and ANCOVA 
# is performed on the posttest scores using the participants' pretest 
# scores as the covariate. The means of pretest scores of group 1, 2, and 3 are 
# 17, 17.7, and 17.4, respectively, whereas the adjusted means of groups 1, 2, and 3 
# are 7.5, 12, and 14, respectively. The error variance in ANCOVA is 29 and thus 
# 5.385165 is the error standard deviation, with the sum of squares within groups 
# from an ANOVA on the covariate is 752.5. 

# To obtained the confidence interval for the standardized adjusted
# mean difference between group 1 and 2, using the ANCOVA error standard
# deviation:
ci.sc.ancova(adj.means=c(7.5, 12, 14), s.ancova=5.385165, c.weights=c(1,-1,0), 
n=10, cov.means=c(17, 17.7, 17.4), SSwithin.x=752.5)

# Or, with less error in rounding:
ci.sc.ancova(adj.means=c(7.54, 11.98, 13.98), s.ancova=5.393, c.weights=c(-1,0,1), 
n=10, cov.means=c(17, 17.7, 17.4), SSwithin.x=752.5)

# Now, using the standard deviation from ANOVA (and not ANCOVA as above), we have:
ci.sc.ancova(adj.means=c(7.54, 11.98, 13.98), s.anova=6.294, s.ancova=5.393, c.weights=c(-1,0,1),
n=10, cov.means=c(17, 17.7, 17.4), SSwithin.x=752.5, standardizer= "s.anova", conf.level=.95)
}
\keyword{design }
