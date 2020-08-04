\name{mr.smd}
\alias{mr.smd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Minimum risk point estimation of the population standardized mean difference}
\description{
A function for the sequential estimation of the standardized mean difference with minimum risk. The function implements the ideas of Chattopadhyay and Kelley (submitted, Psychological Methods), which considers study cost and accuracy of the estimated
standardized mean difference simultaniously. This is important to specify that \code{mr.smd.R} was developed under the assumption of normally distributed data with equal sample size and equal cost of sampling per observation for each group.
}
\usage{
mr.smd(A, structural.cost, epsilon, d, n, sampling.cost, pilot = FALSE, m0 = 4, 
gamma = 0.49)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{is the price one is willing to pay in order to have a maximum allowable difference of \eqn{epsilon^2} between the estimate of the standardized mean difference and its corresponding parameter.
}
  \item{structural.cost}{
%%     ~~Describe \code{structural.cost} here~~
}
  \item{epsilon}{
The maximum desired difference between the estimated standardized mean difference and the population value)
}
  \item{d}{the current estimate of the standardized mean difference}
  \item{n}{current sample size \emph{per group} (thus total sample size is \eqn{2n}); requires equal sample size \emph{per group}.}
  \item{sampling.cost}{
The sampling cost to collect an additional observation. For example, if each survey costs 10 dollars to distribute and score, \code{sampling.cost} would be 10 dollars per additional observation.}
  \item{pilot}{
\code{TRUE} or \code{FALSE} based on whether the users is using the function to plan a pilot sample size (TRUE) or if it is being used to assess if the optimization criterion has been satisfied (FALSE)
}
  \item{m0}{the minimum bound on the initial pilot sample size}
  \item{gamma}{A correction factor in which we suggest .49; see the two Chattopadhyay & Kelley articles for more details (ignorable for most users).}
}
\details{
The standardized mean difference is a widely used measure effect size. In this article, we developed a general theory for estimating the population standardized mean difference by minimizing both the mean square error of the estimator and the total sampling cost. This function implements our ideas discussed in Chattopadhyay and Kelley (submitted). See also Kelley and Rausch (2006) for additional information on the standardized mean difference.

}
\value{
\item{Risk}{\emph{Per group} sample size (this simply repeats what was supplied to the function)}
\item{n1}{Sample size for group 1 (echos the input value)}
\item{n1}{Sample size for group 2 (echos the input value)}
\item{d}{Observed value of the standardized mean difference (i.e., \emph{d}; echos the input value)}
\item{Is.Satisfied?}{A \code{TRUE} or \code{FALSE} statement of that evaluates a stopping rule using the risk function to determine if the optimation criterion has been satisfied (based on the goals of the researcher and current information available)}
}
\references{
Chattopadhyay, B., & Kelley, K. (submitted, minor revision requested). Estimating the standardized mean difference with minimum risk: Maximizing accuracy and minimizing cost with sequential estimation. \emph{Psychological Methods}, \emph{X}, X--X.

Chattopadhyay, B., & Kelley, K. (in press). Estimation of the Coefficient of Variation with Minimum Risk: A Sequential Method for Minimizing Sampling Error and Study Cost. \emph{Multivariate Behavioral Research}, \emph{X}, X--X.

Kelley, K., & Rausch, J. R. (2006). Sample size planning for the standardized mean difference: 
Accuracy in Parameter Estimation via narrow confidence intervals. \emph{Psychological Methods, 11}, 363--385.
}

\author{
Ken Kelley (University of Notre Dame; \email{kkelley@nd.edu}) and Bhargab Chattopadhyay (University of Texas - Dallas; \email{bhargab@utdallas.edu})
}
\note{
When \code{pilot=TRUE} the function returns the size of the pilot sample size, \emph{per group}, that should be used (thus, the total sample size is twice the pilot sample size).}

\seealso{\code{\link{ci.smd}}, \code{\link{mr.cv}}}

\examples{
# To obtain pilot sample size in a situation in which A=10000. Note that 'A' is 
# 'structural.cost' divided by the square of 'epsilon'.

# From Chattopadhyay and Kelley (submitted, minor revision requested)
mr.smd(pilot=TRUE, A=10000, sampling.cost=2.4, gamma=.49)

High.SLS <- c(11, 7, 22, 13, 6, 9, 11, 16, 12, 17, 14, 8, 16)
Low.SLS  <- c(3, 6, 10, 8, 14, 5, 12, 10, 6, 8, 13, 5, 9)

mr.smd(d=1.021484, n=13, A=10000, sampling.cost=2.40, gamma=.49)

# Or, using the smd() function:
mr.smd(d=smd(Group.1=High.SLS, Group.2=Low.SLS), n=13, A=10000, sampling.cost=2.40, gamma=.49)

# Here, for this situation, the stopping rule is satisfied:
mr.smd(d=1.00, n=75, A=10000, sampling.cost=2.40, gamma=.49)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{design}% __ONLY ONE__ keyword per line
\keyword{misc}
\keyword{htest}
