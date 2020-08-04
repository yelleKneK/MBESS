\name{mr.cv}
\alias{mr.cv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Minimum risk point estimation of the population coefficient of variation
}
\description{
A function for the sequential estimation of the coefficient of variations with minimum risk. The function implements the ideas of Chattopadhyay and Kelley (in press), which considers study cost and accuracy of the estimated
coefficient of variation simultaneously. 
}
\usage{
mr.cv(data, A, structural.cost, epsilon, sampling.cost, pilot=FALSE, m0=4, gamma=.49, 
verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{the data for which to evalaute the function
  }
\item{A}{\code{structural.cost}/\code{epsilon}^2; this is the structural cost that one is willing to pay in a study to estimate the coefficient of variation divided by the square of the desired difference (between the estimate and the parameter)
}
\item{structural.cost}{this is the the structural cost of what one is willing to pay in a study (see note below). 
}
\item{epsilon}{
The maximum desired difference between the estimated coefficient of variation and the population value)
}
\item{sampling.cost}{
The sampling cost to collect an additional observation. For example, if each survey costs 10 dollars to distribute and score, \code{sampling.cost} would be 10 dollars per additional observation.
}
\item{pilot}{
\code{TRUE} or \code{FALSE} based on whether the users is using the function to plan a pilot sample size (\code{TRUE}) or if it is being used to assess if the optimization criterion has been satisfied (\code{FALSE})
}
\item{m0}{the minimum bound on the initial pilot sample size}
\item{gamma}{A correction factor in which we suggest .49; see the two Chattopadhyay & Kelley articles for more details (ignorable for most users).}
\item{verbose}{If \code{TRUE}, extra information is printed; defaults to \code{FALSE}
}
}
\details{
The value of \code{epsilon} is context specific; the smaller the value the closer the estimated value will tend to be to the population value.
}
\value{
\item{Risk}{The value of the risk function }
\item{N}{The current sample size}
\item{cv}{The current coefficient of variation}
\item{Is.Satisfied?}{A TRUE/FALSE statement of whether or not the risk function has been satisfied. If TRUE then sampling can stop as the stopping rule has been satisfied.}
}
\references{
Chattopadhyay, B., & Kelley, K. (in press). Estimation of the Coefficient of Variation with Minimum Risk: A Sequential Method for Minimizing Sampling Error and Study Cost. \emph{Multivariate Behavioral Research}, \emph{X}, X--X.

Kelley, K. (2007). Sample size planning for the coefficient of variation from the accuracy in parameter estimation approach. \emph{Behavior Research Methods, 39} \emph{4}, 755--766.

}
\author{
Ken Kelley (University of Notre Dame; \email{kkelley@nd.edu}) and Bhargab Chattopadhyay (University of Texas - Dallas; \email{bhargab@utdallas.edu})
}
\note{
When a study's aim is to estimate a parameter accurately, such as the coefficient of variation, the structural costs and the maximum probable error of the estimate (i.e., \eqn{\epsilon}) are combined to form \eqn{A}. When we say ``what the researcher is willing to pay," we literally mean the structural cost (\eqn{c}) the researcher is willing to invest in a study in order to estimate the parameter of interest with the desired degree of accuracy. This value is implicitly included (along with anticipated sampling cost) in grant applications for empirical studies when a certain amount of money is requested to conduct a study. If a researcher is willing to pay more and/or desire a smaller value of \eqn{\epsilon}, \eqn{A} is larger than it would have been. A larger \eqn{A} value will translate into a more expensive study, holding everything else constant. Notice that \eqn{A} is a fixed value in any investigation, as the researcher specifies \eqn{A} directly or by specifying its two components (structural cost and \eqn{\epsilon}) individually. However, what is not fixed but rather evaluated in multiple steps throughout the process is the sampling cost, as it is unknown the necessary sample size in order to accomplish the study's goal of achieving a sufficiently accurate estimate of the coefficient of variation. This is the core of our contributions: minimizing sampling cost, and thereby study cost, by using a sequential procedure that evaluates a stopping rule using the risk function to determine if the optimation criterion has been satisfied (based on the goals of the researcher and current information available). This function implements the ideas of  sampling error and the study costs are considered simultaneously, so that the cost is not higher than necessary for the tolerable sampling error.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{ci.cv}}, \code{\link{cv}}, \code{\link{mr.smd}}}

\examples{
# Determine pilot sample size:
mr.cv(pilot=TRUE, A=400000, sampling.cost=75, gamma=.49)

# Collect data (the size of which is the pilot sample size)
Data <- c(36, 53, 19, 11, 10, 24, 14, 65, 18, 48, 25, 35, 13, 18, 3, 41, 5, 3)

# Use mr.cv() to assess if the criterion for stopping the sequential study has been satisfied:
mr.cv(data=Data, A=400000, sampling.cost=75, gamma=.49)

# Collect another data (m=1 here) and perform another check:
Data <- c(Data, 44)
mr.cv(data=Data, A=400000, sampling.cost=75, gamma=.49)

# Continue adding obervations, checking each time if m=1, until the minimum risk criteria 
# are satisfied:
Data <- c(Data, 26, 13, 39, 2, 3, 26, 22, 8, 15, 12, 22, 5, 21, 23, 40, 18)
mr.cv(data=Data, A=400000, sampling.cost=75, gamma=.49)
}
\keyword{design}% __ONLY ONE__ keyword per line
\keyword{misc}
\keyword{htest}
