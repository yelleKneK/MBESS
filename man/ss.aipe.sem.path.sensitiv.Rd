\name{ss.aipe.sem.path.sensitiv}
\alias{ss.aipe.sem.path.sensitiv}

\title{ a priori Monte Carlo simulation for sample size planning for SEM targeted effects}
\description{
Conduct a priori Monte Carlo simulation to empirically study the effects of (mis)specifications of input information on the calculated sample size. Random data are generated from the true covariance matrix but fit to the proposed model, whereas sample size is calculated based on the input covariance matrix and proposed model. 
}
\usage{
ss.aipe.sem.path.sensitiv(model, est.Sigma, true.Sigma = est.Sigma, 
which.path, desired.width, N=NULL, conf.level = 0.95, assurance = NULL, 
G = 100, ...)
}

\arguments{
  \item{model}{ the model the researcher proposes, may or may not be the true model. This argument should be an RAM (reticular action model; e.g., McArdle & McDonald, 1984) specification of a structural equation model, and should be of class \code{mod}. The model is specified in the same manner as does the \code{\link[sem]{sem}} package; see \code{\link[sem]{sem}} and \code{\link[sem]{specify.model}} for detailed documentation about model specifications in the RAM notation.}
  \item{est.Sigma}{ the covariance matrix used to calculate sample size, may or may not be the true covariance matrix. The row names and column names of \code{est.Sigma} should be the same as the manifest variables in \code{est.model}. }
  \item{true.Sigma}{ the true population covariance matrix, which will be used to generate random data for the simulation study. The row names and column names of \code{est.Sigma} should be the same as the manifest variables in \code{est.model}.}
  \item{which.path}{ the name of the model parameter of interest, and must be in a double quote }
  \item{desired.width}{desired confidence interval width for the model parameter of interest }
  \item{N}{the sample size of random data. If it is \code{NULL}, it will be determined by the sample size planning method}
  \item{conf.level}{ confidence level (i.e., 1- Type I error rate) }
  \item{assurance}{the assurance that the confidence interval obtained in a particular study will be no wider than desired (must be \code{NULL} or a value between 0.50 and 1)  }
  \item{G}{number of replications in the Monte Carlo simulation }
  \item{\dots}{ allows one to potentially include parameter values for inner functions }
}

\details{This function implements the sample size planning methods proposed in Lai and Kelley (2010). It depends on the 
function \code{\link[sem]{sem}} in the \code{sem} package to calculate the expected information matrix, and uses the same notation to specify SEM 
models as does \code{\link[sem]{sem}}. Please refer to \code{\link[sem]{sem}} for more detailed documentation 
about model specifications, the RAM notation, and model fitting techniques. For technical discussion 
on how to obtain the model implied covariance matrix in the RAM notation given model parameters, see McArdle and McDonald (1984).}

\value{
  \item{w }{the \code{G} random confidence interval widths }
  \item{sample.size }{the sample size calculated }
  \item{path.of.interest}{name of the model parameter of interest}
  \item{desired.width}{desired confidence interval width}
  \item{mean.width}{mean of the \code{G} random confidence interval widths}
  \item{median.width}{median of the \code{G} random confidence interval widths}
  \item{quantile.width}{99, 95, 90, 85, 80, 75, 70, and 60 percentiles of the \code{G} random confidence interval widths}
  \item{width.less.than.desired}{the proportion of confidence interval widths narrower than desired}
  \item{Type.I.err.upper}{the upper empirical Type I error rate}
  \item{Type.I.err.lower}{the lower empirical Type I error rate}
  \item{Type.I.err}{total empirical Type I error rate}
  \item{conf.level}{confidence level}
  \item{rep}{successful replications}
}
\references{ 
Fox, J. (2006). Structural equation modeling with the sem package in R. \emph{Structural Equation Modeling, 13}, 465--486.

Lai, K., & Kelley, K. (in press). Accuracy in parameter estimation for targeted effects in structural equation modeling: Sample size planning for narrow confidence intervals. \emph{Psychological Methods}.

McArdle, J. J., & McDonald, R. P. (1984). Some algebraic properties of the reticular action model. \emph{British Journal of Mathematical and Statistical Psychology, 37}, 234--251.  }
\author{ Keke Lai (University of California -- Merced) and Ken Kelley \email{kkelley@nd.edu} }

\note{Sometimes the simulation stops in the middle of fitting the model to the random data. The reason is that \code{\link{nlm}}, the 
function \code{\link[sem]{sem}} calls to fit the model, fails to converge. We suggest using the \code{\link[base]{try}} function in simulation so that
the simulation can proceed with unsuccessful iterations.}

\seealso{ \code{\link[sem]{sem}}; \code{\link[sem]{specify.model}}; \code{\link{theta.2.Sigma.theta}}; \code{\link{ss.aipe.sem.path}} }

\examples{
\dontrun{
# Suppose the model of interest is Model 2 of the simulation study in 
# Lai and Kelley (2010), and the goal is to obtain a 95\% confidence 
# interval for 'beta21' no wider than 0.3.

library(sem)

# specify a model object in the RAM notation
model.2<-specifyModel()
xi1 -> y1, lambda1, 1
xi1 -> y2, NA, 1
xi1 -> y3, lambda2, 1
xi1 -> y4, lambda3, 0.3
eta1 -> y4, lambda4, 1
eta1 -> y5, NA, 1
eta1 -> y6, lambda5, 1
eta1 -> y7, lambda6, 0.3
eta2 -> y6, lambda7, 0.3
eta2 -> y7, lambda8, 1
eta2 -> y8, NA, 1
eta2 -> y9, lambda9, 1
xi1 -> eta1, gamma11, 0.6
eta1 -> eta2, beta21, 0.6 
xi1 <-> xi1, phi11, 0.49
eta1 <-> eta1, psi11, 0.3136
eta2 <-> eta2, psi22, 0.3136
y1 <-> y1, delta1, 0.51
y2 <-> y2, delta2, 0.51
y3 <-> y3, delta3, 0.51
y4 <-> y4, delta4, 0.2895
y5 <-> y5, delta5, 0.51
y6 <-> y6, delta6, 0.2895
y7 <-> y7, delta7, 0.2895
y8 <-> y8, delta8, 0.51
y9 <-> y9, delta9, 0.51


# to inspect the specified model
model.2

# one way to specify the population covariance matrix is to
# first specify path coefficients and then calcualte the 
# model-implied covariance matrix
theta <- c(1, 1, 0.3, 1,1, 0.3, 0.3, 1, 1, 0.6, 0.6,
0.49, 0.3136, 0.3136, 0.51, 0.51, 0.51, 0.2895, 0.51, 0.2895, 0.2895, 0.51, 0.51)

names(theta) <- c("lambda1","lambda2","lambda3",
"lambda4","lambda5","lambda6","lambda7","lambda8","lambda9",
"gamma11", "beta21",
"phi11", "psi11", "psi22", 
"delta1","delta2","delta3","delta4","delta5","delta6","delta7",
"delta8","delta9")

res<-theta.2.Sigma.theta(model=model.2, theta=theta, 
latent.vars=c("xi1", "eta1","eta2"))

Sigma.theta <- res$Sigma.theta
# thus 'Sigma.theta' is the input covariance matrix for sample size planning procedure.

# the necessary sample size can be calculated as follows.
# ss.aipe.sem.path(model=model.2, Sigma=Sigma.theta, 
# desired.width=0.3, which.path="beta21")

# to verify the sample size calculated
# ss.aipe.sem.path.sensitiv(est.model=model.2, est.Sigma=Sigma.theta, 
# which.path="beta21", desired.width=0.3, G = 300)

# suppose the true covariance matrix ('var(X)' below) is in fact 
# a point close to 'Sigma.theta':

# X<-mvrnorm(n=1000, mu=rep(0,9), Sigma=Sigma.pop)
# var(X)
# ss.aipe.sem.path.sensitiv(est.model=model.2, est.Sigma=Sigma.theta, 
# true.Sigma=var(X), which.path="beta21", desired.width=0.3, G=300)
}
}
\keyword{ design }
\keyword{ multivariate }
