\name{theta.2.Sigma.theta}
\alias{theta.2.Sigma.theta}

\title{Compute the model-implied covariance matrix of an SEM model }
\description{Obtain the model-implied covariance matrix of manifest variables given a structural equation model and its model parameters
}
\usage{
theta.2.Sigma.theta(model, theta, latent.vars)
}

\arguments{
  \item{model}{ an RAM (reticular action model; e.g., McArdle & McDonald, 1984) specification of a structural equation model, and should be of class \code{mod}. The model is specified in the same manner as does the \code{\link[sem]{sem}} package; see \code{\link[sem]{sem}} and \code{\link[sem]{specify.model}} for detailed documentations about model specifications in the RAM notation. }
  \item{theta}{ a vector containing the model parameters. The names of the elements in \code{theta} must be the same as the names of the model parameters specified in \code{model}. }
  \item{latent.vars}{ a vector containing the names of the latent variables }
}
\details{
Part of the codes in this function are adapted from the function \code{\link[sem]{sem}} in the \code{sem} R package (Fox, 2006). This function uses the same notation to specify SEM models as does \code{\link[sem]{sem}}. Please refer to \code{\link[sem]{sem}} and the example below for more detailed documentation about model specification and the RAM notation. For technical discussion on how to obtain the model implied covariance matrix in the RAM notation given model parameters, see McArdle and McDonald (1984).}
\value{
  \item{ram }{RAM matrix, including any rows generated for covariances among fixed exogenous variables; column 5 includes computed start values.}
  \item{t}{number of model parameters (i.e., the length of \code{theta})}
  \item{m}{total number of variables (i.e., manifest variables plus latent variables)}
  \item{n}{number of observed variables}
  \item{all.vars}{the names of all variables (i.e., manifest plus latent)}
  \item{obs.vars}{the names of observed variables}
  \item{latent.vars}{the names of latent variables}
  \item{pars}{the names of model parameters}
  \item{P}{the \emph{P} matrix in RAM notation}
  \item{A}{the \emph{A} matrix in RAM notation}
  \item{Sigma.theta}{the model implied covariance matrix}

}
\references{ 
Fox, J. (2006). Structural equation modeling with the sem package in R. \emph{Structural Equation Modeling, 13}, 465--486.

Lai, K., & Kelley, K. (in press). Accuracy in parameter estimation for targeted effects in structural equation modeling: Sample size planning for narrow confidence intervals. \emph{Psychological Methods}.

McArdle, J. J., & McDonald, R. P. (1984). Some algebraic properties of the reticular action model. \emph{British Journal of Mathematical and Statistical Psychology, 37}, 234--251. 
}

\author{Keke Lai (University of California--Merced) }
%\note{ }

\seealso{\code{\link[sem]{sem}}; \code{\link[sem]{specify.model}}  }
\examples{
\dontrun{
# to obtain the model implied covariance matrix of Model 2 in the simulation 
# study in Lai and Kelley (2010), one can use the present function in the 
# following manner.

library(sem)

# specify a model object in the RAM notation
model.2<-specify.model()
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
}
}

\keyword{ multivariate}
