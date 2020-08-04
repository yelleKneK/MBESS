\name{Sigma.2.SigmaStar}
\alias{Sigma.2.SigmaStar}

\title{Construct a covariance matrix with specified error of approximation}
\description{
This function implements Cudeck & Browne's (1992) method to construct a covariance matrix in the structural equation modeling (SEM) context. Given an SEM model and its model parameters, a covariance matrix is obtained so that (a) the population discrepancy due to approximation equals a certain specified value; and (b) the population model parameter vector is the minimizer of the discrepancy function. 
}
\usage{
Sigma.2.SigmaStar(model, model.par, latent.var, discrep, ML = TRUE)
}

\arguments{
  \item{model}{ an RAM (reticular action model; e.g., McArdle & McDonald, 1984) specification of a structural equation model, and should be of class \code{mod}. The model is specified in the same manner as does the \code{\link[sem]{sem}} package; see \code{\link[sem]{sem}} and \code{\link[sem]{specify.model}} for detailed documentations about model specifications in the RAM notation. }

  \item{model.par}{ a vector containing the model parameters. The names of the elements in \code{theta} must be the same as the names of the model parameters specified in \code{model}. }
  \item{latent.var}{ a vector containing the names of the latent variables }
  \item{discrep}{ the desired discrepancy function minimum value }
  \item{ML}{ the discrepancy function to be used, if \code{ML=TRUE} then the discrepancy function is based on normal theory maximum likelihood}
}

\details{
This function constructs a covariance matrix \eqn{ \Sigma^{*} } such that \eqn{ \Sigma^{*} = \Sigma( \theta ) + E }, where \eqn{ \Sigma(\theta)} is the population model-implied covariance matrix, and \eqn{E} is a matrix containing the errors due to approximation. The matrix \eqn{E} is chosen so that the discrepancy function \eqn{F( \Sigma^{*}, \Sigma (\theta) ) } has the specified discrepancy value. 

This function uses the same notation to specify SEM models as does \code{\link[sem]{sem}}. Please refer to \code{\link[sem]{sem}} for more detailed documentation about model specification and the RAM notation. For technical discussion on how to obtain the model implied covariance matrix in the RAM notation given model parameters, see McArdle and McDonald (1984).
}

\value{  
  \item{Sigma.star}{the population covariance matrix of manifest variables}
  \item{Sigma_theta}{the population model-implied covariance matrix}
  \item{E}{the matrix containing the population errors of approximation,
   i.e., \code{Sigma.star} - \code{Sigma_theta}}
}

\references{ 
Cudeck, R., & Browne, M. W. (1992). Constructing a covariance matrix that yields a specified minimizer and a specified minimum discrepancy function value. \emph{Psychometrika, 57}, 357--369. 

Fox, J. (2006). Structural equation modeling with the sem package in R. \emph{Structural Equation Modeling, 13}, 465--486.

McArdle, J. J., & McDonald, R. P. (1984). Some algebraic properties of the reticular action model. \emph{British Journal of Mathematical and Statistical Psychology, 37}, 234--251.
}

\author{Keke Lai (University of California-Merced) }
%\note{ }

\seealso{\code{\link[sem]{sem}}; \code{\link[sem]{specify.model}}; \code{\link{theta.2.Sigma.theta}} }

\examples{
\dontrun{
library(sem)

###############
## EXAMPLE 1; a CFA model with three latent variables and nine indicators.
###############

# To specify the model
model.cfa<-specify.model()
xi1 -> x1, lambda1, 0.6
xi1 -> x2, lambda2, 0.7
xi1 -> x3, lambda3, 0.8
xi2 -> x4, lambda4, 0.65
xi2 -> x5, lambda5, 0.75
xi2 -> x6, lambda6, 0.85
xi3 -> x7, lambda7, 0.5
xi3 -> x8, lambda8, 0.7
xi3 -> x9, lambda9, 0.9
xi1 <-> xi1, NA, 1
xi2 <-> xi2, NA, 1
xi3 <-> xi3, NA, 1
xi1 <-> xi2, phi21, 0.5
xi1 <-> xi3, phi31, 0.4
xi2 <-> xi3, phi32, 0.6
x1 <-> x1, delta11, 0.36
x2 <-> x2, delta22, 0.5
x3 <-> x3, delta33, 0.9
x4 <-> x4, delta44, 0.4
x5 <-> x5, delta55, 0.5
x6 <-> x6, delta66, 0.6
x7 <-> x7, delta77, 0.6
x8 <-> x8, delta88, 0.7
x9 <-> x9, delta99, 0.7


# To specify model parameters
theta <- c(0.6, 0.7, 0.8,
0.65, 0.75, 0.85,
0.5, 0.7, 0.9,
0.5, 0.4, 0.6,
0.8, 0.6, 0.5,
0.6, 0.5, 0.4,
0.7, 0.7, 0.6)

names(theta) <- c("lambda1", "lambda2", "lambda3", 
"lambda4","lambda5", "lambda6", 
"lambda7", "lambda8","lambda9",
"phi21", "phi31", "phi32", 
"delta11", "delta22","delta33",
"delta44", "delta55","delta66",
"delta77", "delta88","delta99")

res.matrix <- Sigma.2.SigmaStar(model=model.cfa, model.par=theta, 
latent.var=c("xi1", "xi2", "xi3"), discrep=0.06)

# res.matrix

# To verify the returned covariance matrix; the model chi-square
# should be equal to (N-1) times the specified discrepancy value.
# Also the "point estimates" of model parameters should be 
# equal to the specified model parameters

# res.sem<-sem(model.cfa, res.matrix$Sigma.star, 1001)
# summary(res.sem)

# To construct a covariance matrix so that the model has
# a desired population RMSEA value, one can transform the RMSEA
# value to the discrepancy value

res.matrix <- Sigma.2.SigmaStar(model=model.cfa, model.par=theta, 
latent.var=c("xi1", "xi2", "xi3"), discrep=0.075*0.075*24)

# To verify the population RMSEA value
# res.sem<-sem(model.cfa, res.matrix$Sigma.star, 1000000)
# summary(res.sem)

###############
## EXAMPLE 2; an SEM model with five latent variables
###############

model.5f <- specify.model()
eta1 -> y4, NA, 1
eta1 -> y5, lambda5, NA
eta2 -> y1, NA, 1
eta2 -> y2, lambda2, NA 
eta2 -> y3, lambda3, NA
xi1 -> x1, NA, 1
xi1 -> x2, lambda6, NA 
xi1 -> x3, lambda7, NA
xi2 -> x4, NA, 1
xi2 -> x5, lambda8, NA 
xi3 -> x6, NA, 1
xi3 -> x7, lambda9, NA 
xi3 -> x8, lambda10, NA
xi1 -> eta1, gamma11, NA
xi2 -> eta1, gamma12, NA
xi3 -> eta1, gamma13, NA
xi3 -> eta2, gamma23, NA
eta1 -> eta2, beta21, NA
xi1 <-> xi2, phi21, NA
xi1 <-> xi3, phi31, NA
xi3 <-> xi2, phi32, NA
xi1 <-> xi1, phi11, NA
xi2 <-> xi2, phi22, NA
xi3 <-> xi3, phi33, NA
eta1 <-> eta1, psi11, NA
eta2 <-> eta2, psi22, NA
y1 <-> y1, eplison11, NA
y2 <-> y2, eplison22, NA
y3 <-> y3, eplison33, NA
y4 <-> y4, eplison44, NA
y5 <-> y5, eplison55, NA
x1 <-> x1, delta11, NA
x2 <-> x2, delta22, NA
x3 <-> x3, delta33, NA
x4 <-> x4, delta44, NA
x5 <-> x5, delta55, NA
x6 <-> x6, delta66, NA
x7 <-> x7, delta77, NA
x8 <-> x8, delta88, NA


theta <- c(0.84, 0.8, 0.9, 
1.26, 0.75, 1.43, 1.58, 0.83, 
0.4, 0.98, 0.52, 0.6,0.47, 
0.12, 0.14, 0.07,
0.44, 0.22, 0.25, 
0.3, 0.47, 
0.37, 0.5, 0.4, 0.4, 0.58, 
0.56,0.3, 0.6, 0.77, 0.54, 0.75, 0.37, 0.6)

names(theta) <- c(
"lambda5","lambda2","lambda3",
"lambda6","lambda7","lambda8","lambda9","lambda10" , 
"gamma11",  "gamma12","gamma13" ,  "gamma23" ,  "beta21",
"phi21","phi31", "phi32", 
"phi11","phi22",  "phi33",     
"psi11" ,    "psi22"   ,  
"eplison11","eplison22" ,"eplison33", "eplison44" ,"eplison55", 
  "delta11"  , "delta22" ,  "delta33" ,  "delta44" ,  "delta55" ,  "delta66",  
"delta77" , "delta88")

# To construct a covariance matrix so that the model has 
# a population RMSEA of 0.08

res.matrix <- Sigma.2.SigmaStar(model=model.5f, model.par=theta, 
latent.var=c("xi1", "xi2", "xi3", "eta1","eta2"), discrep=0.08*0.08*57)

# To verify
# res.sem<- sem(model.5f, res.matrix$Sigma.star, 1000000)
# summary(res.sem)
}
}

\keyword{ multivariate }