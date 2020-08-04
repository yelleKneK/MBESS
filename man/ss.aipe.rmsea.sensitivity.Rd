\name{ss.aipe.rmsea.sensitivity}
\Rdversion{1.1}
\alias{ss.aipe.rmsea.sensitivity}
\title{a priori Monte Carlo simulation for sample size planning for RMSEA in SEM
}
\description{Conduct a priori Monte Carlo simulation to empirically study the effects of (mis)specifications of input information on the calculated sample size. The sample size is planned so that the expected width of a confidence interval for the population RMSEA is no larger than desired. Random data are generated from the true covariance matrix but fit to the proposed model, whereas sample size is calculated based on the input covariance matrix and proposed model.
}

\usage{
ss.aipe.rmsea.sensitivity(width, model, Sigma, N=NULL, 
conf.level=0.95, G=200, save.file="sim.results.txt", ...)
}

\arguments{
  \item{width}{desired confidence interval width for the model parameter of interest }
  \item{model}{ the model the researcher proposes, may or may not be the true model. This argument should be an RAM (reticular action model; e.g., McArdle & McDonald, 1984) specification of a structural equation model, and should be of class \code{mod}. The model is specified in the same manner as does the \code{\link[sem]{sem}} package; see \code{\link[sem]{sem}} and \code{\link[sem]{specify.model}} for detailed documentation about model specifications in the RAM notation.}
 
  \item{Sigma}{the true population covariance matrix, which will be used to generate random data for the simulation study. The row names and column names of \code{Sigma} should be the same as the manifest variables in \code{model}.
}

  \item{N}{if \code{N} is specified, random sample of the specified \code{N} size will be generated. Otherwise the sample size is calculated with the sample size planning method with the goal that the expected width of a confidence interval for population RMSEA is no larger than desired.
}
  \item{conf.level}{ confidence level (i.e., 1- Type I error rate)}
  \item{G}{number of replications in the Monte Carlo simulation }
  \item{save.file}{the name of the file that simulation results will be saved to}
  \item{\dots}{  allows one to potentially include parameter values for inner functions}

}
\details{This function implements the sample size planning methods proposed in Kelley and Lai (2010). It depends on the 
function \code{\link[sem]{sem}} in the \code{sem} package to fit the proposed model to random data, and uses the same notation to specify SEM 
models as does \code{\link[sem]{sem}}. Please refer to \code{\link[sem]{sem}} for more detailed documentation 
about model specifications, the RAM notation, and model fitting techniques. For technical discussion 
on how to obtain the model implied covariance matrix in the RAM notation given model parameters, see McArdle and McDonald (1984)
}
\value{
\item{successful.replication}{the number of successful replications}
\item{w}{the \code{G} random confidence interval widths }
\item{RMSEA.hat}{the \code{G} estimated RMSEA values based on the \code{G} random samples }
\item{sample.size}{the sample size calculated }
\item{df}{degrees of freedom of the proposed model}
\item{RMSEA.pop}{the input RMSEA value that is used to calculated the necessary sample size}
\item{desired.width}{desired confidence interval width}
\item{mean.width}{mean of the random confidence interval widths}
\item{median.width}{median of the random confidence interval widths}
\item{assurance}{the proportion of confidence interval widths narrower than desired} 
\item{quantile.width}{99, 97, 95, 90, 80, 70, and 60 percentiles of the random confidence interval widths}
\item{alpha.upper}{the upper empirical Type I error rate} 
\item{alpha.lower}{the lower empirical Type I error rate}
\item{alpha}{total empirical Type I error rate}  
\item{conf.level}{confidence level}
\item{sim.results.txt}{a text file that saves the simulation results; it updates after each replication. 'sim.results.txt' is the default file name}
}

\references{
Cudeck, R., & Browne, M. W. (1992). Constructing a covariance matrix that yields a specified minimizer and a specified minimum discrepancy function value. \emph{Psychometrika, 57}, 357--369. 

Fox, J. (2006). Structural equation modeling with the sem package in R. \emph{Structural Equation Modeling, 13}, 465--486.

Kelley, K., & Lai, K. (2010). Accuracy in parameter estimation for the root mean square of approximation: Sample size planning for narrow confidence intervals. \emph{Manuscript under review}.

McArdle, J. J., & McDonald, R. P. (1984). Some algebraic properties of the reticular action model. \emph{British Journal of Mathematical and Statistical Psychology, 37}, 234--251.
}
\author{
Keke Lai (University of California -- Merced) and Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})
}
\note{Sometimes this function jumps out of the loop before it finishes the simulation. The reason is because the 
 \code{\link[sem]{sem}} function that this function calls to fit the model fails to converge when searching for
 maximum likelihood estimates of model parameters. Since the results in previous replications are saved, the
 user can start this function again, and specify the number of replications (i.e., \code{G}) to be the desired
 total number of replications minus the number of previous successful replications.
}


\seealso{
\code{\link[sem]{sem}}; \code{\link[sem]{specify.model}}; \code{\link{ss.aipe.rmsea}}; \code{\link{theta.2.Sigma.theta}}; \code{\link{Sigma.2.SigmaStar}}
}

\examples{
\dontrun{
#########################
EXAMPLE 1 
######################### 
# To replicate the simulation in the first panel, second column of 
# Table 2 (i.e., population RMSEA=0.0268, df=23, desired width=0.02) 
# in Lai and Kelley (2010), the following steps can be used.

## STEP 1: Obtain the (correct) population covariance matrix implied by Model 2
# This requires the model and its population model parameter values.
library(MASS)
library(sem)

# Specify Model 2 in the RAM notation
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


# To inspect the specified model
model.2

# Specify model parameter values
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
# Then 'Sigma.theta' is the (true) population covariance matrix

## STEP 2: Create a misspecified model
# The following model is misspecified in the same way as did Lai and Kelley (2010)
# with the goal to obtain a relatively small population RMSEA

model.2.mis<-specifyModel()
xi1 -> y1, lambda1, 1
xi1 -> y2, NA, 1
xi1 -> y3, lambda2, 1
xi1 -> y4, lambda3, 0.3
eta1 -> y4, lambda4, 1
eta1 -> y5, NA, 1
eta1 -> y6, lambda5, 0.96
eta2 -> y6, lambda7, 0.33
eta2 -> y7, lambda8, 1.33
eta2 -> y8, NA, 1
eta2 -> y9, lambda9, 1
xi1 -> eta1, gamma11, 0.6
eta1 -> eta2, beta21, 0.65 
xi1 <-> xi1, phi11, 0.49
eta1 <-> eta1, psi11, 0.3136
eta2 <-> eta2, psi22, 0.23
y1 <-> y1, delta1, 0.51
y2 <-> y2, delta2, 0.51
y3 <-> y3, delta3, 0.51
y4 <-> y4, delta4, 0.2895
y5 <-> y5, delta5, 0.51
y6 <-> y6, delta6, 0.29
y7 <-> y7, delta7, 0.22
y8 <-> y8, delta8, 0.56
y9 <-> y9, delta9, 0.56


# To verify the population RMSEA of this misspecified model
fit<-sem(ram=model.2.mis, S=Sigma.theta, N=1000000)
summary(fit)$RMSEA

## STEP 3: Conduct the simulation
# The number of replications is set to a very small value just to demonstrate
# and save time. Real simulation studies require a larger number (e.g., 500, 1,000)

ss.aipe.rmsea.sensitivity(width=0.02, model=model.2.mis, Sigma=Sigma.theta, G=10)

## STEP 3+: In cases where this function stops before it finishes the simulation
# Suppose it stops at the 7th replication. The text 
# file "results_ss.aipe.rmsea.sensitivity.txt" saves the results in all 
# previous replications; in this case it contains 6 replications since
# the simulation stopped at the 7th. The user can start this function again and specify
# 'G' to 4 (i.e., 10-6). New results will be appended to previous ones in the same file.

ss.aipe.rmsea.sensitivity(width=0.02, model=model.2.mis, Sigma=Sigma.theta, G=4)

########################################
EXAMPLE 2
########################################
# In addition to create a misspecified model by changing the model
# parameters in the true model as does Example 1, a misspecified
# model can also be created with the Cudeck-Browne (1992) procedure. 
# This procedure is implemented in the 'Sigma.2.SigmaStar( )' function in
# the MBESS package. Please refer to the help file of 'Sigma.2.SigmaStar( )'
# for detailed documentation.

## STEP 1: Specify the model
# This model is the same as the model in the first step of Example 1, but the
# model-implied population covariance matrix is no longer the true population 
# covariance matrix. The true population covariance matrix will be generated
# in Step 2 with the Cudeck-Browne procedure.
library(MASS)
library(sem)

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



theta <- c(1, 1, 0.3, 1,1, 0.3, 0.3, 1, 1, 0.6, 0.6,
0.49, 0.3136, 0.3136, 0.51, 0.51, 0.51, 0.2895, 0.51, 0.2895, 0.2895, 0.51, 0.51)

names(theta) <- c("lambda1","lambda2","lambda3",
"lambda4","lambda5","lambda6","lambda7","lambda8","lambda9",
"gamma11", "beta21",
"phi11", "psi11", "psi22", 
"delta1","delta2","delta3","delta4","delta5","delta6","delta7",
"delta8","delta9")

## STEP 2: Create the true population covariance matrix, so that (a) the model fits
# to this covariance matrix with specified discrepancy; (b) the population model
# parameters (the object 'theta') is the minimizer in fitting the model to the true 
# population covariance matrix.

# Since the desired RMSEA is 0.0268 and the df is 22, the MLE discrepancy value
# is specified to be 22*0.0268*0.0268, given the definition of RMSEA.

res <- Sigma.2.SigmaStar(model=model.2, model.par=theta, 
latent.var=c("xi1", "eta1", "eta2"), discrep=22*0.0268*0.0268)

Sigma.theta.star <- res$Sigma.star

# To verify that the population RMSEA is 0.0268
res2 <- sem(ram=model.2, S=Sigma.theta.star, N=1000000)
summary(res2)$RMSEA

## STEP 3: Conduct the simulation
# Note although Examples 1 and 2 have the same population RMSEA, the
# model df and true population covariance matrix are different. Example 1
# uses 'model.2.mis' and 'Sigma.theta', whereas Example 2 uses 'model.2'
# and 'Sigma.theta.star'. Since the df is different, it requires a different sample
# size to achieve the same desired confidence interval width. 
ss.aipe.rmsea.sensitivity(width=0.02, model=model.2, Sigma=Sigma.theta.star, G=10)
}
}

\keyword{ design }
\keyword{ multivariate }