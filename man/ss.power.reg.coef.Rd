\name{ss.power.reg.coef}
\alias{ss.power.reg.coef}
\title{sample size for a targeted regression coefficient}
\description{
Determine the necessary sample size for a targeted regression coefficient or determine the degree of power given a specified sample size}
\usage{
ss.power.reg.coef(Rho2.Y_X = NULL, Rho2.Y_X.without.j = NULL, p = NULL, 
desired.power = 0.85, alpha.level = 0.05, Directional = FALSE, 
beta.j = NULL, sigma.X = NULL, sigma.Y = NULL, Rho2.j_X.without.j = NULL,
RHO.XX = NULL, Rho.YX = NULL, which.predictor = NULL, Cohen.f2 = NULL, 
Specified.N=NULL, Print.Progress = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Rho2.Y_X}{population squared multiple correlation coefficient predicting the dependent variable (i.e., \emph{Y}) from the \code{p} predictor variables (i.e., the \emph{X} variables)}
  \item{Rho2.Y_X.without.j}{population squared multiple correlation coefficient predicting the dependent variable (i.e., \emph{Y}) from the \code{p}-1 predictor variables, where the one not used is the predictor of interest}
  \item{p}{number of predictor variables}
  \item{desired.power}{desired degree of statistical power for the test of targeted regression coefficient}
  \item{alpha.level}{Type I error rate}
  \item{Directional}{whether or not a direction or a nondirectional test is to be used (usually \code{directional=FALSE})}
  \item{beta.j}{population value of the regression coefficient for the predictor of interest}
  \item{sigma.X}{population standard deviation for the predictor variable of interest}
  \item{sigma.Y}{population standard deviation for the outcome variable}
  \item{Rho2.j_X.without.j}{population squared multiple correlation coefficient predicting the predictor variable of interest from the remaining p-1 predictor variables}
  \item{RHO.XX}{population correlation matrix for the \code{p} predictor variables}
  \item{Rho.YX}{population vector of correlation coefficient between the \code{p} predictor variables and the criterion variable}
  \item{Cohen.f2}{Cohen's (1988) definition for an effect size for a targeted regression coefficient: \code{(Rho2.Y_X-Rho2.Y_X.without.j)/(1-Rho2.Y_X)}}
  \item{which.predictor}{identifies the predictor of interest when \code{RHO.XX} and \code{Rho.YX} are specified}
  \item{Specified.N}{sample size for which power should be evaluated}
  \item{Print.Progress}{if the progress of the iterative procedure is printed to the screen as the iterations are occurring}
}
\details{
Determines the necessary sample size given a desired level of statistical power. Alternatively, determines the statistical power for a given a specified sample size.
There are a number of ways that the specification regarding the size of the regression coefficient can be entered. The most basic, and often the simplest, is to specify \code{Rho2.Y_X} and \code{Rho2.Y_X.without.j}. See the examples section 
for several options.
}

\value{
  \item{Sample.Size}{either the necessary sample size or the specified sample size, depending if one is interested in determining the necessary sample size given a desired degree of statistical power or if one is interested in the determining the value of statistical power given a specified sample size, respectively}
  \item{Actual.Power}{Actual power of the situation described}
  \item{Noncentral.t.Parm}{value of the noncentral distribution for the appropriate \emph{t}-distribution}
   \item{Effect.Size.NC.t}{effect size for the noncentral \emph{t}-distribution; this is the square root of \code{Cohen.f2}, because \code{Cohen.f2} is the effect size using an \emph{F}-distribution}
}
\references{
Cohen, J. (1988). \emph{Statistical power analysis for the behavioral sciences} (2nd ed.). Hillsdale, NJ: Erlbaum.

Kelley, K. & Maxwell, S. E. (2008). Sample Size Planning with applications to multiple regression: Power and accuracy for omnibus and targeted effects. In P. Alasuuta, J. Brannen, & L. Bickman (Eds.), \emph{The Sage handbook of social research methods} (pp. 166--192). Newbury Park, CA: Sage.

Maxwell, S. E. (2000). Sample size for multiple regression. \emph{Psychological Methods, 4}, 434--458.

}
\author{Ken Kelley (University of Notre Dame; \email{KKelley@ND.Edu})}
\seealso{\code{ss.aipe.reg.coef}, \code{ss.power.R2}, \code{conf.limits.ncf}}
\examples{
Cor.Mat <- rbind(
c(1.00, 0.53,  0.58, 0.60, 0.46, 0.66),
c(0.53, 1.00,  0.35, 0.07, 0.14, 0.43),
c(0.58, 0.35,  1.00, 0.18, 0.29, 0.50),
c(0.60, 0.07,  0.18, 1.00, 0.30, 0.26),
c(0.46, 0.14,  0.29, 0.30, 1.00, 0.30),
c(0.66, 0.43,  0.50, 0.26, 0.30, 1.00))

RHO.XX <- Cor.Mat[2:6,2:6]
Rho.YX <- Cor.Mat[1,2:6]

# Method 1
# ss.power.reg.coef(Rho2.Y_X=0.7826786, Rho2.Y_X.without.j=0.7363697, p=5,
# alpha.level=.05, Directional=FALSE, desired.power=.80)

# Method 2
# ss.power.reg.coef(alpha.level=.05, RHO.XX=RHO.XX, Rho.YX=Rho.YX, 
# which.predictor=5,
# Directional=FALSE, desired.power=.80)

# Method 3
# Here, beta.j is the standardized regression coefficient. Had beta.j 
# been the unstandardized regression coefficient, sigma.X and sigma.Y 
# would have been the standard deviation for the 
# X variable of interest and Y, respectively.
# ss.power.reg.coef(Rho2.Y_X=0.7826786, Rho2.j_X.without.j=0.3652136, 
# beta.j=0.2700964,
# p=5, alpha.level=.05,  sigma.X=1, sigma.Y=1, Directional=FALSE, 
# desired.power=.80)

# Method 4
# ss.power.reg.coef(alpha.level=.05, Cohen.f2=0.2130898, p=5, 
# Directional=FALSE,
# desired.power=.80)

# Power given a specified N and squared multiple correlation coefficients.
# ss.power.reg.coef(Rho2.Y_X=0.7826786, Rho2.Y_X.without.j=0.7363697, 
# Specified.N=25,
# p=5, alpha.level=.05, Directional=FALSE)

# Power given a specified N and effect size.
# ss.power.reg.coef(alpha.level=.05, Cohen.f2=0.2130898, p=5, Specified.N=25,
# Directional=FALSE)

# Reproducing Maxwell's (2000, p. 445) Example
Cor.Mat.Maxwell <- rbind(
c(1.00, 0.35,  0.20, 0.20, 0.20, 0.20),
c(0.35, 1.00,  0.40, 0.40, 0.40, 0.40),
c(0.20, 0.40,  1.00, 0.45, 0.45, 0.45),
c(0.20, 0.40,  0.45, 1.00, 0.45, 0.45),
c(0.20, 0.40,  0.45, 0.45, 1.00, 0.45),
c(0.20, 0.40,  0.45, 0.45, 0.45, 1.00))

RHO.XX.Maxwell <- Cor.Mat.Maxwell[2:6,2:6]
Rho.YX.Maxwell <- Cor.Mat.Maxwell[1,2:6]
R2.Maxwell <- Rho.YX.Maxwell%*%solve(RHO.XX.Maxwell)%*%Rho.YX.Maxwell

RHO.XX.Maxwell.no.1 <- Cor.Mat.Maxwell[3:6,3:6]
Rho.YX.Maxwell.no.1 <- Cor.Mat.Maxwell[1,3:6]
R2.Maxwell.no.1 <- 
Rho.YX.Maxwell.no.1%*%solve(RHO.XX.Maxwell.no.1)%*%Rho.YX.Maxwell.no.1


# Note that Maxwell arrives at N=113, whereas this procedure arrives at 111.
# This seems to be the case becuase of rounding error in calculations 
# in Cohen (1988)'s tables. The present procedure is correct and contains no 
# rounding error
# in the application of the method.
# ss.power.reg.coef(Rho2.Y_X=R2.Maxwell, 
# Rho2.Y_X.without.j=R2.Maxwell.no.1, p=5,
# alpha.level=.05, Directional=FALSE, desired.power=.80)
}
\keyword{design}
