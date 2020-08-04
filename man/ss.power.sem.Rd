\name{ss.power.sem}
\alias{ss.power.sem}

\title{
Sample size planning for structural equation modeling from the power analysis perspective
}
\description{
Calculate the necessary sample size for an SEM study, so as to have enough power to reject the 
null hypothesis that (a) the model has perfect fit, or (b) the difference in fit between two nested models equal some
specified amount.
}
\usage{
ss.power.sem(F.ML = NULL, df = NULL, RMSEA.null = NULL, RMSEA.true = NULL, 
F.full = NULL, F.res = NULL, RMSEA.full = NULL, RMSEA.res = NULL, 
df.full = NULL, df.res = NULL, alpha = 0.05, power = 0.8)
}
\arguments{
  \item{F.ML}{
The true maximum likelihood fit function value in the population for the model of interest. Leave this argument NULL
if you are doing nested model significance tests. 
}
  \item{df}{
The degrees of freedom of the model of interest. Leave this argument NULL
if you are doing nested model significance tests. 
}
  \item{RMSEA.null}{
The model's population RMSEA under the null hypothesis. Leave this argument NULL
if you are doing nested model significance tests. 
}
  \item{RMSEA.true}{
The model's population RMSEA under the alternative hypothesis. This should be the model's 
true population RMSEA value. Leave this argument NULL
if you are doing nested model significance tests. 
}
  \item{F.full}{
The maximum likelihood fit function value for the full model.
}
  \item{F.res}{
The maximum likelihood fit function value for the restricted model.
}
  \item{RMSEA.full}{
The population RMSEA value for the full model.
}
  \item{RMSEA.res}{
The population RMSEA value for the restricted model.
}
  \item{df.full}{
The degrees of freedom for the full model. 
}
  \item{df.res}{
The degrees of freedom for the restricted model.
}
  \item{alpha}{
The Type I error rate.
}
  \item{power}{
The desired power.
}
}
%\details{}
%\value{}
%\references{ }
\author{
Keke Lai (University of California - Merced)
}
%\note{}



%\seealso{}
%\examples{ }

\keyword{ design }
\keyword{ multivariate }
