\name{var.ete}
\alias{var.ete}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The Variance of the Estimated Treatment Effect at Selected Covariate Values in a Two-group ANCOVA.
}
\description{
Calculate the variance or an estimated variance of the estimated treatment effect at selected covariate values assuming heterogeneity of regression and a random covariate in a two-group ANCOVA. 
}
\usage{
var.ete(sigma2, sigmaz2, n1, n2, beta1, beta2, muz = 0, c = 0, type = "sample", 
covariate.value = "sample.mean")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sigma2}{
Variance of the residual errors if 'type = population' and sample variance of the residual errors if 'type = sample'}
  \item{sigmaz2}{
Variance of the random covariate if 'type = population' and sample variance of the random covariate if 'type = sample'}
  \item{n1}{
Sample size of group 1}
  \item{n2}{
Sample size of group 2}
  \item{beta1}{
Slope of the random covariate for group 1 if 'type = population' and estimated slope of the random covariate for group 1 if 'type = sample'}
  \item{beta2}{
Slope of the random covariate for group 2 if 'type = population' and estimated slope of the random covariate for group 2 if 'type = sample'}
  \item{muz}{
Population mean of the random covariate if 'type = population' and sample mean of the random covariate if 'type = sample'}
  \item{c}{
Fixed value where the treatment effect is assessed}
  \item{type}{
The type of variance formula: 'population' refers to the variance of the estimated treatment effect using population slopes and variances; 'sample'refers to an unbiased estimate of the variance using sample slopes and variances}
  \item{covariate.value}{
The covariate value is chosen at the sample grand mean if 'covariate.value = sample.mean', at the sample grand mean plus or minus one sample standard deviation if 'covariate.value = SD', and at a fixed value if 'covariate.value = fixed'}
}

\value{
The function yields the variance of the estimated treatment effect for the specified input values. 
}
\references{
Maxwell, S. E., Delaney, H. D., & Kelley, K. (2018).  \emph{Designing experiments and analyzing data:  A model comparison perspective}.  New York:  Routledge.

Li, L., McLouth, C. J., and Delaney, H. D. (submitted). Analysis of Covariance with Heterogeneity of Regression and a Random Covariate: The Variance of the Estimated Treatment Effect at Selected Covariate Values.
}
\author{
Li Li (University of New Mexico; \email{llis@unm.edu})
}

\examples{
# Pygmalion in the Classroom: Teacher Expectation and Pupils' Intellectual Development. 
# This dataset has been used to illustrate heterogeneity of regression 
# by Maxwell, Delaney, and Kelley (2018).
nA <- 64
nB <- 246
muz <- 0
sigma2 <- 175.3251
sigmaz2 <- 348.9099
betaA <- 0.96895
betaB <- 0.77799
var.ete(sigma2=sigma2, sigmaz2=sigmaz2, n1=nA, n2=nB, beta1=betaA, beta2=betaB, 
type="sample", covariate.value = "sample.mean")
var.ete(sigma2=sigma2, sigmaz2=sigmaz2, n1=nA, n2=nB, beta1=betaA, beta2=betaB, 
type="sample", covariate.value = "SD")
var.ete(sigma2=sigma2, sigmaz2=sigmaz2, n1=nA, n2=nB, beta1=betaA, beta2=betaB, 
c = 4.2631, muz=muz, type="sample",covariate.value = "fixed")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
\keyword{design}% __ONLY ONE__ keyword per line
