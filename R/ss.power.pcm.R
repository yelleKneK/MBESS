ss.power.pcm <- function (beta, tau, level.1.variance, frequency, duration, desired.power = NULL, 
N = NULL, alpha.level = 0.05, standardized = TRUE, directional = FALSE) 
{

# Function currently requires that power and sample size is for the linear coefficient of change only
p <- 1
# if(p==0) K.p <- 1
# if(p==1) K.p <- 1/12
# if(p==2) K.p <- 1/720
# if(p==3) K.p <- 1/100800
K.p <- 1/12

# left column, lower, part of p. 392
M <- frequency * duration + 1

if(standardized == FALSE) 
{
beta.standardized <- beta/sqrt(tau) # Equation 13
beta.unstandardized <- beta
}

if(standardized == TRUE)
{
beta.standardized <- beta
beta.unstandardized <- beta*sqrt(tau)
}


if(is.null(N))
{
N.i <- 4
Dif <- 1
while (Dif > 0)
{
N.i <- N.i + 2
V <- (level.1.variance*(frequency^(2*p))*(factorial(M-p-1)))/(K.p*factorial(M+p)) # Can be generalized beyond straight-line model
var.beta <- 4 * (tau + V)/N.i # Equation 10 
reliability <- tau/(tau + V) # Equation 15

if(standardized == FALSE) 
{
Lambda <- (beta^2)/var.beta # Equation 12
}

if(standardized == TRUE)
{
Lambda <- (N.i * beta^2 * reliability)/4 # Equation 14
}

if(directional == FALSE)
{
CV.for.test.of.Null <- qt((1 - alpha.level/2), df = (N.i - 2), lower.tail = TRUE, log.p = FALSE)
}
if (directional == TRUE)
{
CV.for.test.of.Null <- qt((1 - alpha.level), df = (N.i - 2), lower.tail = TRUE, log.p = FALSE)
}
Actual.Power <- 1 - pt(CV.for.test.of.Null, df = (N.i-2), ncp = sqrt(Lambda), lower.tail = TRUE, log.p = FALSE)
Dif <- desired.power - Actual.Power
}
}


if(is.null(desired.power))
{
V <- (level.1.variance*(frequency^(2*p))*(factorial(M-p-1)))/(K.p*factorial(M+p)) # Can be generalized beyond straight-line model

var.beta <- 4 * (tau + V)/N
reliability <- tau/(tau + V) # Equation 15
if (standardized == FALSE)
{
Lambda <- (beta^2)/var.beta # Equation 12
}
if(standardized == TRUE)
{
Lambda <- (N * beta^2 * reliability)/4
}
if (directional == FALSE)
{
CV.for.test.of.Null <- qt((1 - alpha.level/2), df = (N - 2), lower.tail = TRUE, log.p = FALSE)
}
if (directional == TRUE)
{
CV.for.test.of.Null <- qt((1 - alpha.level), df = (N - 2), lower.tail = TRUE, log.p = FALSE)
}

Actual.Power <- 1 - pt(CV.for.test.of.Null, df = (N - 2), ncp = sqrt(Lambda), lower.tail = TRUE, log.p = FALSE)
N.i <- N
}

return(
list(
Design.features = list(Necessary.SS.Control = N.i/2, Necessary.SS.Treatment = N.i/2, Total.SS = N.i, Actual.Power = Actual.Power, Frequency = frequency, Duration = duration, Total.Measurement.Occasions = M), 
Parameters = list(Unstandardized.Regression.Coefficient = beta.unstandardized, Standardized.Regression.Coefficient = beta.standardized, Level.1.error.variance = level.1.variance, true.variance.of.slopes = tau, 
error.variance.of.slopes = V, Reliability = reliability, Noncentral.t.parameter = sqrt(Lambda)))
)
}


