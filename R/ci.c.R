ci.c <- function(means=NULL, s.anova=NULL, c.weights=NULL, n=NULL, N=NULL,
Psi=NULL, conf.level=.95, alpha.lower=NULL, alpha.upper=NULL, df.error=NULL, ...)
{


if(length(n)==1) 
{
print("You specified only one group sample size, the assumption used here is that all gruops have the same sample size.")
n <- rep(n, length(means))
}

# Checks
##########################################
if(is.null(c.weights)) stop("You must specify the vector of contrast weights ('c.weights').")
if(length(n)!=length(c.weights)) stop("The lengths of 'n' and 'c.weights' differ, which should not be the case.")
    
if(identical(sum(c.weights > 1), 0)) stop("Please use fractional values to specify the contrast weights (i.e., no contrast weights above 1.")
if(identical(sum(c.weights < -1), 0)) stop("Please use fractional values to specify the contrast weights (i.e., no contrast weights less than -1.")
if(!identical(round(sum(c.weights), 4), 0)) stop("The sum of the contrast weights ('c.weights') should equal zero.")
    
if(is.null(s.anova)) stop("You must specify the standard deviation of the errors (i.e., the square root of the error variance).")
if(is.null(means)) stop("You must specify the vector of means ('means').")
if(is.null(conf.level)) conf.level <- .95
if(is.null(N) && is.null(n)) stop("You must specify the either total sample size ('N'), or sample sizes per group('n').")
if(is.null(N) && !is.null(n)) N <- sum(n)
if(is.null(df.error)) df.2 <- N - length(means)
##########################################

part.of.se <- sqrt(sum((c.weights^2)/n))

if(!is.null(Psi))
{
if(!is.null(means)) stop("Since the contrast effect (denoted 'Psi') was specified, you should not specify the vector of means ('means').")
if(is.null(s.anova)) stop("You must specify the standard deviation of the errors (i.e., the square root of the error variance).")
if(is.null(n)) stop("You must specify the vector per group/level sample size ('n').")
if(is.null(c.weights)) stop("You must specify the vector of contrast weights ('c.weights').")
}


if(!is.null(means))
{
Psi <- sum(c.weights*means)
}

if(is.null(alpha.lower) & is.null(alpha.upper))
{
alpha.lower <- (1-conf.level)/2
alpha.upper <- (1-conf.level)/2
}

CV.up <- qt(1-alpha.upper, df=df.2, ncp = 0, lower.tail = TRUE, log.p = FALSE)
CV.low <- qt(alpha.lower, df=df.2, ncp = 0, lower.tail = TRUE, log.p = FALSE) 

Result <- list(Lower.Conf.Limit.Contrast = Psi + CV.low*part.of.se*s.anova,
Contrast = Psi, Upper.Conf.Limit.Contrast = Psi + CV.up*part.of.se*s.anova)
        
# print(paste("The", 1 - (alpha.lower + alpha.upper), "confidence limits for the contrast are given as:"))
return(Result)
}
