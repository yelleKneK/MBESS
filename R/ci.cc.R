ci.cc <- function(r, n, conf.level=.95, alpha.lower=NULL, alpha.upper=NULL)
{
if (!is.null(conf.level))
{
if(conf.level >= 1 | conf.level <= 0) stop("You have not properly specified 'conf.level'", call. = FALSE)
if(!is.null(alpha.lower)) stop("You specified both 'conf.level' and 'alpha.lower', specify confidence level using only one approach.", call. = FALSE)
if(!is.null(alpha.upper)) stop("You specified both 'conf.level' and 'alpha.upper', specify confidence level using only one approach.", call. = FALSE)
alpha.lower <- alpha.upper <- (1 - conf.level)/2
}

if(is.null(conf.level))
{
if (is.null(alpha.lower) & is.null(alpha.upper)) stop("You need to specify either 'conf.level', or 'alpha.lower' and 'alpha.upper'.", call. = FALSE)
if (alpha.lower > 0.5 | alpha.lower < 0) stop("You have not properly specified 'alpha.lower' correctly.", call. = FALSE)
if(alpha.upper > 0.5 | alpha.upper < 0) stop("You have not properly specified 'alpha.upper' correctly.", call. = FALSE)
}

CV.Lower <- qnorm(1-alpha.lower)
CV.Upper <- qnorm(1-alpha.upper)

Z <- transform_r.Z(r)
SE.Z <- sqrt(1/(n-3))

CI.Lower_Zeta <- Z - CV.Lower*SE.Z
CI.Upper_Zeta <- Z + CV.Upper*SE.Z

if(alpha.lower > 0) CI.Lower_rho <- transform_Z.r(CI.Lower_Zeta)
if(alpha.upper > 0) CI.Upper_rho <- transform_Z.r(CI.Upper_Zeta)

if(alpha.lower == 0) CI.Lower_rho <- -1
if(alpha.upper == 0) CI.Upper_rho <- 1

return(list(Lower.Limit=CI.Lower_rho, Estimated.Correlation=r, Upper.Limit=CI.Upper_rho))
}