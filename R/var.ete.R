var.ete <- function(sigma2, sigmaz2, n1, n2, beta1, beta2, muz=0, c=0, type="sample", covariate.value="sample.mean"){
  N <- n1+n2
  
  if(covariate.value == "sample.mean"){
  if(type == "population"){
    Variance <- sigma2*(1/n1+1/n2+n2/(N*n1*(n1-3))+n1/(N*n2*(n2-3))) + sigmaz2*(beta1-beta2)^2/N
  }
  if(type == "sample"){
    Variance <- sigma2*(1/n1+1/n2+n2/(N*n1*(n1-3))+n1/(N*n2*(n2-3))) + sigmaz2*(beta1-beta2)^2/N-
      sigma2/N*(1/(N-1)*((N-3)/(n1-3)+(N-3)/(n2-3)))
  }}
  
  if(covariate.value == "SD"){
    if(type == "population"){
      C1 <- 1/N+1-2/(N-1)*exp(2*(lgamma(N/2)-lgamma((N-1)/2)))
      C0 <- 1/n1+1/n2+n2/(N*n1*(n1-3))+n1/(N*n2*(n2-3))+1/(N-1)*((N-3)/(n1-3)+(N-3)/(n2-3))
      Variance <- sigma2*C0+sigmaz2*(beta1-beta2)^2*C1
    }
    if(type == "sample"){
      C1 <- 1/N+1-2/(N-1)*exp(2*(lgamma(N/2)-lgamma((N-1)/2)))
      C0 <- 1/n1+1/n2+n2/(N*n1*(n1-3))+n1/(N*n2*(n2-3))+1/(N-1)*((N-3)/(n1-3)+(N-3)/(n2-3))
      Variance <- sigma2*C0+sigmaz2*(beta1-beta2)^2*C1-C1*sigma2/(N-1)*((N-3)/(n1-3)+(N-3)/(n2-3))
    }
    
  }
  
  if(covariate.value == "fixed"){
    if(type == "population"){
      C1 <- (n1-2)/(n1*(n1-3))+  (n2-2)/(n2*(n2-3))
      C2 <- 1/(n1-3)+1/(n2-3)
      Variance = sigma2*(C1+(c-muz)^2*C2/sigmaz2)
    }
    if(type == "sample"){
      C1 <- (n1-2)/(n1*(n1-3))+  (n2-2)/(n2*(n2-3))
      C2 <- 1/(n1-3)+1/(n2-3)
      Variance <- sigma2*(C1+(c-muz)^2*C2/sigmaz2*(N-3)/(N-1)-C2/N)
    }
  }
return(Variance)
}

