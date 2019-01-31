#Optimization for the Global Covariance Matrix
opt.sigma <- function(nu, lambda, mu, sigprior, gamma=NULL, covar=NULL) {  

  # find the covariance
  if (is.null(mu)) {  # not provided - calculate on the fly
    for (i in 1:nrow(lambda)) {
      for (j in 1:ncol(lambda)) {
        lambda[i, j] <- lambda[i, j] - as.numeric(covar[i,] %*% gamma[,j])
      }
    }
    covariance <- crossprod(lambda)
    # TODO: crossprod probably follows a more complicated logic for dimname transfer
    names <- dimnames(gamma)[[2]]
    if (!is.null(names)) dimnames(covariance) <- list(names, names)
    
  } else {
    if (ncol(mu)==1) {
      covariance <- crossprod(sweep(lambda, 2, STATS=as.numeric(mu), FUN="-"))
    } else {
      covariance <- crossprod(lambda-t(mu)) 
    }    
  }
  
  sigma <- (covariance + nu)/nrow(lambda) #add to estimation variance
  sigma <- diag(diag(sigma),nrow=nrow(nu))*sigprior + (1-sigprior)*sigma #weight by the prior
  return(sigma)
}


