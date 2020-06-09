####################
# Helper functions #
####################

## Creates an n_i x n_i correlation matrix with AR(1) structure
#' @keywords internal
AR1.mat = function(ti, rho) {
  AR1.cov = matrix(0, nrow=length(ti), ncol=length(ti))
  # If time points are not equispaced
  if(length(unique(diff(ti)))!=1){
    AR1.cov = rho^abs(outer(ti,ti,"-"))
    return(AR1.cov)
  } else {
    # if time points are equispaced
    AR1.cov = rho^abs(outer(seq(1:length(ti)),seq(1:length(ti)),"-"))
    return(AR1.cov)
  }
}

## Creates an n_i x n_i correlation matrix with compound symmetry (CS) structure
#' @keywords internal
CS.mat = function(ni, rho) {
  CS.cov = diag(ni)
  CS.cov[lower.tri(CS.cov)] = CS.cov[upper.tri(CS.cov)] = rho
  return(CS.cov)
}

## Fast inversion of AR1 matrix
#' @keywords internal
inv.AR1.mat = function(A,ti){
  # If time points are not equispaced
  if(length(unique(diff(ti)))!=1){
    return(solve(A))
  } else {
   # if time points are equispaced, use this fast inversion trick 
   if(dim(A)[1]==1){
      return(1/A)
    } else {
      rho = A[1,2]
      n.A = dim(A)[1]
      A.inv = matrix(0, nrow=n.A, ncol=n.A)
      A.inv[1,1]=A.inv[n.A,n.A] = 1/(1-rho^2)
      
      for(j in 2:(n.A-1))
        A.inv[j,j] = (1+rho^2)/(1-rho^2)
      for(j in 2:n.A)
        A.inv[j,j-1] = A.inv[j-1,j] = -rho/(1-rho^2)
      return(A.inv)
    }
  }
}

## Fast inversion of a CS matrix
#' @keywords internal
inv.CS.mat = function(A) {
  if(dim(A)[1]==1){
    return(1/A)
  } else {
    rho = A[1,2]
    n.A = dim(A)[1]
    A.inv = (1/(1-rho))*diag(n.A) - (rho/((1-rho)*(1+rho*(n.A-1))))*matrix(1,nrow=n.A,ncol=n.A)
    return(A.inv)
  }
}

## Prior density Psi. No need for normalizing constant C_k as it cancels out
#' @keywords internal
Psi = function(gamma, lambda) {
  m = length(gamma)
  dens = lambda^m * exp(-lambda*sqrt(sum(gamma^2)))
  
  return(dens)
}

## pStar function
#' @keywords internal
pStar = function(gamma, lambda1, lambda0, theta) {
  Psi1 = Psi(gamma=gamma, lambda=lambda1)
  Psi0 = Psi(gamma=gamma, lambda=lambda0)
  
  ## if a coefficient is really large then both these will 
  ## numerically be zero because R can't handle such small numbers
  if ((theta*Psi1) == 0 & (1 - theta)*Psi0 == 0) {
    p = 1
  } else {
    p = (theta*Psi1) / (theta*Psi1 + (1 - theta)*Psi0)
  }
  return(p)
}

## Lambda star function
#' @keywords internal
lambdaStar = function(gamma, lambda1, lambda0, theta) {
  p = pStar(gamma = gamma, lambda1 = lambda1,
            lambda0 = lambda0, theta = theta)
  
  l = lambda1*p + lambda0*(1 - p)
  return(l)
}

## Compute square root of a matrix
#' @keywords internal
sqrt.mat = function(X){
  ## eigendecomposition of X
  E = eigen(X)
  Q = E$vectors
  ## compute square root of X
  sqrt.X = Q %*% diag(sqrt(E$values)) %*% t(Q)
  
  return(sqrt.X)
}

## Compute log-density of SSGL penalty up to normalizing constant
#' @keywords internal
logSSGL = function(gamma, lambda1, lambda0, theta, d){
  
  norm.term = sqrt(sum(gamma)^2)
  term1 = (1-theta)*lambda0^d*exp(-lambda0*norm.term)
  term2 = theta*lambda1^d*exp(-lambda1*norm.term)
  
  return(log(term1+term2))
}

## Multiplies two conformal matrices
#' @keywords internal
mat.mult = function(A, B){
  return(A%*%B)
}

## EM algorithm for AR(1) or CS strctures with fixed rho
#' @keywords internal
nvcssl_EM_auto = function(rho, cov.structure, y.i, U.i, t.i, n.i,
                          n.tot, n.sub, d, p, a, b, c0, d0, groups,
                          gamma.init, theta.init, sigma2.init, 
                          lambda0, lambda1, tol){
  
  ## Initialize the following values
  difference = 100*tol
  counter = 0
  pstar.k = rep(0,p)                        # to hold pstar_k's
  lambdastar.k = rep(0,p)                   # to hold lambdastar_k's
  R.i = vector("list", n.sub)               # to hold within-subject correlation matrices R_1, ..., R_n
  R.i.inv = vector("list", n.sub)           # to hold R_1^{-1}, ..., R_n^{-1}
  R.i.inv.sqrt = vector("list", n.sub)      # to hold R_1^{-1/2}, ..., R_n^{-1/2}
  y.tilde = rep(0, n.tot)                   # to hold y.tilde = R^{-1/2}%*%Y
  y.tilde.i = vector("list",n.sub)          # to hold subvectors Y.tilde_1, ..., Y.tilde_n
  U.tilde = matrix(0, nrow=n.tot, ncol=d*p) # to hold U.tilde = R^{-1/2}%*%Y
  U.tilde.i = vector("list",n.sub)          # to hold submatrices U.tilde_1, ..., U.tilde_n 
  
    
  ## Initialize R_i's
  if(cov.structure=="AR1"){
    R.i = mapply(AR1.mat, ti=t.i, rho=rho, SIMPLIFY=FALSE)
  }
  if(cov.structure=="CS"){
    R.i = mapply(CS.mat, ni=n.i, rho=rho, SIMPLIFY=FALSE)
  }
  
  ## Initialize R_i^{-1}'s
  if(cov.structure=="AR1"){
    R.i.inv = mapply(inv.AR1.mat, A=R.i, ti=t.i, SIMPLIFY=FALSE)
  }
  if(cov.structure=="CS"){
    R.i.inv = lapply(R.i, inv.CS.mat)
  }
  ## Initialize R_i^{-1/2}'s
  R.i.inv.sqrt = lapply(R.i.inv, sqrt.mat)

  ## Initialize y.tilde_i's
  y.tilde.i = mapply(mat.mult, R.i.inv.sqrt, y.i, SIMPLIFY=FALSE)
  ## Combine
  y.tilde = as.matrix(unlist(y.tilde.i))
  
  ## Center y.tilde to get rid of intercept
  y.tilde = scale(y.tilde, center=TRUE, scale=FALSE)
  
  ## Initialize U.tilde_i's
  U.tilde.i = mapply(mat.mult, R.i.inv.sqrt, U.i, SIMPLIFY=FALSE)
  ## Combine
  U.tilde = plyr::rbind.fill.matrix(U.tilde.i)
  
  ## For EM algorithm
  gamma = gamma.init  
  gamma.hat.list = vector("list",p) # to hold gamma_1, ..., gamma_p
  theta = theta.init
  sigma2 = sigma2.init
  
  ### Update the parameters
  while( (difference > tol) & (counter < 100) ){
    
    counter = counter+1
    
    ## Keeps track of old gamma
    gamma.old = gamma
    
    ##############
    ##############
    ### E-step ###
    ##############
    ##############
    for(k in 1:p){
      ## Which groups are active
      active = which(groups == k)
      ## Update pStar
      pstar.k[k] = pStar(gamma.old[active], lambda1, lambda0, theta)
      # Update lambda.k.star
      lambdastar.k[k] = lambda1*pstar.k[k] + lambda0*(1-pstar.k[k])
    }
    
    ############## 
    ##############
    ### M-step ###
    ##############
    ##############
    
    ## Update theta
    theta = (a-1 + sum(pstar.k))/(a+b+p-2)
    
    ## Update gamma
    ## Note that grpreg solves is (1/2n)*||Y-U%*%gamma||_2^2 + pen(gamma)
    ## so we have to divide the multiplier by 1/n
    gamma.mod = grpreg::grpreg(U.tilde, y.tilde, group=groups, lambda=1, 
                               group.multiplier=(sigma2*lambdastar.k)/n.tot)
    gamma = gamma.mod$beta[-1] 
    
    ## Update sigma2
    sigma2 = (d0 + sum((y.tilde-U.tilde%*%as.matrix(gamma))^2))/(n.tot+c0+2)

    ## Update diff
    difference = sum((gamma-gamma.old)^2)
  }
  
  ## Store theta.hat, gamma.hat
  theta.hat = theta
  gamma.hat = gamma
  gamma.hat.list = split(gamma.hat,as.numeric(gl(length(gamma.hat),d,length(gamma.hat))))
  sigma2.hat = sigma2
  
  ## Calculate the log-posterior
  term1 = -(n.tot/2)*log(sigma2.hat) - 0.5*as.double(Reduce("+",lapply(lapply(R.i, Matrix::det), log))) - sum((y.tilde-mat.mult(U.tilde,as.matrix(gamma.hat)))^2)/(2*sigma2.hat) 
  term2 = sum(mapply(logSSGL, gamma.hat.list, lambda1, lambda0, theta, d))
  term3 = (a-1)*log(theta.hat) + (b-1)*log(1-theta.hat)
  
  log.lik = term1+term2+term3
  
  NVC.EM.output <- list(gamma.hat = gamma.hat,
                        theta.hat = theta.hat,
                        sigma2.hat = sigma2.hat,
                        log.lik = log.lik,
                        U.tilde = U.tilde,
                        y.tilde = y.tilde) 
  
  # Return list
  return(NVC.EM.output)
}

## Accepts vector rho as an input
#' @keywords internal
nvcssl_EM_auto.vec = Vectorize(nvcssl_EM_auto, c("rho"))

## EM algorithm for iid errors
## Note the fractional power<1 is for the robustified NVC-SSL model
#' @keywords internal
nvcssl_EM_ind = function(y, U, n.tot, d, p, a, b, c0, d0, groups,
              gamma.init, theta.init, sigma2.init, lambda0, lambda1, tol,
              frac.power=1, update.sigma2=TRUE){
  
  ## Initialize the following values
  difference = 100*tol
  counter = 0
  pstar.k = rep(0,p)              # to hold pstar_k's
  lambdastar.k = rep(0,p)         # to hold lambdastar_k's
  
  ## For EM algorithm
  gamma = gamma.init  
  gamma.hat.list = vector("list",p) # to hold gamma_1, ..., gamma_p
  theta = theta.init
  sigma2 = sigma2.init
  
  ### Update the parameters
  while( (difference > tol) & (counter < 100) ){
    
    counter = counter+1
    
    ## Keeps track of old gamma
    gamma.old = gamma
    
    ##############
    ##############
    ### E-step ###
    ##############
    ##############
    for(k in 1:p){
      ## Which groups are active
      active = which(groups == k)
      ## Update pStar
      pstar.k[k] = pStar(gamma.old[active], lambda1, lambda0, theta)
      # Update lambda.k.star
      lambdastar.k[k] = lambda1*pstar.k[k] + lambda0*(1-pstar.k[k])
    }
    
    ############## 
    ##############
    ### M-step ###
    ##############
    ##############
    
    ## Update theta
    theta = (a-1 + sum(pstar.k))/(a+b+p-2)
    
    ## Update gamma
    ## Note that grpreg solves is (1/2n)*||Y-U%*%gamma||_2^2 + pen(gamma),
    ## so we have to divide multiplier by n.tot
    gamma.mod = grpreg::grpreg(U, y, group=groups, lambda=1, 
                       group.multiplier=(sigma2*lambdastar.k)/(n.tot*frac.power))
    gamma = gamma.mod$beta[-1]
    
    ## Update sigma2 if update.sigma2 is set to TRUE
    if(update.sigma2==TRUE){
      sigma2 = (d0 + sum((y-U%*%as.matrix(gamma))^2))/(n.tot+c0+2)
    } else if(update.sigma2==FALSE){
        sigma2 = sigma2
    }
    ## Update diff
    difference = sum((gamma-gamma.old)^2)
  }
  
  ## Store theta.hat, gamma.hat
  theta.hat = theta
  gamma.hat = gamma
  gamma.hat.list = split(gamma.hat,as.numeric(gl(length(gamma.hat),d,length(gamma.hat))))
  sigma2.hat = sigma2
  
  NVC.EM.output <- list(gamma.hat = gamma.hat,
                        theta.hat = theta.hat,
                        sigma2.hat = sigma2.hat) 
  
  # Return list
  return(NVC.EM.output)
}

## EM algorithm for completely unstructured NVC-SSL model
#' @keywords internal
nvcssl_EM_unstructured = function(y.i, U.i, n.i, n.tot, n.sub, d, p, a, b, groups,
                                  gamma.init, theta.init, Sigma.i.init,
                                  lambda0, lambda1, tol){
  
  ## Initialize the following values
  difference = 100*tol
  counter = 0
  pstar.k = rep(0,p)                        # to hold pstar_k's
  lambdastar.k = rep(0,p)                   # to hold lambdastar_k's
  
  ## Initialize within-subject correlation matrices
  R.i = vector("list", n.sub)               # to hold within-subject correlation matrices R_1, ..., R_n
  R.i.inv = vector("list", n.sub)           # to hold R_1^{-1}, ..., R_n^{-1}
  R.i.inv.sqrt = vector("list", n.sub)      # to hold R_1^{-1/2}, ..., R_n^{-1/2}
  ## Initialize y.bar and U.bar
  y.bar = rep(0, n.tot)                     # to hold y.bar = R^{-1/2}%*%Y
  U.bar = matrix(0, nrow=n.tot, ncol=d*p)   # to hold U.bar = R^{-1/2}%*%Y
  
  ## For EM algorithm
  gamma = gamma.init  
  gamma.hat.list = vector("list",p) # to hold gamma_1, ..., gamma_p
  theta = theta.init
  Sigma.i = Sigma.i.init
  
  ### Update the parameters
  while( (difference > tol) & (counter < 100) ){
    
    counter = counter+1
    
    ## Keeps track of old gamma
    gamma.old = gamma
    
    ##############
    ##############
    ### E-step ###
    ##############
    ##############
    for(k in 1:p){
      ## Which groups are active
      active = which(groups == k)
      ## Update pStar
      pstar.k[k] = pStar(gamma.old[active], lambda1, lambda0, theta)
      # Update lambda.k.star
      lambdastar.k[k] = lambda1*pstar.k[k] + lambda0*(1-pstar.k[k])
    }
    
    ############## 
    ##############
    ### M-step ###
    ##############
    ##############
    
    ## Update theta
    theta = (a-1 + sum(pstar.k))/(a+b+p-2)
    
    ## Update U.bar and Y.bar
    R.i.inv = lapply(Sigma.i, solve)
    R.i.inv.sqrt = lapply(R.i.inv, sqrt.mat)
    U.bar.i = mapply(mat.mult, R.i.inv.sqrt, U.i, SIMPLIFY=FALSE)
    U.bar = plyr::rbind.fill.matrix(U.bar.i)
    y.bar.i = mapply(mat.mult, R.i.inv.sqrt, y.i, SIMPLIFY=FALSE)
    y.bar = as.matrix(unlist(y.bar.i))
    y.bar = scale(y.bar, center=TRUE, scale=FALSE)
    
    ## Update gamma
    ## Note that grpreg solves is (1/2n)*||Y.bar-U.bar%*%gamma||_2^2 + pen(gamma),
    ## so we have to divide multiplier by n.tot
    gamma.mod = grpreg::grpreg(U.bar, y.bar, group=groups, lambda=1, 
                       group.multiplier=lambdastar.k/n.tot)
    gamma = gamma.mod$beta[-1]
    
    ## Update Sigma.i's
    for(r in 1:n.sub){
      Sigma.i[[r]] = (diag(n.i[r])+(y.i[[r]] - U.i[[r]]%*%gamma)%*%t(y.i[[r]]-U.i[[r]]%*%gamma))/(2*n.i[r]+4)
      # to force positive-definiteness
      Sigma.i[[r]] = Matrix::nearPD(Sigma.i[[r]])$mat 
    }
    
    ## Update diff
    difference = sum((gamma-gamma.old)^2)
  }
  
  ## Store theta.hat, gamma.hat
  theta.hat = theta
  gamma.hat = gamma
  gamma.hat.list = split(gamma.hat,as.numeric(gl(length(gamma.hat),d,length(gamma.hat))))
  Sigma.i.hat = Sigma.i
    
  NVC.EM.output <- list(gamma.hat = gamma.hat,
                        theta.hat = theta.hat,
                        Sigma.i.hat = Sigma.i.hat,
                        U.bar = U.bar,
                        y.bar = y.bar) 
  
  # Return list
  return(NVC.EM.output)
  
}

## Obtain log-likelihood for robustified NVC for a given rho
#' @keywords internal
loglik.robustified.NVC = function(rho, working.cov, y.i, U.i, t.i, n.i, 
                                          n.tot, n.sub, d, p, gamma.init){

  ## Initialize within-subject correlation matrices
  R.i = vector("list", n.sub)               # to hold within-subject correlation matrices R_1, ..., R_n
  R.i.inv = vector("list", n.sub)           # to hold R_1^{-1}, ..., R_n^{-1}
  R.i.inv.sqrt = vector("list", n.sub)      # to hold R_1^{-1/2}, ..., R_n^{-1/2}
  y.tilde = rep(0, n.tot)                   # to hold y.tilde = R^{-1/2}%*%Y
  U.tilde = matrix(0, nrow=n.tot, ncol=d*p) # to hold U.tilde = R^{-1/2}%*%Y
  
  if(working.cov=="AR1"){
    R.i = mapply(AR1.mat, ti=t.i, rho=rho, SIMPLIFY=FALSE)
  }
  if(working.cov=="CS"){
    R.i = mapply(CS.mat, ni=n.i, rho=rho, SIMPLIFY=FALSE)
  }
  
  ## Initialize R_i^{-1}'s and R_i^{-1/2}'s
  if(working.cov=="AR1"){
    R.i.inv = mapply(inv.AR1.mat, A=R.i, ti=t.i, SIMPLIFY=FALSE)
  }
  if(working.cov=="CS"){
    R.i.inv = lapply(R.i, inv.CS.mat)
  }

  R.i.inv.sqrt = lapply(R.i.inv, sqrt.mat)
  
  y.tilde = rep(0, n.tot)                   # y.tilde = R^{-1/2}%*%Y
  U.tilde = matrix(0, nrow=n.tot, ncol=d*p) # U.tilde = $R^{-1/2}%*%Y
  y.tilde.i = vector("list",n.sub)          # to hold subvectors Y.tilde_1, ..., Y.tilde_n
  U.tilde.i = vector("list",n.sub)          # to hold submatrices U.tilde_1, ..., U.tilde_n 
  
  ## Initialize y.tilde_i's
  y.tilde.i = mapply(mat.mult, R.i.inv.sqrt, y.i, SIMPLIFY=FALSE)
  ## Combine
  y.tilde = as.matrix(unlist(y.tilde.i))
  
  ## Initialize U.tilde.i_i's
  U.tilde.i = mapply(mat.mult, R.i.inv.sqrt, U.i, SIMPLIFY=FALSE)
  ## Combine
  U.tilde = plyr::rbind.fill.matrix(U.tilde.i)
  
  ## sigma2.hat is (||R^{-1/2}(rho) (Y-U%*%gamma.init) ||_2^2)/N
  sigma2.hat = sum((y.tilde-U.tilde%*%as.matrix(gamma.init))^2)/n.tot 
  
  ## Compute the log-likelihood
  loglik = -(n.tot/2)*log(sigma2.hat) - 0.5*as.double(0.5*Reduce("+",lapply(lapply(R.i, Matrix::det), log))) - n.tot/2
  
  loglik.output = list(sigma2.hat = sigma2.hat, 
                       loglik = loglik, 
                       U.tilde = U.tilde,
                       y.tilde = y.tilde)
  
  ## Return list
  return(loglik.output)
}

## Accepts vector rho as an input
#' @keywords internal
loglik.robustified.NVC.vec = Vectorize(loglik.robustified.NVC, c("rho"))