###################################
###################################
## FUNCTION FOR IMPLEMENTING THE ##
## ROBUSTIFIED NVC-SSL MODEL.    ##
###################################
###################################

# INPUTS:
# y = Nx1 vector of response observations y_11, ..., y_1n_1, ..., y_n1, ..., y_nn_n
# t = Nx1 vector of observation times t_11, ..., t_1n_1, ..., t_n1, ..., t_nn_n
# id = Nx1 vector of labels (different label corresponds to different subjects). Labels should be (1,...,1, 2,...,2, ..., n,....,n)
# X = Nxp design matrix with columns [X_1, ..., X_p], where the kth column
#     is (x_{1k}(t_11), ..., x_{1k}(t_1n_1), ..., x_{nk}(t_n1), ..., x_{nk}(t_nn_n))'
# df = number of basis functions to use
# working.cov = working covariance structure. Can specify AR(1), compound symmetry (CS),
#               or independent
# frac.power = fractional power between (0,1). Default is 0.5.
# lambda0 = a ladder of increasing values for the hyperparameters. Default is {1,...,25}
# lambda1 = a fixed small value for the hyperparameter. Defualt is 1.
# a = hyperparameter in B(a,b) prior on theta. Default is 1.
# b = hyperparameter in B(a,b) prior on theta. Default is p.
# tol = convergence criteria. Default is 1e-6.

# OUTPUT:
# t.ordered = all N time points in order from smallest to largest. Needed for plotting 
# beta.hat = Nxp matrix of the function estimates. The jth column is the smooth function 
#            beta_k(t_ij) evaluated at observation t_ij
# gamma.hat = estimate of basis coefficients. Needed for prediction
# intercept = intercept beta_0 (for prediction)
# classifications = px1 vector of indicator variables. "1" indicates covariate is significant, "0" indicates insignificant.
# AICc = AIC with correction for small sample size. For model selection of fractional power
#        and/or the degrees of freedom. 

robustified_NVC_SSL = function(y, t, id, X, df=8, working.cov=c("AR1","CS","ind"), frac.power=0.5,
                               lambda0=c(seq(from=5,to=30,by=5),seq(from=40,to=100,by=10)), lambda1=1, 
                               a=1, b=ncol(X), tol = 1e-6, print.iteration = TRUE) {
  
  ##################
  ##################
  ### PRE-CHECKS ###
  ##################
  ##################
  
  ## Check that dimensions are conformal
  if( length(t) != length(y) | length(y) != nrow(X) | length(t) != nrow(X) )
    stop("Non-conformal dimensions for y, t, and X.")
  ## Check that degrees of freedom is >=3
  if( df <= 2 )
    stop("Please enter a positive integer for degrees of freedom greater than or equal to 3.")
  ## Check that fractional power is between 0 and 1
  if( frac.power <= 0 || frac.power >= 1)
    stop("Please enter a fractional power greater than 0 and less than 1.")
  ## Check that the ladder is increasing and that all relevant hyperparameters are positive
  if ( !all.equal(cummax(lambda0), lambda0))
    stop("Please enter a strictly increasing ladder of positive values for the spike hyperparameters.")
  if ( !(all(lambda0 > 0) & (lambda1 > 0) & (a > 0) & (b > 0)) )
    stop("Please make sure that all hyperparameters are strictly positive.")
  if (lambda1 > min(lambda0))
    stop("Please make sure that lambda1 is smaller than lambda0.")
  
  ###################
  ###################
  ### PRE-PROCESS ###
  ###################
  ###################

  ## Coercion
  working.cov <- match.arg(working.cov)
  
  ## Make df an integer if not already
  df = as.integer(df)
  ## Make 'id' a factor vector if not already done
  id = factor(id)
  ## Make X a matrix if not already done
  X = as.matrix(X)
  ## Center y
  y.mean = mean(y)
  y = y-y.mean
  
  ## Extract total # of observations, # of subjects, # of observations per subject, # of covariates, and # of basis functions
  n.tot = length(y)           # total number of observations
  n.sub = length(unique(id))  # number of subjects
  n.i = plyr::count(id)$freq        # nx1 vector of per-subject observations (n_1, ..., n_n)'
  p = dim(X)[2]               # number of covariates
  d = df                      # dimension of the vectors of basis coefficients
  groups = rep(1:p, each=d)   # for groups of basis coefficients
  
  ## Extract the subvectors y_i and t.i
  y.i = vector("list", n.sub)
  t.i = vector("list", n.sub)
  for(r in 1:n.sub){
    y.i[[r]] = as.matrix(y[which(id==r)])
    t.i[[r]] = t[which(id==r)]
  }
  
  ## Create Nxd spline matrix with d degrees of freedom
  B = splines::bs(t, df=d, intercept=TRUE)
  
  ## Matrix for storing individual B(t_ij) matrices
  B.t = matrix(0, nrow=p, ncol=p*d)  
  
  ## Create the U matrix. U is an Nxdp matrix
  U = matrix(0, nrow=n.tot, ncol=d*p) 
  for(r in 1:n.tot){
    for(j in 1:p){
      B.t[j,(d*(j-1)+1):(d*j)] = B[r,]
    }
    U[r,] = X[r,] %*% B.t
  }
  ## Free memory
  rm(B.t)
  ## Column means for U
  U.colmeans = colMeans(U)
  
  ## Create the U^i matrices
  U.i = vector("list", n.sub)
  for(r in 1:n.sub){
    U.i[[r]] = U[which(id==r),]
  }
  
  ## Obtain initial estimate of gamma
  opt.lambda = grpreg::cv.grpreg(U, y, group=groups, penalty="grMCP", nfolds=10)$lambda.min
  gamma.init = grpreg::grpreg(U, y, group=groups, penalty="grMCP", lambda=opt.lambda)$beta[-1]
  
  ## Obtain initial estimates of (rho, sigma2) if working AR(1) or CS
  rho.vec = seq(from=0, to=0.9, by=0.1)
  if(working.cov == "AR1" || working.cov == "CS"){
    ## Evaluate log-likelihood for grid of rho's
    loglik_output = loglik.robustified.NVC.vec(rho.vec, working.cov, y.i, U.i, t.i, n.i, 
                                             n.tot, n.sub, d, p, gamma.init)
  
    ## Find the rho that maixmizes the log-likelihood
    max.loglik = loglik_output[,1]$loglik
    opt.rho.index = 1
    
    for(m in 2:length(rho.vec)){
      if(loglik_output[,m]$loglik > max.loglik){
        max.loglik = loglik_output[,m]$loglik
        opt.rho.index = m
      }
    }
    ## Initialize sigma2.hat, U.check, Y.check
    ## degrees of freedom adjustment for # of nonzero elements in gamma.init
    sigma2.hat = (loglik_output[,m]$sigma2.hat)*(n.tot/(n.tot-length(which(gamma.init!=0))))
    U.check = (1/sqrt(sigma2.hat))*loglik_output[,m]$U.tilde
    y.check = (1/sqrt(sigma2.hat))*loglik_output[,m]$y.tilde
  
  } else if(working.cov=="ind"){
    U.check = U
    y.check = y
  }
  ## Center y.check
  y.check = y.check-mean(y.check)
  
  
  ## Create vectors to keep track of iterations to convergence and parameter estimates 
  L = length(lambda0)
  lambda0 = c(lambda0[1],lambda0)
  
  # Initialize values for gamma and theta
  gamma.values = matrix(0, nrow=d*p, ncol=L+1)  # to hold the estimated modes of gamma (original scale)
  gamma.hat.list = vector("list",p)             # to hold gamma_1, ..., gamma_p
  theta.values = rep(0,L+1)                     # to hold the estimated modes of theta
  theta.values[1] = 0.5
  
  
  #################################
  #################################
  ### EM algorithm with dynamic ### 
  ### posterior exploration     ###
  #################################
  #################################
  
  for(s in 2:(L+1)){
    
    ## To output iteration
    if(print.iteration==TRUE){
      cat("lambda0 = ", lambda0[s], "\n")
    }
    
    ## Initialize gamma and theta for the NVC-SSL model
    gamma.init = gamma.values[,s-1]
    theta.init = theta.values[s-1]
    
    ## Run EM algorithm
    EM_output = nvcssl_EM_ind(y=y.check, U=U.check, n.tot, d=d, p=p, a=a, b=b, 
                              groups=groups, gamma.init=gamma.init, 
                              theta.init=theta.init, sigma2.init=1,
                              lambda0=lambda0[s], lambda1=lambda1, 
                              tol=tol, frac.power=frac.power, 
                              update.sigma2=FALSE)
  
    ## Store optimal gamma, theta, sigma2
    gamma.hat = EM_output$gamma.hat
    theta.hat = EM_output$theta.hat
    sigma2.hat = EM_output$sigma2.hat
    
    ## Update tables for the next pass of the dynamic posterior exploration algorithm
    gamma.values[,s] = gamma.hat
    theta.values[s] = theta.hat
  }
  
  ## Store final values
  gamma.hat = gamma.values[,L+1]
  theta.hat = theta.values[L+1]

  ## Intercept estimate
  intercept = y.mean - as.double(crossprod(U.colmeans,gamma.hat))
  
  ## Store gamma_1, ..., gamma_p in a list
  gamma.hat.list = split(gamma.hat,as.numeric(gl(length(gamma.hat),d,length(gamma.hat))))
  
  ## Calculate the estimated beta_k(t)'s and store classifications
  beta.hat = matrix(0, nrow=n.tot, ncol=p)
  beta.ind.hat = rep(0,n.tot)
  classifications = rep(0, p)
  
  for(k in 1:p){
    ## Update beta.hat
    beta.ind.hat = mat.mult(B, as.matrix(gamma.hat.list[[k]]))
    beta.hat[,k] = beta.ind.hat[order(t)]
    
    ## Update classifications
    if(!identical(beta.hat[,k],rep(0,n.tot)))
      classifications[k] = 1    
  } 
  
  ## Order the times
  t.ordered = t[order(t)]
  
  ## Calculate the AICc
  non.zero = which(gamma.hat!=0)
  s.hat = sum(classifications)
  
  rss = sum((y.check-U.check[,non.zero]%*%gamma.hat[non.zero])^2)
  AICc = log(rss/n.tot)+1+(2*(s.hat-1))/(n.tot-s.hat-2)
  
  
  #####################
  #####################
  ### Return a list ###
  #####################
  #####################
  robustified.NVC.SSL.output <- list(t.ordered = t.ordered,
                                     beta.hat = beta.hat,
                                     gamma.hat = gamma.hat,
                                     intercept = intercept,
                                     classifications = classifications,
                                     AICc = AICc) 
  # Return list
  return(robustified.NVC.SSL.output)
}