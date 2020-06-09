###################################
###################################
## FUNCTION FOR IMPLEMENTING THE ##
## NVC-SSL MODEL.                ##
###################################
###################################

# INPUTS:
# y = Nx1 vector of response observations y_11, ..., y_1n_1, ..., y_n1, ..., y_nn_n
# t = Nx1 vector of observation times t_11, ..., t_1n_1, ..., t_n1, ..., t_nn_n
# id = Nx1 vector of labels (different label corresponds to different subjects). Labels should be (1,...,1, 2,...,2, ..., n,....,n)
# X = Nxp design matrix with columns [X_1, ..., X_p], where the kth column
#     is (x_{1k}(t_11), ..., x_{1k}(t_1n_1), ..., x_{nk}(t_n1), ..., x_{nk}(t_nn_n))'
# df = number of basis functions to use
# cov.structure = covariance structure: can specify first-order autoregressive (AR1),
#                 compound symmetry (CS), iid (ind), or unstructured
# lambda0 = a ladder of increasing values for the hyperparameters. Default is {1,...,25}
# lambda1 = a fixed small value for the hyperparameter. Default is 1
# rho = correlation parameter
# a = hyperparameter in B(a,b) prior on theta. Default is 1.
# b = hyperparameter in B(a,b) prior on theta. Default is p.
# c0 = hyperparameter in IG(c_0/2,d_0/2) prior on sigma2. Default is 1.
# d0 = hyperparameter in IG(c_0/2,d_0/2) prior on sigma2. Default is 1.
# rho.support = support for hyperparameter rho. Default is (0, 0.1, ..., 0.9).
#               Ignored if "ind" or "unstructured" is specified as covariance structure.
# tol = convergence criteria. Default is 1e-6.
# print.iteration = flag for printing the current value of spike hyperparameter lambda0 in 
#                    the lambda0 ladder. Default is TRUE.

# OUTPUT:
# t.ordered = all N time points in order from smallest to largest. Needed for plotting 
# beta.hat = Nxp matrix of the function estimates. The kth column is the function estimate 
#            beta_k(t_ij) evaluated at all N time observations t_ij, i=1,...,n, j=1,...,n_i.
# gamma.hat = dfxp estimated basis coefficients (for prediction)
# intercept = intercept beta_0 (for prediction)
# classifications = px1 vector of indicator variables. "1" indicates covariate is significant, "0" indicates insignificant.
# AICc = AIC with correction for small sample size. For model selection of degrees of freedom d

NVC_SSL = function(y, t, id, X, df=8, cov.structure=c("AR1","CS","ind","unstructured"),
                  lambda0=c(seq(from=5,to=30,by=5),seq(from=40,to=100,by=10)), lambda1=1, 
                  a=1, b=ncol(X), c0=1, d0=1, rho.support=seq(0,0.9,by=0.1), 
                  tol = 1e-6, print.iteration=TRUE) {
  
  
  ##################
  ##################
  ### PRE-CHECKS ###
  ##################
  ##################
  
  ## Check that dimensions are conformal
  if( length(t) != length(y) | length(y) != dim(X)[1] | length(t) != dim(X)[1] )
    stop("Non-conformal dimensions for y, t, and X.")
  ## Check that degrees of freedom is >=3
  if( df <= 2 )
    stop("Please enter a positive integer greater than or equal to three for degrees of freedom.")
  ## Check that the ladder is increasing and that all relevant hyperparameters are positive
  if ( !all.equal(cummax(lambda0), lambda0))
    stop("Please enter a strictly increasing ladder of positive values for the spike hyperparameters.")
  if ( !(all(lambda0 > 0) & (lambda1 > 0) & (a > 0) & (b > 0) & (c0 > 0) & (d0 > 0)) )
    stop("Please make sure that all hyperparameters are strictly positive.")
  if (lambda1 > min(lambda0))
    stop("Please make sure that lambda1 is smaller than lambda0.")
  if (min(rho.support) < 0 || max(rho.support) >= 1)
    stop("Please ensure that the support of rho is values greater than or equal to zero and strictly less than one.")

  ###################
  ###################
  ### PRE-PROCESS ###
  ###################
  ###################

  ## Coercion
  cov.structure <- match.arg(cov.structure)
  
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
  
  ## Create vectors to keep track of iterations to convergence and parameter estimates 
  L = length(lambda0)
  lambda0 = c(lambda0)
  
  # Initialize values for gamma and theta
  gamma.hat = rep(0, d*p)             # initial gamma
  gamma.hat.list = vector("list",p)   # to hold gamma_1, ..., gamma_p
  theta.hat = 0.5                     # initial theta
  
  ## if covariance structure is not completely unstructured, initialize sigma2
  if(cov.structure != "unstructured"){

    ## Initialize sigma2
    df.sigma2 <- 3
    sigquant <- 0.9
    sigest <- stats::sd(y)
    qchi <- stats::qchisq(1 - sigquant, df.sigma2)
    ncp <- sigest^2 * qchi / df.sigma2
    # Initial overestimate based on Chipman et al. (2010)
    init.overestimate <- df.sigma2 * ncp / (df.sigma2 + 2)
    
    if(sigest >= 1){
      sigma2.hat = init.overestimate
    } else {
      sigma2.hat = 1
    }
  }
  
  ## If covariance structure is AR(1) or CS, we need the autocorrelation parameter
  if(cov.structure=="AR1" || cov.structure=="CS"){
    # to hold output for different values of rho
    rho.vec = rho.support
  }
  
  ## If covariance structure is unstructured
  if(cov.structure=="unstructured"){
    # initial Sigma_i's, i=1,...,n
    Sigma.i.hat = vector("list", n.sub)
    
    for(r in 1:n.sub)
      Sigma.i.hat[[r]] = diag(n.i[r])
  }
 
  #################################
  #################################
  ### EM algorithm with dynamic ### 
  ### posterior exploration     ###
  #################################
  #################################

  for(s in 1:L){
    
    ## To output iteration
    if(print.iteration==TRUE){
      cat("lambda0 = ", lambda0[s], "\n")
    }
    
    ## Initialize gamma and theta for the NVC-SSL model
    gamma.init = gamma.hat
    theta.init = theta.hat
    
    ## If covariance structure is not unstructured, initialize sigma2
    if (cov.structure != "unstructured"){
      sigma2.init = sigma2.hat
    }
    ## If covariance structure is unstructured
    if (cov.structure == "unstructured"){
      Sigma.i.init = vector("list", n.sub)
      Sigma.i.init = Sigma.i.hat
    }
    
    ## Run EM algorithm
    if(cov.structure=="CS"){
        EM_output = nvcssl_EM_auto.vec(rho=rho.vec, cov.structure="CS", y.i, U.i, t.i, n.i,
                                  n.tot, n.sub, d, p, a, b, c0, d0, groups,
                                  gamma.init, theta.init, sigma2.init, 
                                  lambda0[s], lambda1, tol)
    } else if(cov.structure=="ind"){
        EM_output = nvcssl_EM_ind(y, U, n.tot, d, p, a, b, c0, d0, groups,
                                  gamma.init, theta.init, sigma2.init,
                                  lambda0[s], lambda1, tol)
    } else if(cov.structure=="unstructured"){
        EM_output = nvcssl_EM_unstructured(y.i, U.i, n.i, n.tot, n.sub, d, p, a, b, 
                                           groups, gamma.init, theta.init, Sigma.i.init, 
                                           lambda0[s], lambda1, tol)
    } else {
      ## AR(1) model is the default
      EM_output = nvcssl_EM_auto.vec(rho=rho.vec, cov.structure="AR1", y.i, U.i, t.i, n.i,
                                     n.tot, n.sub, d, p, a, b, c0, d0, groups,
                                     gamma.init, theta.init, sigma2.init, 
                                     lambda0[s], lambda1, tol)
    }
    
    if(cov.structure == "AR1" || cov.structure == "CS"){
      ## Compute optimal rho
      max.log.lik = EM_output[,1]$log.lik
      opt.rho.index = 1
    
      for(m in 2:length(rho.vec)){
        if(EM_output[,m]$log.lik > max.log.lik){
          max.log.lik = EM_output[,m]$log.lik
          opt.rho.index = m
        }
      }
      
      ## Store rho.hat, gamma.hat, theta.hat, sigma2.hat for next pass
      ## of dynamic posterior exploration
      rho.hat = rho.vec[opt.rho.index]

      gamma.hat = EM_output[,opt.rho.index]$gamma.hat
      theta.hat = EM_output[,opt.rho.index]$theta.hat
      sigma2.hat = EM_output[,opt.rho.index]$sigma2.hat

    } else if(cov.structure=="ind"){
      ## Store optimal gamma, theta, sigma2
      gamma.hat = EM_output$gamma.hat
      theta.hat = EM_output$theta.hat
      sigma2.hat = EM_output$sigma2.hat
  
    } else if(cov.structure=="unstructured"){
      ## Store optimal 
      gamma.hat = EM_output$gamma.hat
      theta.hat = EM_output$theta.hat
      Sigma.i.hat = EM_output$Sigma.i.hat
    } 
  }
  
  ## Intercept estimate
  intercept = y.mean - as.double(crossprod(U.colmeans,gamma.hat))
  
  ## Store U.tilde for AR(1) or CS
  if (cov.structure == "AR1" || cov.structure == "CS"){
    U.tilde = EM_output[,opt.rho.index]$U.tilde
    y.tilde = EM_output[,opt.rho.index]$y.tilde
  }
  if (cov.structure=="unstructured"){
    U.bar = EM_output$U.bar
    y.bar = EM_output$y.bar
  }
  
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
  
  if(cov.structure == "AR1" || cov.structure == "CS"){
    rss = sum((y.tilde-U.tilde[,non.zero]%*%gamma.hat[non.zero])^2)
    AICc = log(rss/n.tot)+1+(2*(s.hat-1))/(n.tot-s.hat-2)
  
  } else if(cov.structure == "ind"){
    rss = sum((y-U[,non.zero]%*%gamma.hat[non.zero])^2)
    AICc = log(rss/n.tot)+1+(2*(s.hat-1))/(n.tot-s.hat-2)
  
  } else if(cov.structure == "unstructured"){
    rss = sum((y.bar-U.bar[,non.zero]%*%gamma.hat[non.zero])^2)
    AICc = log(rss/n.tot)+1+(2*(s.hat-1))/(n.tot-s.hat-2)
  }
    
  #####################
  #####################
  ### Return a list ###
  #####################
  #####################
  NVC.SSL.output <- list(t.ordered = t.ordered,
                         beta.hat = beta.hat,
                         gamma.hat = gamma.hat,
                         intercept = intercept,
                         classifications = classifications,
                         AICc = AICc) 
  # Return list
  return(NVC.SSL.output)
}