###############################
###############################
## FUNCTION FOR IMPLEMENTING ##
## PENALIZED NVC MODELS WITH ##
## THE GROUP LASSO           ##
###############################
###############################

# INPUTS:
# y = Nx1 vector of response observations y_11, ..., y_1n_1, ..., y_n1, ..., y_nn_n
# t = Nx1 vector of observation times t_11, ..., t_1n_1, ..., t_n1, ..., t_nn_n
# X = Nxp design matrix with columns [X_1, ..., X_p], where the kth column
#     is (x_{1k}(t_11), ..., x_{1k}(t_1n_1), ..., x_{nk}(t_n1), ..., x_{nk}(t_nn_n))'
# df = number of basis functions to use
# penalty = which penalty to use (group lasso, group SCAD, or group MCP). 
#           The default is group MCP.
# lambda = optional sequence of tuning parameters. If not provided, the program chooses
#          the grid automatically. The optimal lambda is chosen via cross-validation.

# OUTPUT:
# t.ordered = all N time points in order from smallest to largest. Needed for plotting 
# beta.hat = Nxp matrix of the function estimates. The jth column is the smooth function 
#            beta_k(t_ij) evaluated at observation t_ij
# gamma.hat = estimated basis coefficients (for prediction)
# intercept = intercept beta_0 (for prediction)
# classifications = px1 vector of indicator variables. "1" indicates covariate is significant, "0" indicates insignificant.
# AICc = AIC with correction for small sample size. For tuning the degrees of freedom

NVC_frequentist = function(y, t, X, df=8, penalty=c("gLasso","gSCAD","gMCP"),
                           lambda = NULL) {
  
  ##################
  ##################
  ### PRE-CHECKS ###
  ##################
  ##################
  
  ## Check that dimensions are conformal
  if( length(t) != length(y) | length(y) != dim(X)[1] | length(t) != dim(X)[1] )
    stop("Non-conformal dimensions for y, t, and X.")
  if( df <= 2 )
    stop("Please enter an integer greater than or equal to 3 for degrees of freedom.")
  
  ###################
  ###################
  ### PRE-PROCESS ###
  ###################
  ###################

  ## Coercion
  penalty <- match.arg(penalty)
  
  ## Make X a matrix if it is a data already done
  X = as.matrix(X)
  ## Make df an integer if not already
  df = as.integer(df)
  ## Make 'id' a factor vector if not already done
  id = factor(id)
  ## Center y
  y.mean = mean(y)
  y = y-y.mean

  ## Extract total # of observations, # of subjects, # of observations per subject, # of covariates, and # of basis functions
  n.tot = length(y)           # total number of observations
  n.sub = length(unique(id))  # number of subjects
  p = dim(X)[2]               # number of covariates
  d = df                      # dimension of the vectors of basis coefficients
  groups = rep(1:p, each=d)   # for groups of basis coefficients
  gamma.hat = rep(0,d*p)      # to hold estimated basis coefficients
  gamma.hat.list = vector("list", p)     
  
  # Create Nxd spline matrix with d degrees of freedom
  B = splines::bs(t, df=d, intercept=TRUE)
  
  # Matrix for storing individual B(t_ij) matrices
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
  
  ##############################
  ##############################
  ## Fit the penalized method ##
  ##############################
  ##############################
  if(penalty=="gLasso"){

    if(is.null(lambda)){
      opt.lambda = grpreg::cv.grpreg(U, y, group=groups, penalty="grLasso", nfolds=10)$lambda.min
      gamma.hat = grpreg::grpreg(U, y, group=groups, penalty="grLasso", lambda=opt.lambda)$beta[-1]
    } else if(!is.null(lambda)){
      opt.lambda = grpreg::cv.grpreg(U,y, groups=groups, lambda=lambda, penalty="grLasso")$lambda.min
      gamma.hat = grpreg::grpreg(U, y, group=groups, penalty="grLasso", lambda=opt.lambda)$beta[-1]
    }
    gamma.hat.list = split(gamma.hat,as.numeric(gl(length(gamma.hat),d,length(gamma.hat))))
  
  } else if(penalty=="gSCAD"){
  
    if(is.null(lambda)){
      opt.lambda = grpreg::cv.grpreg(U, y, group=groups, penalty="grSCAD", nfolds=10)$lambda.min
      gamma.hat = grpreg::grpreg(U, y, group=groups, penalty="grSCAD", lambda=opt.lambda)$beta[-1]
    } else if(!is.null(lambda)){
      opt.lambda = grpreg::cv.grpreg(U,y, group=groups, lambda=lambda, penalty="grSCAD")$lambda.min
      gamma.hat = grpreg::grpreg(U, y, group=groups, penalty="grSCAD", lambda=opt.lambda)$beta[-1]
    }
    gamma.hat.list = split(gamma.hat,as.numeric(gl(length(gamma.hat),d,length(gamma.hat))))
    
  } else {
    ## Default is MCP
    if(is.null(lambda)){
      opt.lambda = grpreg::cv.grpreg(U, y, group=groups, penalty="grMCP", nfolds=10)$lambda.min
      gamma.hat = grpreg::grpreg(U, y, group=groups, penalty="grMCP", lambda=opt.lambda)$beta[-1]
    } else if(!is.null(lambda)){
      opt.lambda = grpreg::cv.grpreg(U,y, group=groups, lambda=lambda, penalty="grMCP")$lambda.min
      gamma.hat = grpreg::grpreg(U, y, group=groups, penalty="grMCP", lambda=opt.lambda)$beta[-1]
    }
    gamma.hat.list = split(gamma.hat,as.numeric(gl(length(gamma.hat),d,length(gamma.hat))))
  }
  
  ## Intercept estimate
  intercept = y.mean - as.double(crossprod(U.colmeans,gamma.hat))
  
  ## Order the times
  t.ordered = t[order(t)]
  
  ## Calculate the estimated beta_k(t)'s and store classifications
  beta.hat = matrix(0, nrow=n.tot, ncol=p)
  beta.ind.hat = rep(0,n.tot)
  classifications = rep(0, p)
  
  for(k in 1:p){
      ## Update beta.hat
      beta.ind.hat = B %*% gamma.hat.list[[k]]
      beta.hat[,k] = beta.ind.hat[order(t)]
      
      ## Update classifications
      if(!identical(beta.hat[,k],rep(0,n.tot)))
        classifications[k] = 1    
  } 
  
  ## Calculate AICc
  non.zero = which(gamma.hat!=0)
  s.hat = sum(classifications)
  rss = sum((y-U[,non.zero]%*%gamma.hat[non.zero])^2)
  AICc = log(rss/n.tot)+1+(2*(s.hat-1))/(n.tot-s.hat-2)
  
  #####################
  #####################
  ### Return a list ###
  #####################
  #####################
  NVC.frequentist.output <- list(t.ordered = t.ordered,
                         beta.hat = beta.hat,
                         gamma.hat = gamma.hat,
                         intercept = intercept,
                         classifications = classifications,
                         AICc = AICc) 
  # Return list
  return(NVC.frequentist.output)
}