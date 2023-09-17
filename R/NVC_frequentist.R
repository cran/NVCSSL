############################################
## Implements frequentist NVC models with ##
## different choices of penalty function  ##
############################################

# INPUTS:
# y = Nx1 vector of response observations y_11, ..., y_1m_1, ..., y_n1, ..., y_nm_n
# t = Nx1 vector of observation times t_11, ..., t_1m_1, ..., t_n1, ..., t_nm_n
# X = Nxp design matrix with columns [X_1, ..., X_p], where the kth column
#     is (x_{1k}(t_11), ..., x_{1k}(t_1m_1), ..., x_{nk}(t_n1), ..., x_{nk}(t_nm_n))
# n_basis = number of basis functions to use. Default is n_basis=8.
# penalty = which penalty to use (group lasso, group SCAD, or group MCP). 
# lambda = tuning parameter. If lambda is not specified, then the
#          program automatically chooses a grid of lambdas. 
# include_intercept = whether or not to include intercept function beta_0(t).
#                     Default is TRUE.

# OUTPUT:
# t_ordered = all N time points ordered from smallest to largest. Needed for plotting
# classifications = px1 vector of indicator variables, where "1" indicates that the covariate is 
#                   selected and "0" indicates that it is not selected.
#                   These classifications are determined by the optimal lambda chosen from BIC.
#                   Note that this vector does not include an intercept function.
# beta_hat = Nxp matrix of the varying coefficient function estimates for the optimal lambda 
#            chosen from BIC. The kth column in the matrix is the kth estimated function 
#            at the observation times in t_ordered
# beta0_hat =  intercept function estimate at the observation times in t_ordered for the
#              optimal lambda chosen from BIC. This is not returned if include_intercept = FALSE.
# gamma_hat = estimated basis coefficients (needed for prediction) for the optimal lambda.
# lambda_min = lambda which minimizes the BIC. If only one value was pass for lambda, then
#              this just returns that lambda.
# lambda_all = grid of all L lambdas.  Note that since the objective function is scaled
#              by 1/N for the penalized frequentist methods, the lambda grid that is chosen
#              automatically by the program consists of smaller values than the lambda grid 
#              for the NVC-SSL model. 
# BIC_all = Lx1 vector of BIC values corresponding to all L lambdas in lambda_all.
#           The lth entry corresponds to the lth entry in lambda_all.
# beta_est_all_lambda = list of length L of the estimated varying coefficients corresponding to
#                       all L lambdas in lambda_all. The lth entry corresponds to the lth 
#                       entry in lambda_all.
# beta0_est_all_lambda = NxL matrix of estimated intercept functions corresponding to all
#                        L lambdas in lambda_all. The lth column corresponds to the lth
#                        entry in lambda_all.
# gamma_est_all_lambda = dpxL matrix of estimated basis coefficients corresponding to all
#                        the lambdas in lambda_all. The lth column corresponds to the lth
#                        entry in lambda_all.
# classifications_all_lambda = pxL matrix of classifications corresponding to all the 
#                              lambdas in lambda_all. The lth column corresponds to the
#                              lth entry in lambda_all.
# iters_to_converge = number of iterations it took for the group ascent algorithm to converge.

NVC_frequentist = function(y, t, X, n_basis=8, penalty=c("gLASSO","gSCAD","gMCP"),
                           lambda=NULL, include_intercept=TRUE) {
  
  ##################
  ##################
  ### PRE-CHECKS ###
  ##################
  ##################
  
  ## Check that dimensions are conformal
  if( length(t) != length(y) | length(y) != dim(X)[1] | length(t) != dim(X)[1] )
    stop("Non-conformal dimensions for y, t, and X.")
  if( n_basis <= 2 )
    stop("Please enter an integer greater than or equal to 3 for degrees of freedom.")
  
  ###################
  ###################
  ### PRE-PROCESS ###
  ###################
  ###################

  ## Coercion
  penalty <- match.arg(penalty)
  if(penalty=="gLASSO"){
    penalty="grLasso"
  } else if(penalty=="gSCAD"){ 
    penalty="grSCAD"
  } else if(penalty=="gMCP"){
    penalty="grMCP"
  }
  
  ## Make df an integer if not already
  df = as.integer(n_basis)
  
  ## Extract total # of observations, # of subjects, # of observations per subject, # of covariates, and # of basis functions
  n_tot = length(y)           # total number of observations
  d = df                      # dimension of the vectors of basis coefficients
  
  ## If include_intercept=TRUE, augment X with a column of ones
  if(include_intercept==TRUE){
    X = cbind(rep(1,n_tot), X)
  }
  ## Make X a matrix
  X = as.matrix(X)
  p = dim(X)[2]               # number of columns of X
  groups = rep(1:p, each=d)   # for groups of basis coefficients

  # Create Nxd spline matrix with d degrees of freedom
  B = splines::bs(t, df=d, intercept=TRUE)
  
  # Matrix for storing individual B(t_ij) matrices
  B_t = matrix(0, nrow=p, ncol=p*d)
  
  ## Create the U matrix. U is an Nxdp matrix
  U = matrix(0, nrow=n_tot, ncol=d*p) 
  for(r in 1:n_tot){
    for(j in 1:p){
      B_t[j,(d*(j-1)+1):(d*j)] = B[r,]
    }
    U[r,] = X[r,] %*% B_t
  }
  ## Free memory
  rm(B_t)
  
  ## Don't want an intercept, so we center U and y. 
  ## Note that centering does not change the estimates of the basis coefficients.
  U = scale(U, center=TRUE, scale=FALSE)
  y = scale(y, center=TRUE, scale=FALSE)
  
  ##############################
  ##############################
  ## Fit the penalized method ##
  ##############################
  ##############################
  if(is.null(lambda)){
    
    ## Fit model
    grp_mod = grpreg::grpreg(U, y, group=groups, penalty=penalty)
    mod_select = grpreg::select(grp_mod, crit="BIC", df="active")
  
    ## Store results  
    lambda = grp_mod$lambda
    BIC = mod_select$IC
    gamma_mat = grp_mod$beta[-1, ]
    iter_to_converge = grp_mod$iter
  
  } else if(!is.null(lambda)){
    
    ## Fit model
    grp_mod = grpreg(U, y, group=groups, penalty=penalty, lambda=lambda)
        
    if(length(lambda)>1){
      ## Store results
      gamma_mat = grp_mod$beta[-1, ]
      iter_to_converge = grp_mod$iter
      mod_select = grpreg::select(grp_mod, crit="BIC", df="active")
      ## Store results  
      BIC = mod_select$IC
    } else {
      gamma_mat = grp_mod$beta[-1]
      iter_to_converge = grp_mod$iter
      num_nonzero = length(which(gamma_mat!=0))
      BIC = grp_mod$loss + log(n_tot)*num_nonzero
    }
  }
  
  ## lambda which minimizes the BIC
  min_index = as.double(which.min(BIC))
  lambda_min = lambda[min_index]
  
  ## Order the times
  t_ordered = t[order(t)]
  
  ## For storing the results
  beta_list = vector(mode='list', length=length(lambda))
  classifications_mat = matrix(0, nrow=p, ncol=length(lambda))
    
  beta_hat = matrix(0, nrow=n_tot, ncol=p)
  beta_ind_hat = rep(0, n_tot)
  classifications = rep(0, p)
    
  for(l in 1:length(lambda)){
      
    ## Calculate the estimated beta_k(t)'s and store classifications
    if(length(lambda)>1){
      gamma_hat_list = split(gamma_mat[,l], as.numeric(gl(length(gamma_mat[,l]),d,length(gamma_mat[,l]))))
    } else {
      gamma_hat_list = split(gamma_mat, as.numeric(gl(length(gamma_mat),d,length(gamma_mat))))
    }
    
    for(k in 1:p){
      ## Update beta_hat
      beta_ind_hat = B %*% gamma_hat_list[[k]]
      beta_hat[,k] = beta_ind_hat[order(t)]
    
      ## Update classifications
      if(!identical(beta_hat[,k],rep(0,n_tot)))
        classifications[k] = 1    
      
      ## Store results
      beta_list[[l]] = beta_hat
      if(length(lambda)==1){
        classifications_mat = classifications
      } else {
      classifications_mat[,l] = classifications
      }
    }
  }    

  
  if(include_intercept==FALSE){
    
    beta_hat = beta_list[[min_index]]
    
    if(length(lambda)>1){
      gamma_hat = gamma_mat[,min_index]
      classifications = classifications_mat[,min_index] 
    } else {
      gamma_hat = gamma_mat[min_index]
      classifications = classifications_mat
    }
    
    NVC_frequentist_output <- list(t_ordered = t_ordered,
                                   classifications = classifications,
                                   beta_hat = beta_hat,
                                   gamma_hat = gamma_hat,
                                   lambda_min = lambda_min,
                                   lambda_all = lambda,
                                   BIC_all = BIC,
                                   beta_est_all_lambda = beta_list,
                                   gamma_est_all_lambda = gamma_mat,
                                   classifications_all_lambda = classifications_mat,
                                   iter_to_converge = iter_to_converge)
  
  } else if(include_intercept==TRUE) {
    
      ## Separate the intercept function from the other functions
      beta0_mat = matrix(0, nrow=n_tot, ncol=length(lambda))
      
      for(l in 1:length(lambda)){
        temp = beta_list[[l]]
        beta0_mat[,l] = temp[,1]
        beta_hat = temp[,-1]
        beta_list[[l]] = beta_hat
      }
      beta_hat = beta_list[[min_index]]
 
      
      if(length(lambda)>1){
        gamma_hat = gamma_mat[,min_index]
        beta0_hat = beta0_mat[,min_index]
        classifications_mat = classifications_mat[-1,]
        classifications = classifications_mat[,min_index]
      } else {
        gamma_hat = gamma_mat
        beta0_hat = as.double(beta0_mat)
        classifications_mat = classifications_mat[-1]
        classifications = classifications_mat
      }
    
      NVC_frequentist_output <- list(t_ordered = t_ordered,
                                     classifications = classifications,
                                     beta_hat = beta_hat,
                                     beta0_hat = beta0_hat,
                                     gamma_hat = gamma_hat,
                                     lambda_min = lambda_min,
                                     lambda_all = lambda,
                                     BIC_all = BIC,
                                     beta_est_all_lambda = beta_list,
                                     beta0_est_all_lambda = beta0_mat,
                                     gamma_est_all_lambda = gamma_mat,
                                     classifications_all_lambda = classifications_mat,
                                     iter_to_converge = iter_to_converge)
  }
  
  #####################
  #####################
  ### Return a list ###
  #####################
  #####################
  # Return list
  return(NVC_frequentist_output)
  
}