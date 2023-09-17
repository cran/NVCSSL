#####################################
## GIBBS SAMPLING FUNCTION FOR     ##
## IMPLEMENTING THE NVC-SSL MODEL. ##                
#####################################

# INPUTS:
# y = Nx1 vector of response observations y_11, ..., y_1m_1, ..., y_n1, ..., y_nm_n
# t = Nx1 vector of observation times t_11, ..., t_1m_1, ..., t_n1, ..., t_nm_n
# id = Nx1 vector of labels, where each unique label corresponds to one of the subjects
# X = Nxp design matrix with columns [X_1, ..., X_p], where the kth column
#     is (x_{1k}(t_11), ..., x_{1k}(t_1m_1), ..., x_{nk}(t_n1), ..., x_{nk}(t_nm_n))
# n_basis = number of basis functions to use. Default is n_basis=8.
# lambda0 = grid of spike hyperparameters. Default is (300, 290, ..., 20, 10) 
# lambda1 = slab hyperparameter. Default is 1.
# a = hyperparameter in B(a,b) prior on theta. Default is 1.
# b = hyperparameter in B(a,b) prior on theta. Default is p.
# c0 = hyperparameter in Inverse-Gamma(c_0/2,d_0/2) prior on sigma2. Default is 1.
# d0 = hyperparameter in Inverse-Gamma(c_0/2,d_0/2) prior on sigma2. Default is 1.
# nu = degrees of freedom for Inverse-Wishart prior on Omega. Default is n_basis+2.
# Phi = scale matrix in the Inverse-Wishart prior on Omega. Default is the identity matrix.
# include_intercept = whether or not to include intercept function beta_0(t).
#                     Default is TRUE.
# tol = convergence criteria. Default is 1e-6.
# max_iter = maximum number of iterations to run ECM algorithm. Default is 100.
# return_CI = whether or not to return the 95% pointwise credible bands. Set return_CI=TRUE
#             if credible bands are desired.
# approx_MCMC = whether or not to run the approximate MCMC algorithm instead of the exact
#               MCMC algorithm. If approx_MCMC=TRUE, then the approximate MCMC algorithm is 
#               used. Otherwise. the exact MCMC algorithm is used. This argument is ignored
#               if return_CI=FALSE.
# n_samples = number of MCMC samples to save for posterior inference. Default is 1500.
# burn = number of initial MCMC samples to discard during the warmup pereiod. Default is 500.
# print_iter = Prints the progress of the algorithms. Default is TRUE.

# OUTPUT:
# t_ordered = all N time points ordered from smallest to largest. Needed for plotting
# classifications = px1 vector of indicator variables, where "1" indicates that the covariate is 
#                   selected and "0" indicates that it is not selected.
#                   These classifications are determined by the optimal lambda0 chosen from BIC.
#                   Note that this vector does not include an intercept function.
# beta_hat = Nxp matrix of the varying coefficient function MAP estimates for the optimal lambda0 
#            chosen from BIC. The kth column in the matrix is the kth estimated function 
#            at the observation times in t_ordered
# beta0_hat =  intercept function MAP estimate at the observation times in t_ordered for the
#              optimal lambda0 chosen from BIC. This is not returned if include_intercept = FALSE.
# gamma_hat = MAP estimates for basis coefficients (needed for prediction) using the optimal lambda0.
# beta_post_mean = Nxp matrix of the varying coefficient function posterior mean estimates.
#                  The kth column in the matrix is the kth posterior mean function estimate at
#                  the observation times in t_ordered. This is not returned if return_CI=FALSE.
# beta_CI_lower = Nxp matrix of the lower endpoints of the 95% pointwise posterior credible interval (CI) for
#                 the varying coefficient functions. The kth column in the matrix is the lower
#                 endpoint for the kth varying coefficient's CI at the observation times in
#                 t_ordered. This is not returned if return_CI=FALSE.
# beta_CI_upper = Nxp matrix of the upper endpoints of the 95% pointwise posterior credible interval (CI) for
#                 the varying coefficient functions. The kth column in the matrix is the upper
#                 endpoint for the kth varying coefficient's CI at the observation times in
#                 t_ordered. This is not returned if return_CI=FALSE.
# beta0_post_mean = intercept function posterior mean estimate. This is only returned if return_CI=TRUE
#                   and include_intercept=TRUE.
# beta0_CI_lower = lower endpoints of the 95% pointwise posterior credible intervals for the intercept function.
#                 This is only returned if return_CI=TRUE and include_intercept=TRUE.
# beta0_CI_upper = upper endpoints of the 95% pointwise posterior credible intervals for the intercept function.
#                 This is only returned if return_CI=TRUE and include_intercept=TRUE.
# gamma_post_mean = estimated posterior means for basis coefficients. This is not returned if return_CI=FALSE.
# gamma_CI_lower = lower endpoints of the 95% posterior credible intervals for the basis coefficients.
# gamma_CI_upper = upper endpoints of the 95% posterior credible intervals for the basis coefficients.
# post_incl = estimated posterior inclusion probabilities for each of the varying coefficients
# lambda0_min = lambda0 which minimizes the BIC for the MAP estimator. If only one value was pass 
#               for lambda0, then this just returns that lambda.
# lambda0_all = grid of all L lambda0's.
# BIC_all = Lx1 vector of BIC values corresponding to all L lambda0's in lambda0_all.
#           The lth entry corresponds to the lth entry in lambda0_all.
# beta_est_all_lambda0 = list of length L of the MAP estimates for the varying coefficients 
#                        corresponding to all L lambda0's in lambda0_all. The lth entry corresponds 
#                        to the lth entry in lambda0_all.
# beta0_est_all_lambda0 = NxL matrix of MAP estimates for the intercept function corresponding 
#                         to all L lambda0's in lambda0_all. The lth column corresponds to the lth
#                        entry in lambda0_all.
# gamma_est_all_lambda0 = dpxL matrix of MAP estimates for the basis coefficients corresponding to all
#                        the lambda0's in lambda0_all. The lth column corresponds to the lth
#                        entry in lambda0_all.
# classifications_all_lambda0 = pxL matrix of classifications corresponding to all the 
#                              lambda0's in lambda0_all. The lth column corresponds to the
#                              lth entry in lambda0_all.
# ECM_iters_to_converge = Lx1 vector of the number of iterations it took for the ECM algorithm 
#                         to converge for each lambda0 in lambda0_all. The lth entry corresponds
#                         to the lth entry in lambda0_all.
# ECM_runtimes = Lx1 vector of the number of seconds it took for the ECM algorithm to converge
#                for each lambda0 in lambda0_all. The lth entry corresponds to the lth entry
#                in lambda0_all.
# gibbs_runtime = number of minutes it took for the Gibbs sampling algorithm to run total
#                 number of MCMC iterations given in gibbs_iters
# gibbs_iters = total number of MCMC iterations


NVC_SSL = function(y, t, id, X, n_basis=8, 
                   lambda0=seq(from=300,to=10,by=-10), lambda1=1, 
                   a=1, b=ncol(X), c0=1, d0=1, nu=n_basis+2, Phi=diag(n_basis),
                   include_intercept=TRUE, tol=1e-6, max_iter=100, 
                   return_CI=FALSE, approx_MCMC=FALSE,
                   n_samples=1500, burn=500, print_iter=TRUE) {
  
  ################
  ## PRE-CHECKS ##
  ################
  
  ## Check that dimensions are conformal
  if( length(t) != length(y) | length(y) != dim(X)[1] | length(t) != dim(X)[1] | length(t) != length(id) )
    stop("Non-conformal dimensions for y, t, id, and X.")
  ## Check that number of basis functions is >=3
  if( n_basis <= 2 )
    stop("Please enter a positive integer greater than or equal to three for number of basis functions.")
  ## Check that all relevant hyperparameters are positive
  if ( !(all(lambda0 > 0) & (lambda1 > 0) & (a > 0) & (b > 0) & (c0 > 0) & (d0 > 0) 
         & (nu > n_basis-1) & all(eigen(Phi)$values > 0)) )
    stop("Inappropriate values for some hyperparameters.")
  if (any(lambda0<=lambda1))
    stop("Please make sure that lambda1 is smaller than lambda0.")
  if(tol<=0)
    stop("Please set a tolerance greater than 0.")
  if(max_iter<10)
    stop("Please set maximum number of iterations to at least 10.")
  if( (n_samples<=0) || (burn<0))
    stop("Please enter a valid number of samples and burnin for the Gibbs sampler.")
  
  ###################
  ## PREPROCESSING ##
  ###################
  
  ## Make id consecutive integers
  id = as.integer(id)
  id = cumsum(c(0, diff(id)) != 0) + 1
  ## Sort data according to id
  order_id = order(id)
  id = id[order(id)]
  y = y[order_id]
  t = t[order_id]
  X = X[order_id, ]
  ## Make id a factor
  id = as.factor(id)
  
  ## Make n_basis and max_iter an integer if not already
  n_basis = as.integer(n_basis)
  max_iter = as.integer(max_iter)
  n_samples = as.integer(n_samples)
  burn = as.integer(burn)
  
  ##################################
  # Construct appropriate matrices #
  # and vectors                    #
  ##################################
  
  ## Extract total # of observations, # of subjects, # of observations per subject, # of covariates, and # of basis functions
  n_tot = length(y)           # total number of observations
  n_i = plyr::count(id)$freq        # number of within-subject observations per subject
  ## Make X a matrix if not already done
  ## If include_intercept=TRUE, augment X with a column of ones
  if(include_intercept==TRUE){
    X = cbind(rep(1,n_tot), X)
  }
  X = as.matrix(X)
  n_sub = length(unique(id))  # number of subjects
  p = dim(X)[2]               # number of covariates
  d = n_basis                  # dimension of the vectors of basis coefficients
  
  ## Create Nxd spline matrix with d degrees of freedom
  B = splines::bs(t, df=d, intercept=TRUE)
  
  ## Matrix for storing individual B(t_ij) matrices
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
  ## Create the y_i subvectors and U_i submatrices
  y_i = vector("list", n_sub)
  U_i = vector("list", n_sub)
  
  for(r in 1:n_sub){
    y_i[[r]] = matrix(y[which(id==r)])
    U_i[[r]] = U[which(id==r),]
  }
  
  ## Create the Z_i matrices
  Z = splines::bs(t, df=d, intercept=FALSE)
  Z_i = vector("list", n_sub)
  for(r in 1:n_sub){
    Z_i[[r]] = Z[which(id==r),]
  }
  ## Make Z the direct sum of the Z_i's
  Z = dae::mat.dirsum(Z_i)
  
  ## Time-saving
  UtU = crossprod(U,U)
  Uty = crossprod(U,y)
  UtZ = crossprod(U,Z)
  
  Zit_Zi = vector("list", n_sub)
  Zit_yi = vector("list", n_sub)
  Zit_Ui = vector("list", n_sub)
  
  for(r in 1:n_sub){
    Zit_Zi[[r]] = crossprod(Z_i[[r]], Z_i[[r]])
    Zit_yi[[r]] = crossprod(Z_i[[r]], y_i[[r]])
    Zit_Ui[[r]] = crossprod(Z_i[[r]], U_i[[r]])
  }
  
  ###################
  ## ECM algorithm ##
  ###################
  
  ## Sort lambda0
  lambda0 = sort(lambda0, decreasing = TRUE)
  # lambda0 = sort(lambda0)
  ## Initialize gamma_old
  gamma_old = rep(0, d*p)
  ## Initialize theta_old
  theta_old = 0.5

  ## For storing gamma, beta, beta0, classifications, BIC, iterations to converge
  gamma_mat = matrix(0, nrow=length(gamma_old), ncol=length(lambda0))
  beta_list = vector(mode='list', length=length(lambda0))
  if(include_intercept==TRUE){
    beta0_mat = matrix(0, nrow=n_tot, ncol=length(lambda0))
    classifications_mat = matrix(0, nrow=p-1, ncol=length(lambda0))
  } else{
    classifications_mat = matrix(0, nrow=p, ncol=length(lambda0))
  }
  BIC_vec = rep(0, length(lambda0))
  iters_to_converge_vec = rep(0, length(lambda0))
  ECM_runtimes_vec = rep(0, length(lambda0))
    
  ## Begin algorithm
  for(l in 1:length(lambda0)){
    ## Print iteration number
    if(print_iter==TRUE)
      cat("lambda0 =", lambda0[l], "\n")
    
    ## Initialize
    gamma_init = gamma_old
    theta_init = theta_old
    lambda0_current = lambda0[l]

    start_time = Sys.time()
    
    ## ECM algorithm to update gamma_old
    EM_mod = NVC_SSL_EM(y=y, t=t, Z=Z, U=U, B=B, Zit_Zi=Zit_Zi, Z_i=Z_i, y_i=y_i, U_i=U_i,  
                        lambda0=lambda0_current, lambda1=lambda1, 
                        n_basis=d, n_cov=p, n_sub=n_sub, n_tot=n_tot, 
                        a=a, b=b, c0=c0, d0=d0, nu=nu, Psi=Phi,
                        tol=tol, max_iter=max_iter, include_intercept=include_intercept,
                        gamma_init=gamma_init, theta_init=theta_init)
    
    end_time = Sys.time()
    
    ## Store ECM runtime in seconds
    ECM_runtimes_vec[l] <- as.double(difftime(end_time, start_time, units="secs"))
    
    ## Update gamma_old
    gamma_old = EM_mod$gamma_hat
    ## Update theta_old
    theta_old = EM_mod$theta_hat
    
    ## Store results
    gamma_mat[,l] = gamma_old
    beta_list[[l]] = EM_mod$beta_hat
    if(include_intercept==TRUE){
      beta0_mat[,l] = EM_mod$beta0_hat
    }
    classifications_mat[,l] = EM_mod$classifications
    BIC_vec[l] = EM_mod$BIC
    iters_to_converge_vec[l] = EM_mod$iters_to_converge
  }
  
  ## lambda which minimizes the BIC
  min_index = as.double(which.min(BIC_vec))
  lambda0_min = lambda0[min_index]
  
  ## Best model according to BIC
  beta_hat = beta_list[[min_index]]
  if(include_intercept==TRUE){
    beta0_hat = beta0_mat[,min_index]
  }
  gamma_hat = gamma_mat[,min_index]
  classifications = classifications_mat[,min_index]  
  
  ## If no posterior credible intervals are desired
  if(return_CI==FALSE){
    
    if(include_intercept==TRUE){
      NVC_SSL_output <- list(t_ordered = EM_mod$t_ordered,
                             classifications = classifications,
                             beta_hat = beta_hat,
                             beta0_hat = beta0_hat,
                             gamma_hat = gamma_hat,
                             lambda0_min = lambda0_min,
                             lambda0_all = lambda0,
                             BIC_all = BIC_vec,
                             beta_est_all_lambda0 = beta_list,
                             beta0_est_all_lambda0 = beta0_mat,
                             gamma_est_all_lambda0 = gamma_mat,
                             classifications_all_lambda0 = classifications_mat,
                             ECM_iters_to_converge = iters_to_converge_vec,
                             ECM_runtimes = ECM_runtimes_vec)
      return(NVC_SSL_output)
      
    } else if(include_intercept==FALSE){
      NVC_SSL_output <- list(t_ordered = EM_mod$t_ordered,
                             classifications = classifications,
                             beta_hat = beta_hat,
                             gamma_hat = gamma_hat,
                             lambda0_min = lambda0_min,
                             lambda0_all = lambda0,
                             BIC_all = BIC_vec,
                             beta_est_all_lambda0 = beta_list,
                             gamma_est_all_lambda0 = gamma_mat,
                             classifications_all_lambda0 = classifications_mat,
                             ECM_iters_to_converge = iters_to_converge_vec,
                             ECM_runtimes = ECM_runtimes_vec)
      
      return(NVC_SSL_output)
    }
  }
  
  ## If posterior credible intervals are desired, then 
  ## we will also run the Gibbs sampler
  if(return_CI==TRUE){
    
    ## Run Gibbs sampler initialized at the estimate which minimizes the BIC
    
    gamma_init = gamma_hat
    theta_init = theta_old
    
    start_time = Sys.time()
    
    gibbs_mod = NVC_SSL_gibbs(y=y, t=t, U=U, B=B, Z=Z, UtU=UtU, Uty=Uty, UtZ=UtZ, 
                              Zit_Zi=Zit_Zi, Zit_yi=Zit_yi, Zit_Ui=Zit_Ui,
                              lambda0=5, lambda1=1, 
                              d=d, p=p, n_sub=n_sub, n_tot=n_tot,
                              a=a, b=b, c0=c0, d0=d0, nu=nu, Psi=Phi,
                              n_samples=n_samples, burn=burn, thres=0.5, 
                              include_intercept=include_intercept, print_iter=print_iter, 
                              gamma_init=gamma_init, theta_init=theta_init, approx_MCMC=approx_MCMC)

    end_time = Sys.time()
    
    gibbs_runtime = as.double(difftime(end_time, start_time, units="mins"))
    
    if(include_intercept==TRUE){
      
      NVC_SSL_output <- list(t_ordered = EM_mod$t_ordered,
                             classifications = classifications,
                             beta_hat = beta_hat,
                             beta0_hat = beta0_hat,
                             gamma_hat = gamma_hat,
                             beta_post_mean = gibbs_mod$beta_hat,
                             beta_CI_lower = gibbs_mod$beta_CI_lower,
                             beta_CI_upper = gibbs_mod$beta_CI_upper,
                             beta0_post_mean = gibbs_mod$beta0_hat,
                             beta0_CI_lower = gibbs_mod$beta0_CI_lower,
                             beta0_CI_upper = gibbs_mod$beta0_CI_upper,
                             gamma_post_mean = gibbs_mod$gamma_hat,
                             gamma_CI_lower = gibbs_mod$gamma_CI_lower,
                             gamma_CI_upper = gibbs_mod$gamma_CI_upper,
                             post_incl = gibbs_mod$post_incl,
                             lambda0_min = lambda0_min,
                             lambda0_all = lambda0,
                             BIC_all = BIC_vec,
                             beta_est_all_lambda0 = beta_list,
                             beta0_est_all_lambda0 = beta0_mat,
                             gamma_est_all_lambda0 = gamma_mat,
                             classifications_all_lambda0 = classifications_mat,
                             ECM_iters_to_converge = iters_to_converge_vec,
                             ECM_runtimes = ECM_runtimes_vec,
                             gibbs_runtime = gibbs_runtime,
                             gibbs_iters = gibbs_mod$gibbs_iters)
      
      return(NVC_SSL_output)
      
    } else if(include_intercept==FALSE){
      NVC_SSL_output <- list(t_ordered = EM_mod$t_ordered,
                             classifications = classifications,
                             beta_hat = beta_hat,
                             gamma_hat = gamma_hat,
                             beta_post_mean = gibbs_mod$beta_hat,
                             beta_CI_lower = gibbs_mod$beta_CI_lower,
                             beta_CI_upper = gibbs_mod$beta_CI_upper,
                             gamma_post_mean = gibbs_mod$gamma_hat,
                             gamma_CI_lower = gibbs_mod$gamma_CI_lower,
                             gamma_CI_upper = gibbs_mod$gamma_CI_upper,
                             post_incl = gibbs_mod$post_incl,
                             lambda0_min = lambda0_min,
                             lambda0_all = lambda0,
                             BIC_all = BIC_vec,
                             beta_est_all_lambda0 = beta_list,
                             gamma_est_all_lambda0 = gamma_mat,
                             classifications_all_lambda0 = classifications_mat,
                             ECM_iters_to_converge = iters_to_converge_vec,
                             ECM_runtimes = ECM_runtimes_vec,
                             gibbs_runtime = gibbs_runtime,
                             gibbs_iters = gibbs_mod$gibbs_iters)
      
      return(NVC_SSL_output)
    }
  }
}