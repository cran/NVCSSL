#############################################
## ECM algorithm for finding MAP estimator ##
#############################################
#' @keywords internal
NVC_SSL_EM = function(y, t, Z, U, B, Zit_Zi, Z_i, y_i, U_i,
                      lambda0, lambda1, n_basis, n_cov, n_sub, n_tot,
                      a, b, c0, d0, nu, Psi,
                      tol, max_iter, include_intercept,
                      gamma_init=NULL, theta_init=NULL) {
  
  
  ## Initialize values
  d = n_basis
  p = n_cov
  eta = rep(0, d*n_sub)            # to hold eta (vector)
  # Initialize gamma vector
  if(is.null(gamma_init)){
    gamma_vec = rep(0, d*p)
  } else if(!is.null(gamma_init)){
    gamma_vec = gamma_init
  }
  ## Initialize theta
  if(is.null(theta_init)){
    theta = 0.5
  } else if(!is.null(theta_init)){
    theta = theta_init
  }
  groups = rep(1:p, each=d)        # for groups of basis coefficients
  ## Initialize Omega
  Omega = diag(nrow=d)             
  ## Initialize sigma2
  df_sigma2 <- 3
  sig_quant <- 0.9
  sig_est <- stats::sd(y)
  qchi <- stats::qchisq(1 - sig_quant, df_sigma2)
  ncp <- sig_est^2 * qchi / df_sigma2
  # Initial overestimate based on Chipman et al. (2010)
  init_overestimate <- (df_sigma2 * ncp)/(df_sigma2 + 2)
  if(sig_est >= 1){
    sigma2 = init_overestimate
  } else {
    sigma2 = 1
  }
  
  ## Initialize the following values
  difference = 100*tol
  counter = 0
  pstar = rep(0,p)            # to hold pstar_k's
  lambdastar = rep(0,p)       # to hold lambdastar_k's
  
  
  ###############################
  ## EM algorithm for updating ## 
  ## the parameters            ##
  ###############################
  while( (difference > tol) & (counter < max_iter) ){
    
    counter = counter+1
    
    ## Keeps track of old gamma
    gamma_old = gamma_vec
    
    ############
    ## E-step ##
    ############
    
    for(k in 1:p){
      ## Which groups are active
      active = which(groups == k)
      ## Update pstar[k]
      pi1 = theta*lambda1^d*exp(-lambda1*norm(gamma_old[active], type="2"))
      pi2 = (1-theta)*lambda0^d*exp(-lambda0*norm(gamma_old[active], type="2"))
      if(pi1 == 0 & pi2 == 0){
        pstar[k] = 1
      } else {
        pstar[k] = pi1/(pi1+pi2)
      }
      # Update lambdastar[k]
      lambdastar[k] = lambda1*pstar[k] + lambda0*(1-pstar[k])
    }
    
    #############
    ## CM-step ##
    #############
    
    ## Update theta
    theta = (a-1 + sum(pstar))/(a+b+p-2)
    
    ## Update eta
    for(r in 1:n_sub){
      
      ## Update Bz_i
      Bz_i = Matrix::chol2inv(Matrix::chol(Zit_Zi[[r]]+sigma2*Matrix::chol2inv(Matrix::chol(Omega))))%*%t(Z_i[[r]])
      
      ## Update eta_i
      eta[(d*(r-1)+1):(d*r)] = Bz_i%*%(y_i[[r]]-U_i[[r]]%*%gamma_old)
    }
    
    ## Update gamma
    ## Note that grpreg solves is (1/2n)*||Y.bar-U.bar%*%gamma||_2^2 + pen(gamma),
    ## so we have to divide multiplier by n_tot
    y_tilde = y-Z%*%eta
    gamma_mod = grpreg::grpreg(U, y_tilde, group=groups, lambda=1, 
                               group.multiplier=sigma2*lambdastar/n_tot)
    gamma_vec = gamma_mod$beta[-1] 
    
    ## The indices which have lambdastar approximately lambda0 should also be 
    # thresholded to zero.
    # Give gamma_vec some time to 'warm up', so we only do this after 5 iterations
    # have been run.
    if( (lambda0>=5) & (counter > 5) ){
      lambdastar_inds = which(abs(lambdastar-lambda0)<1e-6)
      if(length(lambdastar_inds)>0){
        gamma_vec[which(groups%in%lambdastar_inds)] = 0
      }
    }
    
    ## Update Omega
    eta_etat = matrix(0, nrow=d, ncol=d)
    for(r in 1:n_sub){
      etai_etait= eta[(d*(r-1)+1):(d*r)]%*%t(eta[(d*(r-1)+1):(d*r)])
      eta_etat = eta_etat + etai_etait
    }
    Omega = (Psi+eta_etat)/(n_sub+nu+d+1)
    
    ## Update sigma2
    resid = y-U%*%gamma_vec-Z%*%eta
    sigma2 = as.double((crossprod(resid,resid)+d0)/(n_tot+c0+2))
    
    ## Update diff
    difference = norm(gamma_vec-gamma_old, type="2")^2
  }
  ## End of EM algorithm
  
  
  ## Calculate the estimated beta_k(t)'s and store classifications
  beta_hat = matrix(0, nrow=n_tot, ncol=p)
  beta_ind_hat = rep(0, n_tot)
  classifications = rep(0, p)
  
  for(k in 1:p){
    ## Update beta_hat
    beta_ind_hat = B %*% gamma_vec[((k-1)*d+1):(k*d)]
    beta_hat[,k] = beta_ind_hat[order(t)]
    
    ## Update classifications
    if(!identical(beta_hat[,k], rep(0,n_tot)))
      classifications[k] = 1    
  } 
  
  ## Order the times
  t_ordered = t[order(t)]
  
  ## Calculate the BIC
  mu = U%*%gamma_vec
  Sigma = sigma2*diag(n_tot) + Z%*%kronecker(diag(n_sub), Omega)%*%t(Z)
  
  logLik = mvtnorm::dmvnorm(y, mean=mu, sigma=Sigma, log=TRUE)
  num_nonzero = sum(classifications)*d
  BIC = -2*logLik + log(n_tot)*num_nonzero
  
  ## Return list
  if(include_intercept==FALSE){
    NVC_EM_output <- list(t_ordered = t_ordered,
                          beta_hat = beta_hat,
                          gamma_hat = gamma_vec,
                          theta_hat = theta,
                          classifications = classifications,
                          BIC = BIC,
                          iters_to_converge = counter)
    
  } else if(include_intercept==TRUE){
    ## Extract the intercept function
    beta0_hat = beta_hat[,1]
    beta_hat = beta_hat[,-1]
    classifications = classifications[-1]
    
    NVC_EM_output <- list(t_ordered = t_ordered,
                          beta_hat = beta_hat,
                          beta0_hat = beta0_hat,
                          gamma_hat = gamma_vec,
                          theta_hat = theta,
                          classifications = classifications,
                          BIC = BIC,
                          iters_to_converge = counter)
  }
  
  return(NVC_EM_output)
}


######################################################
## Gibbs sampling algorithm for posterior inference ##
######################################################
#' @keywords internal
NVC_SSL_gibbs = function(y, t, U, B, Z, UtU, Uty, UtZ, Zit_Zi, Zit_yi, Zit_Ui,
                         lambda0, lambda1, d, p, n_sub, n_tot,
                         a, b, c0, d0, nu, Psi,
                         n_samples, burn, thres,
                         include_intercept, print_iter, 
                         gamma_init=NULL, theta_init=NULL, approx_MCMC=FALSE) {
  
  
  eta = rep(0, d*n_sub)            # to hold eta (vector)
  # to hold gamma (vector)
  if(is.null(gamma_init)){
    gamma_vec = rep(0, d*p)
  } else if(!is.null(gamma_init)){
    gamma_vec = gamma_init
  }
  # to hold theta
  if(is.null(theta_init)){
    theta = 0.5
  } else if(!is.null(theta_init)){
    theta = theta_init
  }
  groups = rep(1:p, each=d)        # for groups of basis coefficients
  tau = rep(0, p)                  # to hold tau (vector)
  xi = rep(1,p)                    # to hold xi (vector)
  lambda_star = rep(0,d)           # to hold lambda_star (vector)
  Omega = diag(nrow=d)             # to hold Omega
  sigma2 = 1                       # to hold sigma2
  
  ## To save the gamma and tau samples
  n_iterations = n_samples+burn
  tau_samples = matrix(0, n_iterations, p)
  gamma_samples = matrix(0, n_iterations, d*p)
  
  ###########################
  # Start the Gibbs sampler #
  ###########################
  j = 0
  while (j < n_iterations){
    
    ## Print iteration number
    j = j + 1
    if ( (print_iter==TRUE) & (j %% 100 == 0) ) {
      cat("Gibbs sampling iteration:", j, "\n")
    }
    
    ################
    ## Sample eta ##
    ################
    for(r in 1:n_sub){
      ridge = Zit_Zi[[r]]/sigma2 +Matrix::chol2inv(Matrix::chol(Omega))
      Sigma_i = Matrix::chol2inv(Matrix::chol(ridge))
      mu_i = Sigma_i %*% (Zit_yi[[r]]-Zit_Ui[[r]]%*%gamma_vec)/sigma2
      eta[(d*(r-1)+1):(d*r)] = MASS::mvrnorm(1, mu_i, Sigma_i)                   
    }
    
    ##################
    ## Sample Omega ##
    ##################
    Omega_df = nu+n_sub
    eta_etat = matrix(0, nrow=d, ncol=d)
    for(r in 1:n_sub){
      etai_etait= eta[(d*(r-1)+1):(d*r)]%*%t(eta[(d*(r-1)+1):(d*r)])
      eta_etat = eta_etat + etai_etait
    }
    Omega_scale = Psi+eta_etat
    Omega = MCMCpack::riwish(Omega_df, Omega_scale)
    
    ###################
    ## Sample sigma2 ##
    ###################
    sigma2_shape = (n_tot+c0)/2
    sigma2_rate = (sum((y-U%*%gamma_vec-Z%*%eta)^2)+d0)/2
    sigma2 = 1/stats::rgamma(1, sigma2_shape, sigma2_rate)
    
    #######################
    ## Sample tau and xi ##
    #######################
    
    for(k in 1:p){
      pi1 = theta*lambda1^(d+1)*exp(-lambda1^2*xi[k]/2)
      pi2 = (1-theta)*lambda0^(d+1)*exp(-lambda0^2*xi[k]/2)
      if(pi1 == 0 & pi2 == 0){
        prop = 1
      } else {
        prop = pi1/(pi1+pi2)
      }
      # Sample tau[k]
      tau[k] = stats::rbinom(1, 1, prop)
      # set lambda[k]
      lambda_star[k] = tau[k]*lambda1 + (1-tau[k])*lambda0
      # sample xi[k]
      xi[k] = GIGrvg::rgig(1, lambda=0.5, 
                           chi=sum(gamma_vec[which(groups==k)]^2), 
                           psi=lambda_star[k]^2)
    }
    ## Store tau samples
    tau_samples[j, ] = tau
    
    ##################
    ## Sample theta ##
    ##################
    sum_tau = sum(tau)
    theta = stats::rbeta(1, a+sum_tau, b+p-sum_tau)
    
    ################## 
    ## Sample gamma ##
    ##################
    xi_diag = as.vector(sapply(xi, function(x) rep(x,d)))
    
    if (d*p <= n_tot){
      ridge_gamma = UtU/sigma2 + diag(1/xi_diag)
      inv_ridge_gamma = Matrix::chol2inv(Matrix::chol(ridge_gamma))
      mu_gamma = (inv_ridge_gamma %*% (Uty -UtZ%*%eta))/sigma2  
      gamma_vec = MASS::mvrnorm(1, mu_gamma, inv_ridge_gamma)
    } else if (d*p > n_tot){
      
      if(approx_MCMC==FALSE){
        # Use the exact Bhattacharya et al. (2016) algorithm
        sigma = sqrt(sigma2)
        little_u = stats::rnorm(d*p, mean=rep(0,d*p), sd=sqrt(xi_diag))
        delta = stats::rnorm(n_tot, mean=0, sd=1)
        v = (U/sigma)%*%little_u+delta
        v_star = (y-Z%*%eta)/sigma-v
        Phi = (U/sigma2)%*%(xi_diag*t(U))+diag(n_tot)
        w = Matrix::chol2inv(Matrix::chol(Phi))%*%v_star
        gamma_vec = little_u + (xi_diag*t(U/sigma))%*%w
      } else {
        # Use modified Johndrow et al. (2020) algorithm
        sigma = sqrt(sigma2)
        tau_nonzero = which(tau!=0)
        tau_zero = which(tau==0)
        little_u = stats::rnorm(d*p, mean=rep(0,d*p), sd=sqrt(xi_diag))
        delta = stats::rnorm(n_tot, mean=0, sd=1)
        v = (U/sigma)%*%little_u+delta
        v_star = (y-Z%*%eta)/sigma-v
        
        if(length(tau_nonzero)!=0){
          ## Low-rank approximation
          S = which(groups %in% tau_nonzero)
          S_C = which(groups %in% tau_zero)
          U_S = U[,S]
          approx_ridge = diag(1/xi_diag[S]) + crossprod(U_S, U_S)/sigma2
          Phi_inv = diag(n_tot) - U_S%*%Matrix::chol2inv(Matrix::chol(approx_ridge))%*%t(U_S/sigma2)
          w = Phi_inv%*%v_star
          new_xi_diag = xi_diag
          new_xi_diag[S_C] = 0 
          gamma_vec = little_u + (new_xi_diag*t(U/sigma))%*%w
        } else {
          w = Matrix::chol2inv(Matrix::chol((U/sigma2)%*%(xi_diag*t(U))+diag(n_tot)))%*%v_star
          gamma_vec = little_u + (xi_diag*t(U/sigma))%*%w
        }
      }
    }
    ## Save gamma samples
    gamma_samples[j, ] = gamma_vec
  }
  ## End of Gibbs sampler  
  
  #####################################
  ## Extract the posterior summaries ##
  #####################################
  ## Discard burn-in
  tau_samples = tau_samples[(burn+1):n_iterations, ]
  all_gamma_samples = gamma_samples
  gamma_samples = gamma_samples[(burn+1):n_iterations, ]
  
  ## Estimated posterior inclusion probabilities
  post_incl = colMeans(tau_samples)
  if(include_intercept==TRUE) {
    post_incl = post_incl[-1]
  }
  
  ## Classifications according to the thresholding rule
  classifications = rep(0, p)
  classifications[which(post_incl>=thres)] = 1
  if(include_intercept==TRUE){
    classifications = classifications[-1]
  }
  
  ## Extract the posterior mean, 2.5th quantile, and 97.5th quantile for gamma
  gamma_hat = colMeans(gamma_samples)
  gamma_intervals = apply(gamma_samples, 2, function(x) stats::quantile(x, prob=c(.025,.975)))
  gamma_CI_lower = gamma_intervals[1,]
  gamma_CI_upper = gamma_intervals[2,]
  
  ## Order the times
  t_ordered = t[order(t)]
  # For updating the beta_k(t)'s
  gamma_hat_list = split(gamma_hat, as.numeric(gl(length(gamma_hat),d,length(gamma_hat))))
  gamma_CI_lower_list = split(gamma_CI_lower, as.numeric(gl(length(gamma_CI_lower),d,length(gamma_CI_lower))))
  gamma_CI_upper_list = split(gamma_CI_upper, as.numeric(gl(length(gamma_CI_upper),d,length(gamma_CI_upper))))
  
  ## Calculate the estimated beta_k(t)'s and store classifications
  beta_ind_hat = rep(0, n_tot)
  beta_hat = matrix(0, nrow=n_tot, ncol=p)
  beta_CI_lower = matrix(0, nrow=n_tot, ncol=p)
  beta_CI_upper = matrix(0, nrow=n_tot, ncol=p)
  
  for(k in 1:p){
    ## Update beta_hat
    beta_ind_hat = B %*% gamma_hat_list[[k]]
    beta_hat[,k] = beta_ind_hat[order(t)]
    
    ## Update beta_CI_lower
    beta_ind_hat = B %*% gamma_CI_lower_list[[k]]
    beta_CI_lower[,k] = beta_ind_hat[order(t)]
    
    ## Update beta_CI_upper
    beta_ind_hat = B %*% gamma_CI_upper_list[[k]]
    beta_CI_upper[,k] = beta_ind_hat[order(t)]
  }
  
  if(include_intercept==FALSE){
    NVC_gibbs_output <- list(t_ordered = t_ordered,
                             beta_hat = beta_hat,
                             gamma_hat = gamma_hat,
                             classifications = classifications,
                             post_incl = post_incl,
                             beta_CI_lower = beta_CI_lower,
                             beta_CI_upper = beta_CI_upper,
                             gamma_CI_lower = gamma_CI_lower,
                             gamma_CI_upper = gamma_CI_upper,
                             all_gamma_samples = all_gamma_samples,
                             gibbs_iters = n_iterations)
    
  } else if(include_intercept==TRUE){ 
    
    ## Extract the intercept function
    beta0_hat = beta_hat[,1]
    beta_hat = beta_hat[,-1]
    beta0_CI_lower = beta_CI_lower[,1]
    beta_CI_lower = beta_CI_lower[,-1]
    beta0_CI_upper = beta_CI_upper[,1]
    beta_CI_upper = beta_CI_upper[,-1]
    
    NVC_gibbs_output <- list(t_ordered = t_ordered,
                             beta_hat = beta_hat,
                             beta0_hat = beta0_hat,
                             gamma_hat = gamma_hat,
                             classifications = classifications,
                             post_incl = post_incl,
                             beta_CI_lower = beta_CI_lower,
                             beta_CI_upper = beta_CI_upper,
                             beta0_CI_lower = beta0_CI_lower,
                             beta0_CI_upper = beta0_CI_upper,
                             gamma_CI_lower = gamma_CI_lower,
                             gamma_CI_upper = gamma_CI_upper,
                             all_gamma_samples = all_gamma_samples,
                             gibbs_iters = n_iterations)
  }
  
  ## Return a list
  return(NVC_gibbs_output)
}