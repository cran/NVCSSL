###############################################
## Gives predictions for new subjects based  ##
## on an estimated varying coefficient model ##
###############################################


# INPUTS:
# NVC_mod = an object with a fitted NVC model
# t_new = vector of new observation times
# id_new = vector of new labels, where a label corresponds to one of the new subjects
# X_new = new design matrix with columns [X_1, ..., X_p], where the kth column
#         corresponds to the kth covariate

# OUTPUT:
# id = subject ID
# time = within-subject observation time
# y_pred = predicted response, where the (ij)th entry is the predicted y for the jth
#          observation time of subject i

NVC_predict = function(NVC_mod, t_new, id_new, X_new) {
  
  ## Need gamma_hat
  gamma_hat = NVC_mod$gamma_hat
  
  ## Get n_tot, p, and d=degrees of freedom
  n_tot = dim(X_new)[1]
  
  ## Check if there is an intercept function
  if(is.null(NVC_mod$beta0_hat)){
    p = dim(X_new)[2]
    d = length(gamma_hat)/p

    ## Make X_new a matrix if not already done
    X_new = as.matrix(X_new)
  } else {
    ## If there is an intercept included
    p = dim(X_new)[2]+1
    d = length(gamma_hat)/p
    
    ## Add a column of ones to X_new and make X_new a matrix
    X_new = cbind(rep(1,n_tot), X_new)
    X_new = as.matrix(X_new)
  }
  
  ############################
  ## Construct new U matrix ##
  ############################

  ## First create Nxd spline matrix with d degrees of freedom
  B = splines::bs(t_new, df=d, intercept=TRUE)
  
  ## Matrix for storing individual B(t_ij) matrices
  B_t = matrix(0, nrow=p, ncol=d*p)
  
  ## Create the U matrix. U.new is an Nxdp matrix
  U_new = matrix(0, nrow=n_tot, ncol=d*p) 
  for(r in 1:n_tot){
    for(j in 1:p){
      B_t[j,(d*(j-1)+1):(d*j)] = B[r,]
    }
    U_new[r,] = X_new[r,] %*% B_t
  }
  ## Free memory
  rm(B_t)
  
  #########################
  ## Compute predictions ##
  #########################
  
  y_pred = U_new %*% gamma_hat

  results = list(id=id_new, time = t_new, y_pred = y_pred)
  
  ## Return results
  return(results)
}