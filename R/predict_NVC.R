#############################
#############################
## FUNCTION FOR PREDICTION ##  
## USING THE NVC-SSL MODEL ##
#############################
#############################

# INPUTS:
# NVC.mod = a list with an NVC fitted model 
# t.new = Nx1 vector of new observation times t_11, ..., t_1n_1, ..., t_n1, ..., t_nn_n
# id.new = Nx1 vector of new labels (different label corresponds to different subjects). Labels should be (1,...,1, 2,...,2, ..., n,....,n)
# X.new = new Nxp design matrix with columns [X_1, ..., X_p], where the kth column
#     is (x_{1k}(t_11), ..., x_{1k}(t_1n_1), ..., x_{nk}(t_n1), ..., x_{nk}(t_nn_n))'

# OUTPUT:
# y.hat = Nx1 vector of predicted values

predict_NVC = function(NVC.mod, t.new, id.new, X.new) {
  
  ## Need gamma.hat and intercept
  gamma.hat = NVC.mod$gamma.hat
  intercept = NVC.mod$intercept
  
  ## Get n.tot, p, and d=degrees of freedom
  n.tot = dim(X.new)[1]
  p = dim(X.new)[2]
  d = length(gamma.hat)/p
  ## Make 'id' a factor vector if not already done
  id = factor(id.new)
  ## Make X a matrix if not already done
  X.new = as.matrix(X.new)
  
  ############################
  ############################
  ## Construct new U matrix ##
  ############################
  ############################
  
  ## First create Nxd spline matrix with d degrees of freedom
  B = splines::bs(t.new, df=d, intercept=TRUE)
  
  ## Matrix for storing individual B(t_ij) matrices
  B.t = matrix(0, nrow=p, ncol=d*p)
  
  ## Create the U matrix. U.new is an Nxdp matrix
  U.new = matrix(0, nrow=n.tot, ncol=d*p) 
  for(r in 1:n.tot){
    for(j in 1:p){
      B.t[j,(d*(j-1)+1):(d*j)] = B[r,]
    }
    U.new[r,] = X.new[r,] %*% B.t
  }
  ## Free memory
  rm(B.t)
  
  #########################
  #########################
  ## Compute predictions ##
  #########################
  #########################
  
  y.hat = U.new %*% gamma.hat
  
  results = cbind(id.new, y.hat)
  colnames(results) = c("id","y.hat")

  ## Return results
  return(results)
}