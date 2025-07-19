\name{NVC_SSL}
\alias{NVC_SSL}

\title{Nonparametric Varying Coefficient Spike-and-Slab Lasso (NVC-SSL)}

\description{This function implements the Nonparametric Varying Coefficient Spike-and-Slab Lasso (NVC-SSL) model of Bai et al. (2023) for high-dimensional NVC models. The function returns the MAP estimator for the varying coefficients \eqn{\beta_k(t), k = 1, ..., p,} obtained from the ECM algorithm described in Bai et al. (2023). The BIC criterion is used to select the optimal spike hyperparameter \eqn{lambda_0}.  

If the user specifies \code{return_CI=TRUE}, then this function will also return the 95 percent pointwise posterior credible intervals for the varying coefficients \eqn{\beta_k(t), k = 1, ..., p,} obtained from Gibbs sampling. If the number of covariates \eqn{p} is large, then the user can additionally use the approximate MCMC algorithm introduced in Bai et al. (2023) (\code{approx_MCMC=TRUE}) which is much faster than the exact Gibbs sampler and gives higher simultaneous coverage.

Finally, this function returns the number of iterations and the runtime for the ECM algorithms and MCMC algorithms which can be used for benchmarking and timing comparisons.}

\usage{
NVC_SSL(y, t, id, X, n_basis=8, 
        lambda0=seq(from=300,to=10,by=-10), lambda1=1, 
        a=1, b=ncol(X), c0=1, d0=1, nu=n_basis+2, Phi=diag(n_basis),
        include_intercept=TRUE, tol=1e-6, max_iter=100, 
        return_CI=FALSE, approx_MCMC=FALSE,
        n_samples=1500, burn=500, print_iter=TRUE)
}

\arguments{
\item{y}{\eqn{N \times 1} vector of response observations \eqn{y_{11}, ..., y_{1m_1}, ..., y_{n1}, ..., y_{nm_n}}}
\item{t}{\eqn{N \times 1} vector of observation times \eqn{t_{11}, ..., t_{1m_1}, ..., t_{n1}, ..., t_{nm_n}}}
\item{id}{\eqn{N \times 1} vector of labels, where each unique label corresponds to one of the subjects}
\item{X}{\eqn{N \times p} design matrix with columns \eqn{[X_1, ..., X_p]}, where the \eqn{k}th column contains the entries \eqn{x_{ik}(t_{ij})}'s}
\item{n_basis}{number of basis functions to use. Default is \code{n_basis=8}.}
\item{lambda0}{grid of spike hyperparameters. Default is to tune \code{lambda0} from the grid of decreasing values \eqn{(300, 290, ..., 20, 10)}.} 
\item{lambda1}{slab hyperparameter. Default is \code{lambda1=1}.}
\item{a}{hyperparameter in \eqn{B(a,b)} prior on mixing proportion \eqn{\theta}. Default is \eqn{a=1}.}
\item{b}{hyperparameter in \eqn{B(a,b)} prior on mixing proportion \eqn{\theta}. Default is \eqn{b=p}.}
\item{c0}{hyperparameter in Inverse-Gamma\eqn{(c_0/2,d_0/2)} prior on measurement error variance \eqn{\sigma^2}. Default is \eqn{c_0=1}.}
\item{d0}{hyperparameter in Inverse-Gamma\eqn{(c_0/2,d_0/2)} prior on measurement error variance \eqn{\sigma^2}. Default is \eqn{d_0=1}.}
\item{nu}{degrees of freedom for Inverse-Wishart prior on \eqn{\Omega}. Default is \code{n_basis+2}.}
\item{Phi}{scale matrix in the Inverse-Wishart prior on \eqn{\Omega}. Default is the identity matrix.}
\item{include_intercept}{Boolean variable for whether or not to include an intercept function \eqn{\beta_0(t)} in the estimation. Default is \code{include_intercept=TRUE}.}
\item{tol}{convergence criteria for the ECM algorithm. Default is \code{tol=1e-6}.}
\item{max_iter}{maximum number of iterations to run ECM algorithm. Default is \code{max_iter=100}.}
\item{return_CI}{Boolean variable for whether or not to return the 95 percent pointwise credible bands. Set \code{return_CI=TRUE} if posterior credible bands are desired.}
\item{approx_MCMC}{Boolean variable for whether or not to run the approximate MCMC algorithm instead of the exact MCMC algorithm. If \code{approx_MCMC=TRUE}, then an approximate MCMC algorithm isused. Otherwise, if \code{approx_MCMC=FALSE}, the exact MCMC algorithm is used. This argument is ignored if \code{return_CI=FALSE}.}
\item{n_samples}{number of MCMC samples to save for posterior inference. The default is to save \code{n_samples=1500}. This is ignored if \code{return_CI=FALSE}.}
\item{burn}{number of initial MCMC samples to discard during the warm-up period. Default is \code{burn=500}. This is ignored if \code{return_CI=FALSE}.}
\item{print_iter}{Boolean variable for whether or not to print the progress of the algorithms. Default is \code{print_iter=TRUE}.}
}

\value{The function returns a list containing the following components:
\item{t_ordered}{all \eqn{N} time points ordered from smallest to largest. Needed for plotting}
\item{classifications}{\eqn{p \times 1} vector of indicator variables, where \code{"1"} indicates that the covariate is selected and \code{"0"} indicates that it is not selected. These classifications are determined by the optimal \code{lambda0} chosen from BIC. Note that this vector does not include an intercept function.}
\item{beta_hat}{\eqn{N \times p} matrix of the MAP estimates for varying coefficient functions \eqn{\beta_k(t), k = 1, ..., p}, using the optimal \code{lambda0} chosen from BIC. The \eqn{k}th column in the matrix is the \eqn{k}th estimated function at the observation times in \code{t_ordered}.}
\item{beta0_hat}{MAP estimate of the intercept function \eqn{\beta_0(t)} at the observation times in \code{t_ordered} for the optimal \code{lambda0} chosen from BIC. This is not returned if \code{include_intercept = FALSE}.}
\item{gamma_hat}{MAP estimates of the basis coefficients (needed for prediction) for the optimal \code{lambda0}.}
\item{beta_post_mean}{\eqn{N \times p} matrix of the posterior mean estimates of the varying coefficient functions. The \eqn{k}th column in the matrix is the \eqn{k}th posterior mean estimate for \eqn{\beta_k(t)} at the observation times in \code{t_ordered}. This is not returned if \code{return_CI=FALSE}.}
\item{beta_CI_lower}{\eqn{N \times p} matrix of the lower endpoints of the 95 percent pointwise posterior credible interval (CI) for the varying coefficient functions. The \eqn{k}th column in the matrix is the lower endpoint for the CI of \eqn{\beta_k(t)} at the observation times in \code{t_ordered}. This is not returned if \code{return_CI=FALSE}.}
\item{beta_CI_upper}{\eqn{N \times p} matrix of the upper endpoints of the 95 percent pointwise posterior credible interval (CI) for the varying coefficient functions. The \eqn{k}th column in the matrix is the upper endpoint for the CI of \eqn{\beta_k(t)} at the observation times in \code{t_ordered}. This is not returned if \code{return_CI=FALSE}.}
\item{beta0_post_mean}{Posterior mean estimate of the intercept function \eqn{\beta_0(t)} at the observation times in \code{t_ordered}. This is not returned if \code{return_CI=FALSE}.} 
\item{beta0_CI_lower}{Lower endpoints of the 95 percent pointwise posterior credible intervals (CIs) for the intercept function \eqn{\beta_0(t)} at the observation times in \code{t_ordered}. This is not returned if \code{return_CI=FALSE}.} 
\item{beta0_CI_upper}{Upper endpoints of the 95 percent pointwise posterior credible intervals (CIs) for the intercept function \eqn{\beta_0(t)} at the observation times in \code{t_ordered}. This is not returned if \code{return_CI=FALSE}.} 
\item{gamma_post_mean}{Posterior mean estimates of all the basis coefficients. This is not returned if \code{return_CI=FALSE}.}
\item{gamma_CI_lower}{Lower endpoints of the 95 percent posterior credible intervals for the basis coefficients. This is not returned if \code{return_CI=FALSE}.}
\item{gamma_CI_upper}{Upper endpoints of the 95 percent posterior credible intervals for the basis coefficients. This is not returned if \code{return_CI=FALSE}.}
\item{post_incl}{\eqn{p \times 1} vector of estimated posterior inclusion probabilities (PIPs) for each of the varying coefficients. The \eqn{k}th entry in \code{post_incl} is the PIP for \eqn{\beta_k}. This is not returned if \code{return_CI=FALSE}.}
\item{lambda0_min}{the individual \code{lambda0} which minimizes the BIC. If only one value was originally passed for \code{lambda0}, then this just returns that \code{lambda0}.}
\item{lambda0_all}{grid of all \eqn{L} regularization parameters in \code{lambda0}.  Note that since the objective function is scaled by \eqn{1/N} for the penalized frequentist methods in the \code{NVC_frequentist} function, the \code{lambda_all} grid that is chosen automatically by \code{NVC_frequentist} typically consists of smaller values than the default values in the \code{lambda0_all} grid for \code{NVC_SSL}.} 
\item{BIC_all}{\eqn{L \times 1} vector of BIC values corresponding to all \eqn{L} entries in \code{lambda0_all}. The \eqn{l}th entry corresponds to the \eqn{l}th entry in \code{lambda0_all}.}
\item{beta_est_all_lambda0}{list of length \eqn{L} of the estimated varying coefficients \eqn{\beta_k(t), k = 1, ..., p}, corresponding to all \eqn{L} lambdas in \code{lambda0_all}. The \eqn{l}th entry corresponds to the \eqn{l}th entry in \code{lambda0_all}.}
\item{beta0_est_all_lambda0}{\eqn{N \times L} matrix of estimated intercept function \eqn{\beta_0(t)} corresponding to all \eqn{L} entries in \code{lambda0_all}. The \eqn{l}th column corresponds to the \eqn{l}th entry in \code{lambda0_all}. This is not returned if \code{include_intercept=FALSE}.}
\item{gamma_est_all_lambda0}{\eqn{dp \times L} matrix of estimated basis coefficients corresponding to all entries in \code{lambda0_all}. The \eqn{l}th column corresponds to the \eqn{l}th entry in \code{lambda0_all}.}
\item{classifications_all_lambda0}{\eqn{p \times L} matrix of classifications corresponding to all the entries in \code{lambda0_all}. The \eqn{l}th column corresponds to the \eqn{l}th entry in \code{lambda0_all}.}
\item{ECM_iters_to_converge}{number of iterations it took for the ECM algorithm to converge for each entry in \code{lambda0_all}. The \eqn{l}th entry corresponds to the \eqn{l}th entry in \code{lambda0_all}.}
\item{ECM_runtimes}{\eqn{L \times 1} vector of the number of seconds it took for the ECM algorithm to converge for each entry in \code{lambda0_all}. The \eqn{l}th entry corresponds to the \eqn{l}th entry in \code{lambda0_all}.}
\item{gibbs_runtime}{number of minutes it took for the Gibbs sampling algorithm to run for the total number of MCMC iterations given in \code{gibbs_iters}}
\item{gibbs_iters}{total number of MCMC iterations run for posterior inference}
}

\references{
Bai, R., Boland, M. R., and Chen, Y. (2023). "Scalable high-dimensional Bayesian varying coefficient models with unknown within-subject covariance." \emph{Journala of Machine Learning Research}, \bold{24}:1-49.

Bai, R., Moran, G. E., Antonelli, J. L., Chen, Y., and Boland, M.R. (2022). "Spike-and-slab group lassos for grouped regression and sparse generalized additive models." \emph{Journal of the American Statistical Association}, \bold{117}:184-197.
}

\examples{
## Load data
data(SimulatedData)
attach(SimulatedData)
y = SimulatedData$y
t = SimulatedData$t
id = SimulatedData$id
X = SimulatedData[,4:103]

## Fit NVC-SSL model. Default implementation uses a grid of 30 lambdas.
## Below illustration uses just two well-chosen lambdas

NVC_SSL_mod = NVC_SSL(y, t, id, X, lambda0=c(60,50))

## NOTE: Should use default, which will search for lambda0 from a bigger grid
# NVC_SSL_mod = NVC_SSL(y, t, id, X)

## Classifications. First 6 varying coefficients are selected as nonzero
NVC_SSL_mod$classifications

## Optimal lambda chosen from BIC
NVC_SSL_mod$lambda0_min

## Plot first estimated varying coefficient function
t_ordered = NVC_SSL_mod$t_ordered
beta_hat= NVC_SSL_mod$beta_hat
plot(t_ordered, beta_hat[,1], lwd=3, type='l', col='blue',
     xlab="Time", ylim = c(-12,12), ylab=expression(beta[1]))

## Plot third estimated varying coefficient function
plot(t_ordered, beta_hat[,3], lwd=3, type='l', col='blue',
     xlab="Time", ylim = c(-4,2), ylab=expression(beta[3]))

## Plot fifth estimated varying coefficient function
plot(t_ordered, beta_hat[,5], lwd=3, type='l', col='blue',
     xlab="Time", ylim = c(0,15), ylab=expression(beta[5]))

\donttest{
## If you want credible intervals, then set return_CI=TRUE to also run Gibbs sampler.
## Below, we run a total of 1000 MCMC iterations, discarding the first 500 as burnin
## and keeping the final 500 samples for inference.

NVC_SSL_mod_2 = NVC_SSL(y, t, id, X, return_CI=TRUE, approx_MCMC=FALSE, 
                        n_samples=500, burn=500)

## Note that NVC_SSL() always computes a MAP estimator first and then
## initializes the Gibbs sampler with the MAP estimator.

## Plot third varying coefficient function and its credible bands
t_ordered = NVC_SSL_mod_2$t_ordered
beta_MAP = NVC_SSL_mod_2$beta_hat
beta_mean = NVC_SSL_mod_2$beta_post_mean
beta_CI_lower = NVC_SSL_mod_2$beta_CI_lower
beta_CI_upper = NVC_SSL_mod_2$beta_CI_upper

plot(t_ordered, beta_MAP[,3], lwd=3, type='l', col='blue', xlab="Time", ylim=c(-5,3), lty=1, 
     ylab=expression(beta[3]), cex.lab=1.5)
lines(t_ordered, beta_mean[,3], lwd=3, type='l', col='red', lty=4)
lines(t_ordered, beta_CI_lower[,3], lwd=4, type='l', col='purple', lty=3)
lines(t_ordered, beta_CI_upper[,3], lwd=4, type='l', col='purple', lty=3)
legend("bottomleft", c("MAP", "Mean", "95 percent CI"), lty=c(1,4,3), lwd=c(2,2,3), 
       col=c("blue","red","purple"), inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n")

## Plot fifth varying coefficient function and its credible bands
plot(t_ordered, beta_MAP[,5], lwd=3, type='l', col='blue', xlab="Time", ylim=c(-1,14), lty=1, 
     ylab=expression(beta[5]), cex.lab=1.5)
lines(t_ordered, beta_mean[,5], lwd=3, type='l', col='red', lty=4)
lines(t_ordered, beta_CI_lower[,5], lwd=4, type='l', col='purple', lty=3)
lines(t_ordered, beta_CI_upper[,5], lwd=4, type='l', col='purple', lty=3)
legend("bottomleft", c("MAP", "Mean", "95 percent CI"), lty=c(1,4,3), lwd=c(2,2,3), 
       col=c("blue","red","purple"), inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n")
}
}