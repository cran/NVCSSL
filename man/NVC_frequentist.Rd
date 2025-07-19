\name{NVC_frequentist}
\alias{NVC_frequentist}

\title{Fits frequentist penalized nonparametric varying coefficient (NVC) models}

\description{This function implements frequentist penalized nonparametric varying coefficient (NVC) models. It supports the following penalty functions: the group lasso penalty of Yuan and Lin (2006), the group minimax concave penalty (MCP) of Breheny and Huang (2015), and the group smoothly clipped absolute deviation (SCAD) penalty of Breheny and Huang (2015). This function solves a penalized regression problem of the form,

\deqn{ argmax_{\gamma} \frac{1}{N} \ell(\gamma) + pen_{\lambda}(\gamma),}
 
 where \eqn{N} is the total number of observations, \eqn{\ell(\gamma)} is the loss function, and \eqn{pen_{\lambda}(\cdot)} is a penalty function with regularization parameter \eqn{\lambda > 0}. Since the objective function is rescaled by \eqn{1/N}, the penalty \eqn{\lambda} is typically smaller than the spike hyperparameter \eqn{\lambda_0} used by the \code{NVC_SSL} function. The BIC criterion is used to select the optimal tuning parameter \eqn{\lambda}.
}

\usage{
NVC_frequentist(y, t, X, n_basis=8, penalty=c("gLASSO","gSCAD","gMCP"),
                lambda=NULL, include_intercept=TRUE) 
}

\arguments{
\item{y}{\eqn{N \times 1} vector of response observations \eqn{y_{11}, ..., y_{1m_1}, ..., y_{n1}, ..., y_{nm_n}}}
\item{t}{\eqn{N \times 1} vector of observation times \eqn{t_{11}, ..., t_{1m_1}, ..., t_{n1}, ..., t_{nm_n}}}
\item{X}{\eqn{N \times p} design matrix with columns \eqn{[X_1, ..., X_p]}, where the \eqn{k}th column contains the entries \eqn{x_{ik}(t_{ij})}'s}
\item{n_basis}{number of basis functions to use. Default is \code{n_basis=8}.}
\item{penalty}{string specifying which penalty function to use. Specify \code{"gLASSO"} for group lasso, \code{"gSCAD"} for group SCAD, or \code{"gMCP"} for group MCP.}
\item{lambda}{grid of tuning parameters. If \code{lambda} is not specified (i.e. \code{lambda=NULL}), then the program automatically chooses a grid for \code{lambda}. Note that since the objective function is scaled by \eqn{1/N}, the automatically chosen grid for \code{lambda} typically consists of smaller values than the default grid for \code{lambda0} used by the function \code{NVC_SSL}.}
\item{include_intercept}{Boolean variable for whether or not to include an intercept function \eqn{\beta_0(t)} in the estimation. Default is \code{include_intercept=TRUE}.}
}

\value{The function returns a list containing the following components:
\item{t_ordered}{all \eqn{N} time points ordered from smallest to largest. Needed for plotting.}
\item{classifications}{\eqn{p \times 1} vector of indicator variables, where \code{"1"} indicates that the covariate is selected and \code{"0"} indicates that it is not selected. These classifications are determined by the optimal \code{lambda} chosen from BIC. Note that this vector does not include an intercept function.}
\item{beta_hat}{\eqn{N \times p} matrix of the estimates for varying coefficient functions \eqn{\beta_k(t), k = 1, ..., p}, using the optimal \code{lambda} chosen from BIC. The \eqn{k}th column in the matrix is the \eqn{k}th estimated function at the observation times in \code{t_ordered}.}
\item{beta0_hat}{estmate of the intercept function \eqn{\beta_0(t)} at the observation times in \code{t_ordered} for the optimal \code{lambda} chosen from BIC. This is not returned if \code{include_intercept = FALSE}.}
\item{gamma_hat}{estimated basis coefficients (needed for prediction) for the optimal \code{lambda}.}
\item{lambda_min}{the individual \code{lambda} which minimizes the BIC. If only one value was originally passed for \code{lambda}, then this just returns that \code{lambda}.}
\item{lambda0_all}{grid of all \eqn{L} regularization parameters in \code{lambda}.  Note that since the objective function is scaled by \eqn{1/N} for the penalized frequentist methods in the \code{NVC_frequentist} function, the \code{lambda_all} grid that is chosen automatically by \code{NVC_frequentist} typically consists of smaller values than the default values in the \code{lambda0_all} grid for \code{NVC_SSL}.} 
\item{BIC_all}{\eqn{L \times 1} vector of BIC values corresponding to all \eqn{L} entries in \code{lambda_all}. The \eqn{l}th entry corresponds to the \eqn{l}th entry in \code{lambda_all}.}
\item{beta_est_all_lambda}{list of length \eqn{L} of the estimated varying coefficients \eqn{\beta_k(t), k = 1, ..., p}, corresponding to all \eqn{L} lambdas in \code{lambda_all}. The \eqn{l}th entry corresponds to the \eqn{l}th entry in \code{lambda_all}.}
\item{beta0_est_all_lambda}{\eqn{N \times L} matrix of estimated intercept function \eqn{\beta_0(t)} corresponding to all \eqn{L} entries in \code{lambda_all}. The \eqn{l}th column corresponds to the \eqn{l}th entry in \code{lambda_all}. This is not returned if \code{include_intercept=FALSE}.}
\item{gamma_est_all_lambda}{\eqn{dp \times L} matrix of estimated basis coefficients corresponding to all entries in \code{lambda_all}. The \eqn{l}th column corresponds to the \eqn{l}th entry in \code{lambda_all}.}
\item{classifications_all_lambda}{\eqn{p \times L} matrix of classifications corresponding to all the entries in \code{lambda_all}. The \eqn{l}th column corresponds to the \eqn{l}th entry in \code{lambda_all}.}
\item{iters_to_converge}{number of iterations it took for the group ascent algorithm to converge for each entry in \code{lambda_all}. The \eqn{l}th entry corresponds to the \eqn{l}th entry in \code{lambda_all}.}
}

\references{
Bai, R., Boland, M. R., and Chen, Y. (2023). "Scalable high-dimensional Bayesian varying coefficient models with unknown within-subject covariance." \emph{Journal of Machine Learning Research}, \bold{24}:1-49.

Breheny, P. and Huang, J. (2015). "Group descent algorithms for nonconvex penalized linear and logistic regression models with grouped predictors." \emph{Statistics and Computing}, \bold{25}:173-187.
  
Wei, F., Huang, J., and Li, H. (2011). "Variable selection and estimation in high-dimensional varying coefficient models." \emph{Statistica Sinica}, \bold{21}:1515-1540.  

Yuan, M. and Lin, Y. (2006). "Model selection and estimation in regression with grouped variables." \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{68}:49-67.
  
}

\examples{
## Load data
data(SimulatedData)
attach(SimulatedData)
y = SimulatedData$y
t = SimulatedData$t
id = SimulatedData$id
X = SimulatedData[,4:103]

## Fit frequentist penalized NVC model with the SCAD penalty. 
## Can set penalty as "gLASSO", "gSCAD", or "gMCP". 
## No need to specify an 'id' argument when using NVC_frequentist() function

NVC_gSCAD_mod = NVC_frequentist(y, t, X, penalty="gSCAD")

## Classifications. First varying coefficients are selected as nonzero
NVC_gSCAD_mod$classifications

## Optimal lambda chosen from BIC
NVC_gSCAD_mod$lambda_min

## Plot first estimated varying coefficient function
t_ordered = NVC_gSCAD_mod$t_ordered
beta_hat= NVC_gSCAD_mod$beta_hat
plot(t_ordered, beta_hat[,1], lwd=3, type='l', col='blue',
     xlab="Time", ylim = c(-12,12), ylab=expression(beta[1]))

## Plot third estimated varying coefficient function
plot(t_ordered, beta_hat[,3], lwd=3, type='l', col='blue',
     xlab="Time", ylim = c(-4,2), ylab=expression(beta[3]))

## Plot fifth estimated varying coefficient function
plot(t_ordered, beta_hat[,5], lwd=3, type='l', col='blue',
     xlab="Time", ylim = c(0,15), ylab=expression(beta[5]))
}