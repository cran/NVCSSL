\name{NVC_SSL}
\alias{NVC_SSL}

\title{Nonparametric Varying Coefficient Spike-and-Slab Lasso (NVC-SSL)}

\description{This function implements the nonparametric varying coefficient spike-and-slab lasso (NVC-SSL) model of Bai et al. (2020) for repeated measures data. The NVC-SSL model aims to estimate the functions of time, \eqn{\beta_1(t), \ldots, \beta_p(t)}, in the model,
\deqn{y_i(t_{ij}) = \sum_{k=1}^{p} x_{ik}(t_{ij}) \beta_k(t_{ij}) + \varepsilon_{ij},}
where the subscript \eqn{ij} refers to the \eqn{j}th observation for the \eqn{i}th subject.

The NVC-SSL model can model within-subject temporal correlations, and the data does \emph{not} need to be regularly spaced. In particular, we have \eqn{n} subjects with \eqn{n_i, i=1, \ldots, n} within-study observations each, for a total of \eqn{N=\sum_{i=1}^{n} n_i} observations. 
}

\usage{NVC_SSL(y, t, id, X, df=8, cov.structure=c("AR1","CS","ind","unstructured"),
               lambda0=c(seq(from=5,to=30,by=5),seq(from=40,to=100,by=10)), lambda1=1, 
               a=1, b=ncol(X), c0=1, d0=1, rho.support=seq(0,0.9,by=0.1), 
               tol = 1e-6, print.iteration=TRUE)}
\arguments{
  \item{y}{\eqn{N \times 1} vector of all observed responses, \eqn{(y_{11}, \ldots, y_{1 n_1}, \ldots, y_{n1}, \ldots, y_{n n_n})'}.}
  \item{t}{\eqn{N \times 1} vector of all time points, \eqn{(t_{11}, \ldots, t_{1 n_1}, \ldots, t_{n1}, \ldots, t_{n n_n})'}.}
  \item{id}{\eqn{N \times 1} vector of subject identifiers (\eqn{1} through \eqn{n}).}
  \item{X}{\eqn{N \times p} design matrix of possibly time-varying covariates. The \eqn{k}th column of \eqn{X} is \eqn{(x_{1k}(t_{11}), \ldots, x_{1k}(t_{1 n_1}), \ldots, x_{nk}(t_{n1}), \ldots, x_{nk}(t_{n n_n}))'}.}
  \item{df}{number of basis functions to use in each basis expansion. Default is \code{df=8}, but the user may specify degrees of freedom as any integer greater than or equal to 3. }
  \item{cov.structure}{within-subject covariance structure with which to fit the NVC-SSL model. The options are: \code{"AR1"} for first-order autoregressive, \code{"CS"} for compound symmetry, \code{"ind"} for independent errors (i.e. no within-subject temporal correlations), and \code{"unstructured"} for completely unstructured. The default is \code{"AR1"}. Note that \code{"unstructured"} is generally not recommended. Instead, we recommend using the \code{robustified_NVC_SSL} function when the covariance structure may be misspecified.}
  \item{lambda0}{ladder of spike hyperparameters in the spike-and-slab group lasso (SSGL) prior for dynamic posterior exploration. The user may specify either a scalar or a vector. The default is \code{lambda0=c(5,10,15,20,25,30,40,50,60,70,80,90,100)}. }
  \item{lambda1}{slab hyperparameter in the SSGL prior, which is held fixed. The default is \code{lambda1=1}.}
  \item{a}{shape hyperparameter for the \eqn{B(a,b)} prior on the mixing proportion \eqn{\theta}. Default is \code{a=1}.}
  \item{b}{shape hyperparameter for the \eqn{B(a,b)} prior on the mixing proportion \eqn{\theta}. Default is \code{b=ncol(X)}.}
  \item{c0}{shape hyperparameter for the \eqn{IG(c_0/2,d_0/2)} prior on the unknown global variance \eqn{\sigma^2}. Default is \code{c0=1}. Ignored if the covariance structure is unstructured, i.e. \code{cov.structure="unstructured"}. }
  \item{d0}{rate hyperparameter for the \eqn{IG(c_0/2,d_0/2)} prior on the unknown global variance \eqn{\sigma^2}. Default is \code{d0=1}. Ignored if the covariance structure is unstructured, i.e. \code{cov.structure="unstructured"}. }
  \item{rho.support}{the support for the discrete uniform prior on the autocorrelation parameter \eqn{\rho}. Default is \code{rho.support=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)}, but the user may specify a finer grid (as long as the atoms are greater than or equal to 0 and strictly less than 1). This is ignored if \code{cov.structure="ind"} or \code{cov.structure="unstructured"}.}
\item{tol}{convergence threshold for the EM algorithm. Default is \code{tol=1e-6}.}
\item{print.iteration}{flag for printing the current value of \code{lambda0} in the ladder of spike hyperparameters. Default is \code{TRUE}.}
}

\value{The function returns a list containing the following components:
\item{t.ordered}{all \eqn{N} time points in order from smallest to largest. Needed for plotting.}
\item{beta.hat}{\eqn{N \times p} matrix of the function estimates. The \eqn{k}th column is the function estimate \eqn{\beta_k(t_{ij})} evaluated at the \eqn{N} time observations \eqn{t_{ij}, i=1,...,n, j=1,...,n_i}.}
\item{gamma.hat}{\eqn{df \times p} estimated basis coefficients \eqn{\gamma}. Needed for prediction.}
\item{intercept}{Estimate of an intercept (or grand mean). Needed for prediction.}
\item{classifications}{\eqn{p \times 1} vector of binary variables for the \eqn{p} covariates. "1" indicates that the covariate was selected, and "0" indicates that it was not selected.}
\item{AICc}{AIC with correction for small sample size. This value can be used for tuning the degrees of freedom \code{df}. The user should pick the model which obtains the lowest \eqn{AIC_c}.}
}

\references{ 
Bai, R., Boland, M. R., and Chen, Y. (2020). "Fast algorithms and theory for high-dimensional Bayesian varying coefficient models." \emph{arXiv preprint arXiv:1907.06477}.
}

\examples{

#####################################
# Example on the synthetic data set #
#####################################
data(SimulatedData)
attach(SimulatedData)
y = SimulatedData$y
t = SimulatedData$t
id = SimulatedData$id
X = SimulatedData[,4:103]

## Training set. First 20 subjects
training.indices = which(id <= 20)
y.train = y[training.indices]
t.train = t[training.indices]
id.train = id[training.indices]
X.train = X[training.indices,]

## Test set. Last 10 subjects
test.indices = which(id > 20)
y.test = y[test.indices]
t.test = t[test.indices]
id.test = id[test.indices]
X.test = X[test.indices, ]

## Set seed
set.seed(123)

###############################
## NVC-SSL with fixed lambda0 #
###############################
# We use with 6 degrees of freedom and FIXED lambda0=15,
# so we are estimating 600 total unknown basis coefficients

nvcssl.mod = NVC_SSL(y=y.train, t=t.train, id=id.train, X=X.train, 
                     df=6, cov.structure="AR1", lambda0=15)

## The first 6 functionals were classified as nonzero
nvcssl.mod$classifications 

## AIC with correction for sample size. Can use this to tune df
nvcssl.mod$AICc

## Make plot of the 6 active functionals

oldpar <- par(no.readonly = TRUE)

t.ordered = nvcssl.mod$t.ordered
beta.hat = nvcssl.mod$beta.hat

par(mfrow=c(3,2), mar=c(4,4,4,4))
plot(t.ordered, beta.hat[,1], lwd=3, type='l', col='blue', 
     xlab="T", ylim=c(-10,10), ylab=expression(beta[1]))
plot(t.ordered, beta.hat[,2], lwd=3, type='l', col='blue', 
     xlab="T", ylim=c(-8,8), ylab=expression(beta[2]))
plot(t.ordered, beta.hat[,3], lwd=3, type='l', col='blue', 
     xlab="T", ylim=c(-7,4), ylab=expression(beta[3]))
plot(t.ordered, beta.hat[,4], lwd=3, type='l', col='blue', 
     xlab="T", ylim=c(-6,4), ylab=expression(beta[4]))
plot(t.ordered, beta.hat[,5], lwd=3, type='l', col='blue', 
     xlab="T", ylim=c(0,15), ylab=expression(beta[5]))
plot(t.ordered, beta.hat[,6], lwd=3, type='l', col='blue', 
     xlab="T", ylim=c(-5,0), ylab=expression(beta[6]))

par(oldpar)

## Prediction on test set
y.preds = predict_NVC(nvcssl.mod, t.new=t.test, id.new=id.test, X.new=X.test)

## MSPE
nvcssl.mod.mspe = mean((y.test-y.preds[,2])^2)
nvcssl.mod.mspe


\donttest{
##############################
# NVC-SSL model with dynamic #
# posterior exploration      #
##############################

## Now we fit the NVC-SSL mod with default 8 degrees of freedom 
## (800 unknown basis coefficients), AR(1) covariance structure,
## and dynamic posterior exploration for spike hyperparameter lambda0.

## Dynamic posterior exploration will generally produce a better fit
## than arbitrarily fixing lambda0.

## Note if the standard deviation of y is very small, the algorithm
## may be initially slow for values of lambda0 that are too small. 

## In this case, consider increasing the lambda0 ladder to begin 
## at a larger value.

data(SimulatedData)
attach(SimulatedData)
y = SimulatedData$y
t = SimulatedData$t
id = SimulatedData$id
X = SimulatedData[,4:103]

## Fit NVC-SSL model with dynamic posterior exploration
nvcssl.mod = NVC_SSL(y=y, t=t, id=id, X=X, cov.structure="AR1")

## Look at which covariates were selected
nvcssl.mod$classifications
}

}